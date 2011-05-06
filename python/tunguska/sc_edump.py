
import sys, calendar, time, os, logging
from os.path import join as pjoin
from pyrocko import trace, model

logger = logging.getLogger('tunguska.sc_edump')

try:
    
    import seiscomp3.Client
    import seiscomp3.Communication
    import seiscomp3.DataModel
    from seiscomp3.Core import TimeSpan
            
    def get_nsl(x):
        return x.network, x.station, x.location
            

    def gmctime(t):
        return time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(t))
    

    def pyrocko_trace(network_code, station, location, stream):
        return trace.Trace(  network=network_code,
                                station=station.code(),
                                location=location.code(),
                                channel=stream.code(),
                                deltat=stream.sampleRateDenominator()/stream.sampleRateNumerator(),
                                tmin=0.,tmax=0.)
            
    def pyrocko_station(network_code, station, location):
        return  model.Station(network=network_code,
                                station=station.code(),
                                location=location.code(),
                                lat=location.latitude(),
                                lon=location.longitude(),
                                elevation=location.elevation())
    
    
    def sctime(instant):
        return seiscomp3.Core.Time( *time.gmtime(instant)[:6] )
    
    def fill_stream_tmpl( template, network, station, location, channel, **kwargs ):
    
        d = dict( network=network,
                station=station,
                location=location,
                channel=channel )
                
        d.update(kwargs)
        return template % d
    
    def ensuredirs(dst):
        d,x = os.path.split(dst)
        dirs = []
        while d and not os.path.exists(d):
            dirs.append(d)
            d,x = os.path.split(d)
            
        dirs.reverse()
        
        for d in dirs:
            if not os.path.exists(d):
                os.mkdir(d)
    
    class ResponseUnavailable(Exception):
        pass
    
    class EventDumper(seiscomp3.Client.Application):
        
        def __init__(self, argv,
                trace_selector = lambda tr: True,
                station_selector = lambda sta: True,
                station_weeder = lambda stations: stations,
                redundant_channel_priorities = None, # use something like this to activate:
                                                     # [ (('BH1', 'BH2'), ('BHN', 'BHE')),
                                                     #   (('HH1', 'HH2'), ('HHN', 'HHE')) ]
                skip_responseless = False,
                time_range = (-100,100),
                event_path = 'events/%(event_name)s',
                events = []):
            
            seiscomp3.Client.Application.__init__(self, len(argv), argv)
                    
            subscriptions = []
            
            self.setMessagingEnabled(True)
            self.setDatabaseEnabled(True, True)
            self.setLoadConfigModuleEnabled(True)
            self.setLoadInventoryEnabled(True)
            self.setRecordStreamEnabled(True)
            
            st_nslc = '%(network)s_%(station)s_%(location)s_%(channel)s'
            
            self._trace_selector = trace_selector
            self._station_selector = station_selector
            self._station_weeder = station_weeder
            self._redundant_channel_priorities = redundant_channel_priorities
            self._timewindow = time_range
            self._event_path_template = event_path
            self._skip_responseless = skip_responseless
            
            self._mseed_fn_template = 'raw-%s.mseed' % st_nslc
            self._component_fn_template = 'component-%s.txt' % st_nslc
            self._polezero_fn_template = 'polezero-%s.txt' % st_nslc
            self._station_fn_template = 'stations.txt'
            self._fseed_fn_template = 'dataless.seed'
            self._event_fn_template = 'event.txt'
            
            self.enableTimer(1)
            self._events = events
            
        def handleTimeout(self):
            for event in self._events:
                self.dumpEvent(event)
                
            self.enableTimer(0)
            self.quit()
        
        def get_response(self, network_code, station, location, stream):
            
            inventory = seiscomp3.Client.Inventory_Instance().inventory()
            
            try:
                sensor = inventory.findSensor(stream.sensor())
                datalogger = inventory.findDatalogger(stream.datalogger())
                
                if sensor is None or datalogger is None:
                    raise ResponseUnavailable( 'response information incomplete for stream %s.%s.%s.%s' % 
                    (network_code, station.code(), location.code(), stream.code()) )
                
                if sensor.response().startswith('ResponsePAZ'):
                    resppaz = inventory.findResponsePAZ(sensor.response())
                else:
                    raise ResponseUnavailable('response information for stream %s.%s.%s.%s is in an unexpected format: "%s"' % 
                    (network_code, station.code(), location.code(), stream.code(),
                    seismometer.response()[:4]) )
                    
                return datalogger, resppaz
            
            except KeyError:
                raise ResponseUnavailable( 'response information incomplete for stream %s.%s.%s.%s' % 
                    (network_code, station.code(), location.code(), stream.code()) )
    
        def _redundant_channel_weeder(self, nslcs):
            
            if self._redundant_channel_priorities is None: return []
            
            # n=network,s=station,l=location,c=channel
            # channels by station
            by_nsl = {}
            for nslc in nslcs:
                nsl = nslc[:3]
                if nsl not in by_nsl:
                    by_nsl[nsl] = []
                
                by_nsl[nsl].append( nslc )
            
            # figure out what items to remove
            to_delete = []
            for ((h1, h2), (l1, l2)) in self._redundant_channel_priorities:
                for nsl, nslcs in by_nsl.iteritems():
                    channels = [ nslc[3] for nslc in nslcs ]
                    if h1 in channels and h2 in channels and l1 in channels and l2 in channels:
                        to_delete.append(nslc[:3] + (l1,))
                        to_delete.append(nslc[:3] + (l2,))
            
            return to_delete
                    
        def iterNSLS(self, anno):
            '''Iterate over possibly available streams.
            
            Yields tuples like
            
                (network_code, station, location, stream)
                
            '''
            
            inventory = seiscomp3.Client.Inventory_Instance()
            cm = self.configModule()
            
            for i in xrange( cm.configStationCount() ):
                cs = cm.configStation(i)
                if cs.enabled():
                    station = inventory.getStation(cs.networkCode(), cs.stationCode(), anno)
                    if station:
                        for ilocation in xrange(station.sensorLocationCount()):
                            location = station.sensorLocation(ilocation)
                            
                            for istream in xrange(location.streamCount()):
                                stream = location.stream(istream)
                                yield (cs.networkCode(), station, location, stream)
                                
        def dumpEvent(self, event):
                                    
            phase_counter = {}
            timewindows = {}
            streams = {}
            responses = {}
            used_stations = {}
            pyr_stations = {}
            
            # gather 
            for x in self.iterNSLS(sctime(event.time)):
                network_code, station, location, stream = x
                nslc = network_code, station.code(), location.code(), stream.code()
            
                tr = pyrocko_trace(network_code, station, location, stream)
                sta = pyrocko_station(network_code, station, location)
                sta.set_event_relative_data(event)
                
                if not self._trace_selector(tr):
                    continue
                
                if not self._station_selector(sta):
                    continue
                
                tmin = self._timewindow[0]
                tmax = self._timewindow[1]
                
                nsl = nslc[:3]
                
                try:
                    responses[nslc] = self.get_response(network_code, station, location, stream)
                except ResponseUnavailable, e:
                    logger.info('No response information available for stream %s.%s.%s.%s' % nslc)
                    
                    if self._skip_responseless:
                        continue
                
                # we need timewindows per station, not per component because of
                # problem with SeedLink
                if nslc not in timewindows:
                    timewindows[nsl] = (tmin, tmax)
                else:
                    ptmin, ptmax = timewindows[nsl]
                    timewindows[nsl] = (min(tmin,ptmin), max(tmax,ptmax))
                 
                streams[nslc] = stream
                used_stations[nsl] = (station, location)
                pyr_stations[nsl] = sta

            # reduce number of stations
            pyr_stations_wanted = self._station_weeder(pyr_stations.values())
            
            stations_wanted_nsl = [ get_nsl(s) for s in pyr_stations_wanted ]
            
            stream_to_timewindow = {}
            for nslc in streams.keys():
                if nslc[:3] in stations_wanted_nsl:
                    stream_to_timewindow[nslc] = timewindows[nslc[:3]]
                else:
                    streams.pop(nslc)
                    responses.pop(nslc, None)
            
            
            for nslc in self._redundant_channel_weeder(streams.keys()):
                stream_to_timewindow.pop(nslc)
                streams.pop(nslc)
                responses.pop(nslc, None)
                    
            for nsl in used_stations.keys():
                if nsl not in stations_wanted_nsl:
                    used_stations.pop(nsl)
                                
            self.dumpEventInfo(event)
            self.dumpStreams(stream_to_timewindow, event.name, sctime(event.time))
            self.dumpComponents(streams, event.name)
            self.dumpPoleZeros(responses, event.name)
            self.dumpStations(used_stations, event.name)
        
        def dumpEventInfo(self, event):
            ofpath_tmpl = pjoin(self._event_path_template,  self._event_fn_template)
            ofpath = ofpath_tmpl % dict(event_name=event.name)
            ensuredirs(ofpath)
            f = open(ofpath, 'w')
            f.write('name = %s\n' % event.name)
            f.write('time = %s\n' % gmctime(event.time) )
            f.write('latitude = %15.8e\n' % event.lat)
            f.write('longitude = %15.8e\n' % event.lon)
            f.write('depth = %15.8e\n' % (event.depth))
            f.write('magnitude = %15.8e\n' % event.magnitude)
            if event.region:
                f.write('region = %s\n' % event.region)
            f.close()
    
        def dumpComponents(self, streams, eventid):
            for (net, sta, loc, cha), stream in streams.iteritems():
                ofpath_tmpl = pjoin(self._event_path_template,  self._component_fn_template)
                ofpath = fill_stream_tmpl( ofpath_tmpl, net, sta, loc, cha, event_name=eventid)
                ensuredirs(ofpath)
                f = open(ofpath, 'w')
                f.write('%e %e %e\n' % (stream.gain(), stream.azimuth(), stream.dip()))
                f.close()
                
        def dumpPoleZeros(self, responses, eventid):
            for (net, sta, loc, cha), (datalogger, resp)  in responses.iteritems():
                ofpath_tmpl = pjoin(self._event_path_template,  self._polezero_fn_template)
                ofpath = fill_stream_tmpl( ofpath_tmpl, net, sta, loc, cha, event_name=eventid)
                ensuredirs(ofpath)
                f = open(ofpath, 'w')
                f.write('ZEROS %i\n' % resp.numberOfZeros())
                for v in resp.zeros().content():
                    f.write('%e %e\n' % (v.real, v.imag))
                f.write('POLES %i\n' % resp.numberOfPoles())
                for v in resp.poles().content():
                    f.write('%e %e\n' % (v.real, v.imag))
                f.write('CONSTANT %e\n' % (resp.gain()*resp.normalizationFactor()*datalogger.gain()))
                f.close()
                
        def dumpStations(self, stations, eventid):
            ofpath_tmpl = pjoin(self._event_path_template,  self._station_fn_template)
            ofpath = ofpath_tmpl % dict(event_name=eventid)
            ensuredirs(ofpath)
            f = open(ofpath, 'w')
            f.write('# net_code.sta_code.loc_code latitude longitude elevation depth\n')
            for (net,sta,loc), (station,location) in stations.iteritems():
                nsl = '.'.join((net,sta,loc))
                f.write('%-10s %15.8e %15.8e %15.8e %15.8e\n' % (nsl, 
                                                                location.latitude(),
                                                                location.longitude(),
                                                                location.elevation(),
                                                                0.0))
            f.close()
                        
        def dumpStreams(self, stream_to_timewindow, eventid, event_time):
            streamurl = self.recordStreamURL()
            
            all = list(stream_to_timewindow.iteritems())
            all.sort()
            nblock = 100
            i = 0

            while i < len(all):
                rs = seiscomp3.IO.RecordStream.Open(streamurl)
                if not rs:
                    logger.error("could not open recordstream '%s'" % streamurl)
                
                now = sctime(time.time())
                for (net,sta,loc,cha), (tmin,tmax) in all[i:i+nblock]:
                        rs.addStream(net,sta,loc,cha, 
                                    event_time+TimeSpan(tmin),
                                    event_time+TimeSpan(tmax))
                i += nblock
                
                recordInput = seiscomp3.IO.RecordInput(rs, seiscomp3.Core.Array.DOUBLE, seiscomp3.Core.Record.SAVE_RAW)
                if not recordInput:
                    logger.error("could not get recordinput from recordstream '%s'" % streamurl)
                    return
                
                files = {}
                for record in recordInput:
                    net, sta, loc, cha = record.streamID().split('.')
                    if (net,sta,loc,cha) in stream_to_timewindow:
                        ofpath_tmpl = pjoin(self._event_path_template,  self._mseed_fn_template)
                        ofpath = fill_stream_tmpl( ofpath_tmpl, net, sta, loc, cha, event_name=eventid)
                        
                        ensuredirs(ofpath)
                        
                        if not ofpath in files:
                            files[ofpath] = open(ofpath, 'wb')
                            
                        rawdata = record.raw().str()
        
                        files[ofpath].write(rawdata)
                    
        
                for f in files.values():
                    f.close()

except ImportError:
    #sys.stderr.write('no seiscomp3 available\n')
    pass
