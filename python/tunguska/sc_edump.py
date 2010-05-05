
import sys, calendar, time, os, logging
from os.path import join as pjoin
from pyrocko import trace, model

try:
    
    import seiscomp3.Client
    import seiscomp3.Communication
    import seiscomp3.DataModel
    from seiscomp3.Core import TimeSpan

    def gmctime(t):
        return time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(t))
    

    def pyrocko_trace(network_code, station, stream, component):
        return trace.Trace(  network=network_code,
                                station=station.code(),
                                location=stream.locCode(),
                                channel=stream.code()+component.code(),
                                deltat=1./stream.sampleRate(),
                                tmin=0.,tmax=0.)
            
    def pyrocko_station(network_code, station, stream):
        return  model.Station(network=network_code,
                                station=station.code(),
                                location=stream.locCode(),
                                lat=station.latitude(),
                                lon=station.longitude(),
                                elevation=station.elevation())
    
    class SimpleInventory:
        def __init__(self, inventory):
            
            for key in ('seismometer', 'respPaz', 'respFir', 'datalogger'):
                d = {}
                
                count = getattr(inventory, key+'Count')
                get = getattr(inventory, key)
                for i in xrange(count()):
                    obj = get(i)
                    if obj.name() in d:
                        logging.warn( 'multiple entries with %s.name == %s' % (key, obj.name()) )
                        
                    d[obj.name()] = obj
                setattr(self, key, d)
                
    
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
                time_range = (-100,100),
                event_path = 'events/%(event_name)s',
                events = []):
            
            seiscomp3.Client.Application.__init__(self, len(argv), argv)
                    
            subscriptions = []
            
            self.setMessagingEnabled(True)
            self.setDatabaseEnabled(True, True)
            self.setLoadConfigModuleEnabled(True)
            self.setLoadInventoryEnabled(True)
            
            st_nslc = '%(network)s_%(station)s_%(location)s_%(channel)s'
            
            self._trace_selector = trace_selector
            self._station_selector = station_selector
            self._station_weeder = station_weeder
            self._timewindow = time_range
            self._event_path_template = event_path
            
            self._mseed_fn_template = 'raw-%s.mseed' % st_nslc
            self._component_fn_template = 'component-%s.txt' % st_nslc
            self._polezero_fn_template = 'polezero-%s.txt' % st_nslc
            self._station_fn_template = 'stations.txt'
            self._fseed_fn_template = 'dataless.seed'
            self._event_fn_template = 'event.txt'
            
            self._simple = None
            self.enableTimer(1)
            self._events = events
            
            
        def handleTimeout(self):
            for event in self._events:
                self.dumpEvent(event)
                
            self.enableTimer(0)
            self.quit()
        
        def get_response(self, network_code, station, stream):
            inventory = seiscomp3.Client.Inventory_Instance().inventory()
            
            if self._simple is None:
                self._simple = SimpleInventory(inventory)
            
            simple = self._simple
            
            try:
                seismometer = simple.seismometer[stream.seismometer()]
                datalogger = simple.datalogger[stream.datalogger()]
                if seismometer.response()[:4] == 'paz:':
                    resppaz = simple.respPaz[seismometer.response()[4:]]
                else:
                    raise ResponseUnavailable('response information for stream %s.%s.%s.%s is in an unexpected format: "%s"' % 
                    (network_code, station.code(), stream.locCode(), stream.code(),
                    seismometer.response()[:4]) )
                    
                return datalogger, resppaz
            
            except KeyError:
                raise ResponseUnavailable( 'response information incomplete for stream %s.%s.%s.%s' % 
                    (network_code, station.code(), stream.locCode(), stream.code()) )
    
        def iterNSSC(self, anno):
            '''Iterate over possibly available streams and components.
            
            Yields tuples like
            
                (network_code, station, stream, component)
                
            '''
            
            inventory = seiscomp3.Client.Inventory_Instance()
            
            cm = self.configModule()
            
            for i in xrange( cm.configStationCount() ):
                cs = cm.configStation(i)
                if cs.enabled():
                    station = inventory.getStation(cs.networkCode(), cs.stationCode(), anno)
                    if station:
                        for istream in xrange(station.seisStreamCount()):
                            
                            stream = station.seisStream(istream)
                            for icomp in xrange(stream.componentCount()):
                                component = stream.component(icomp)
                                yield (cs.networkCode(), station, stream, component)
        

                
        def dumpEvent(self, event):
                                    
            phase_counter = {}
            timewindows = {}
            components = {}
            responses = {}
            used_stations = {}
            
            # preselection 
            pyr_stations = []
            for x in self.iterNSSC(sctime(event.time)):
                network_code, station, stream, component = x
                
                tr = pyrocko_trace(network_code, station, stream, component)
                sta = pyrocko_station(network_code, station, stream)
                sta.set_event_relative_data(event)
                
                if not self._trace_selector(tr):
                    continue
                
                if not self._station_selector(sta):
                    continue
                
                pyr_stations.append(sta)
            
            def get_nsl(x):
                return x.network, x.station, x.location
            
            pyr_stations_wanted = self._station_weeder(pyr_stations)
            stations_wanted_nsl = [ get_nsl(s) for s in pyr_stations_wanted ]
            
            for x in self.iterNSSC(sctime(event.time)):
                network_code, station, stream, component = x
            
                tr = pyrocko_trace(network_code, station, stream, component)
                sta = pyrocko_station(network_code, station, stream)
                sta.set_event_relative_data(event)
                
                if not self._trace_selector(tr):
                    continue
                
                if not self._station_selector(sta):
                    continue
                
                if get_nsl(sta) not in stations_wanted_nsl:
                    continue
                
                tmin = self._timewindow[0]
                tmax = self._timewindow[1]
                
                nslcc = network_code, station.code(), stream.locCode(), stream.code(), component.code()
                nslc = nslcc[:4]
                
                try:
                    responses[nslcc] = self.get_response(network_code, station, stream)
                except ResponseUnavailable, e:
                    pass
                    
                components[nslcc] = component
                
                # we need timewindows per stream, not per component because of
                # problem with SeedLink
                if nslc not in timewindows:
                    timewindows[nslc] = (tmin, tmax)
                else:
                    ptmin, ptmax = timewindows[nslc]
                    timewindows[nslc] = (min(tmin,ptmin), max(tmax,ptmax))
                
                # treat locations as different stations
                used_stations[network_code, station.code(), stream.locCode()] = station
                    
            stream_to_timewindow = {}
            for nslcc in components.keys():
                stream_to_timewindow[nslcc] = timewindows[nslcc[:4]]
            
            
            self.dumpEventInfo(event)
            self.dumpStreams(stream_to_timewindow, event.name, sctime(event.time))
            self.dumpComponents(components, event.name)
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
            for (net, sta, loc, cha, com), component in streams.iteritems():
                ofpath_tmpl = pjoin(self._event_path_template,  self._component_fn_template)
                ofpath = fill_stream_tmpl( ofpath_tmpl, net, sta, loc, cha+com, event_name=eventid)
                ensuredirs(ofpath)
                f = open(ofpath, 'w')
                f.write('%e %e %e\n' % (component.gain(), component.azimuth(), component.dip()))
                f.close()
                
        def dumpPoleZeros(self, responses, eventid):
            for (net, sta, loc, cha, com), (datalogger, resp)  in responses.iteritems():
                ofpath_tmpl = pjoin(self._event_path_template,  self._polezero_fn_template)
                ofpath = fill_stream_tmpl( ofpath_tmpl, net, sta, loc, cha+com, event_name=eventid)
                ensuredirs(ofpath)
                f = open(ofpath, 'w')
                f.write('ZEROS %i\n' % resp.nzeros())
                for v in resp.zeros().content():
                    f.write('%e %e\n' % (v.real, v.imag))
                f.write('POLES %i\n' % resp.npoles())
                for v in resp.poles().content():
                    f.write('%e %e\n' % (v.real, v.imag))
                f.write('CONSTANT %e\n' % (resp.gain()*resp.normFac()*datalogger.gain()))
                f.close()
                
        def dumpStations(self, stations, eventid):
            ofpath_tmpl = pjoin(self._event_path_template,  self._station_fn_template)
            ofpath = ofpath_tmpl % dict(event_name=eventid)
            ensuredirs(ofpath)
            f = open(ofpath, 'w')
            f.write('# net_code.sta_code.loc_code latitude longitude elevation depth\n')
            for (net,sta,loc), station in stations.iteritems():
                nsl = '.'.join((net,sta,loc))
                f.write('%-10s %15.8e %15.8e %15.8e %15.8e\n' % (nsl, 
                                                                station.latitude(),
                                                                station.longitude(),
                                                                station.elevation(),
                                                                station.depth()))
            f.close()
                        
        def dumpStreams(self, stream_to_timewindow, eventid, event_time):
            #streamurl = self.recordStreamURL()
            mh = self.messagingHost()
            streamurl = 'combined://%s:18000;%s:18001' % (mh,mh)
            
            all = list(stream_to_timewindow.iteritems())
            nblock = 100
            i = 0
            while i < len(all):
                rs = seiscomp3.IO.RecordStream.Open(streamurl)
                if not rs:
                    logging.error("could not open recordstream '%s'" % streamurl)
                
                now = sctime(time.time())
                
                added = {}
                for (net,sta,loc,cha,com), (tmin,tmax) in all[i:i+nblock]:
                    if (net,sta,loc,cha) not in added:
                        rs.addStream(net,sta,loc,cha+'?', 
                                    event_time+TimeSpan(tmin),
                                    event_time+TimeSpan(tmax))
                        added[net,sta,loc,cha] = True
                
                i += nblock
                
                recordInput = seiscomp3.IO.RecordInput(rs, seiscomp3.Core.Array.DOUBLE, seiscomp3.Core.Record.SAVE_RAW)
                if not recordInput:
                    logging.error("could not get recordinput from recordstream '%s'" % streamurl)
                    return
                
                files = {}
                for record in recordInput:
                    net, sta, loc, chacom = record.streamID().split('.')
                    cha = chacom[:2]
                    com = chacom[-1]
                    if (net,sta,loc,cha,com) in stream_to_timewindow:
                        ofpath_tmpl = pjoin(self._event_path_template,  self._mseed_fn_template)
                        ofpath = fill_stream_tmpl( ofpath_tmpl, net, sta, loc, chacom, event_name=eventid)
                        
                        ensuredirs(ofpath)
                        
                        if not ofpath in files:
                            files[ofpath] = open(ofpath, 'wb')
                            
                        rawdata = record.raw().str()
        
                        files[ofpath].write(rawdata)
                    
        
                for f in files.values():
                    f.close()

except ImportError:
    sys.stderr.write('no seiscomp3 available\n')
