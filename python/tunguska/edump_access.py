from pyrocko import pile, trace, eventdata, util, model, pz
import os, calendar, time
from os.path import join as pjoin

from pyrocko.eventdata import FileNotFound

class BadEventFile(Exception):
    pass

st_nslc = '%(network)s_%(station)s_%(location)s_%(channel)s'

class EventDumpAccess(eventdata.EventDataAccess):
    
    def __init__(self, dirpath):
        eventdata.EventDataAccess.__init__(self)
        
        self._dirpath = dirpath
        self._pile = pile.Pile()
        self._pile.add_files(util.select_files([self._dirpath], regex=r'raw-[^/]+\.mseed$'))
        
    def get_restitution(self, tr, allowed_methods):
        
        if 'polezero' in allowed_methods:
            try:
                zeros, poles, constant = self._get_polezero(tr)
                zeros.append(0.0j) # for displacement
                return trace.PoleZeroResponse( poles, zeros, 1./constant )
            except FileNotFound:
                pass
        
        if 'integration' in allowed_methods:
            try:
                cha = self.get_channel(tr)
                if cha is None:
                    raise eventdata.NoRestitution('No gain information available')
                
                return trace.IntegrationResponse(1./cha.gain)
            
            except FileNotFound, e:
                raise eventdata.NoRestitution(e)
        
        raise eventdata.NoRestitution('no working restitution method available')
        
    def _get_channel_description_from_file(self, nslc):
        fnt = pjoin(self._dirpath, 'component-%s.txt' % st_nslc)
        fn = fnt % dict(zip( ('network', 'station', 'location', 'channel'), nslc))
        if os.path.exists(fn):
            f = open(fn,'r')
            gain, azimuth, dip = [ float(x) for x in f.read().split() ]
            f.close()
            
            return model.Channel(nslc[3], azimuth, dip, gain=gain)
        else:
            raise FileNotFound(fn)
    
    def _get_polezero(self, tr):
        fnt = pjoin(self._dirpath, 'polezero-%s.txt' % st_nslc)
        fn = tr.fill_template(fnt)
        if os.path.exists(fn):
            return  pz.read_sac_zpk(fn)
        else:
            raise FileNotFound(fn)
    
    def _get_stations_from_file(self):
        fn = pjoin(self._dirpath, 'stations.txt')
        f = open(fn, 'r')
        stations = []
        for line in f:
            if line.strip().startswith('#'): continue
            toks = line.split()
            if len(toks) != 5: continue
            net,sta,loc = toks[0].split('.')
            lat, lon, elevation, depth  =  [ float(x) for x in toks[1:] ]
            
            station = model.Station(net, sta, loc, lat, lon, elevation)
            stations.append(station)
            
        f.close()
        
        return stations
            
    def _get_events_from_file(self):
        fn = pjoin(self._dirpath, 'event.txt')
        ev = model.Event(load=fn)
        return [ev]
        