from pyrocko import pile, trace, eventdata, util, model, pz
import os, calendar, time
from os.path import join as pjoin

class BadEventFile(Exception):
    pass

class FileNotFound(Exception):
    
    def __init__(self, s):
        self.s = s
        
    def __str__(self):
        return 'File not found: %s' % self.s

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
                gain = self._get_gain(tr)
                return trace.IntegrationResponse(1./gain)
            
            except FileNotFound, e:
                raise eventdata.NoRestitution(e)
        
        raise eventdata.NoRestitution('no working restitution method available')
        
    def _get_gain(self, tr):
        fnt = pjoin(self._dirpath, 'gain-%s.txt' % st_nslc)
        fn = tr.fill_template(fnt)
        if os.path.exists(fn):
            f = open(fn,'r')
            gain = float(f.read())
            f.close()
            return gain
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
        f = open(fn, 'r')
        inp = {}
        for line in f:
            k,v = line.split( '=')
            k = k.strip()
            v = v.strip()
            inp[k] = v
            
        required = 'time', 'latitude', 'longitude'
        for k in required:
            if k not in inp:
                raise BadEventFile('key "%s" missing in file "%s"' % (k, fn))     
        
        ev = model.Event()
        ev.time = calendar.timegm( time.strptime( inp['time'][:19], '%Y-%m-%d %H:%M:%S' ))
        ev.lat = float(inp['latitude'])
        ev.lon = float(inp['longitude'])
        if 'depth' in inp:
            ev.depth = float(inp['depth'])*1000.
        if 'publicID' in inp:
            ev.name = inp['publicID']
        if 'magnitude' in inp:
            ev.magnitude = float(inp['magnitude'])
            
        return [ev]
        