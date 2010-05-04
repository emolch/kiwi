import pymseed

class Receiver:
    
    def __init__(self, lat=0.0, lon=0.0, depth=0.0, components=None, name=None, from_string=None ):
        
        # treat as immutable, although it would be possible to change the attributes
        
        if not from_string:
            if components is None: components = 'ned'
            self.lat = lat
            self.lon = lon
            self.depth = depth
            self.components = components
            self.name = name
        else:
            toks = from_string.split()
            if len(toks) >= 4:
                self.lat = float(toks[0])
                self.lon = float(toks[1])
                self.depth = float(toks[2])
                if components is None:
                    self.components = toks[3]
                else:
                    comps = ''
                    for c in components:
                        if c in toks[3]:
                            comps += c
                    self.components = comps
                    
            if len(toks) == 5:
                self.name = toks[4]
        
        self.cumulative_shift = 0.
        
        # used and set by Seismosizer, when attached to it
        # treat as read-only
        self.enabled = True
        self.distance_m = None
        self.distance_deg = None
        self.azimuth = None
        self.proc_id = None
        self.misfits = [0.]*len(self.components)
        self.misfit_norm_factors = [0.]*len(self.components)
        # filled by seismosizer.get_receivers_copy()
        self.ref_seismograms        = [None]*len(self.components)
        self.syn_seismograms        = [None]*len(self.components)
        self.ref_spectra            = [None]*len(self.components)
        self.syn_spectra            = [None]*len(self.components)
        self.comp_ind = dict( [ (c,i) for (i,c) in enumerate(self.components) ] )
       
    def set_distazi(self, distance_deg, distance_m, azimuth):
        # should only be called by Seismosizer
        self.distance_deg = distance_deg
        self.distance_m   = distance_m
        self.azimuth      = azimuth
    
    def get_location(self):
        toks = self.name.split('.')
        if len(toks) == 3:
            return toks[2]
        else:
            return ''
    
    def get_station(self):
        toks = self.name.split('.')
        if len(toks) == 3:
            return toks[1]
        
        # compatibility with old preprocessing scripts
        else:
            return toks[0]
       
    
    def get_network(self):
        toks = self.name.split('.')
        if len(toks) == 3:
            return toks[0]
            
        # compatibility with old preprocessing scripts
        elif len(toks) == 2:
            return toks[1]
        else:
            return ''
        
    def __str__(self):
        s = ' '.join( (str(self.lat), str(self.lon), str(self.depth), self.components) )
        return s
    
    def save_traces_mseed(self, filename_tmpl='%{whichset}s_%(network)s_%(station)s_%(location)s_%(channel)s.mseed' ):
        
        station, network = self.get_station(), self.get_network()
        for icomp, comp in enumerate(self.components):
            channel = comp
            for (whichset, sgram) in zip(('references', 'synthetics'), 
                             (self.ref_seismograms[icomp], self.syn_seismograms[icomp])):
                if sgram and len(sgram[0]) > 1:
                    starttime = sgram[0][0]
                    endtime = sgram[0][-1]
                    deltat = (endtime-starttime)/(len(sgram[0])-1)
                    data = sgram[1]
                    location = whichset
                    trace = (network, station, location, channel, starttime*pymseed.HPTMODULUS, endtime*pymseed.HPTMODULUS, 1.0/deltat, data)
                    fn = filename_tmpl % { 'whichset': whichset,
                                           'network': network,
                                           'station': station,
                                           'location': location,
                                           'channel': channel }
                                           
                    pymseed.store_traces([trace], fn)

    def get_misfit(self, component):
        return self.misfits[self.comp_ind[component]]
    
    def get_misfit_norm_factor(self, component):
        return self.misfit_norm_factors[self.comp_ind[component]]

    def get_misfit_and_norm_factor(self, component):
        i = self.comp_ind[component]
        return self.misfits[i], self.misfit_norm_factors[i]

def load_table( filename, components=None ):
    receivers = []
    
    file = open(filename, "r")
    irec = 0
    for line in file:        
        if line.lstrip().startswith('#') or line.strip() == '': continue
        if components is not None:
            r = Receiver( from_string=line, components=components[irec%len(components)] )
        else:
            r = Receiver( from_string=line )
        receivers.append(r)
        irec += 1
    
    file.close()
    
    return receivers

    