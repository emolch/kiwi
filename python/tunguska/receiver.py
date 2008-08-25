

class Receiver:
    
    def __init__(self, lat=0.0, lon=0.0, components=None, name=None, from_string=None ):
        
        # treat as immutable, although it would be possible to change the attributes
        
        if not from_string:
            if components is None: components = 'ned'
            self.lat = lat
            self.lon = lon
            self.components = components
            self.name = name
        else:
            toks = from_string.split()
            if len(toks) >= 3:
                self.lat = float(toks[0])
                self.lon = float(toks[1])
                if components is None:
                    self.components = toks[2]
                else:
                    comps = ''
                    for c in components:
                        if c in toks[2]:
                            comps += c
                    self.components = comps
                    
            if len(toks) == 4:
                self.name = toks[3]
        
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
       
    def set_distazi(self, distance_deg, distance_m, azimuth):
        # should only be called by Seismosizer
        self.distance_deg = distance_deg
        self.distance_m   = distance_m
        self.azimuth      = azimuth
    
    def __str__(self):
        s = ' '.join( (str(self.lat), str(self.lon), self.components) )
        return s
    
    

def load_table( filename, components ):
    receivers = []
    
    file = open(filename, "r")
    irec = 0
    for line in file:        
        if line.lstrip().startswith('#') or line.strip() == '': continue
        r = Receiver( from_string=line, components=components[irec%len(components)] )
        receivers.append(r)
        irec += 1
    
    file.close()
    
    return receivers

    