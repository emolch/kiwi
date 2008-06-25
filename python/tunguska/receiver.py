

class Receiver:
    
    def __init__(self, lat=0.0, lon=0.0, components='ned', name=None, from_string=None ):
        
        # treat as immutable, although it would be possible to change the attributes
        
        if not from_string:
            self.lat = lat
            self.lon = lon
            self.components = components
            self.name = name
        else:
            toks = from_string.split()
            if len(toks) >= 3:
                self.lat = float(toks[0])
                self.lon = float(toks[1])
                self.components = toks[2]
            if len(toks) == 4:
                self.name = toks[3]
        
        # used and set by Seismosizer, when attached to it
        # treat as read-only
        self.distance_m = None
        self.distance_deg = None
        self.azimuth = None
        self.proc_id = None
        self.misfits = [0.]*len(self.components)
        self.misfit_norm_factors = [0.]*len(self.components)
       
    def set_distazi(self, distance_deg, distance_m, azimuth):
        # should only be called by Seismosizer
        self.distance_deg = distance_deg
        self.distance_m   = distance_m
        self.azimuth      = azimuth
    
    def __str__(self):
        return ' '.join( (str(self.lat), str(self.lon), self.components) )
    
    

def load_table( filename ):
    receivers = []
    
    file = open(filename, "r")
    
    for line in file:        
        if line.lstrip().startswith('#') or line.strip() == '': continue
        r = Receiver( from_string=line )
        receivers.append(r)
    
    file.close()
    
    return receivers

    