import util

import sys, os

class Phase:
    def __init__(self,name,filename=None):
    
        self.name = name
        self.filename = filename
        if filename is None:
            filename = os.path.join(util.invearthquake_aux_dir(), 'phases', name)
        
        f = open(filename,'r')
        self.ref_points = []
        for line in f:
            distance, time = [float(x) for x in line.split()]
            self.ref_points.append( (distance, time) )
        f.close()
        
    def __call__(self, distance):
        distance = float(distance)
        for (low,high) in zip( self.ref_points[0:-1],self.ref_points[1:len(self.ref_points)]):
            if low[0] <= distance <= high[0]:
                return low[1] + (distance-low[0])/(high[0]-low[0])*(high[1]-low[1])
        return None
    
    def __repr__(self):
        s = "Phase(name='%s'" % self.name
        if self.filename is not None:
            s += ", filename='%s'" % self.filename
        s += ')'
        return s
        
            
class Timing:
    '''"Intelligent" Phase, e.g. "S or Sn, whatever is available minus 10 s".'''
    
    def __init__(self, phases, offset=0.):
        '''usage:
           Timing('S')
           Timing(('S','Sn'),+5)
           Timing(Phase('R'),-10)'''
        
        if isinstance(phases, str):
            phases = [ phases ]
        self.phases = []
        for p in phases:
            if isinstance(p, str):
                self.phases.append( Phase(p) )
            else:
                self.phases.append( p )
                
        self.offset = float(offset)
        
    def __call__(self, distance):
        '''Returns timing at given distance.'''
        for phase in self.phases:
            t = phase(distance)
            if not t is None:
                return t+self.offset

    def __repr__(self):
        s = 'Timing(phases=['
        s += ', '.join([repr(p) for p in self.phases])
        s += '], offset='+str(self.offset)+')'
        return s
        
class Taper:
    def __init__(self, timings=None,
                       phases=None, offsets=None, amplitude=1.):
        '''usage:
           t = Taper( timings=(Timing('P',-10),Timing('P', 0), ...)
           t = Taper( phases=('S','Sn'), offsets=(-10,0,40,50) )'''
        
        if phases and offsets:
            timings = [ Timing(phases,offset) for offset in offsets ]
            
        assert(len(timings) == 4)
        self.timings = timings
        self.amplitude = amplitude
        
    def __call__(self, distance):
        '''Returns representation which can be used for Seismosizer.do_set_misfit_taper().'''
        return ( self.timings[0](distance), 0.,
                 self.timings[1](distance), self.amplitude,
                 self.timings[2](distance), self.amplitude,
                 self.timings[3](distance), 0. )

    def __repr__(self):
        s = 'Taper(timings=[\n    '
        s += ',\n    '.join([repr(t) for t in self.timings])
        s += '\n])'
        return s
        
if __name__ == '__main__':
    p = Phase('P')
    pn = Phase('Pn')
    for i in range(15):
        print p(i*1000000.), pn(i*1000000.)

    sany = ('S','Sn')
    timings = [ Timing(sany,offset) for offset in [-10., 0., 40., 50. ] ] 
    print Taper(timings)(1000000.)
    t = Taper(phases=('S','Sn'), offsets=(-10,0, 40,50))
    print t(1000000.)
    print repr( t )
    