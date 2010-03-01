import util
from bisect import bisect
import sys, os

class OutOfBounds(Exception):
    pass

class PLF:
    '''Nestable piecewise linear function'''
    
    def __init__(self, xdata, ydata):
        self.xdata = xdata
        self.ydata = ydata
        
    def __call__(self, *args):
        x = args[0]
        y0, y1, frac = self.ip(x)
        if isinstance(y0,PLF):
            y0 = y0(*args[1:])
        if isinstance(y1,PLF):
            y1 = y1(*args[1:])
        return y0 + frac*(y1-y0)
        
    def ip(self,x):
        xdata = self.xdata
        ydata = self.ydata
        if x < xdata[0]: raise OutOfBounds()
        if x > xdata[-1]: raise OutOfBounds()
        i = bisect(xdata, x)
        i = max(1, i)
        i = min(len(xdata)-1, i)
        frac = (x-xdata[i-1])/(xdata[i]-xdata[i-1])
        return ydata[i-1], ydata[i], frac
        
class Phase:
    def __init__(self,name,filename=None):
    
        self.name = name
        self.filename = filename
        if filename is None:
            if os.path.isfile(name+'.phase'):
                filename = name+'.phase'
            else:
                filename = os.path.join(util.kiwi_aux_dir(), 'phases', name)
            
        
        f = open(filename,'r')
        self.ref_points = []
        dists = {}
        distances, depths, times = [], [], []
        have_seen = {}
        have_depth = False
        for line in f:
            toks = line.split()
            dist = float(toks[0])
            if len(toks) == 3:
                depth = float(toks[1])
                have_depth = True
            else:
                depth = 10000.
            
            if (dist,depth) not in have_seen:
                times.append(float(toks[-1]))
                distances.append(dist)
                depths.append(depth)
            
            have_seen[(dist,depth)] = True
            
        f.close()
        if have_depth:
            dists = {}
            for di, de, ti in zip(distances, depths, times):
                if not di in dists:
                    dists[di] = ([], [])
                dists[di][0].append(de)
                dists[di][1].append(ti)
            distances1 = []
            depth_plfs = []
            for (di, (des,tis)) in sorted(dists.items()):
                distances1.append(di)
                depth_plfs.append( PLF(des,tis) )
            self.lookup = PLF(distances1, depth_plfs)
        else:
            self.lookup = PLF(distances, times)
            
        self.have_depth = have_depth
                
    def __call__(self, distance, depth=10000.):
        distance = float(distance)
        depth = float(depth)
        try:
            return self.lookup(distance, depth)
        except OutOfBounds:
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
        
    def __call__(self, distance, depth=10000.):
        '''Returns timing at given distance.'''
        for phase in self.phases:
            t = phase(distance, depth)
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
        
    def __call__(self, distance, depth=10000.):
        '''Returns representation which can be used for Seismosizer.do_set_misfit_taper().'''
        return ( self.timings[0](distance,depth), 0.,
                 self.timings[1](distance,depth), self.amplitude,
                 self.timings[2](distance,depth), self.amplitude,
                 self.timings[3](distance,depth), 0. )

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
    