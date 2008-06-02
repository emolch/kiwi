
import config
from subprocess import Popen, PIPE
from util import gform
import os
import shutil
import numpy as num
from scipy.io import read_array
import tempfile

pjoin = os.path.join

class Gfdb:
    def __init__(self, gfdbpath):
   
        self.path = gfdbpath
        gfdb_infos_str = {}
        cmd = [ config.gfdb_info_prog, gfdbpath ]
        
        self.string = Popen( cmd, stdout=PIPE ).communicate()[0]
        
        for line in self.string.splitlines():
            k,v = line.split('=')
            gfdb_infos_str[k] = v.strip()
        
        for k in [ 'dt', 'dx', 'dz', 'firstx', 'firstz' ]:
            setattr(self, k, float( gfdb_infos_str[k] ))
            
        for k in [ 'nchunks', 'nx', 'nz', 'ng' ]:
            setattr(self, k, int( gfdb_infos_str[k] ))
            
        self.extractor = None
        self.tempdir = None
        self.tempfilebase = None
    
    
    def get_traces_slow( self, x, z ):
        
        if not self.extractor:
            self.extractor = Popen( [config.gfdb_extract_prog, self.path], stdin=PIPE, stdout=PIPE )
            self.tempdir = tempfile.mkdtemp("","gfdb_extract-")
            self.tempfilebase = pjoin( self.tempdir, 'trace' )
        
        
        fns = []
        for ig in range(self.ng):
            fn = '%s-%i.table' % (self.tempfilebase, ig)
            self.extractor.stdin.write("%f %f %i '%s'\n" % (x,z,ig+1,fn))
            self.extractor.stdin.flush()
            answer = self.extractor.stdout.readline()
            if answer.strip() == 'ok':
                fns.append(fn)
            else:
                fns.append(None)
        
        
        traces = []
        for fn in fns:
            if fn:
                f = open(fn,'r')
                tabX = []
                for line in f:
                    x,y = [float(x) for x in line.split()]
                    tabX.append([x,y])
                
                tab = num.array( tabX, dtype=num.float ).transpose()
                f.close()
                
                if tab.ndim == 2:
                    time = tab[0].copy()
                    data = tab[1].copy()
                else:
                    time = num.array([tab[0].copy()])
                    data = num.array([tab[1].copy()])
                
                traces.append( (time,data) )
            else:
                traces.append( None )
        
        
        return traces
            
    def __del__(self):
        if self.extractor:
            self.extractor.stdin.close()
            self.extractor.stdout.close()
            self.extractor.wait()
            self.extractor = None
        if self.tempdir:
            shutil.rmtree(self.tempdir)
            self.tempdir = None
            self.tempfilebase = None
           
    def __str__(self):
        return '''
GFDB: %s
    dt     [s]: %s
    dx     [m]: %s
    dz     [m]: %s
    firstx [m]: %s
    firstz [m]: %s
    nx        : %6i
    nz        : %6i
    ng        : %6i
'''.strip() % tuple([self.path] + [ gform(x,5) for x in (self.dt, self.dx, self.dz,
                                    self.firstx, self.firstz) ] + [ self.nx, self.nz, self.ng ])
    

if __name__ == '__main__':
    import sys
    g = Gfdb(sys.argv[1])
    print g
    print g.get_traces_slow( 10000.,10000. )
    
    
    