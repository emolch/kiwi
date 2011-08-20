
import config
from subprocess import Popen, PIPE
from util import gform
import math
import os
import shutil
import numpy as num
import scipy
import tempfile
import sys
from pyrocko import io, util
pjoin = os.path.join

class Gfdb:
    def __init__(self, gfdbpath):
   
        self.path = gfdbpath
        self.extractor = None
        self.builder = None
        self.tempdir = None
        self.tempfilebase = None
        
        gfdb_infos_str = {}
        cmd = [ config.gfdb_info_prog, gfdbpath ]
        
        gfdb_info_process = Popen( cmd, stdout=PIPE)
        self.string = gfdb_info_process.communicate()[0]
        if gfdb_info_process.poll():
            sys.exit('fatal: could not get get gfdb information')
        
        for line in self.string.splitlines():
            k,v = line.split('=')
            gfdb_infos_str[k] = v.strip()
        
        for k in [ 'dt', 'dx', 'dz', 'firstx', 'firstz' ]:
            setattr(self, k, float( gfdb_infos_str[k] ))
            
        for k in [ 'nchunks', 'nx', 'nz', 'ng' ]:
            setattr(self, k, int( gfdb_infos_str[k] ))    
    
    def dump_traces(self, x, z, format='table'):
        if not self.extractor:
            self.extractor = Popen( [config.gfdb_extract_prog, self.path], stdin=PIPE, stdout=PIPE, close_fds=True)
            self.tempdir = tempfile.mkdtemp("","gfdb_extract-")
            self.tempfilebase = pjoin( self.tempdir, 'trace' )
        
        
        fns = []
        for ig in range(self.ng):
            fn = '%s-%i.%s' % (self.tempfilebase, ig, format)
            self.extractor.stdin.write("%f %f %i '%s'\n" % (x,z,ig+1,fn))
            self.extractor.stdin.flush()
            answer = self.extractor.stdout.readline()
            if answer.strip() == 'ok':
                fns.append(fn)
            else:
                fns.append(None)
        
        return fns
    
    def get_traces_slow( self, x, z ):
        fns = self.dump_traces(x,z)
        traces = []
        for fn in fns:
            if fn:
                tab = num.loadtxt(fn, dtype=num.float).transpose()
                
                if tab.ndim == 2:
                    time = tab[0]
                    data = tab[1]
                else:
                    time = num.array([tab[0].copy()])
                    data = num.array([tab[1].copy()])
                
                traces.append( (time,data) )
            else:
                traces.append( None )
        
        
        return traces
    
    def get_traces_pyrocko(self, x,z):
        fns = self.dump_traces(x,z, format='mseed')
        traces = []
        for igm, fn in enumerate(fns):
            ig = igm+1
            if fn:
                for tr in io.load(fn):
                    sx = util.base36encode(int(round(x)))
                    sz = util.base36encode(int(round(z)))
                    tr.meta = {'x':x, 'z':z, 'ig':ig}
                    tr.set_codes(network=sz, station=sx, channel='%i' % ig)
                    traces.append(tr)
                    
        return traces
                
    def put_traces_slow( self, x,z, traces):
        assert len(traces) == self.ng
        if not self.builder:
            self.builder = Popen( [config.gfdb_build_prog, self.path], stdin=PIPE, stdout=PIPE, close_fds=True)
            self.tempdir = tempfile.mkdtemp('','gfdb_build-')
            self.tempfilebase = pjoin(self.tempdir, 'trace')
            
        fns = []
        
        for ig, xx in enumerate(traces):
            if xx is not None:
                (time, data)  = xx
                tab = num.array((time,data)).transpose()
                if len(tab) > 0:
                    fn = '%s-%i.table' % (self.tempfilebase, ig)
                    num.savetxt(fn, tab)
                    self.builder.stdin.write("%f %f %i '%s'\n" % (x,z,ig+1,fn))
                    fns.append(fn)
        
        self.builder.stdin.flush()

        for fn in fns:
            answer = self.builder.stdout.readline()
            assert answer.strip() == fn
                
            
    def terminate(self):
        if self.extractor:
            self.extractor.stdin.close()
            self.extractor.stdout.close()
            self.extractor.wait()
            self.extractor = None
        
        if self.builder:
            self.builder.stdin.close()
            self.builder.stdout.close()
            self.builder.wait()
            self.builder = None
            
        if self.tempdir:
            shutil.rmtree(self.tempdir)
            self.tempdir = None
            self.tempfilebase = None
    
    def __del__(self):
        self.terminate()
           
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
    

    
