
from StringIO import StringIO
import numpy as num
import logging, shutil, os, copy, sys, math, datetime, time
import gmtpy

from tempfile import mkdtemp
from subprocess import Popen, PIPE
from os.path import join as pjoin
from pyrocko import trace, io, util
from tunguska.phase import Phase
from tunguska.gfdb import Gfdb

km = 1000.

# how to call the programs
program_bins = {
    'poel': 'poel',
    'gfdb_build': 'gfdb_build',
}

poel_components = 'uz ur ut ezz err ett ezr ert etz tr p vz vr vt'.split()

def str_float_vals(vals):
    return ' '.join( [ '%g' % val for val in vals ] )

def str_int_vals(vals):
    return ' '.join( [ '%i' % val for val in vals ] )

def str_str_vals(vals):
    return ' '.join( [ "'%s'" % val for val in vals ] )

def str_complex_vals(vals):
    return ', '.join( ['(%g, %g)' % (val.real, val.imag) for val in vals ] )

def gfdb_redeploy(in_db_path, out_db_path):
    
    in_db = Gfdb( in_db_path )
    out_db = Gfdb( out_db_path )
    
    cmd = [str(x) for x in ['gfdb_redeploy', in_db_path, out_db_path ]]
    p = Popen( cmd, stdin=PIPE )
    
    for ix in xrange(in_db.nx):
        x = in_db.firstx + ix * in_db.dx
        if not (out_db.firstx <= x and x <= out_db.firstx + (out_db.nx-1)*out_db.dx): continue
        
        for iz in xrange(in_db.nz):
            z = in_db.firstz + iz * in_db.dz
            if not (out_db.firstz <= z and z <= out_db.firstz + (out_db.nz-1)*out_db.dz): continue
    
            logging.info('redeploy distance: %10g km, depth: %10g km\n' % (x/1000., z/1000.))
    
            p.stdin.write( '%g %g\n' % ( x, z ) )
            
        p.stdin.flush()
    
    p.stdin.close()
    p.wait()
    if p.returncode != 0:
        sys.exit('gfdb_phaser: gfdb_redeploy returned an error')


class PoelSourceFunction:

    def __init__(self):
        self.data = [ [ 0., 0. ],
                      [ 0., 1. ] ]

    def __str__(self):
        return '\n'.join( [ '%i %s' % (i, str_float_vals(row)) for (i, row) in enumerate(self.data) ] )

class PoelLayeredModel:
    
    def __init__(self):
        self.data = None
        
    def set_model_from_string(self, s ):
        self.data = num.loadtxt(StringIO(s))
        if self.data.ndim == 1:
            self.data = self.data[num.newaxis,:]

    def set_model(self, depth, mu, nu, nu_u, b, d):
        self.data = num.zeros((len(depth),6), dtype=num.float)
        self.data[:,0] = depth
        self.data[:,1] = mu 
        self.data[:,2] = nu
        self.data[:,3] = nu_u
        self.data[:,4] = b
        self.data[:,5] = d
        
    def get_depth(self):
        return self.data[:,0]
    
    def get_mu(self):
        return self.data[:,1]
    
    def get_nu(self):
        return self.data[:,2]
    
    def get_nu_u(self):
        return self.data[:,3]
    
    def get_b(self):
        return self.data[:,4]
    
    def get_d(self):
        return self.data[:,5]
    
    def get_nlines(self):
        return self.data.shape[0]

    def __str__(self):
        assert self.data is not None 
        srows = []
        for i, row in enumerate(self.data):
            srows.append( '%i %s #' % (i+1, str_float_vals(row)) )
        
        return '\n'.join(srows)


class PoelConfig:
    
    def __init__(self):
        
        # The following are default values which may be overridden at different 
        # steps. The meaning of the markers is as follows:
        #
        #    X:  These are maintained and overridden inside PoelGFDBBuilder;
        #        they should be adjusted when running a plain PoelRunner.
        #
        #    M:  Can probably be modified without breaking the functionality of
        #        this wrapper module.
       
        self.s_start_depth = 50.0  # X
        self.s_end_depth = 50.0    # X
        self.s_radius = 1.0        # M
        self.source_function = PoelSourceFunction() # M
        self.receiver_depth = 0.0  # M
        self.sw_equidistant = 1    # X
        self.no_distances = 10     # X 
        self.distances = [10.0, 100.] # X
        self.t_window = 20.0       # M
        self.no_t_samples = 120    # X
        self.accuracy = 0.025      # M
        self.t_files =  [ x+'.t' for x in poel_components ]  # M
        self.sw_t_files = [ 1 for x in self.t_files ] # M
        self.isurfcon = 1 # M
        self.model = PoelLayeredModel()  # M
        self.model.set_model_from_string('''
   0.00    0.4E+09   0.2   0.4    0.75  5.00
 200.00    0.4E+09   0.2   0.4    0.75  5.00
 '''.lstrip())
    
    def copy(self):
        return copy.deepcopy(self)
       
    def get_output_filenames(self, rundir):
        return [ pjoin(rundir, fn) for fn in self.t_files ]

    def __str__(self):
        
        d = self.__dict__.copy()

        if not self.sw_equidistant:
            d['no_distances'] = len(self.distances)
        d['str_distances'] = str_float_vals(self.distances)
       
        d['sw_t_files_1_3'] = ' '.join([ '%i' % i for i in self.sw_t_files[0:3] ])
        d['t_files_1_3'] = ' '.join([ "'%s'" % s for s in self.t_files[0:3] ])
        d['sw_t_files_4_10'] = ' '.join([ '%i' % i for i in self.sw_t_files[3:10] ])
        d['t_files_4_10'] = ' '.join([ "'%s'" % s for s in self.t_files[3:10] ])
        d['sw_t_files_11_14'] = ' '.join([ '%i' % i for i in self.sw_t_files[10:14] ])
        d['t_files_11_14'] = ' '.join([ "'%s'" % s for s in self.t_files[10:14] ])

        d['no_model_lines'] = self.model.get_nlines()

        template = '''
# This is the input file of FORTRAN77 program "poel06" for modeling
# coupled deformation-diffusion processes based on a multi-layered (half-
# or full-space) poroelastic media induced by an injection (pump) of
# from a borehole or by a (point) reservoir loading.
#
# by R. Wang,
# GeoForschungsZentrum Potsdam
# e-mail: wang@gfz-potsdam.de
# phone 0049 331 2881209
# fax 0049 331 2881204
#
# Last modified: Potsdam, Nov, 2006
#
##############################################################
##                                                          ##
## Cylindrical coordinates (Z positive downwards!) are used ##
## If not others specified, SI Unit System is used overall! ##
##                                                          ##
## Tilt is positive when the upper end of a borehole tilt-  ##
## meter body moves away from the pumping well.             ##
##                                                          ##
##############################################################
#
#	SOURCE PARAMETERS
#	=================
# 1. start and end depth [m] of the injection (pump) screen, and well radius [m];
# 2. number of input data lines for the source time function;
# 3. listing of the source time series.
#-------------------------------------------------------------------------------
  %(s_start_depth)g %(s_end_depth)g  %(s_radius)g                 |dble: s_start_depth, s_end_depth, s_radius;
#-------------------------------------------------------------------------------
 2
#-------------------------------------------------------------------------------
# no    time    source_function
# [-]   [hour]  [m^3/h for injection(+)/pump(-) rate]
#-------------------------------------------------------------------------------
  %(source_function)s
#-------------------------------------------------------------------------------
#
#	RECEIVER PARAMETERS
#	===================
# 1. observation depth [m]
# 2. switch for equidistant trace steping (1/0 = yes/no)
# 3. number of distance samples [m] (<= nrmax defined in "peglobal.h")
# 4. if equidistant, start trace distance, end trace distance; else list of
#    trace distances (all >= 0 and ordered from small to large!)
# 5. length of time window [h], number of time samples
#-------------------------------------------------------------------------------
 %(receiver_depth)g              |dble: r_depth;
 %(sw_equidistant)i              |int: sw_equidistant;
 %(no_distances)i                |int: no_distances;
 %(str_distances)s               |dble: d_1,d_n; or d_1,d_2, ...;
 %(t_window)s %(no_t_samples)i   |dble: t_window; int: no_t_samples;
#-------------------------------------------------------------------------------
#
#	WAVENUMBER INTEGRATION PARAMETERS
#	=================================
# 1. relative accuracy (0.01 for 1%% error) for numerical wavenumber integration;
#-------------------------------------------------------------------------------
 %(accuracy)s                           |dble: accuracy;
#-------------------------------------------------------------------------------
#
#	OUTPUTS A: DISPLACEMENT
#	=======================
# 1. select the 3 displacement time series (1/0 = yes/no)
# 2. file names of these 3 time series
#-------------------------------------------------------------------------------
 %(sw_t_files_1_3)s                                        |int: sw_t_files(1-3);
 %(t_files_1_3)s                                   |char: t_files(1-3);
#-------------------------------------------------------------------------------
#
#	OUTPUTS B: STRAIN TENSOR & TILT
#	===============================
# 1. select strain time series (1/0 = yes/no): Ezz, Err, Ett, Ezr, Ert, Etz
#    (6 tensor components) and Tr (= -dur/dz, the radial component of the
#    vertical tilt). Note Tt can be derived from Etz and Ut
# 2. file names of these 7 time series
#-------------------------------------------------------------------------------
 %(sw_t_files_4_10)s      |int: sw_t_files(4-10);
 %(t_files_4_10)s |char: t_files(4-10);
#-------------------------------------------------------------------------------
#
#	OUTPUTS C: PORE PRESSURE & DARCY VELOCITY
#	=========================================
# 1. select pore pressure and Darcy velocity time series (1/0 = yes/no):
#    P (excess pore pressure), Vz, Vr, Vt (3 Darcy velocity components)
# 2. file names of these 4 time series
#-------------------------------------------------------------------------------
 %(sw_t_files_11_14)s                              |int: sw_t_files(11-14);
 %(t_files_11_14)s                         |char: t_files(11-14);
#-------------------------------------------------------------------------------
#
#	GLOBAL MODEL PARAMETERS
#	=======================
# 1. switch for surface conditions:
#    0 = without free surface (whole space),
#    1 = unconfined free surface (p = 0),
#    2 = confined free surface (dp/dz = 0).
# 2. number of data lines of the layered model (<= lmax as defined in
#    "peglobal.h") (see Note below)
#-------------------------------------------------------------------------------
 %(isurfcon)i                   |int: isurfcon
 %(no_model_lines)i             |int: no_model_lines;
#-------------------------------------------------------------------------------
#
#	MULTILAYERED MODEL PARAMETERS
#	=============================
#
# no depth[m] mu[Pa]    nu    nu_u   B     D[m^2/s]   Explanations
#-------------------------------------------------------------------------------
%(model)s
#--------------------------end of all inputs------------------------------------

Note for the model input format and the step-function approximation for models
with depth gradient:

The surface and the upper boundary of the half-space as well as the interfaces
at which the poroelastic parameters are continuous, are all defined by a single
data line; All other interfaces, at which the poroelastic parameters are
discontinuous, are all defined by two data lines (upper-side and lower-side
values). This input format would also be needed for a graphic plot of the
layered model. Layers which have different parameter values at top and bottom,
will be treated as layers with a constant gradient, and will be discretised to a
number of homogeneous sublayers. Errors due to the discretisation are limited
within about 5%% (changeable, see peglobal.h).
'''.lstrip()

        return template % d
        
class PoelError(Exception):
    pass
        
class PoelRunner:
    
    def __init__(self, tmp=None):
        self.tempdir = mkdtemp(prefix='poelrun', dir=tmp)
        self.program = program_bins['poel']
        self.config = None
    
    
    def run(self, config):
        self.config = config
        
        input_fn = pjoin(self.tempdir, 'input')
                
        f = open(input_fn, 'w')
        poel_input = str(config) 
        
        logging.debug('===== begin poel input =====\n%s===== end poel input =====' % poel_input)
        
        f.write( poel_input )
        f.close()
        program = self.program
        
        old_wd = os.getcwd()
        os.chdir(self.tempdir)
        
        try:
            proc = Popen(program, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        except OSError:
            os.chdir(old_wd)
            raise PoelError('could not start poel: "%s"' % program)
        
        (poel_output, poel_error) = proc.communicate('input\n')
       
        logging.debug('===== begin poel output =====\n%s===== end poel output =====' % poel_output)
        if poel_error:
            logging.error('===== begin poel error =====\n%s===== end poel error =====' % poel_error)

        errmess = []
        if proc.returncode != 0:
            errmess.append('poel had a non-zero exit state: %i' % proc.returncode)
        if poel_error:
            errmess.append('poel emitted something via stderr')
        if poel_output.lower().find('error') != -1:
            errmess.append("the string 'error' appeared in poel output")

        if errmess:
            os.chdir(old_wd)
            raise PoelError('''===== begin poel input =====\n%s===== end poel input =====
===== begin poel output =====\n%s===== end poel output =====
===== begin poel error =====\n%s===== end poel error =====
%s
poel has been invoked as "%s"''' % (poel_input, poel_output, poel_error, '\n'.join(errmess), program))
        
        self.poel_output = poel_output
        self.poel_error = poel_error
        
        os.chdir(old_wd)
        
    def get_traces(self):
       
        if self.config.sw_equidistant == 1:
            nx = self.config.no_distances
            xmin, xmax = self.config.distances
            if nx > 1:
                dx = (xmax-xmin)/(nx-1)
            else:
                dx = 1.0
            distances = [ xmin + ix*dx for ix in xrange(nx) ]
        else:
            distances = self.config.distances
        
        fns = self.config.get_output_filenames(self.tempdir)
        traces = []
        for comp, fn in zip(poel_components, fns):
            data = num.loadtxt(fn, skiprows=1, dtype=num.float)
            nsamples, ntraces = data.shape
            ntraces -= 1
            tmin = data[0,0]
            deltat = (data[-1,0] - data[0,0])/(nsamples-1)
            for itrace in xrange(ntraces):
                x = distances[itrace]
                tr = trace.Trace( '', '%i' % itrace, 'c', comp,
                        tmin=tmin, deltat=deltat, ydata=data[:,itrace+1],
                        meta={'itrace': itrace, 'x':x})
                
                traces.append(tr)
        
        return traces
        
    def __del__(self):
        shutil.rmtree(self.tempdir)
    
class GFTrace(trace.Trace):
    def __init__(self, x, z, ig, trace):
        self.x = x
        self.z = z
        self.ig = ig
        self.trace = trace
    
    def save(self, fn):
        io.save([self.trace], fn)
        
    def cmpkey(self):
        return self.x, self.z, self.ig
    
class GFDBBuilder:
    def __init__(self, partial_db_path, output_db_path, gfdb_config, block_nx, extra_traces_dir=None, tmp=None):
        self.partial_db_path = partial_db_path
        self.output_db_path = output_db_path
        self.gfdb_config = gfdb_config
        self.block_nx = block_nx
        self.gfdb_build_program = program_bins['gfdb_build']
        self.tmp = tmp
        self.extra_traces_dir = extra_traces_dir
        
    def combine(self):
        
        config = self.gfdb_config.copy()
        self._run_gfdb_build(self.output_db_path, config=config)
        
        c = self.gfdb_config
        for iz in xrange(c['nz']):
            z = c['firstz'] + iz*c['dz']
            partial_db_path = self.partial_db_path % { 'depth': z }
            
            gfdb_redeploy(partial_db_path, self.output_db_path)
            
    def _run_gfdb_build(self, database, config=None, traces=None):
        args = []
        if config is not None:
            args = [ str(config[k]) for k in 'nchunks nx nz ng dt dx dz firstx firstz'.split() ]
            
        tempdir = mkdtemp(prefix='gfdbbuilder', dir=self.tmp)

        ignore = open('/dev/null','w')
        proc = Popen([self.gfdb_build_program, database]+args, stdin=PIPE, stdout=ignore)
        if traces is not None:
            for itr, tr in enumerate(traces):
                fn = pjoin(tempdir, '%i.mseed' % itr)
                tr.save(fn)
                proc.stdin.write("%e %e %i '%s'\n" % (tr.x, tr.z, tr.ig, fn))
            
        proc.stdin.close()
        proc.wait()
        ignore.close()
        
        shutil.rmtree(tempdir)
        
        if proc.returncode != 0:
            raise Exception('running gfdb_build failed.')
        
    def create_partial_db(self, z):
        config = self.gfdb_config.copy()
        config['nchunks'] = 1
        config['firstz'] = z
        config['nz'] = 1
        partial_db_path = self.partial_db_path % { 'depth': z }
        self._run_gfdb_build(partial_db_path, config=config)
        
    def fill_partial_db(self, z, traces):
        partial_db_path = self.partial_db_path % { 'depth': z }
        self._run_gfdb_build(partial_db_path, traces=traces)
        
    def extra_save_traces(self, traces):
        if self.extra_traces_dir is None: return
        
        for tr in traces:
            fn = pjoin(self.extra_traces_dir, 'x%08i.z%08i.g%02i.mseed' % (int(tr.x),int(tr.z),tr.ig) )
            tr.save(fn)
    
    def all_depths(self):
        c = self.gfdb_config
        depths = []
        for iz in xrange(c['nz']):
            z = c['firstz'] + iz*c['dz']
            depths.append(z)
        
        return depths
    
    def work_all(self):
        map(self.work_depth, self.all_depths())
    
    def work_depth(self, z):
        self.create_partial_db(z)
        gfdb_config = self.gfdb_config
        
        logging.info( 'Starting work for source depth at %g km' % (z/km) )
        
        for ix in xrange(0, gfdb_config['nx'], self.block_nx):

            this_block_nx = min(gfdb_config['nx']-ix, self.block_nx)
            block_firstx = ix*gfdb_config['dx'] + gfdb_config['firstx']
            block_lastx = (ix+this_block_nx-1)*gfdb_config['dx'] + gfdb_config['firstx']
            tstart = time.time()
            logging.info( 'Processing block (depth %g km, distances %g km - %g km)' % (z/km, block_firstx/km, block_lastx/km) )
        
            traces = self.work_block(block_firstx, block_lastx, this_block_nx, z)
            self.post_processing(traces)
            traces.sort(lambda a,b: cmp(a.cmpkey(), b.cmpkey()))
            self.fill_partial_db( z, traces)
            self.extra_save_traces(traces)
            
            logging.info( 'Done with block (depth %g km, distances %g km - %g km)' % (z/km, block_firstx/km, block_lastx/km) )
            logging.info( 'Block processing time: %s' % str(datetime.timedelta(seconds=time.time()-tstart)))
            
        logging.info( 'Done with work for source depth at %g km' % (z/km) )
            
    def work_block(self, firstx, lastx, nx, z):
        # dummy implementation, should be overloaded in subclass
        if nx > 1:
            dx = (lastx-firstx)/(nx-1)
        else:
            dx = 1.0
        traces = []
        for ix in xrange(nx):
            x = firstx + ix*dx
            for i in xrange(poel_components):
                ig = i+1
                data = num.random.random(100)
                tmin = 0.
                deltat = self.gfdb_config['dt']
                rawtrace = trace.Trace('', '%i' % ix, '', '%i' % ig,
                                           tmin=tmin, deltat=deltat, ydata=data)
                                                              
                tr = GFTrace(  x,z,ig,  rawtrace)
                traces.append(tr)
                
        return traces        
    
    def post_processing(self, traces):
        return traces
    
    def __del__(self):
        import shutil

class PoelGFDBBuilder(GFDBBuilder):
    def __init__(self, partial_db_path, output_db_path, gfdb_config, block_nx, poel_config, extra_traces_dir=None, tag='', tmp=None ):
        GFDBBuilder.__init__(self, partial_db_path, output_db_path, gfdb_config, block_nx, extra_traces_dir=extra_traces_dir, tmp=tmp)
        self.poel_config = poel_config
        self.tag = tag
        
        assert gfdb_config['ng'] == len(poel_components)
        
    def work_block(self, firstx, lastx, nx, z):
        traces = []
        have_gfs = False
        runner = PoelRunner(tmp=self.tmp)
            
        conf = self.poel_config.copy()
        
        conf.s_start_depth = z
        conf.s_end_depth = z
        conf.sw_equidistant = 1
        conf.distances = [ firstx, lastx ]
        conf.no_distances = nx
        conf.no_t_samples = int(round(conf.t_window / self.gfdb_config['dt']))+1
        conf.t_window = (conf.no_t_samples - 1) * self.gfdb_config['dt']
        
        runner.run(conf)
        
        comp2ig = dict( [ (c, ig+1) for (ig, c) in enumerate(poel_components) ] )

        rawtraces = runner.get_traces()
        for tr in rawtraces:
            
            x = tr.meta['x']
            ig = comp2ig[tr.channel]
            
            ix = int(round((x-self.gfdb_config['firstx'])/self.gfdb_config['dx']))
            iz = int(round((z-self.gfdb_config['firstz'])/self.gfdb_config['dz']))
            nxdig = int(math.log(self.gfdb_config['nx'])/math.log(10.))+1
            nzdig = int(math.log(self.gfdb_config['nz'])/math.log(10.))+1
            xtemp = '%%0%ii' % nxdig
            ztemp = '%%0%ii' % nzdig
            tr.set_codes(network=ztemp % iz, station=xtemp % ix, location=self.tag, channel='%02i' % ig)
            if ix < self.gfdb_config['nx']:
                traces.append( GFTrace(x,z,ig, tr) )
            
        return traces
   
