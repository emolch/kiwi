
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
from pyrocko.moment_tensor import MomentTensor, symmat6


km = 1000.

# number of processes to run in parallel
nworkers = 1

# how to call the programs
program_bins = {
    'qseis': 'qseis',
    'gfdb_build': 'gfdb_build',
}

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

class QSeisLayeredModel:
    
    def __init__(self):
        self.data = None
        
    def set_model_from_string(self, s, units=('standard','ugly')[0]):
        self.data = num.loadtxt(StringIO(s))
        if self.data.ndim == 1:
            self.data = self.data[num.newaxis,:]
            
        if units == 'ugly':
            self.data[:,0] *= 1000.
            self.data[:,1] *= 1000.
            self.data[:,2] *= 1000.
            self.data[:,3] *= 1000.
            
    def set_model(self, depth, vp,vs, density, qp, qs):
        
        self.data = num.zeros((len(depth),6), dtype=num.float)
        self.data[:,0] = depth
        self.data[:,1] = vp
        self.data[:,2] = vs
        self.data[:,3] = density
        self.data[:,4] = qp
        self.data[:,5] = qs
        
    def get_depth(self):
        return self.data[:,0]
    
    def get_vp(self):
        return self.data[:,1]
    
    def get_vs(self):
        return self.data[:,2]
    
    def get_density(self):
        return self.data[:,3]
    
    def get_qp(self):
        return self.data[:,4]
    
    def get_qs(self):
        return self.data[:,5]
    
    def __str__(self):
        if self.data is None:
            return '0'
        
        srows = []
        for i, row in enumerate(self.data):
            row_ugly = row[0]/1000., row[1]/1000., row[2]/1000., row[3]/1000., row[4], row[5]
            srows.append( '%i %s' % (i+1, str_float_vals(row_ugly)) )
        
        return ('%i\n' % self.data.shape[0]) + '\n'.join(srows)


class QSeisConfig:
    
    def __init__(self):
        
        # The following are default values which may be overridden at different 
        # steps. The meaning of the markers is as follows:
        #
        #    X:  These are maintained and overridden inside QSeisGFDBBuilder;
        #        they should be adjusted when running a plain QSeisRunner.
        #
        #    O:  These are configured by QSeisConfig.autoconf_modelling(),
        #        based on gfdb_config values and some tuning variables.
        # 
        #    M:  Can probably be modified without breaking the functionality of
        #        this wrapper module.
         
        self.source_depth_km = 10.0       # X
        
        self.receiver_depth_km = 0.0      # M
        self.sw_equidistant = 1           # X      # 1 is eqdist, 0 if not
        self.sw_d_unit = 1                # X      # 1 is km,  0 is deg
        self.no_distances = 100           # X
        self.distances_km = [100., 600.]  # X      # first, last, if equidist
        self.t_start = -20.               # O
        self.t_window = 1024./2           # O
        self.no_t_samples = 1024          # O
        self.sw_t_reduce = 1              # O      # 1 is vel [km/s], 0 is 
                                                   # slowness [s/deg]
        self.t_reduce = 12.               # O
        
        self.sw_algorithm = 0             # O      # 0 full wave-field, 1 or 2 
                                                   # slowness window
        self.slw = 0.01, 0.02, 0.5, 0.6   # O      # cos taper
        self.sample_rate = 2.5            # M      # wavenumber integration 
                                                   # sampling rate in units of
                                                   # spacial nyquist freq
        self.supp_factor = 0.01           # M      # suppress temporal aliasing
                                                   # factor
        
        self.isurf = 0                    # M
        self.sw_path_filter = 0           # M
        self.shallow_depth_limit = 560.0  # M
        self.no_of_depth_ranges = 0       # M
        
        self.wavelet_duration = 4.0       # M      # in samples
        self.sw_wavelet = 2               # M      # 1=pulse, 2=heaviside
        
        self.norm_factor = 1.0            # M
        self.filter_no_roots = 0          # M
        self.roots = []                   # M
        self.filter_no_poles = 0          # M
        self.poles = []                   # M
        
        
        self.gf_sw_source_types = \
            1,1,1,1,0,0                   # X       # explosion, strike-slip, 
                                                    # dip-slip, clvd, 
                                                    # single_f_v, single_f_h
        self.gf_filenames = \
            'ex', 'ss', 'ds', \
            'cl', 'fz', 'fh'              # M
        
        
        self.source_type = 1              # X      # source:  1=mt   File
        self.source_vals = [ 1.,1.,1.,0.,0.,0. ]   # X      # Mxx Myy Mzz Mxy Myz Mzx  
        self.seismogram_filename = 'seis' # M      # length limits may apply due
                                                   # to fixed length strings in
                                                   # Qseis.
        
        self.sw_irregular_station_azimuths = 0       # X
        self.station_azimuths = [ 0.0 ]              # X
        
        self.sw_flat_earth_transform = 0             # M
        self.gradient_resolutions = 0.25, 0.25, 5.0  # M   # (vp, vs, density) 
                                                           # given in [%]
        
        self.layered_model = QSeisLayeredModel()     # M
        self.receiver_model = QSeisLayeredModel()    # M
    
    
    def autoconf_modelling(self, gfdb_config, 
        length_factor=1.0,
        tlead_in=0., tlead_out=0., 
        slowness_window_factors=(0.005,0.01, 2.,4.),
        allow_time_reduction=True):
        
        '''Setup QSeis configuration for creating a Kiwi GF database.
        
        This method considers the extent and sampling intervals of a GFDB
        configuration and the minimum and maximum seismic velocities in the
        earth model and tries to set reasonable values for some of the basic
        QSeis configuration parameters.

        The setup is guided by a few tuning parameters:
        
            length_factor -- increase automatically estimated length by factor.
            tlead_in -- seconds of padding to be added at the beginning.
            tlead_out -- seconds of padding to be added at the end.
            slowness_window_factors -- factors used to create the slowness 
                window the first two are relative to 1/vmax the latter are 
                relative to 1/vmin.
            allow_time_reduction -- whether to allow time reduction to shorten
                the required time windows in the modelling.
                
        The following QSeis configuration parameters are set:
        
            t_start 
            t_window 
            no_t_samples 
            sw_t_reduce
            t_reduce
            sw_algorithm
            slw
        '''
        
        
        xmax = gfdb_config['firstx'] + (gfdb_config['nx']-1)*gfdb_config['dx']
        xmin = gfdb_config['firstx']
        
        vmin = self.layered_model.get_vs().min()
        vmax = self.layered_model.get_vp().max()
        
        if allow_time_reduction:
            vred = vmax
        else:
            vred = None
        
        tmin = xmin/vmax - tlead_in
        tmax = xmax/vmin*length_factor + tlead_out
        
        if vred is not None:
            tmin_red = xmin/vmax - xmin/vred - tlead_in
            tmax_red = xmax/vmin*length_factor - xmax/vred + tlead_out
        else:
            tmin_red = tmin
            tmax_red = tmax
            
        nsamples_phys = (tmax_red-tmin_red)/gfdb_config['dt']
        nsamples = 2**(int(num.log(nsamples_phys)/num.log(2))+1)
        
        slowness_window = (1./vmax*slowness_window_factors[0], 1./vmax*slowness_window_factors[1], 
                           1./vmin*slowness_window_factors[2], 1./vmin*slowness_window_factors[3])
        
        logging.info( 'Min distance: %g km' % (xmin/km) )
        logging.info( 'Max distance: %g km' % (xmax/km) )
        logging.info( 'Min velocity: %g km/s' % (vmin/km) )
        logging.info( 'Max velocity: %g km/s' % (vmax/km) )
        logging.info( 'Time range (without time reduction): [%g, %g] s' % (tmin, tmax) )
        logging.info( 'Time range (with time reduction): [%g, %g] s' % (tmin_red, tmax_red) )
        logging.info( 'Number of samples needed (physical): %i' % nsamples_phys )
        logging.info( 'Number of samples needed (computational): %i' % nsamples )
        logging.info( 'Slowness window: [%g, %g, %g, %g] s/km' % tuple( [ _s*km for _s in slowness_window ]) )
        
        self.t_start = tmin_red
        self.t_window = (nsamples-1)*gfdb_config['dt']
        self.no_t_samples = nsamples
        self.sw_t_reduce = 1                   # 1 is vel [km/s], 0 is slowness [s/deg]
        if vred is not None:
            self.t_reduce = vred/km
        else:
            self.t_reduce = 0
        
        self.sw_algorithm = 0                  # 0 full wave-field, 1 or 2 slowness window
        self.slw = tuple( [ _s*km for _s in slowness_window ]) # cos taper
    
    def copy(self):
        return copy.deepcopy(self)
    
    def get_seismogram_filenames_zrt(self, rundir):
        fn = self.seismogram_filename
        return (pjoin(rundir, fn + '.tz'),
                pjoin(rundir, fn + '.tr'),
                pjoin(rundir, fn + '.tt'))
        
    def __str__(self):
        
        d = self.__dict__.copy()
        if not self.sw_equidistant:
            d['no_distances'] = len(self.distances_km)
        d['str_distances'] = str_float_vals(self.distances_km)
        d['str_slw'] = str_float_vals(self.slw)
        if self.roots:
            d['str_roots'] = '\n'+str_complex_vals(self.roots)
        else:
            d['str_roots'] = '\n#'
        if self.poles:
            d['str_poles'] = '\n'+str_complex_vals(self.poles)
        else:
            d['str_poles'] = '\n#'
        d['str_gf_sw_source_types'] = str_int_vals(self.gf_sw_source_types)
        d['str_gf_filenames'] = str_str_vals(self.gf_filenames)
        d['str_source_vals'] = str_float_vals(self.source_vals)
        d['str_station_azimuths'] = str_float_vals(self.station_azimuths)
        d['str_gradient_resolutions'] = str_float_vals(self.gradient_resolutions)
        
        template = '''
# source_depth_km
%(source_depth_km)g
# 
# receiver_depth_km
%(receiver_depth_km)g
# sw_equidistant sw_d_unit
%(sw_equidistant)i %(sw_d_unit)i
# no_distances
%(no_distances)i
%(str_distances)s
# t_start t_window no_t_samples
%(t_start)g %(t_window)g %(no_t_samples)i
# sw_t_reduce t_reduce
%(sw_t_reduce)i %(t_reduce)g
#
# sw_algorithm
%(sw_algorithm)i
# slowness_window
%(str_slw)s
# sl_sample_rate
%(sample_rate)g
# supp_factor
%(supp_factor)g
#
# isurf
%(isurf)i
# sw_path_filter shallow_depth_limit
%(sw_path_filter)i %(shallow_depth_limit)g
# no_of_depth_ranges
%(no_of_depth_ranges)i
#
# wavelet_duration sw_wavelet
%(wavelet_duration)g %(sw_wavelet)i
#
# norm_factor
%(norm_factor)g
# roots
%(filter_no_roots)i%(str_roots)s
# poles
%(filter_no_poles)i%(str_poles)s
#
# gf_sw_source_types
%(str_gf_sw_source_types)s
%(str_gf_filenames)s
#
# source_type source_vals seismogram_filename
%(source_type)i %(str_source_vals)s '%(seismogram_filename)s'
# sw_irregular_station_azimuths
%(sw_irregular_station_azimuths)i
%(str_station_azimuths)s
#
# sw_flat_earth_transform
%(sw_flat_earth_transform)i
# gradient_resolutions
%(str_gradient_resolutions)s
#
%(layered_model)s
%(receiver_model)s
'''.lstrip()

        return template % d
        
class QSeisError(Exception):
    pass
        
class QSeisRunner:
    
    def __init__(self, tmp=None):
        self.tempdir = mkdtemp(prefix='qseisrun', dir=tmp)
        self.program = program_bins['qseis']
        self.config = None
    
    
    def run(self, config):
        self.config = config
        
        input_fn = pjoin(self.tempdir, 'input')
                
        f = open(input_fn, 'w')
        
        qseis_input = str(config) % { 'tempdir': self.tempdir }
        
        logging.debug('===== begin qseis input =====\n%s===== end qseis input =====' % qseis_input)
        
        f.write( qseis_input )
        f.close()
        program = self.program
        
        old_wd = os.getcwd()
        os.chdir(self.tempdir)
        
        try:
            proc = Popen(program, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        except OSError:
            os.chdir(old_wd)
            raise QSeisError('could not start qseis: "%s"' % program)
        
        (qseis_output, qseis_error) = proc.communicate('input\n')
       
        logging.debug('===== begin qseis output =====\n%s===== end qseis output =====' % qseis_output)
        logging.debug('===== begin qseis error =====\n%s===== end qseis error =====' % qseis_error)

        if proc.returncode != 0 or qseis_error or qseis_output.lower().find('error'):
            os.chdir(old_wd)
            raise QSeisError('''===== begin qseis input =====\n%s===== end qseis input =====
===== begin qseis output =====\n%s===== end qseis output =====
===== begin qseis error =====\n%s===== end qseis error =====
non-zero exit status from qseis or the string 'error' has been found in qseis output (qseis has been invoked as "%s")''' % (qseis_input, qseis_output, qseis_error, program))
        
        self.qseis_output = qseis_output
        self.qseis_error = qseis_error
        
        os.chdir(old_wd)
        
    def get_traces(self):
        assert self.config.sw_d_unit == 1, 'can only handle distances given in km'
        assert self.config.sw_t_reduce == 1, 'can only handle t_reduce given in km/s'
       
        if self.config.sw_equidistant == 1:
            nx = self.config.no_distances
            xmin, xmax = self.config.distances_km
            xmin *= km
            xmax *= km
            if nx > 1:
                dx = (xmax-xmin)/(nx-1)
            else:
                dx = 1.0
            distances = [ xmin + ix*dx for ix in xrange(nx) ]
        else:
            distances = [ x*km for x in self.config.distances_km ]
        
        vred = self.config.t_reduce*km
        if vred == 0.0:
            vred = None
        
        fns = self.config.get_seismogram_filenames_zrt(self.tempdir)
        traces = []
        for comp, fnt in zip(('z', 'r', 't'), fns):
            fn = fnt % { 'tempdir': self.tempdir }
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
                if vred is not None:
                    tr.shift(x/vred)
                
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
            for ig in xrange(1,9):
                
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
       


class QSeisGFDBBuilder(GFDBBuilder):
    def __init__(self, partial_db_path, output_db_path, gfdb_config, block_nx, qseis_config, cutting=None, extra_traces_dir=None, tag='', tmp=None ):
        GFDBBuilder.__init__(self, partial_db_path, output_db_path, gfdb_config, block_nx, extra_traces_dir=extra_traces_dir, tmp=tmp)
        self.qseis_config = qseis_config
        self.tag = tag
        self.cutting = cutting
        self.gfmapping = [
            (MomentTensor( m=symmat6(1,0,0,1,0,0) ), {'r': (1, +1), 't': (4, +1), 'z': (6, +1) }),
            (MomentTensor( m=symmat6(0,0,0,0,1,1) ), {'r': (2, +1), 't': (5, +1), 'z': (7, +1) }),
            (MomentTensor( m=symmat6(0,0,1,0,0,0) ), {'r': (3, +1),               'z': (8, +1) }),
        ]
        
        if gfdb_config['ng'] == 10:
            self.gfmapping.append(
                (MomentTensor( m=symmat6(0,1,0,0,0,0) ), {'r': (9, +1),               'z': (10, +1) }),
            )
        
    def work_block(self, firstx, lastx, nx, z):
        traces = []
        have_gfs = False
        runner = QSeisRunner(tmp=self.tmp)
        for mt, gfmap in self.gfmapping:
            
            conf = self.qseis_config.copy()
            
            if not have_gfs:
                conf.gf_sw_source_types = 1,1,1,1,0,0
            else:
                conf.gf_sw_source_types = 0,0,0,0,0,0

            m = mt.m()
            
            # Qseis wants this ordering: Mxx Myy Mzz Mxy Myz Mzx
            conf.source_type = 1
            conf.source_vals = [m[0,0], m[1,1], m[2,2], m[0,1], m[1,2], m[0,2]]
            conf.source_depth_km = z/km
            conf.sw_equidistant = 0                # 1 is eqdist, 0 if not
            conf.sw_d_unit = 1                     # 1 is km,  0 is deg
            conf.no_distances = nx+1
            conf.distances_km = [firstx/km, lastx/km]       # first, last, if equidist
            
            distances = num.linspace(firstx, lastx, nx).tolist()
            distances_wanted = distances
            distances_km = [ dist/km for dist in distances_wanted ]
            
            onebeyond = self.gfdb_config['firstx'] + self.gfdb_config['dx'] * self.gfdb_config['nx']
            distances_km.append(onebeyond/km)
            conf.distances_km = distances_km
            conf.sw_irregular_station_azimuths = 0
            conf.station_azimuths = [ 0.0 ]
            
            runner.run(conf)
            have_gfs = True
            
            rawtraces = runner.get_traces()
            for tr in rawtraces:
                
                if tr.channel not in gfmap: continue
                (ig, factor) = gfmap[tr.channel]
                x = tr.meta['x']
                if factor != 1.0:
                    tr.ydata *= factor
                
                if self.cutting is not None:
                    tmin_cut = self.cutting[0](x,z)
                    tmax_cut = self.cutting[1](x,z)
                    tr.chop(tmin_cut, tmax_cut, inplace=True)
                
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
   
