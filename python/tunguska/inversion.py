
import receiver
import gfdb
import seismosizer
import config
import source as source_model
import util
import plotting
import moment_tensor
import orthodrome
import gridsearch
import filtering
from util import gform

import shutil
import os
import sys
import re
import math
import datetime, time, calendar
import numpy as num
from scipy.optimize import fmin_l_bfgs_b
import progressbar
import cPickle as pickle
import copy
import logging
from subprocess import call
import subprocess

from os.path import join as pjoin

def backticks(command):
    return subprocess.Popen(command, stdout=subprocess.PIPE).communicate()[0]

def convert_graph(in_filename, out_filename):
    if backticks(['which', 'pdftoppm']).strip():
        inter_filename = out_filename.replace('.png','.oversized')
        call(['pdftoppm', in_filename, inter_filename])
        for candidate in ('-000001.ppm', '-1.ppm'):   # version-dependent grrr.
            fn = inter_filename+candidate
            if os.path.exists(fn):
                call(['convert', fn, '-resize', '50%', out_filename])
                break
    else:
        logging.info("Using convert; install pdftoppm from poppler for better quality PNG output.")
        call(['convert', in_filename, out_filename])

def all(x):
    for e in x:
        if not e: return False
    return True

def ra(a):
    return a.min(), a.max()

def grow(r,*args):
    for a in args:
        if r[0] is None: 
            r[0] = a[0]
        else:
            r[0] = min(r[0],a[0])
        
        if r[1] is None:
            r[1] = a[1]
        else:
            r[1] = max(r[1],a[1])

def nonzero_range(a):
    nz, = a[1].nonzero()
    return a[0][nz.min()], a[0][nz.max()]
    
def d2u(s):
    if '_' in s: raise Exception('uuups, found underscore in param name where not expected: %s' % s)
    return s.replace('-','_')

def u2d(s):
    if '-' in s: raise Exception('uuups, found dash in param name where not expected: %s' % s)
    return s.replace('_','-')

    
def grid_defi( param, oldval, descr ):
    mi, ma, inc = [float(x) for x in descr[:3]]
    mode = 'absolute'
    if len(descr) == 4: mode = descr[3]
    if mode == 'exp':
        x = mi
        vals = []
        while x <= ma:
            vals.append(x)
            x *= inc
        return param, num.array(vals)
        
    elif mode == 'symexpinc':
        x = mi
        cinc = mi
        vals = [ 0. ]
        while x <= ma:
            vals.append(x)
            vals.append(-x)
            x += cinc
            cinc *= inc
        return param, num.array(sorted(vals))
            
    else:    
        if mode == 'mult':
            mi *= oldval
            ma *= oldval
            inc *= oldval
        elif mode == 'add':
            mi += oldval
            ma += oldval
        return param, gridsearch.mimainc_to_gvals(mi, ma, inc)
    




def standard_setup( datadir,
                    gfdb_path,
                    components,
                    effective_dt=1,
                    spacial_undersampling = [ 1, 1 ],
                    hosts = ['localhost'],
                    balance_method = '123321',
                    crustal_thickness_limit = None,
                    constraining_planes = None,
                    shifts = None,
                    blacklist = None,
                    xblacklist = None,
                    verbose = False,
    
                    local_interpolation = 'bilinear',
                    source_origin_file = 'source-origin.table',
                    receivers_file = 'receivers.table',
                    ref_seismogram_stem = 'reference',
                    ref_seismogram_format = 'mseed', **kwargs):
                    
    '''Start seismosizers and setup source location, Greens functions, 
       receivers, and reference seismograms.'''

    source_origin_file      = pjoin(datadir, source_origin_file)
    ref_seismogram_stem     = pjoin(datadir, ref_seismogram_stem)
    receivers_file          = pjoin(datadir, receivers_file)

    # setup database
    database = gfdb.Gfdb(gfdb_path)
    seis = seismosizer.Seismosizer(hosts, balance_method)
    if verbose: seis.set_verbose('T')
    seis.set_database(database)
    seis.set_effective_dt(effective_dt)
    seis.set_local_interpolation(local_interpolation)
    seis.set_spacial_undersampling(*spacial_undersampling)
    
    # setup source origin
    f = open( source_origin_file, 'r' )
    (slat, slon, stime) = [ float(x) for x in f.read().split() ]
    f.close()
    seis.set_source_location(slat, slon, stime)
    if crustal_thickness_limit is not None:
        seis.set_source_crustal_thickness_limit( crustal_thickness_limit )
    
    if constraining_planes is not None:
        values = []
        for plane in constraining_planes:
            for vect in plane:
                values.extend(vect)
                
        seis.set_source_constraints( *values )
    
    # setup receivers
    receivers = receiver.load_table(receivers_file, set_components=components)
    
    if len(receivers) == 0: sys.exit('no receivers')
    
    seis.set_receivers(receivers)
    
    seis.set_ref_seismograms( ref_seismogram_stem, ref_seismogram_format )
    if blacklist:
        seis.blacklist_receivers( blacklist )
    if xblacklist:
        seis.xblacklist_receivers( xblacklist )
        
    # apply reference seismograms shifts
    if shifts is not None:
        seis.shift_ref_seismograms( shifts )
    
    return seis
    
standard_setup.required = set(('datadir', 'gfdb_path', 'components')) 
standard_setup.optional = set(('effective_dt', 'spacial_undersampling', 'hosts', 'balance_method',
                               'crustal_thickness_limit', 'constraining_planes', 'shifts', 'local_interpolation', 
                               'source_origin_file', 'receivers_file',
                               'ref_seismogram_stem', 'ref_seismogram_format', 'blacklist', 'xblacklist', 'verbose'))

def gen_dweights( seis, base_source, datadir,
                                 ref_seismogram_stem = 'reference',
                                 ref_seismogram_format = 'mseed', **kwargs):
    
    
    base_source = copy.deepcopy(base_source)
    moment = base_source['moment']
    base_source['moment'] = 0.0
    seis.set_source(base_source)
    seis.set_synthetic_reference()    
    
    strike_grid    = ('strike',  -180., 150., 30.)
    dip_grid       = ('dip',  0., 90., 30.)
    slip_rake_grid = ('slip-rake', -180., 150., 30.)
    
    base_source['moment'] = moment
    sdr_grid = gridsearch.MisfitGrid( base_source, [strike_grid, dip_grid, slip_rake_grid])
    
    sdr_grid.compute( seis )
    
    means = sdr_grid.get_mean_misfits_by_r()
    means /= num.mean(means[means>0.])
    dweights = num.where(means>0., 1./means, 0.)
    # reset reference seismograms
    if ref_seismogram_stem is not None:
        ref_seismogram_stem     = pjoin(datadir, ref_seismogram_stem)
        seis.set_ref_seismograms( ref_seismogram_stem, ref_seismogram_format )
    
    return dweights
    
gen_dweights.required = set(('datadir',))
gen_dweights.optional = set(('ref_seismogram_stem', 'ref_seismogram_format'))

class Step:
    inner_misfit_method_params = set(('inner_norm', 'taper', 'filter', 'nsets', 'depth'))
    outer_misfit_method_params = set(('outer_norm', 'bootstrap_iterations', 'anarchy', 'receiver_weights'))
    
    def __init__(self, workdir, name, dump_processing='filtered'):
        self.baseworkdir = workdir
        self.stepname = name
        self.in_config = None
        self.out_config = None
        self.stepdir = pjoin(self.baseworkdir, self.stepname)
        self.seismosizer = None
        self.required = set(standard_setup.required)
        self.optional = set(standard_setup.optional)
        self.dump_processing = dump_processing
        
    def make_rundir_path(self, run_id):
        return pjoin(self.stepdir, str(run_id))
    
    def make_plotdir_path(self, run_id):
        return pjoin(self.stepdir, str(run_id), 'plots')
        
    def next_available_rundir(self):
        entries = os.listdir(self.stepdir)
        is_int = re.compile('^\d+$')
        int_entries = [ int(x) for x in entries if is_int.search(x) ]
        if int_entries:
            inew = max(int_entries)+1
        else:
            inew = 1
            
        return pjoin(self.stepdir, '%03i' % inew)
    
    def get_config(self):
        if not self.out_config:
            # step has not been run, get config from last successful run
            rundir = self.make_rundir_path('current')
            c = config.Config( pjoin( rundir, 'config-out.pickle' ) )
            return c.get_config()
        else:
            return self.out_config.get_config()
    
    def pre_work(self, start_seismosizer=True):
        assert(self.in_config is not None)
        
        have = set(self.in_config.get_config().keys())
        for k in self.required - have:
            logging.warn('Required parameter missing for step %s: %s' % (self.stepname, k))
            
        #for k in have - (self.optional | self.required):
        #    logging.info('unused parameter in config for step %s: %s' % (self.stepname, k))
            
        
        logging.info('Starting work on step %s' % self.stepname)
        rundir = self.make_rundir_path('incomplete')
        if os.path.exists( rundir ): 
            shutil.rmtree( rundir )
        os.makedirs( rundir )
        self.in_config.dump( pjoin(rundir,'config-in.pickle') )
        self.out_config = config.Config()
        
        self.work_started = time.time()
        
        datadir = self.in_config.get_config()['datadir']
        ref_time_filename = pjoin(datadir,'reference-time.txt')
        self.set_ref_time(ref_time_filename)
        
        if start_seismosizer:
            sconf = self.in_config.get_config(keys=standard_setup.required|standard_setup.optional)
            self.seismosizer = standard_setup( **sconf )
            self.out_config.source_location = self.seismosizer.source_location
        
    def set_ref_time(self, filename):
        f = open(filename, 'r')
        dtstr = f.read().strip().split(None,1)[1]
        f.close()
        format = '%Y/%m/%d %H:%M:%S'
        self.ref_time = calendar.timegm(time.strptime(dtstr, format))
        
        
    def setup_inner_misfit_method(self):
        conf = self.in_config.get_config(keys=Step.inner_misfit_method_params)
        seis = self.seismosizer
        
        tapers_by_set = conf['taper']
        assert(len(tapers_by_set) == conf['nsets'])
        tapers = [ tapers_by_set[i%len(tapers_by_set)] for i in range(len(seis.receivers)) ]
        
        seis.set_taper(tapers, conf['depth'])
        if conf['filter']:
            seis.set_filter(conf['filter'])
            
        seis.set_misfit_method(conf['inner_norm'])
        
    def post_work(self, stop_seismosizer=True):
        
        self.make_alternative_stats()
        self.make_locations()

        if stop_seismosizer:
            self.seismosizer.close()
        
        rundir = self.make_rundir_path('incomplete')
        self.out_config.dump(pjoin(rundir, 'config-out.pickle'))
        if os.path.exists( self.make_rundir_path('current') ):
            shutil.move(self.make_rundir_path('current'), self.next_available_rundir())
        shutil.move(rundir, self.make_rundir_path('current'))
        
        time_wasted = time.time() - self.work_started
        logging.info('Time wasted on step %s: %s' % (self.stepname, datetime.timedelta(seconds=round(time_wasted)) ))
        
        logging.info('Done with work on step %s' % self.stepname)
        
    def make_alternative_stats(self):
        oc = self.out_config
        c = self.out_config.__dict__
        
        if 'moment_stats' in c:
            oc.magnitude_stats = oc.moment_stats.converted('magnitude', moment_tensor.moment_to_magnitude)
            
        if ('north_shift_stats' in c) and ('east_shift_stats' in c):
            
            rlat, rlon, rtime = self.seismosizer.source_location
            cnorth, ceast =  oc.north_shift_stats.best, oc.east_shift_stats.best
            clat, clon = orthodrome.ne_to_latlon( rlat, rlon, cnorth, ceast )
            
            approx_lat = lambda north: clat + (north-cnorth) * 180. / (config.earthradius*math.pi)
            approx_lon = lambda east: clon + (east-ceast) * 180. / (config.earthradius*math.pi*math.cos(clat*math.pi/180.))
            
            oc.latitude_stats = oc.north_shift_stats.converted( 'latitude', approx_lat )
            oc.longitude_stats = oc.east_shift_stats.converted( 'longitude', approx_lon )
            
    def make_locations(self):
        oc = self.out_config
        icd = self.in_config.get_config()
        ocd = self.out_config.__dict__
        
        
        if self.seismosizer:
            rlat, rlon, rtime = self.seismosizer.source_location
            if 'north_shift' in ocd:
                cnorth = oc.north_shift
            elif 'north_shift' in icd:
                cnorth = icd['north_shift']
            else:
                cnorth = 0.
                
            if 'east_shift' in ocd:
                ceast = oc.east_shift
            elif 'east_shift' in icd:
                ceast = icd['east_shift']
            else:
                ceast = 0.
                
            if 'time' in ocd:
                ctime = oc.time
            elif 'time' in icd:
                ctime = icd['time']
            else:
                ctime = 0.
                
            clat, clon = orthodrome.ne_to_latlon( rlat, rlon, cnorth, ceast )
    
            oc.centroid_latitude = clat
            oc.centroid_longitude = clon
            oc.centroid_time = self.ref_time + ctime
            oc.reference_time = self.ref_time
        

    def snapshot(self, source, ident, mm_conf):
        
        if not self.seismosizer:
            logging.warn('Cannot create snapshot, because no seismosizers are running')
            return
        
        self.seismosizer.set_source( source )
        receivers = self.seismosizer.get_receivers_snapshot(which_processing=self.dump_processing )
        for irec, rec in enumerate(receivers):
            if 'receiver_weights' in mm_conf:
                rec.weight = mm_conf['receiver_weights'][irec]
            else:
                rec.weight = 1.0
        
        self.dump( source, 'snapshot_source_%s' % ident )
        self.dump( receivers, 'snapshot_%s' % ident )
        self.dump( self.seismosizer.get_source_location(), 'snapshot_source_location_%s' % ident )
        
        self.dump( self.seismosizer.get_psm_infos(), 'source_infos_%s' % ident )
        rundir = self.make_rundir_path('incomplete')
        tracesdir = pjoin(rundir, 'snapshot_%s' % ident)
        if os.path.exists(tracesdir):
            shutil.rmtree(tracesdir)
        
        os.mkdir(tracesdir)
        for r in receivers:
            r.save_traces_mseed(pjoin(tracesdir, '%(whichset)s_%(network)s_%(station)s_%(channel)s.mseed'))
        
    def get_snapshot(self, ident, run_id='current'):
        return self.load( ident='snapshot_%s' % ident, run_id=run_id )

    def get_snapshot_source_location(self, ident, run_id='current'):
        return self.load( ident='snapshot_source_location_%s' % ident, run_id=run_id )

    def get_snapshot_source_infos(self, ident, run_id='current'):
        return self.load( ident='source_infos_%s' % ident, run_id=run_id)

    def get_snapshot_source(self, ident, run_id='current'):
        return self.load( ident='snapshot_source_%s' % ident, run_id=run_id)

    def dump(self, object, ident, run_id='incomplete'):
        rundir = self.make_rundir_path(run_id)
        filename = pjoin(rundir, '%s.pickle' % ident)
        f = open(filename, 'w')
        pickle.dump(object, f)
        f.close()
        
    def result(self, string, ident, run_id='incomplete'):
        rundir = self.make_rundir_path(run_id)
        filename = pjoin(rundir, '%s.result' % ident)
        f = open(filename, 'w')
        f.write( '%s\n' % string )
        f.close()

    def load(self, ident, run_id='current'):
        rundir = self.make_rundir_path(run_id)
        filename = pjoin(rundir, '%s.pickle' % ident)
        f = open(filename, 'r')
        object = pickle.load(f)
        f.close()
        return object

    def plot(self, run_id='current'):
         logging.info('Starting plotting for step: %s' % self.stepname)
         
         plotdir = self.make_plotdir_path(run_id)
         if os.path.exists( plotdir ): 
            shutil.rmtree( plotdir )
         os.makedirs( plotdir )
         
         plot_files = self._plot( run_id )
         
         # convert plots to png
         for infn in plot_files:
             inpath = pjoin( plotdir, infn )
             fnbase = os.path.splitext(infn)[0]
             outfn = fnbase+'.png'
             outpath = pjoin( plotdir, outfn )
             convert_graph( inpath, outpath )
             
         logging.info('Done with plotting for step: %s' % self.stepname)

    def _plot(self, rundir_id):
        logging.info('Nothing to do')
        return []

    def work(self, *args, **kwargs):
        pass
    
    
    def ic(self, run_id='current'):
        try:
            c = self.load('config-in', run_id=run_id)
        except:
            c = None
        return c
    
    def oc(self, run_id='current'):
        try:
            c = self.load('config-out', run_id=run_id)
        except:
            c = None
        return c
            
    def ocf(self, fmt, mult=1.0, run_id='current'):
        
        c = self.oc(run_id)
        
        if mult != 1.0:
            c = copy.deepcopy(c)
            for k in c.keys():
                v = c[k]
                if isinstance(v,float) or isinstance(v,int) or isinstance(v,long):
                    c[k] = mult*v
                    
        return fmt % c
        
            
    def gx(self, arg):
        plotdir = self.make_plotdir_path( 'current' )
        location = pjoin(plotdir, arg)
        img_snippet = '<a href="%s"><img style="border:none;" src="%s" /></a>'
        return img_snippet % (location+'.pdf', location+'.png')
    
    def gxi(self, arg, format='png'):
        plotdir = self.make_plotdir_path( 'current' )
        location = pjoin(plotdir, arg)
        return location+('.%s' % format)

    def show_in_config(self, run_id='current'):
        c = self.load('config-in', run_id=run_id)
        for k in sorted(c.keys()):
            print k, '=', repr(c[k])
            
    def show_out_config(self, run_id='current'):
        c = self.load('config-out', run_id=run_id)
        for k in sorted(c.keys()):
            print k, '=', repr(c[k])
            
    def show_active_config(self, run_id='current'):
        c = self.in_config.get_config()
        for k in sorted(c.keys()):
            print k, '=', repr(c[k])


class Informer(Step):
    
    def __init__(self, workdir, name='weightmaker'):
        Step.__init__(self, workdir, name)
        
    
    def work(self, **kwargs):
        self.pre_work(True)        
        seis = self.seismosizer
        conf = self.in_config.get_config()
        
        distances_m = num.array([ r.distance_m for r in seis.receivers ], dtype=num.float)
        imin = num.argmin(distances_m)
        imax = num.argmax(distances_m)
        
        def sx(r):
            return '%10s  %s km   %s deg' % (r.name, gform(r.distance_m/1000.,4), gform(r.distance_deg,3))
        
        nsta = len(distances_m)/conf['nsets']
        
        print 'closest station: %s' % sx(seis.receivers[imin])
        self.out_config.__dict__['closest_station'] = sx(seis.receivers[imin])
        print 'farthest station:  %s' % sx(seis.receivers[imax])
        self.out_config.__dict__['farthest_station'] = sx(seis.receivers[imax])
        print 'number of stations: %i' % nsta
        self.out_config.__dict__['nstations'] = nsta
        
        save = {}
        save['receivers'] = seis.receivers
        save['source_location'] = seis.source_location
        
        self.dump(save, 'source_receivers')
        

        self.post_work(True)
        
        #if nsta < 3:
        #    sys.exit('too few stations')
    

    def _plot(self, run_id='current'):

        saved = self.load('source_receivers', run_id=run_id)
        conf = self.in_config.get_config()
        
        #
        # Station map
        #
        receivers = saved['receivers']
        source_location = saved['source_location']
               
        lats = num.array( [ r.lat for r in receivers ], dtype='float' )
        lons = num.array( [ r.lon for r in receivers ], dtype='float' )
        dists = num.array( [ r.distance_deg for r in receivers ], dtype='float' )
        rnames = [ ' '.join(r.name.split('.')) for r in receivers ]
        slat, slon = source_location[:2]
        
        station_size = [ 0.05 ] * len(receivers)
        station_color = [ (0.,1.)[r.enabled] for r in receivers ]
        
        plotdir = self.make_plotdir_path(run_id)
        source = None
        
        fn = pjoin(plotdir, 'stations.pdf')
        plotting.station_plot( slat, slon, lats, lons, rnames, station_color, station_size, 
                               source, num.amax(dists)*1.05, fn, {}, zexpand=1.02, nsets=conf['nsets'], symbols=('t','t'))
        
        #
        # Region map
        #
        fn = pjoin(plotdir, 'region.pdf')
        plotting.location_map( fn, slat, slon, 5., {}, receivers=(lats, lons, rnames))
        
        
        return [ 'stations.pdf', 'region.pdf' ]
        
        
class WeightMaker(Step):
    
    def __init__(self, workdir, name='weightmaker'):
        Step.__init__(self, workdir, name)
        
        self.required |= Step.inner_misfit_method_params | gen_dweights.required \
                        | set(('depth', 'moment', 'rise_time')) 
        
        self.optional |=  gen_dweights.optional
    
    def work(self, **kwargs):
        self.pre_work(True)        
        self.setup_inner_misfit_method()
        seis = self.seismosizer
        conf = self.in_config.get_config()
        
        sourcetype = 'eikonal'
        base_source = source_model.Source( sourcetype, 
                               {"time": float(conf['time']),
                                "depth": float(conf['depth']),
                                "bord-radius": 0.,
                                "moment": float(conf['moment']),
                                "rise-time": float(conf['rise_time']) } )
        
        self.out_config.receiver_weights = gen_dweights(seis, base_source, **conf )
        self.post_work(True)


class EffectiveDtTester(Step):
    
    def __init__(self, workdir, name='effective_dt_tester'):
        Step.__init__(self, workdir, name)
        
        self.required |= Step.inner_misfit_method_params | gen_dweights.required \
                        | set(('depth', 'moment', 'rise_time')) 
        
        self.optional |=  gen_dweights.optional
    
    def work(self, **kwargs):
        self.pre_work(True)        
        self.setup_inner_misfit_method()
        seis = self.seismosizer
        conf = self.in_config.get_config()
        
        
        sourcetype = 'eikonal'
        base_source = source_model.Source( sourcetype, 
                               {"depth": float(conf['depth']),
                                "bord-radius": float(conf['bord_radius']),
                                "moment": float(conf['moment']),
                                "rise-time": float(conf['rise_time']) } )
                                
        ref_seismogram_stem = 'reference'
        ref_seismogram_format = 'sac'
        
        seis.set_source(base_source)
        seis.set_synthetic_reference()
        
        for i in range(20):
            effdt = i*0.25 + 0.5
            seis.set_effective_dt(effdt)
            seis.make_misfits_for_source( base_source )
            ms = []
            for ireceiver, receiver in enumerate(seis.receivers):
                if receiver.enabled:
                    for icomp, comp in enumerate(receiver.components):
                        ms.append( receiver.misfits[icomp]/receiver.misfit_norm_factors[icomp] )
        
        
        datadir = conf['datadir']
        # reset reference seismograms
        ref_seismogram_stem     = pjoin(datadir, ref_seismogram_stem)
        seis.set_ref_seismograms( ref_seismogram_stem, ref_seismogram_format )
        
        self.post_work(True)


class Shifter(Step):
    def __init__(self, workdir, name='shifter'):
        Step.__init__(self, workdir, name)
        
        self.required |= set(('taper', 'filter', 'autoshift_range', 'autoshift_limit'))
                        
        self.optional |=  set([d2u(p) for p in source_model.param_names('eikonal')])
        
    def work(self, **kwargs):
        self.pre_work(True)
        seis = self.seismosizer
        conf = self.in_config.get_config()
        
        sourcetype = 'eikonal'
        base_source = source_model.Source( sourcetype )
        
        for p in source_model.param_names(sourcetype):
            if d2u(p) in conf:
                base_source[p] = float(conf[d2u(p)])
        
        seis.set_source(base_source)
        
        tapers_by_set = conf['taper']
        tapers = [ tapers_by_set[i%len(tapers_by_set)] for i in range(len(seis.receivers)) ]
        seis.set_taper(tapers, base_source['depth'])
        seis.set_filter(conf['filter'])
        
        shifts = seis.autoshift_ref_seismograms( conf['autoshift_range'] )
        fails = []
        limit = conf['autoshift_limit']
        for i,s in enumerate(shifts):
            if s < limit[0] or limit[1] < s:
                seis.shift_ref_seismograms( [-s], [i+1] )
                shifts[i] = 0.
                fails.append(True)
            else:
                fails.append(False)
           
        save = {}
        save['receivers'] = seis.receivers
        save['source_location'] = seis.source_location
        save['source'] = base_source
        save['shifts'] = shifts
        save['fails'] = fails
        self.dump(save, 'source_receivers')
        
        shifts_debug = []
        for i,s in enumerate(shifts):
            trace = {}
            trace['name'] = seis.receivers[i].name
            trace['shift'] = s
            trace['failed'] = fails[i]
            shifts_debug.append(trace)
        
        self.out_config.shifts = shifts
        self.out_config.shifts_debug = shifts_debug
        self.post_work(True)
        
    def _plot( self, run_id='current' ):
                
        saved = self.load('source_receivers', run_id=run_id)
        conf = self.in_config.get_config()
        #
        # Station plot
        #
        receivers = saved['receivers']
        source_location = saved['source_location']
        source = saved['source']
        shifts = saved['shifts']
        fails = saved['fails']
        
        lats = num.array( [ r.lat for r in receivers ], dtype='float' )
        lons = num.array( [ r.lon for r in receivers ], dtype='float' )
        dists = num.array( [ r.distance_deg for r in receivers ], dtype='float' )
        rnames = [ ' '.join(r.name.split('.')) for r in receivers ]
        slat, slon = source_location[:2]
        station_color = []
        for s,f in zip(shifts,fails):
            if f:
                station_color.append(num.NaN)
            else:
                station_color.append(s)
        
        station_size = [ 1. ] * len(receivers)
        
        plotdir = self.make_plotdir_path(run_id)

        fn = pjoin(plotdir, 'stations.pdf')
        plotting.station_plot( slat, slon, lats, lons, rnames, station_color, station_size, 
                               source, num.amax(dists)*1.05, fn, {}, zexpand=1.02, nsets=conf['nsets'])
               
        return [ 'stations.pdf' ]


class ExtConfigurator(Step):
    
    def __init__(self, workdir, name, generate=('filter', 'constraining_planes', 
                    'bord_radius_range', 'nukl_shift_x_range', 'nukl_shift_y_range' ), size_factor=4000., steps=5.):
        Step.__init__(self, workdir, name)
        self.required |= Step.inner_misfit_method_params | set(('depth', 'rise_time'))
        self.generate = generate
        self.size_factor = size_factor
        self.steps = steps
        
    def work(self, **kwargs):
        self.pre_work(False)
        oc = self.out_config
        ic = self.in_config.get_config()
        
        rise_time = ic['rise_time']
        depth = ic['depth']
        
        if 'filter' in self.generate:
            filter = ic['filter']
            filter.set(2, 1./(rise_time*2./3.))
            filter.set(3, 1./(rise_time*1./2.))
            oc.filter = filter

        maxradius = self.size_factor*rise_time
        if 'bord_radius_range' in self.generate:
            oc.bord_radius_range     = 0., maxradius, maxradius/self.steps
        
        if 'nukl_shift_x_range' in self.generate:
            oc.nukl_shift_x_range    = -maxradius, maxradius, maxradius/self.steps
        
        if 'nukl_shift_y_range' in self.generate:
            oc.nukl_shift_y_range    = -maxradius, maxradius, maxradius/self.steps

        if 'constraining_planes' in self.generate:
            cp = ic['constraining_planes']
    
            oc.constraining_planes   = [((0.,0.,cp[0][0][2]),(0.,0.,-1.)),
                                        ((0.,0.,min(depth*2.,cp[1][0][2])),(0.,0.,+1.))]
                                        
        self.post_work(False)

class ParamTuner(Step):
     
    def __init__(self, workdir, sourcetype='eikonal', params=['time'], name=None, 
            xblacklist_level=None, dump_processing='filtered', ref_source_from=None):
        if name is None: name = '-'.join(params)+'-tuner'
        Step.__init__(self, workdir, name, dump_processing)
        self.sourcetype = sourcetype
        self.params = params
        self.xblacklist_level = xblacklist_level
        self.ref_source_from = ref_source_from
        
        self.required |= Step.outer_misfit_method_params | Step.inner_misfit_method_params \
                        | set([param+'_range' for param in self.params]) \
                        | set(self.params)
                        
        self.optional |= set([d2u(p) for p in source_model.param_names(self.sourcetype)])
        
    def work(self, search=True, forward=True, run_id='current'):
        self.pre_work(search or forward)
        seis = self.seismosizer
        conf = self.in_config.get_config()
        mm_conf = self.in_config.get_config(keys=Step.outer_misfit_method_params)
        
        sourcetype = self.sourcetype
        base_source = source_model.Source( sourcetype )
        
        for p in source_model.param_names(sourcetype):
            if d2u(p) in conf:
                base_source[p] = float(conf[d2u(p)])
                
        if 'plane' in conf and conf['plane'] == 2: 
            strike, dip, slip_rake = float(conf['strike']), float(conf['dip']), float(conf['slip_rake'])
            strike, dip, slip_rake = moment_tensor.other_plane( strike, dip, slip_rake )
            base_source['strike'] = strike
            base_source['dip'] = dip
            base_source['slip-rake'] = slip_rake
            
        if 'plane' in conf:
            for param in 'strike', 'dip', 'slip-rake':
                self.out_config.__dict__['active_'+d2u(param)] = base_source[param]
                
        grid_def = []
        for param in self.params:
            oldval = base_source[u2d(param)]
            descr = conf[param+'_range']            
            grid_def.append(grid_defi(u2d(param),oldval,descr))
            
        if search or forward: self.setup_inner_misfit_method()
        if search:
            if self.ref_source_from is not None:
                step, ident = self.ref_source_from
                ref_source = step.get_snapshot_source(ident)
            else:
                ref_source = None
            
            self.setup_inner_misfit_method()
            finder = gridsearch.MisfitGrid( base_source, param_values=grid_def, ref_source=ref_source)
            finder.compute(seis)
        else:
            finder = self.load(self.stepname, run_id=run_id)
            
        finder.postprocess(**mm_conf)
        self.dump(finder, self.stepname)
        
        stats = finder.stats
        for param in self.params:
            str_result = stats[u2d(param)].str_best_and_confidence()
            logging.info(str_result)
            self.result(str_result, param )
            base_source[u2d(param)] = stats[u2d(param)].best
            self.out_config.__dict__[param] = stats[u2d(param)].best
            self.out_config.__dict__[param+'_stats'] = stats[u2d(param)]
        
        self.out_config.min_misfit = finder.get_best_misfit()
        self.out_config.nstations_total  = finder.nreceivers
        self.out_config.nstations_used  = finder.nreceivers_enabled
        
        logging.info('Misfit = %f, Source = %s' % (finder.get_best_misfit(), str(base_source)))
        mt = base_source.moment_tensor()
        logging.info(str(mt))
        
        self.result(str(mt), 'moment_tensor')
        
        misfit_median = finder.get_median_of_misfits_by_r()
        if self.xblacklist_level is not None:
            ir = 0
            
            if 'xblacklist' in conf:
                xblacklist = set(conf['xblacklist'])
            else:
                xblacklist = set()
                
            for r, mm in zip(seis.receivers, finder.misfits_by_r):
                if mm/misfit_median > self.xblacklist_level:
                    xblacklist.add(ir) 
                    logging.info('Blacklisting:  %i, %s, %g' % (ir+1, r.name, mm/misfit_median))
                    
                ir += 1
            self.out_config.xblacklist = sorted(list(xblacklist))
                    
        if forward:
            self.snapshot( base_source, 'best', mm_conf )
            
        self.post_work(search or forward)
        
    def _plot( self, run_id='current' ):
        plot_files = []
        
        conf = self.in_config.get_config()
        plotdir = self.make_plotdir_path(run_id)
        
        source_infos = self.get_snapshot_source_infos('best', run_id=run_id)
        
        fn = pjoin(plotdir, 'rupture.pdf')
        plotting.rupture_plot(fn, source_infos)
        plot_files.append('rupture.pdf')
        
        finder = self.load(self.stepname, run_id=run_id)
        plot_files.extend(finder.plot( plotdir, conf['nsets'], source_model_infos = source_infos))
        
        return plot_files
    
class EnduringPointSource(Step):
    
    def __init__(self, workdir, name='extension'):
        Step.__init__(self, workdir, name)
        
        self.params = ('rise_time',)
        
        self.required |= Step.outer_misfit_method_params | Step.inner_misfit_method_params \
                        | set([param+'_range' for param in self.params]) \
                        | set(self.params)
                        
        self.optional |= set([d2u(p) for p in source_model.param_names('eikonal')])
    
    def work(self, search=True, forward=True, run_id='current'):
        self.pre_work(search or forward)
        seis = self.seismosizer
        conf = self.in_config.get_config()
        mm_conf = self.in_config.get_config(keys=Step.outer_misfit_method_params)
        
        sourcetype = 'eikonal'
        
        base_source = source_model.Source( sourcetype, {} )
    
        for p in source_model.param_names(sourcetype):
            if d2u(p) in conf:
                base_source[p] = float(conf[d2u(p)])
                
        grid_def = []
        for param in self.params:
            oldval = base_source[u2d(param)]
            descr = conf[param+'_range']            
            grid_def.append(grid_defi(u2d(param),oldval,descr))
            
        if search or forward: self.setup_inner_misfit_method()
        if search:
            self.setup_inner_misfit_method()
            finder = gridsearch.MisfitGrid( base_source, param_values=grid_def )
            finder.compute(seis)
        else:
            finder = self.load(self.stepname, run_id=run_id)
            
        finder.postprocess(**mm_conf)
        self.dump(finder, self.stepname)
        
        stats = finder.stats
        for param in self.params:
            altparam = param
            if param == 'rise_time': altparam = 'duration'
            str_result = stats[u2d(param)].str_best_and_confidence()
            logging.info(str_result)
            base_source[u2d(param)] = stats[u2d(param)].best        
        
        logging.info('Reweighting')
        xweights = num.where(finder.misfits_by_r != 0.,  1./finder.misfits_by_r, 0.)
        
        mm_conf_copy = copy.copy(mm_conf)
        mm_conf_copy['receiver_weights'] = xweights*mm_conf['receiver_weights']
        finder.postprocess(**mm_conf_copy)
        self.dump(finder, self.stepname)
        
        stats = finder.stats
        for param in self.params:
            altparam = param
            if param == 'rise_time': altparam = 'duration'
            str_result = stats[u2d(param)].str_best_and_confidence()
            logging.info(str_result)
            self.result(str_result, altparam )
            base_source[u2d(param)] = stats[u2d(param)].best
            self.out_config.__dict__[altparam] = stats[u2d(param)].best
            self.out_config.__dict__[altparam+'_stats'] = stats[u2d(param)]
        
        self.out_config.receiver_weights = mm_conf_copy['receiver_weights']
        self.out_config.best_point_source = base_source
        
        if forward:
            self.snapshot( base_source, 'best', mm_conf )
            
        self.post_work(search or forward)
        
    def _plot( self, run_id='current' ):
        
        conf = self.in_config.get_config()

        plotdir = self.make_plotdir_path(run_id)
        
        finder = self.load(self.stepname, run_id=run_id)
        return finder.plot( plotdir, conf['nsets'])
                

class TracePlotter(Step):
    
    def __init__(self, workdir, snapshots, name='traceplotter'):
        Step.__init__(self, workdir, name)
        self.snapshots = snapshots
        self.required |= set()
        self.optional |= set()
        
    def work(self, search=True, forward=True, run_id='current'):
        self.pre_work(False)
        
        conf = self.in_config.get_config()
        
        # make private copy of snapshots 
        loaded_snapshots = []
        for step, ident in self.snapshots:
            snapshot = step.get_snapshot(ident, run_id='current')
            loaded_snapshots.append(snapshot)
            self.dump( snapshot, 'snapshot_%s_%s' % (step.stepname, ident) )
            sloc = step.get_snapshot_source_location(ident, run_id='current')
            self.dump( sloc, 'snapshot_source_location_%s_%s' % (step.stepname, ident) )
        
        nsets = conf['nsets']
        set_names = conf['set_names']
        plots = []
        for irec, recs in enumerate(zip(*loaded_snapshots)):
            proto = recs[0]
            name = '????'
            if all([r.name == proto.name for r in recs]):
                name = r.name
                
            if all([len(r.components) == 0 for r in recs]): continue
            
            plots.append( { 'name': name, 'number': irec+1, 'set': set_names[irec%nsets] } )
            
            
        self.out_config.__dict__['plots'] = plots
            
        self.post_work(False)
        
    def _plot(self, run_id='current'):
    
        plotdir = self.make_plotdir_path(run_id)
        
        loaded_snapshots = []
        loaded_snapshots_source_locations = []
        for step, ident in self.snapshots:
            loaded_snapshots.append(self.get_snapshot("%s_%s" % (step.stepname, ident)))
            loaded_snapshots_source_locations.append(self.get_snapshot_source_location("%s_%s" % (step.stepname, ident)))
        
        allfilez = plotting.multi_seismogram_plot2( loaded_snapshots, loaded_snapshots_source_locations,  plotdir )
        
        return [ os.path.basename(fn) for fn in allfilez ]

def snap(val,coords):
    ipos = num.argmin(num.abs(val-coords))
    return ipos, coords[ipos]

class Greeper(Step):
    '''Grid search over gradient searches
    
    aaahrrggg! this has to be rewritten!'''
     
    def __init__(self, workdir, sourcetype='eikonal', params=['time'], name=None):
        if name is None: name = '-'.join(params)+'-greeper'
        Step.__init__(self, workdir, name)
        self.params = params
        self.sourcetype = sourcetype
        
        self.required |= Step.outer_misfit_method_params | Step.inner_misfit_method_params \
                        | set([param+'_range' for param in self.params]) \
                        | set(self.params)
                        
        self.optional |= set([d2u(p) for p in source_model.param_names(sourcetype)])
        
    def work(self, search=True, forward=True, run_id='current'):
        self.pre_work(search or forward)
        seis = self.seismosizer
        conf = self.in_config.get_config()
        mm_conf = self.in_config.get_config(keys=Step.outer_misfit_method_params)
        
        sourcetype = self.sourcetype
        base_source = source_model.Source( sourcetype )
        
        for p in source_model.param_names(sourcetype):
            if d2u(p) in conf:
                base_source[p] = float(conf[d2u(p)])
                
        if 'plane' in conf and conf['plane'] == 2: 
            strike, dip, slip_rake = float(conf['strike']), float(conf['dip']), float(conf['slip_rake'])
            strike, dip, slip_rake = moment_tensor.other_plane( strike, dip, slip_rake )
            base_source['strike'] = strike
            base_source['dip'] = dip
            base_source['slip-rake'] = slip_rake
            
        if 'plane' in conf:
            for param in 'strike', 'dip', 'slip-rake':
                self.out_config.__dict__['active_'+d2u(param)] = base_source[param]
        
        # setup grid of starting points for gradient se
        starter_grid_def = []
        for param in self.params:
            if param+'_start_range' in conf:
                oldval = base_source[u2d(param)]
                descr = conf[param+'_start_range']
                starter_grid_def.append( grid_defi(u2d(param),oldval,descr) )
        
        if starter_grid_def:
            starter_sources = base_source.grid( starter_grid_def )
        else:
            starter_sources = [ base_source ]
        
        grid_def = []
        for param in self.params:
            oldval = base_source[u2d(param)]
            descr = conf[param+'_range']
            grid_def.append(grid_defi(u2d(param),oldval,descr))
        
        dparams = [ u2d(param) for param in self.params ]
        self.nminfunccalls = 0

        norms = [ num.min(vals[1:]-vals[:-1]) for (param, vals) in grid_def ]
        bounds = [ (num.min(vals)/norm,num.max(vals)/norm) for  ((param, vals),norm) in zip(grid_def, norms) ]

        if search or forward: self.setup_inner_misfit_method()
        
        def x_to_source(x):
            source = copy.deepcopy(base_source)
            for i,param in enumerate(self.params):
                source[u2d(param)] = x[i]*norms[i]
            return source
        
        def source_to_x(source):
            x = []
            for param, norm in zip(self.params, norms):
                x.append(source[u2d(param)]/norm)
            return x
        
        
        def minfunc(x):
            self.nminfunccalls += 1
            source = x_to_source(x)
            try:
                logging.debug('Evaluating source: %s', source.pretty_str(dparams) )
                best, misfit = seis.best_source([source], **mm_conf)
                logging.debug('Misfit: %g' % misfit)
            except NoValidSources:
                logging.warn('Gradient search aborted at invalid source:' )
                for str_param in current_base.pretty_str(dparams).splitlines():
                    logging.warn( str_param )
                raise
            return misfit
            
        if search:
            self.setup_inner_misfit_method()
                        
            # fix depth range by trying out different depths
            for iparam, (param, vals) in enumerate(grid_def):
                if param == 'depth':
                    depth_sources = base_source.grid( [ (param,vals) ] )
                    dummy_source, dummy_misfit, failings = seis.best_source(depth_sources, return_failings=True, **mm_conf)
                    ok = []
                    for isource, source in enumerate(depth_sources):
                        if isource not in failings:
                            ok.append(source['depth'])
                    miok = min(ok)
                    maok = max(ok)
                    bounds[iparam] = ((miok+gridsearch.step_at(ok,miok)*0.3)/norms[iparam],
                                      (maok-gridsearch.step_at(ok,maok)*0.3)/norms[iparam])
            
            # grid search over gradient searches
            min_misfit = None
            very_best_source = None
            ngood = 0
            ntotal = 0
            for starter_source in starter_sources:
                ntotal += 1
                # look if starting source is valid
                current_base = copy.deepcopy(starter_source)
                try:
                    current_base, misfit = seis.best_source([current_base], **mm_conf)
                except NoValidSources:
                    logging.warn('Skipping invalid starting source:')
                    for str_param in current_base.pretty_str(dparams).splitlines():
                        logging.warn( str_param )
                    continue
                if min_misfit == None:
                    min_misfit = misfit
                    very_best_source = current_base
                
                x0 = source_to_x(current_base)
                try:
                    x, misfit, d = fmin_l_bfgs_b(minfunc, x0, approx_grad=True, bounds=bounds, epsilon=0.2, factr=1e10)
                except NoValidSources:
                    continue
                current = x_to_source(x)
                if d['warnflag'] != 0: continue
                
                # restart gradient search at current minimum
                
                x0 = source_to_x(current)
                try:
                    x, misfit, d = fmin_l_bfgs_b(minfunc, x0, approx_grad=True, bounds=bounds, epsilon=0.05, factr=1e7)
                except NoValidSources:
                    continue
                current = x_to_source(x)
                if d['warnflag'] != 0: continue

                
                minkind = ''
                if misfit < min_misfit:
                    logging.info('Found possible minimum: Misfit = %f' % misfit)
                    for str_param in current.pretty_str(dparams).splitlines():
                        logging.info( str_param )
                    if 'strike' in self.params:
                        mt = moment_tensor.MomentTensor(strike=current['strike'], dip=current['dip'], rake=current['slip-rake'])
                        fps1, fps2 = mt.str_fault_planes().strip().splitlines()
                        logging.info(fps1)
                        logging.info(fps2)
                    min_misfit = misfit
                    very_best_source = current
                
                ngood += 1
                
            if min_misfit is None:
                raise Exception('No valid starting points found')
        else:
            logging.error('Fixme: inversion.py')
        
        very_best_source_normalform = copy.deepcopy(very_best_source)
        if 'strike' in base_source.keys():
            very_best_source_normalform.disambigue_sdr()
        
        for param, coords in grid_def:
            val = very_best_source[param]
            remark = ''
            if snap(val,coords)[0] in (0,len(coords)-1): remark = ' (?: at end of range)'
            str_result = '%s = %g%s' % (param.title(), val, remark)
            logging.info(str_result)
            self.result(str_result, d2u(param))
            base_source[param] = very_best_source_normalform[param]
            self.out_config.__dict__[d2u(param)] = very_best_source_normalform[param]
            
        str_result = 'Misfit = %g' % min_misfit
        logging.info(str_result)
        self.result(str_result, 'Misfit')
        self.out_config.__dict__['min_misfit'] = min_misfit
        
        logging.info('Total number iterations in gradient searches: %i' % self.nminfunccalls)
        logging.info('Number of gradient searches tried: %i' % ntotal)
        logging.info('Number of successful gradient searches: %i' % ngood)
        self.out_config.__dict__['greeper_nminfunccalls'] = self.nminfunccalls
        self.out_config.__dict__['greeper_ntotal'] = ntotal
        self.out_config.__dict__['greeper_ngood'] = ngood
           
        if forward:
            self.snapshot( base_source, 'best', mm_conf )
            
        self.post_work(search or forward)
    


    