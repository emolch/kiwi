import receiver
import gfdb
import seismosizer
import config
import source as source_model
import util
import plotting
import moment_tensor
import orthodrome
from util import gform

import shutil
import os
import sys
import re
import datetime, time
import numpy as num
import scipy
import scipy.stats
from scipy.optimize import fmin_l_bfgs_b
import progressbar
import cPickle as pickle
import copy
import logging
from subprocess import call
import subprocess
from Cheetah.Template import Template

from os.path import join as pjoin

def backticks(command):
    return subprocess.Popen(command, stdout=subprocess.PIPE).communicate()[0]

def convert_graph(in_filename, out_filename):
    if backticks(['which', 'pdf2png']).strip():
        inter_filename = out_filename.replace('.png','.oversized.png')
        call(['pdf2png', in_filename, inter_filename, '1'])
        call(['convert', inter_filename, '-resize', '66%', out_filename])
    else:
        logging.info("Using convert; install pdf2png from libcairo for better quality PNG output.")
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

def mimainc_to_gvals(mi,ma,inc):
    vmin, vmax, vinc = float(mi), float(ma), float(inc)
    n = int(round((vmax-vmin)/vinc))+1
    vinc = (vmax-vmin)/(n-1)
    return num.array([ vmin+i*vinc for i in xrange(n) ], dtype=num.float)
    
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
    else:    
        if mode == 'mult':
            mi *= oldval
            ma *= oldval
            inc *= oldval
        elif mode == 'add':
            mi += oldval
            ma += oldval
            
        return param, mimainc_to_gvals(mi, ma, inc)
    
    

def step_at(values, value):
    if len(values) <= 1: return 1.
    i = num.clip(num.searchsorted(values, value), 1, len(values)-1)
    return values[i]-values[i-1]

def values_to_bin_edges(values):
    if len(values) == 0: return []
    edges = num.zeros(len(values)+1, dtype=num.float)
    edges[1:-1] = (values[1:]+values[:-1])/2.
    edges[0] = values[0]-step_at(values, values[0])/2.
    edges[-1] = values[-1]+step_at(values, values[-1])/2.
    return edges

def standard_setup( datadir,
                    gfdb_path,
                    components,
                    effective_dt=1,
                    spacial_undersampling = [ 1, 1 ],
                    hosts = ['localhost'],
                    crustal_thickness_limit = None,
                    constraining_planes = None,
                    shifts = None,
                    blacklist = None,
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
    seis = seismosizer.Seismosizer(hosts)
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
    receivers = receiver.load_table(receivers_file, components=components)
    seis.set_receivers(receivers)
    seis.set_ref_seismograms( ref_seismogram_stem, ref_seismogram_format )
    if blacklist:
        seis.blacklist_receivers( blacklist )
    
    # apply reference seismograms shifts
    if shifts is not None:
        seis.shift_ref_seismograms( shifts )
    
    return seis
    
standard_setup.required = set(('datadir', 'gfdb_path', 'components')) 
standard_setup.optional = set(('effective_dt', 'spacial_undersampling', 'hosts', 
                               'crustal_thickness_limit', 'constraining_planes', 'shifts', 'local_interpolation', 
                               'source_origin_file', 'receivers_file',
                               'ref_seismogram_stem', 'ref_seismogram_format', 'blacklist', 'verbose'))

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
    sdr_grid = MisfitGrid( base_source, [strike_grid, dip_grid, slip_rake_grid])
    
    sdr_grid.compute( seis )
    
    means = sdr_grid.get_mean_misfits_by_r()
    means /= num.mean(means)
    dweights = num.where(means>0., 1./means, 0.)
    # reset reference seismograms
    ref_seismogram_stem     = pjoin(datadir, ref_seismogram_stem)
    seis.set_ref_seismograms( ref_seismogram_stem, ref_seismogram_format )
    
    return dweights
    
gen_dweights.required = set(('datadir',))
gen_dweights.optional = set(('ref_seismogram_stem', 'ref_seismogram_format'))

def other_plane( strike, dip, rake ):
    mt = moment_tensor.MomentTensor( strike=strike, dip=dip, rake=rake )
    both_sdr = mt.both_strike_dip_rake()
    w = [ sum( [ abs(x-y) for x,y in zip(both_sdr[i], (strike, dip, rake)) ] ) for i in (0,1) ]
    if w[0]<w[1]:
        return both_sdr[1]
    else:
        return both_sdr[0]


def make_misfits_for_sources(seis, sources, show_progress=False, progress_title='grid search'):
    nsources = len(sources)
    nreceivers = len(seis.receivers)
    ncomponents = max([ len(r.components) for r in seis.receivers ])
    
    # results gathered by (source,receiver,component)
    misfits_by_src = num.zeros( (nsources, nreceivers, ncomponents), dtype=num.float)
    norms_by_src = num.zeros(  (nsources, nreceivers, ncomponents), dtype=num.float)
    
    if show_progress:
        widgets = [progress_title, ' ',
                progressbar.Bar(marker='-',left='[',right=']'), ' ',
                progressbar.Percentage(), ' ',]
        
        pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(sources)).start()
    
    failings = []
    for isource, source in enumerate(sources):
        try:
            seis.make_misfits_for_source( source )
            for ireceiver, receiver in enumerate(seis.receivers):
                if receiver.enabled:
                    for icomp, comp in enumerate(receiver.components):
                        misfits_by_src[isource, ireceiver, icomp] = receiver.misfits[icomp]
                        norms_by_src[isource, ireceiver, icomp] = receiver.misfit_norm_factors[icomp]
                    
        except seismosizer.SeismosizersReturnedErrors:
            failings.append(isource)
        
        if show_progress: pbar.update(isource+1)
                    
    if show_progress: pbar.finish()
    
    return misfits_by_src, norms_by_src, failings

def make_global_misfits(misfits_by_src, norms_by_src, receiver_weights=1., outer_norm='l2norm', anarchy=False, bootstrap=False, **kwargs):
    
    nreceivers = misfits_by_src.shape[1]
    
    if isinstance(receiver_weights, float):
        rweights = receiver_weights
    else:
        rweights = receiver_weights[num.newaxis,:].copy()
    
    if bootstrap:
        bweights = num.zeros(nreceivers, dtype=num.float)
        bweights_x = num.bincount(num.random.randint(0,nreceivers,nreceivers))
        bweights[:len(bweights_x)] = bweights_x[:]
    
    if outer_norm == 'l1norm':
        misfits_by_sr = num.sum(misfits_by_src,2)
        norms_by_sr   = num.sum(norms_by_src,2)
        
        if anarchy:
            xrweights = num.zeros(norms_by_sr.shape, dtype=float)
            xrweights[:,:] = rweights
            xrweights /= num.where( norms_by_sr != 0., norms_by_sr, -1.)
            rweights = num.maximum(xrweights, 0.)
        
        if bootstrap:
            rweights *= bweights
        
        misfits_by_sr *= rweights
        norms_by_sr *= rweights
        
        ms = num.sum( misfits_by_sr, 1 )
        ns = num.sum( norms_by_sr, 1 )
        
        misfits_by_s  = num.where(ns > 0., ms/ns, -1.)
        misfits_by_s = num.where(misfits_by_s<0, num.NaN, misfits_by_s)
        
    elif outer_norm == 'l2norm':
        misfits_by_sr = num.sqrt(num.sum(misfits_by_src**2,2))
        norms_by_sr   = num.sqrt(num.sum(norms_by_src**2,2))
        
        if anarchy:
            rweights /= num.where( norms_by_sr != 0., norms_by_sr, -1.)
            rweights = num.maximum(rweights, 0.)
            
        if bootstrap:
            rweights *= num.sqrt(bweights)
        
        misfits_by_sr *= rweights
        norms_by_sr *= rweights
        
        ms = num.sum( (misfits_by_sr)**2, 1 )
        ns = num.sum( (norms_by_sr)**2, 1 )
        
        misfits_by_s  = num.where(ns > 0., num.sqrt(ms/ns), -1.)
        misfits_by_s = num.where(misfits_by_s<0, num.NaN, misfits_by_s)
    
    else:
        raise Exception('unknown norm method: %s' % outer_norm)
    
    return misfits_by_s, misfits_by_sr

class NoValidSources(Exception):
    pass

def best_source(seis, sources, return_failings=False, **outer_misfit_config):
    misfits_by_src, norms_by_src, failings = make_misfits_for_sources(seis, sources)
    misfits_by_s, misfits_by_sr = make_global_misfits( misfits_by_src, norms_by_src, **outer_misfit_config)
    misfits_by_s = num.where( misfits_by_s > 0, misfits_by_s, num.NaN)
    ibest = num.nanargmin(misfits_by_s)
    if num.isnan(ibest) or num.isnan(misfits_by_s[ibest]): raise NoValidSources()
    if return_failings:
        return sources[ibest], misfits_by_s[ibest], failings
    return sources[ibest], misfits_by_s[ibest]

class Step:
    inner_misfit_method_params = set(('inner_norm', 'taper', 'filter', 'nsets'))
    outer_misfit_method_params = set(('outer_norm', 'bootstrap_iterations', 'anarchy', 'receiver_weights'))
    
    def __init__(self, workdir, name):
        self.baseworkdir = workdir
        self.stepname = name
        self.in_config = None
        self.out_config = None
        self.stepdir = pjoin(self.baseworkdir, self.stepname)
        self.seismosizer = None
        self.required = set(standard_setup.required)
        self.optional = set(standard_setup.optional)
        
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
        
        if start_seismosizer:
            sconf = self.in_config.get_config(keys=standard_setup.required|standard_setup.optional)
            self.seismosizer = standard_setup( **sconf )

    def setup_inner_misfit_method(self):
        conf = self.in_config.get_config(keys=Step.inner_misfit_method_params)
        seis = self.seismosizer
        
        tapers_by_set = conf['taper']
        assert(len(tapers_by_set) == conf['nsets'])
        tapers = [ tapers_by_set[i%len(tapers_by_set)] for i in range(len(seis.receivers)) ] 
        seis.set_taper(tapers)
        if conf['filter']:
            seis.set_filter(conf['filter'])
        seis.set_misfit_method(conf['inner_norm'])
        
    def post_work(self, stop_seismosizer=True):
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
        

    def snapshot(self, source, ident):
        
        if not self.seismosizer:
            logging.warn('Cannot create snapshot, because no seismosizers are running')
            return
        
        self.seismosizer.set_source( source )
        self.dump( self.seismosizer.get_receivers_snapshot(), 'snapshot_%s' % ident )
        self.dump( self.seismosizer.get_psm_infos(), 'source_infos_%s' % ident )
        
    def get_snapshot(self, ident, run_id='current'):
        return self.load( ident='snapshot_%s' % ident, run_id=run_id )

    def get_snapshot_source_infos(self, ident, run_id='current'):
        return self.load( ident='source_infos_%s' % ident, run_id=run_id)

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
            
    def gx(self, arg):
        plotdir = self.make_plotdir_path( 'current' )
        location = pjoin(plotdir, arg)
        img_snippet = '<a href="%s"><img style="border:none;" src="%s" /></a>'
        return img_snippet % (location+'.pdf', location+'.png')
    
    def gxi(self, arg):
        plotdir = self.make_plotdir_path( 'current' )
        location = pjoin(plotdir, arg)
        return location+'.png'

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
        rnames = [ re.sub(r'\..*$', '', r.name) for r in receivers ]
        slat, slon = source_location[:2]
        
        station_size = [ 0.05 ] * len(receivers)
        station_color = [ 1. ] * len(receivers)
        
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
                               {"depth": float(conf['depth']),
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
        seis.set_taper(tapers)
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
        rnames = [ re.sub(r'\..*$', '', r.name) for r in receivers ]
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
        

class ParamTuner(Step):
     
    def __init__(self, workdir, params=['time'], name=None):
        if name is None: name = '-'.join(params)+'-tuner'
        Step.__init__(self, workdir, name)
        self.params = params
        
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
        base_source = source_model.Source( sourcetype )
        
        for p in source_model.param_names(sourcetype):
            if d2u(p) in conf:
                base_source[p] = float(conf[d2u(p)])
                
        if 'plane' in conf and conf['plane'] == 2: 
            strike, dip, slip_rake = float(conf['strike']), float(conf['dip']), float(conf['slip_rake'])
            strike, dip, slip_rake = other_plane( strike, dip, slip_rake )
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
            self.setup_inner_misfit_method()
            finder = MisfitGrid( base_source, param_values=grid_def )
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
        
        logging.info('Misfit = %f, Source = %s' % (finder.get_best_misfit(), str(base_source)))

        if forward:
            self.snapshot( base_source, 'best' )
            
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
    
    
        
class PlaneTuner(Step):
     
    def __init__(self, workdir, name='planetuner'):
        Step.__init__(self, workdir, name)
        
        self.required |= Step.outer_misfit_method_params | Step.inner_misfit_method_params \
                        | set(('rise_time',)) \
                        | set(('time', 'depth', 'moment', 'strike', 'dip', 'slip_rake')) \
                        | set([param+'_range' for param in 'strike', 'dip', 'slip_rake', 'moment', 'rough_moment' ])
        self.optional |= set(('north_shift', 'east_shift'))
         
    def work(self, search=True, forward=True, run_id='current'):
        self.pre_work(search or forward)

        seis = self.seismosizer
        conf = self.in_config.get_config()
        mm_conf = self.in_config.get_config(keys=Step.outer_misfit_method_params)
        
        sourcetype = 'eikonal'
        base_source = source_model.Source( sourcetype, { 'bord-radius': 0. } )
        
        for param in 'time', 'north-shift', 'east-shift', 'depth', 'moment', 'strike', 'dip', 'slip-rake', 'rise-time':
            if d2u(param) in conf:
                base_source[param] = float(conf[d2u(param)])
        
        #
        # Rough Moment:
        #
        for param in 'moment',:
            oldval = base_source[param]
            descr = conf['rough_'+d2u(param)+'_range']
            grid_def = [ grid_defi(param,oldval,descr) ]
            
            if search or forward: self.setup_inner_misfit_method()
            if search:
                self.setup_inner_misfit_method()
                finder = MisfitGrid( base_source, param_values=grid_def )
                finder.compute(seis)
            else:
                finder = self.load('rough_'+param, run_id=run_id)
                
            finder.postprocess(**mm_conf)
            self.dump(finder, 'rough_'+param)
            
            stats = finder.stats[param]
            str_result = stats.str_mean_and_stddev()
            logging.info(str_result)
            base_source[param] = stats.mean
        
        
        #
        # Strike, Dip, Rake, Moment
        #
        sdrparams = ('strike', 'dip', 'slip-rake', 'moment')
        
        grid_def = []
        for param in sdrparams:
            oldval = base_source[param]
            descr = conf[d2u(param)+'_range']            
            grid_def.append(grid_defi(param,oldval,descr))
            
        if search or forward: self.setup_inner_misfit_method()
        if search:
            self.setup_inner_misfit_method()
            sdrfinder = MisfitGrid( base_source, param_values=grid_def )
            sdrfinder.compute(seis)
        else:
            sdrfinder = self.load(self.stepname, run_id=run_id)
                    
        sdrfinder.postprocess(**mm_conf)
        self.dump(sdrfinder, self.stepname)
        
        stats = sdrfinder.stats
        for param in sdrparams:
            str_result = stats[param].str_best_and_confidence()
            logging.info(str_result)
            self.result(str_result, param )
            base_source[param] = stats[param].best
            self.out_config.__dict__[d2u(param)] = stats[param].best
            self.out_config.__dict__[d2u(param)+'_stats'] = stats[param]
        
        if forward:
            self.snapshot( base_source, 'best' )
            
        self.post_work(search or forward)
        
    def _plot( self, run_id='current' ):
        #rundir = self.make_rundir_path(run_id)
        #finder = self.load('rough_moment', run_id=run_id)
        #finder.plot( pjoin(rundir, 'plots_rough_moment') )
        conf = self.in_config.get_config()

        plotdir = self.make_plotdir_path(run_id)
        finder = self.load(self.stepname, run_id=run_id)
        return finder.plot( plotdir, conf['nsets'] )
    
    
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
            finder = MisfitGrid( base_source, param_values=grid_def )
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
            self.snapshot( base_source, 'best' )
            
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
        for step, ident in self.snapshots:
            loaded_snapshots.append(self.get_snapshot("%s_%s" % (step.stepname, ident)))
        
        allfilez = plotting.multi_seismogram_plot( loaded_snapshots, plotdir )
        
        return [ os.path.basename(fn) for fn in allfilez ]

def snap(val,coords):
    ipos = num.argmin(num.abs(val-coords))
    return ipos, coords[ipos]

class Greeper(Step):
    '''Grid search over gradient searches
    
    aaahrrggg! this has to be rewritten!'''
     
    def __init__(self, workdir, params=['time'], name=None):
        if name is None: name = '-'.join(params)+'-greeper'
        Step.__init__(self, workdir, name)
        self.params = params
        
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
        base_source = source_model.Source( sourcetype )
        
        for p in source_model.param_names(sourcetype):
            if d2u(p) in conf:
                base_source[p] = float(conf[d2u(p)])
                
        if 'plane' in conf and conf['plane'] == 2: 
            strike, dip, slip_rake = float(conf['strike']), float(conf['dip']), float(conf['slip_rake'])
            strike, dip, slip_rake = other_plane( strike, dip, slip_rake )
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
                mi,ma,inc = grid_defi(u2d(param),oldval,descr)
                starter_grid_def.append( (u2d(param), gdef_to_gvals(mi,ma,inc)) )
        
        if starter_grid_def:
            starter_sources = base_source.grid( param_values=starter_grid_def )
        else:
            starter_sources = [ base_source ]
        
        grid_def = []
        for param in self.params:
            oldval = base_source[u2d(param)]
            descr = conf[param+'_range']
            grid_def.append(grid_defi(u2d(param),oldval,descr))
        
        dparams = [ u2d(param) for param in self.params ]
        self.nminfunccalls = 0

        # for each dimension in parameter space create a ladder with the grid coords
        grid_coords = []
        for param, mi, ma, inc in grid_def:
            ncoords = int(round((ma-mi)/inc))+1
            grid_coords.append((param, num.linspace(mi, ma, ncoords)))

        norms = [ inc for (param, mi, ma, inc) in grid_def ]
        bounds = [ (mi/norm,ma/norm) for  ((param, mi, ma, inc),norm) in zip(grid_def, norms) ]

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
                best, misfit = best_source(seis, [source], **mm_conf)
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
            for iparam, (param, mi, ma, inc) in enumerate(grid_def):
                if param == 'depth':
                    depth_sources = base_source.grid( [ (param,mi,ma,inc) ] )
                    dummy_source, dummy_misfit, failings = best_source(seis, depth_sources, return_failings=True, **mm_conf)
                    ok = []
                    for isource, source in enumerate(depth_sources):
                        if isource not in failings:
                            ok.append(source['depth'])
                    bounds[iparam] = ((min(ok)+inc*0.3)/norms[iparam], (max(ok)-inc*0.3)/norms[iparam])
            
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
                    current_base, misfit = best_source(seis, [current_base], **mm_conf)
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
        very_best_source_normalform.disambigue_sdr()
        
        for param, coords in grid_coords:
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
            self.snapshot( base_source, 'best' )
            
        self.post_work(search or forward)
    
class MisfitGridStats:
    
    def str_best_and_confidence(self, factor=1., unit =''):
        lw = ''
        uw = ''
        if self.percentile16_warn: lw = ' (?)'
        if self.percentile84_warn: uw = '(?) '
           
        return '%s = %.3g %s  (confidence interval 68%%) = [ %.3g%s, %.3g %s] %s' % \
                (self.paramname.title(), self.best*factor, unit, self.percentile16*factor, lw, self.percentile84*factor, uw, unit)
        
    def str_best(self, factor=1., unit =''):
        
        return '%s = %.3g %s' % (self.paramname.title(), self.best*factor, unit)
        
    def str_mean_and_stddev(self):
        return '%(paramname)s = %(mean)g +- %(std)g' % self.__dict__
    

    
class MisfitGrid:
    '''Brute force grid search minimizer with builtin bootstrapping.'''
    
    def __init__( self, base_source,
                        param_ranges=None,
                        param_values=None,
                        source_constraints=None,
                        ref_source=None):
        
        self.base_source = copy.deepcopy(base_source)
        if ref_source:
            self.ref_source = copy.deepcopy(ref_source)
        else:
            self.ref_source = copy.deepcopy(base_source)
        
        if param_values is not None:
            self.param_values = param_values
        else:
            self.param_values = []
            for param, mi, ma, inc in param_ranges:
                self.param_values.append( (param, gdef_to_gvals(mi,ma,inc)) )
            
        self.sources = self.base_source.grid( self.param_values, 
                                              source_constraints=source_constraints)

        self.sourceparams = [ x[0] for x in self.param_values ]

        # will be set by compute()
        self.misfits_by_src = None
        self.norms_by_src = None
        self.ref_misfits_by_src = None
        self.ref_norms_by_src = None        
        self.receivers = None
        self.nreceivers = None
        
        # will be set by postprocess()
        self.best_source = None
        self.misfits_by_s = None
        self.misfits_by_r = None
        self.ref_misfits_by_r = None
        self.variability_by_r = None
        self.bootstrap_sources = None
        self.stats = None
        
    def compute(self, seis):
        '''Let seismosizer calculate the trace misfits.'''
        
        if len(self.sourceparams) == 1:
            progress_title = 'Grid search param: ' + ', '.join(self.sourceparams)
        else:
            progress_title = 'Grid search params: ' + ', '.join(self.sourceparams)
            
        nreceivers = len(seis.receivers)
        
        # results, gathered by (source,receiver,component)
        misfits_by_src, norms_by_src, failings = make_misfits_for_sources(seis,
                                                                          self.sources,
                                                                          show_progress=config.show_progress,
                                                                          progress_title=progress_title)
              
        # results for reference source, gathered by (source,receiver,component)
        ref_misfits_by_src, ref_norms_by_src, failings = make_misfits_for_sources(seis, [self.ref_source])
        
        self.misfits_by_src = misfits_by_src
        self.norms_by_src = norms_by_src
        self.ref_misfits_by_src = ref_misfits_by_src
        self.ref_norms_by_src = ref_norms_by_src
        
        self.nreceivers = nreceivers
        self.receivers = seis.receivers
        self.source_location = seis.source_location
        
        self.best_source = None
        self.misfits_by_s = None
        self.misfits_by_r = None
        self.ref_misfits_by_r = None
        self.variability_by_r = None
        self.bootstrap_sources = None
        self.stats = None
    
    def postprocess(self, **outer_misfit_config):
        '''Combine trace misfits to global misfits, find best source, make statistics.'''
        
        self.best_source, self.misfits_by_s, self.misfits_by_r, self.variability_by_r = self._best_source(return_misfits_by_r=True, **outer_misfit_config)
        
        self.ref_misfits_by_r = self._ref_misfits_by_r(**outer_misfit_config)
        
        self.bootstrap_sources = self._bootstrap(**outer_misfit_config)
        self.stats = self._stats(self.param_values, self.best_source, self.bootstrap_sources)
    
    def get_mean_misfits_by_r(self):
        '''Get mean raw misfits by receiver, e.g. to auto-create weights.''' 
        mean_misfits_by_r = num.zeros( self.nreceivers, dtype=num.float )
        for irec in range(self.nreceivers):
            ncomps = len(self.receivers[irec].components)
            if ncomps != 0:
                x = num.sum(self.misfits_by_src[:,irec,:]) /( ncomps*len(self.sources))
            else:
                x = -1.0
            mean_misfits_by_r[irec] = x
        return mean_misfits_by_r
        
    def get_best_misfit(self):
        return num.nanmin(self.misfits_by_s)
        
    def _best_source(self, return_misfits_by_r=False, **outer_misfit_config):
        misfits_by_s, misfits_by_sr = make_global_misfits( self.misfits_by_src, self.norms_by_src, **outer_misfit_config)
        ibest = num.nanargmin(misfits_by_s)
        if not return_misfits_by_r:
            return self.sources[ibest], misfits_by_s
        else:
            # misfit variability by receiver
            misfits_varia_by_r = num.std(misfits_by_sr,0)
            return self.sources[ibest], misfits_by_s, misfits_by_sr[ibest,:], misfits_varia_by_r
        
    def _ref_misfits_by_r(self, **outer_misfit_config):
        misfits_by_s, misfits_by_sr = make_global_misfits(self.ref_misfits_by_src, self.ref_norms_by_src, **outer_misfit_config)
        return misfits_by_sr[0,:]
        
    def _bootstrap(self, bootstrap_iterations=1000, **outer_misfit_config):
        bootstrap_sources = []
        if config.show_progress:
            widgets = ['Bootstrapping', ' ',
                    progressbar.Bar(marker='-',left='[',right=']'), ' ',
                    progressbar.Percentage(), ' ',]
            
            pbar = progressbar.ProgressBar(widgets=widgets, maxval=bootstrap_iterations).start()
            
        for i in xrange(bootstrap_iterations):
            bootstrap_sources.append( self._best_source(bootstrap=True, **outer_misfit_config)[0] )
            if config.show_progress: pbar.update(i+1)
            
        if config.show_progress: pbar.finish()
        
        return bootstrap_sources
        
    def _stats(self, param_values, best_source, bootstrap_sources):
        
        results = {}
        for param, gvalues in param_values:
            mi = num.min(gvalues)
            ma = num.max(gvalues)
            param_results = num.zeros(len(bootstrap_sources), dtype=num.float)
            for i, source in enumerate(bootstrap_sources):
                param_results[i] = source[param]
                
            result = MisfitGridStats()
            result.best = best_source[param]
            if len(param_results) == 0: continue
            result.paramname = param
            result.mean = num.mean(param_results)
            result.std = num.std(param_results)
            result.median = num.median(param_results)
            result.percentile16 = scipy.stats.scoreatpercentile(param_results, 16.)
            result.percentile84 = scipy.stats.scoreatpercentile(param_results, 84.)
            result.percentile16 -= step_at( gvalues, result.percentile16 )/2.
            result.percentile84 += step_at( gvalues, result.percentile84 )/2.
            result.percentile16_warn = result.percentile16 < mi
            result.percentile84_warn = result.percentile84 > ma
            result.distribution = param_results
            results[param] = result
            
        return results
        

        
    def plot(self, dirname, nsets, source_model_infos=None):
        
        best_source = self.best_source
        bootstrap_sources = self.bootstrap_sources
        
        
        plot_files = []
        
        for iparam, param in enumerate(self.sourceparams):
            #
            # 1D misfit cross section
            #
            xy = []
            for s,m in zip(self.sources, self.misfits_by_s):
                if not num.isnan(m):
                    xy.append((s[param], m))
            
            xdata, ydata = num.array(xy,dtype=num.float).transpose()
            
            conf = dict( xlabel = param.title(),
                         xunit = self.base_source.sourceinfo(param).unit )
            
            plotting.km_hack(conf)
            fn = 'misfit-%s.pdf' % param
            plotting.misfit_plot_1d( [(xdata,ydata)],
                                     pjoin(dirname, fn),
                                     conf )
            plot_files.append(fn)
           
            #
            # 1D histogram
            #
            gvalues = self.param_values[iparam][1]
            gedges = values_to_bin_edges(gvalues)
            hist, edges = num.histogram(self.stats[param].distribution,
                                        bins=gedges,
                                        new=True)
            
            hist = hist/float(len(bootstrap_sources))
            
            conf = dict( xlabel = param.title(),
                         xunit = self.base_source.sourceinfo(param).unit)
            
            plotting.km_hack(conf)
            fn = 'histogram-%s.pdf' % param
            plotting.histogram_plot_1d( edges, hist, 
                                        pjoin(dirname, fn),
                                        conf )
            
            plot_files.append( fn )
            
        for ixparam, xparam in enumerate(self.sourceparams):
            for iyparam, yparam in enumerate(self.sourceparams):
                if ixparam == iyparam: continue
                
                #
                # 2D misfogram plot
                #
                vmisfits = {}
                maxmisfit = num.nanmax(self.misfits_by_s)
                for (source, misfit) in zip(self.sources, self.misfits_by_s):
                    if num.isnan(misfit): continue
                    vx = source[xparam]
                    vy = source[yparam]
                    vmisfits[(vx,vy)] = min(vmisfits.get((vx,vy), maxmisfit), misfit)
                    
                vcounts = {}
                for source in bootstrap_sources:
                    vx = source[xparam]
                    vy = source[yparam]
                    vcounts[(vx,vy)] = vcounts.get((vx,vy), 0) + 1
                
                lx, ly, lz = [], [], []
                for ((vx,vy),vmisfit) in vmisfits.iteritems():
                    lx.append(vx)
                    ly.append(vy)
                    lz.append(vmisfit)
                    
                lxc, lyc, lc =  [], [], []
                for ((vx,vy),vcount) in vcounts.iteritems():
                    lxc.append(vx)
                    lyc.append(vy)
                    lc.append(vcount)
                
                ax = num.array(lx,dtype=num.float)
                ay = num.array(ly,dtype=num.float)
                az = num.array(lz,dtype=num.float)
                
                
                axc = num.array(lxc,dtype=num.float)
                ayc = num.array(lyc,dtype=num.float)
                ac = num.array(lc,dtype=num.float)
                ac = ac/len(bootstrap_sources)
                
                conf = dict( xlabel = xparam.title(),
                             xunit = self.base_source.sourceinfo(xparam).unit,
                             ylabel = yparam.title(),
                             yunit = self.base_source.sourceinfo(yparam).unit,
                         )
                
                best_loc = (num.array([ best_source[xparam]]), num.array([best_source[yparam]]))
                plotting.km_hack(conf)
                plotting.nukl_hack(conf)
                
                fn = 'misfogram-%s-%s.pdf' % (xparam, yparam)
                plotting.misfogram_plot_2d_gmtpy( [(ax, ay, az), best_loc, (axc, ayc, ac)],
                                        pjoin(dirname, fn),
                                        conf )
                plot_files.append(fn)
        
        #
        # Station plot
        #
        lats = num.array( [ r.lat for r in self.receivers ], dtype='float' )
        lons = num.array( [ r.lon for r in self.receivers ], dtype='float' )
        dists = num.array( [ r.distance_deg for r in self.receivers ], dtype='float' )
        rnames = [ re.sub(r'\..*$', '', r.name) for r in self.receivers ]
        slat, slon = self.source_location[:2]
        station_misfits = self.misfits_by_r-self.ref_misfits_by_r
        #station_varia = self.variability_by_r / num.sum(self.variability_by_r) * len(self.receivers)
        station_varia = self.misfits_by_r / num.sum(self.misfits_by_r) * len(self.receivers)
        plotting.station_plot( slat, slon, lats, lons, rnames, station_misfits, station_varia, 
                              best_source, num.amax(dists)*1.05, pjoin(dirname, 'stations.pdf'), {}, nsets=nsets)
        plot_files.append('stations.pdf')
        
        if source_model_infos:
            delta_lat = 1.5*max(best_source['bord-radius']*2.,50000)/(20000.*1000.)*180
            plotting.location_map(pjoin(dirname, 'location.pdf'), slat, slon, delta_lat, {}, source=best_source,
                                  source_model_infos=source_model_infos, receivers=(lats,lons,rnames))
            plot_files.append('location.pdf')
        
        return plot_files
    

def progress_off(option, opt_str, value, parser):
    config.show_progress = False

def install(src, dst):
    dirs = os.path.split(dst)
    d,x = os.path.split(dst)
    dirs = []
    while d and not os.path.isdir(d):
        dirs.append(d)
        d,x = os.path.split(d)
        
    dirs.reverse()
    
    for d in dirs:
        if not os.path.exists(d):
            os.mkdir(d)
    
    shutil.copy(src, dst)
    
def kiwi_main(steps):
    
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--loglevel', action='store', dest='loglevel', type='choice', choices=('info', 'debug'), default='info')
    parser.add_option('--no-progress', action='callback', callback=progress_off)
    parser.add_option('--no-search', action='store_false', dest='do_search', default=True)
    parser.add_option('--no-forward', action='store_false', dest='do_forward', default=True)
    parser.add_option('--run-id', action='store', dest='run_id', type='string', default='current')
    
    (options, args) = parser.parse_args()
    
    logformat =  '[%(asctime)s] %(levelname)-8s %(message)s'
    dateformat = '%Y-%m-%d %H:%M:%S'
    levels = { 'info' : logging.INFO, 'debug' : logging.DEBUG }
    logging.basicConfig( level   = levels[options.loglevel],
                         format  = logformat,
                         datefmt = dateformat,
                         filename='kiwi.log' )
    console = logging.StreamHandler()
    formatter = logging.Formatter(logformat, datefmt=dateformat)
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    if options.loglevel == 'debug': progress_off
   
    commands = 'work', 'replot', 'show-steps', 'show-in-config', 'show-out-config', 'show-active-config', 'report'
    try:
        command = args.pop(0)
        
        assert(command in  commands)
    except:
        parser.error('available commands are: '+' '.join(commands))
    
    if command == 'show-steps':
        for step in steps:
            print step.stepname
        sys.exit()
    
    stepnames_all = [ step.stepname for step in steps ]
    
    if command == 'report':
        data = {}
        for step in steps: 
            data[step.stepname] = step
        
        t = Template(file='report.html', searchList=[ data ])
        
        page = str(t)
        files = [ x[1] or x[3] for x in re.findall(r'("([^"]+\.(png|pdf))"|\'([^\']+\.(png|pdf))\')', page) ]
        report_dir = 'report'
        
        for file in files:
            install(file, pjoin(report_dir, file))
            
        f = open(pjoin(report_dir,'index.html'),'w')
        f.write( page )
        f.close()
            
        sys.exit()
    
    if args:
        for arg in args:
            if arg != '-' and arg not in stepnames_all:
                parser.error('unknown stepname: %s\n' % arg + 'available stepnames are: '+' '.join(stepnames_all))
                
        stepnumber = dict([ (stepname,i) for i,stepname in enumerate(stepnames_all) ])
        xargs = []
        iarg = 0
        if args[0] == '-': args.insert(0,stepnames_all[0])
        if args[-1] == '-': args.append(stepnames_all[-1])
        while iarg < len(args)-2:
            if args[iarg+1] == '-':
                xargs.extend(stepnames_all[stepnumber[args[iarg]]:stepnumber[args[iarg+2]]+1])
                iarg += 3
            else:
                xargs.append(args[iarg])
                iarg += 1
        while iarg < len(args):
            xargs.append(args[iarg])
            iarg += 1
        stepnames_to_do = xargs
        
    else:
        stepnames_to_do = list(stepnames_all)
    
    for stepname in stepnames_to_do:
        if stepname not in stepnames_all:
            parser.error('unknown stepname: %s\n' % stepname + 'available stepnames are: '+' '.join(stepnames_all))
    
    for step in steps:
        if step.stepname in stepnames_to_do:
            if command == 'work':
                step.work(search=options.do_search, forward=options.do_forward, run_id=options.run_id)
                step.plot(run_id=options.run_id)
                if step.stepname != stepnames_to_do[-1]: logging.info('---')
            if command == 'replot':
                step.plot(run_id=options.run_id)
                if step.stepname != stepnames_to_do[-1]: logging.info('---')
            if command == 'show-in-config':
                step.show_in_config(run_id=options.run_id)
            if command == 'show-out-config':
                step.show_out_config(run_id=options.run_id)
            if command == 'show-active-config':
                step.show_active_config(run_id=options.run_id)
    
    