import receiver
import gfdb
import seismosizer
import config
import source
import util
import plotting
import moment_tensor
from util import gform

import shutil
import os
import sys
import re
import numpy as num
import scipy
import scipy.stats
import progressbar
import cPickle as pickle
import copy
import logging
from subprocess import call
from Cheetah.Template import Template

from os.path import join as pjoin

def convert_graph(in_filename, out_filename):
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
    return s.replace('-','_')

def u2d(s):
    return s.replace('_','-')

def grid_defi( param, oldval, descr ):
    mi, ma, inc = descr[:3]
    mode = 'absolute'
    if len(descr) == 4: mode = descr[3]
    if mode == 'mult':
        mi *= oldval
        ma *= oldval
        inc *= oldval
    elif mode == 'add':
        mi += oldval
        ma += oldval
    return param, mi, ma, inc

def standard_setup( datadir,
                    gfdb_path,
                    components,
                    effective_dt=1,
                    spacial_undersampling = [ 1, 1 ],
                    hosts = ['localhost'],
                    crustal_thickness_limit = None,
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
                               'crustal_thickness_limit', 'shifts', 'local_interpolation', 
                               'source_origin_file', 'receivers_file',
                               'ref_seismogram_stem', 'ref_seismogram_format', 'blacklist', 'verbose'))

def gen_dweights( seis, base_source, datadir,
                                 ref_seismogram_stem = 'reference',
                                 ref_seismogram_format = 'mseed', **kwargs):
    
    
    base_source = copy.copy(base_source)
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
            logging.warn('required parameter missing for step %s: %s' % (self.stepname, k))
            
        for k in have - (self.optional | self.required):
            logging.info('unused parameter in config for step %s: %s' % (self.stepname, k))
            
        
        logging.info('Starting work on step: %s' % self.stepname)
        rundir = self.make_rundir_path('incomplete')
        if os.path.exists( rundir ): 
            shutil.rmtree( rundir )
        os.makedirs( rundir )
        self.in_config.dump( pjoin(rundir,'config-in.pickle') )
        self.out_config = config.Config()
        
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
        logging.info('Done with work on step: %s' % self.stepname)

    def snapshot(self, source, ident):
        
        if not self.seismosizer:
            logging.warn('cannot create snapshot, because no seismosizers are running')
            return
        
        self.seismosizer.set_source( source )
        snapshot = self.seismosizer.get_receivers_snapshot()
        self.dump( snapshot, 'snapshot_%s' % ident )
        
    def get_snapshot(self, ident, run_id='current'):
        return self.load( ident='snapshot_%s' % ident, run_id=run_id )

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
        c = self.load('config-in', run_id=run_id)
        return c
    
    def oc(self, run_id='current'):
        c = self.load('config-out', run_id=run_id)
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
        
        print 'closest station: %s' % sx(seis.receivers[imin])
        print 'farest station:  %s' % sx(seis.receivers[imax])
        
        self.post_work(True)


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
        base_source = source.Source( sourcetype, 
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
        base_source = source.Source( sourcetype, 
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
            print i, sum(ms)
        
        
        datadir = conf['datadir']
        # reset reference seismograms
        ref_seismogram_stem     = pjoin(datadir, ref_seismogram_stem)
        seis.set_ref_seismograms( ref_seismogram_stem, ref_seismogram_format )
        
        self.post_work(True)




class Shifter(Step):
    def __init__(self, workdir, name='shifter'):
        Step.__init__(self, workdir, name)
        
        self.required |= set(('taper', 'filter', 'autoshift_range', 'autoshift_limit'))
                        
        self.optional |=  set([d2u(p) for p in source.param_names('eikonal')])
        
    def work(self, **kwargs):
        self.pre_work(True)
        seis = self.seismosizer
        conf = self.in_config.get_config()
        
        sourcetype = 'eikonal'
        base_source = source.Source( sourcetype )
        
        for p in source.param_names(sourcetype):
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
                        
        self.optional |= set([d2u(p) for p in source.param_names('eikonal')])
        
    def work(self, search=True, forward=True, run_id='current'):
        self.pre_work(search or forward)
        seis = self.seismosizer
        conf = self.in_config.get_config()
        mm_conf = self.in_config.get_config(keys=Step.outer_misfit_method_params)
        
        sourcetype = 'eikonal'
        base_source = source.Source( sourcetype )
        
        for p in source.param_names(sourcetype):
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
            finder = MisfitGrid( base_source, grid_def )
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
        
        if forward:
            self.snapshot( base_source, 'best' )
            
        self.post_work(search or forward)
        
    def _plot( self, run_id='current' ):
        
        conf = self.in_config.get_config()
        plotdir = self.make_plotdir_path(run_id)
        finder = self.load(self.stepname, run_id=run_id)
        plot_files = finder.plot( plotdir, conf['nsets'] )
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
        base_source = source.Source( sourcetype, { 'bord-radius': 0. } )
        
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
                finder = MisfitGrid( base_source, grid_def )
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
            sdrfinder = MisfitGrid( base_source, grid_def )
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
                        
        self.optional |= set([d2u(p) for p in source.param_names('eikonal')])
    
    def work(self, search=True, forward=True, run_id='current'):
        self.pre_work(search or forward)
        seis = self.seismosizer
        conf = self.in_config.get_config()
        mm_conf = self.in_config.get_config(keys=Step.outer_misfit_method_params)
        
        sourcetype = 'eikonal'
        
        base_source = source.Source( sourcetype, {} )
    
        for p in source.param_names(sourcetype):
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
            finder = MisfitGrid( base_source, grid_def )
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
        
        logging.info('reweighting')
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
        
        
class MisfitGridStats:
    
    def str_best_and_confidence(self, factor=1., unit =''):
        lw = ''
        uw = ''
        if self.percentile16_warn: lw = ' (?)'
        if self.percentile84_warn: uw = '(?) '
           
        return '%s = %g %s  (confidence interval 68%%) = [ %g%s, %g %s] %s' % \
                (self.paramname.title(), self.best*factor, unit, self.percentile16*factor, lw, self.percentile84*factor, uw, unit)

    def str_mean_and_stddev(self):
        return '%(paramname)s = %(mean)g +- %(std)g' % self.__dict__
    
class MisfitGrid:
    '''Brute force grid search minimizer with builtin bootstrapping.'''
    
    def __init__( self, base_source,
                        param_ranges,
                        source_constraints=None,
                        ref_source=None):
        
        self.base_source = copy.copy(base_source)
        if ref_source:
            self.ref_source = copy.copy(ref_source)
        else:
            self.ref_source = copy.copy(base_source)
            
        self.param_ranges = param_ranges
        self.sources = self.base_source.grid( self.param_ranges, 
                                              source_constraints=source_constraints)

        self.sourceparams = [ x[0] for x in self.param_ranges ]

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
        
        nsources = len(self.sources)
        nreceivers = len(seis.receivers)
        ncomponents = max([ len(r.components) for r in seis.receivers ])
        
        
        # results gathered by (source,receiver,component)
        misfits_by_src = num.zeros( (nsources, nreceivers, ncomponents), dtype=num.float)
        norms_by_src = num.zeros(  (nsources, nreceivers, ncomponents), dtype=num.float)
        ref_misfits_by_src = num.zeros( (1, nreceivers, ncomponents), dtype=num.float)
        ref_norms_by_src = num.zeros(  (1, nreceivers, ncomponents), dtype=num.float)
        
        if config.show_progress:
            if len(self.sourceparams) == 1:
                title = 'Grid search param: ' + ', '.join(self.sourceparams)
            else:
                title = 'Grid search params: ' + ', '.join(self.sourceparams)
            widgets = [title, ' ',
                    progressbar.Bar(marker='-',left='[',right=']'), ' ',
                    progressbar.Percentage(), ' ',]
            
            pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(self.sources)).start()
        
        failings = []
        for isource, source in enumerate(self.sources):
            try:
                seis.make_misfits_for_source( source )
          
                for ireceiver, receiver in enumerate(seis.receivers):
                    if receiver.enabled:
                        for icomp, comp in enumerate(receiver.components):
                            misfits_by_src[isource, ireceiver, icomp] = receiver.misfits[icomp]
                            norms_by_src[isource, ireceiver, icomp] = receiver.misfit_norm_factors[icomp]
                        
            except seismosizer.SeismosizersReturnedErrors:
                failings.append(isource)
            
            if config.show_progress: pbar.update(isource+1)
                        
        if config.show_progress: pbar.finish()
              
        # reference misfits
        seis.make_misfits_for_source( self.ref_source )
        for ireceiver, receiver in enumerate(seis.receivers):
             for icomp, comp in enumerate(receiver.components):
                 ref_misfits_by_src[0, ireceiver, icomp] = receiver.misfits[icomp]
                 ref_norms_by_src[0, ireceiver, icomp] = receiver.misfit_norm_factors[icomp]
                 
                 
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
        stepsizes = [ x[3] for x in self.param_ranges ]
        self.stats = self._stats(self.param_ranges, self.best_source, self.bootstrap_sources)
    
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
        
    def _best_source(self, return_misfits_by_r=False, **outer_misfit_config):
        misfits_by_s, misfits_by_sr = self._make_global_misfits( **outer_misfit_config)
        ibest = num.argmin(misfits_by_s)
        if not return_misfits_by_r:
            return self.sources[ibest], misfits_by_s
        else:
            # misfit variability by receiver
            misfits_varia_by_r = num.std(misfits_by_sr,0)
            return self.sources[ibest], misfits_by_s, misfits_by_sr[ibest,:], misfits_varia_by_r
        
    def _ref_misfits_by_r(self, **outer_misfit_config):
        misfits_by_s, misfits_by_sr = self._make_global_misfits(process_ref_source=True, **outer_misfit_config)
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
        
    def _stats(self, param_ranges, best_source, bootstrap_sources):
        
        results = {}
        for param, mi, ma, step in param_ranges:
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
            result.percentile16 = scipy.stats.scoreatpercentile(param_results, 16.) - step/2.
            result.percentile84 = scipy.stats.scoreatpercentile(param_results, 84.) + step/2.
            result.percentile16_warn = result.percentile16 < mi
            result.percentile84_warn = result.percentile84 > ma
            result.distribution = param_results
            results[param] = result
            
        return results
        
    def _make_global_misfits(self, process_ref_source=False, receiver_weights=1., outer_norm='l2norm', anarchy=False, bootstrap=False, **kwargs):
        
        if process_ref_source:
            misfits_by_src = self.ref_misfits_by_src
            norms_by_src = self.ref_norms_by_src
        else:
            misfits_by_src = self.misfits_by_src
            norms_by_src = self.norms_by_src
        
        if isinstance(receiver_weights, float):
            rweights = receiver_weights
        else:
            rweights = receiver_weights[num.newaxis,:].copy()
        
        if bootstrap:
            bweights = num.zeros(self.nreceivers, dtype=num.float)
            bweights_x = num.bincount(num.random.randint(0,self.nreceivers,self.nreceivers))
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
            maxm = num.amax(misfits_by_s)
            misfits_by_s = num.where(misfits_by_s<0, maxm, misfits_by_s)
            
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
            maxm = num.amax(misfits_by_s)
            misfits_by_s = num.where(misfits_by_s<0, maxm, misfits_by_s)
        
        else:
            raise Exception('unknown norm method: %s' % outer_norm)
        
        return misfits_by_s, misfits_by_sr
        
    def plot(self, dirname, nsets):
        
        best_source = self.best_source
        bootstrap_sources = self.bootstrap_sources
        
        
        plot_files = []
        
        for iparam, param in enumerate(self.sourceparams):
            
            #
            # 1D misfit cross section
            #
            xdata = num.array([ s[param] for s in self.sources ], dtype=num.float)
            conf = dict( xlabel = param.title(),
                         xunit = self.base_source.sourceinfo(param).unit )
            
            plotting.km_hack(conf)
            fn = 'misfit-%s.pdf' % param
            plotting.misfit_plot_1d( [(xdata, self.misfits_by_s)],
                                     pjoin(dirname, fn),
                                     conf )
            plot_files.append(fn)
            
           
            #
            # 1D histogram
            #
            # mi = num.min(r[param].distribution)
            # ma = num.max(r[param].distribution)
            mi, ma, step = self.param_ranges[iparam][1:]
            n = int(round((ma-mi)/step))+1
            hist, edges = num.histogram(self.stats[param].distribution,
                                        bins=n,
                                        range=(mi-step/2., ma+step/2.),
                                        new=True)
            
            hist = hist/float(len(bootstrap_sources))
            xdata = edges[:-1]+step/2.
            conf = dict( xlabel = param.title(),
                         xunit = self.base_source.sourceinfo(param).unit,
                         xrange = (mi-step/2., ma+step/2.), )
            
            plotting.km_hack(conf)
            fn = 'histogram-%s.pdf' % param
            plotting.histogram_plot_1d( [(xdata, hist)],
                                        pjoin(dirname, fn),
                                        conf )
            
            plot_files.append( fn )
            
        for ixparam, xparam in enumerate(self.sourceparams):
            for iyparam, yparam in enumerate(self.sourceparams):
                if ixparam == iyparam: continue
                
                #
                # 2D misfit plot
                #
                vmisfits = {}
                maxmisfit = num.amax(self.misfits_by_s)
                for (source, misfit) in zip(self.sources, self.misfits_by_s):
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
                fn = 'misfit-%s-%s.pdf' % (xparam, yparam)
                plotting.misfit_plot_2d( [(ax, ay, az)],
                                         pjoin(dirname, fn),
                                         conf )
                plot_files.append(fn)
                
                conf['xrange'] = (num.amin(ax), num.amax(ax))
                conf['yrange'] = (num.amin(ay), num.amax(ay))
                
                fn = 'histogram-%s-%s.pdf' % (xparam, yparam)
                plotting.histogram_plot_2d( [(axc, ayc, ac)],
                                         pjoin(dirname, fn),
                                         conf )
                plot_files.append(fn)
                
                conf['zrange'] = (num.amin(az), num.amax(az))
                fn = 'misfogram-%s-%s.pdf' % (xparam, yparam)
                plotting.misfogram_plot_2d( [(ax, ay, az), best_loc, (axc, ayc, ac)],
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
    
def main(steps):
    
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--loglevel', action='store', dest='loglevel', type='choice', choices=('info', 'debug'), default='info')
    parser.add_option('--no-progress', action='callback', callback=progress_off)
    parser.add_option('--no-search', action='store_false', dest='do_search', default=True)
    parser.add_option('--no-forward', action='store_false', dest='do_forward', default=True)
    parser.add_option('--run-id', action='store', dest='run_id', type='string', default='current')
    
    (options, args) = parser.parse_args()
    
    levels = { 'info' : logging.INFO, 'debug' : logging.DEBUG }
    logging.basicConfig( level = levels[options.loglevel], format='%(relativeCreated)s: %(message)s' )
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
        stepnames_to_do = args
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
            if command == 'replot':
                step.plot(run_id=options.run_id)
            if command == 'show-in-config':
                step.show_in_config(run_id=options.run_id)
            if command == 'show-out-config':
                step.show_out_config(run_id=options.run_id)
            if command == 'show-active-config':
                step.show_active_config(run_id=options.run_id)
    
    