import receiver
import gfdb
import seismosizer
import config
import source
import util
import plotting
import moment_tensor

import shutil
import os
import sys
import re
import numpy as num
import scipy
import scipy.stats
import progressbar
import pickle
import copy
import logging

from os.path import join as pjoin    

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
                    effective_dt=1,
                    spacial_undersampling = [ 1, 1 ],
                    hosts = ['localhost'],
                    crustal_thickness_limit = None,
                    
                    local_interpolation = 'bilinear',
                    source_origin_file = 'source-origin.table',
                    receivers_file = 'receivers.table',
                    ref_seismogram_stem = 'reference',
                    ref_seismogram_format = 'sac', **kwargs):
                    
    '''Start seismosizers and setup source location, Greens functions, 
       receivers, and reference seismograms.'''

    source_origin_file      = pjoin(datadir, source_origin_file)
    ref_seismogram_stem     = pjoin(datadir, ref_seismogram_stem)
    receivers_file          = pjoin(datadir, receivers_file)

    # setup database
    database = gfdb.Gfdb(gfdb_path)
    seis = seismosizer.Seismosizer(hosts)
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
    receivers = receiver.load_table(receivers_file)
    seis.set_receivers(receivers)
    seis.set_ref_seismograms( ref_seismogram_stem, ref_seismogram_format )
    return seis
    
standard_setup.required = set(('datadir', 'gfdb_path')) 
standard_setup.optional = set(('effective_dt', 'spacial_undersampling', 'hosts', 
                               'crustal_thickness_limit', 'local_interpolation', 
                               'source_origin_file', 'receivers_file',
                               'ref_seismogram_stem', 'ref_seismogram_format'))

def gen_dweights( seis, base_source, datadir,
                                 ref_seismogram_stem = 'reference',
                                 ref_seismogram_format = 'sac', **kwargs):
    
    
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
    
    dweights = 1./means
    
    # reset reference seismograms
    ref_seismogram_stem     = pjoin(datadir, ref_seismogram_stem)
    seis.set_ref_seismograms( ref_seismogram_stem, ref_seismogram_format )
    
    return dweights
    
gen_dweights.required = set(('datadir',))
gen_dweights.optional = set(('ref_seismogram_stem', 'ref_seismogram_format'))

class Step:
    
    
    inner_misfit_method_params = set(('inner_norm', 'taper', 'filter'))
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
        seis.set_taper(conf['taper'])
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
        f.write( string )
        f.close()

    def load(self, ident, run_id='current'):
        rundir = self.make_rundir_path(run_id)
        filename = pjoin(rundir, '%s.pickle' % ident)
        f = open(filename, 'r')
        object = pickle.load(f)
        f.close()
        return object

    def pre_plot(self):
         logging.info('Starting plotting for step: %s' % self.stepname)
         
    def post_plot(self):
         logging.info('Done with plotting for step: %s' % self.stepname)

    def work(self, *args, **kwargs):
        pass
    
    def plot(self, *args, **kwargs):
        pass

class WeightMaker(Step):
    
    def __init__(self, workdir, name='weightmaker'):
        Step.__init__(self, workdir, name)
        
        self.required |= Step.inner_misfit_method_params | gen_dweights.required \
                        | set(('depth', 'moment', 'duration')) 
        
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
                                "rise-time": float(conf['duration']) } )
        
        self.out_config.receiver_weights = gen_dweights(seis, base_source, **conf )
        self.post_work(True)

class ParamTuner(Step):
     
    def __init__(self, workdir, params=['time'], name=None):
        if name is None: name = '-'.join(params)+'-tuner'
        Step.__init__(self, workdir, name)
        self.params = params
        
        self.required |= Step.outer_misfit_method_params | Step.inner_misfit_method_params \
                        | set(('duration',)) \
                        | set(('depth', 'moment', 'strike', 'dip', 'slip_rake')) \
                        | set([param+'_range' for param in self.params])
                        
        self.optional |= set(('north_shift', 'east_shift', 'time', 'bord_radius'))
        
    def work(self, search=True, forward=True, run_id='current'):
        self.pre_work(search or forward)
        seis = self.seismosizer
        conf = self.in_config.get_config()
        mm_conf = self.in_config.get_config(keys=Step.outer_misfit_method_params)
        
        sourcetype = 'eikonal'
        base_source = source.Source( sourcetype, { 'bord-radius': 0.,
                                                   'rise-time': float(conf['duration']),
                                                   'moment': float(conf['moment'])  } )
        
        for p in 'time', 'north-shift', 'east-shift', 'depth', 'strike', 'dip', 'slip-rake', 'bord-radius':
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
            str_result = stats[u2d(param)].str_best_and_confidence()
            logging.info(str_result)
            self.result(str_result, param )
            base_source[u2d(param)] = stats[u2d(param)].best
            self.out_config.__dict__[param] = stats[u2d(param)].best
        
        if forward:
            self.snapshot( base_source, 'best' )
            
        self.post_work(search or forward)
        
    def plot( self, run_id='current' ):
        self.pre_plot()
        rundir = self.make_rundir_path(run_id)
                
        finder = self.load(self.stepname, run_id=run_id)
        finder.plot( pjoin(rundir,self.stepname) )
        self.post_plot()
        
class PlaneTuner(Step):
     
    def __init__(self, workdir, name='planetuner'):
        Step.__init__(self, workdir, name)
        
        self.required |= Step.outer_misfit_method_params | Step.inner_misfit_method_params \
                        | set(('duration',)) \
                        | set(('time', 'depth', 'moment', 'strike', 'dip', 'slip_rake')) \
                        | set([param+'_range' for param in 'strike', 'dip', 'slip_rake', 'moment', 'rough_moment' ])
        self.optional |= set(('north_shift', 'east_shift'))
         
    def work(self, search=True, forward=True, run_id='current'):
        self.pre_work(search or forward)

        seis = self.seismosizer
        conf = self.in_config.get_config()
        mm_conf = self.in_config.get_config(keys=Step.outer_misfit_method_params)
        
        sourcetype = 'eikonal'
        base_source = source.Source( sourcetype, { 'bord-radius': 0.,
                                                   'rise-time': float(conf['duration'])  } )
        
        for param in 'time', 'north-shift', 'east-shift', 'depth', 'moment', 'strike', 'dip', 'slip-rake':
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
        
        if forward:
            self.snapshot( base_source, 'best' )
            
        self.post_work(search or forward)
        
    def plot( self, run_id='current' ):
        self.pre_plot()
        rundir = self.make_rundir_path(run_id)
        
        for substep in 'rough_moment', self.stepname:
            finder = self.load(substep, run_id=run_id)
            finder.plot( pjoin(rundir, substep) )
        self.post_plot()
        
class ExtensionFinder(Step):
    
    
    
    def __init__(self, workdir, name='extension'):
        Step.__init__(self, workdir, name)
        
        self.required |= Step.outer_misfit_method_params | Step.inner_misfit_method_params \
                        | set(('strike', 'dip', 'slip_rake')) \
                        | set(('time', 'depth', 'moment', 'rise_time', 'rel_rupture_velocity')) \
                        | set(('maxradius', 'delta', 'rel_rupture_velocity', 'plane')) \
                          
        self.optional |= set(('north_shift', 'east_shift'))
        
    def work(self, search=True, forward=True, run_id='current'):
        self.pre_work(search or forward)
        seis = self.seismosizer
        conf = self.in_config.get_config()
        mm_conf = self.in_config.get_config(keys=Step.outer_misfit_method_params)
        
        strike, dip, slip_rake = float(conf['strike']), float(conf['dip']), float(conf['slip_rake'])
        mt = moment_tensor.MomentTensor( strike=strike, dip=dip, rake=slip_rake )
        both_sdr = mt.both_strike_dip_rake()
        w = [ sum( [x-y for x,y in zip(both_sdr[i], (strike, dip, slip_rake)) ] ) for i in (0,1) ]
        if w[0]<w[1]:
            strike, dip, slip_rake = both_sdr[conf['plane']-1]
        else:
            strike, dip, slip_rake = both_sdr[2-conf['plane']]
                
        sourcetype = 'eikonal'
        
        base_source = source.Source( sourcetype, {  'strike':strike,
                                                    'dip':dip,
                                                    'slip-rake':slip_rake,
                                                    'bord-radius': 0.  } )
        
        for param in 'time', 'depth', 'moment', 'rise-time', 'rel-rupture-velocity', 'north-shift', 'east-shift':
            if d2u(param) in conf:
                base_source[param] = float(conf[d2u(param)])
        
        maxradius = conf['maxradius']
        delta = conf['delta']
        
        radius_grid         = ( 'bord-radius', 0., maxradius, delta ) 
        nukl_shift_x_grid   = ( 'nukl-shift-x', -maxradius, maxradius, delta )
        nukl_shift_y_grid   = ( 'nukl-shift-y', -maxradius, maxradius, delta )
        
        if search or forward: self.setup_inner_misfit_method()
        if search:
            finder = MisfitGrid( base_source, [radius_grid, nukl_shift_x_grid, nukl_shift_y_grid] )
            finder.compute(seis)
        else:
            finder = self.load(self.stepname, run_id=run_id)
        self.dump(finder, self.stepname)
        finder.postprocess(**mm_conf)
        self.dump(finder, self.stepname)
        
        stats = finder.stats
        for param in 'bord-radius', 'nukl-shift-x', 'nukl-shift-y':
            str_result = stats[param].str_best_and_confidence()
            logging.info(str_result)
            self.result(str_result, param )
            base_source[param] = stats[param].best
            self.out_config.__dict__[d2u(param)] = stats[param].best
    
        if forward:
            self.snapshot( base_source, 'best' )
            
        self.post_work(search or forward)        
   
    def plot( self, run_id='current' ):
        self.pre_plot()
        rundir = self.make_rundir_path(run_id)
        
        finder = self.load(self.stepname, run_id=run_id)
        finder.plot( pjoin(rundir,self.stepname) )
        self.post_plot()
    
class TracePlotter(Step):
    
    def __init__(self, workdir, snapshots, name='traceplotter'):
        Step.__init__(self, workdir, name)
        self.snapshots = snapshots
        
        self.required |= set()
                        
        self.optional |= set()
    
    def plot(self, run_id='current'):
    
        rundir = self.make_rundir_path(run_id)
        dirname = pjoin(rundir, 'plots')
        
        if os.path.exists( dirname ): 
            shutil.rmtree( dirname )
        os.makedirs( dirname )
        
        loaded_snapshots = []
        for i, (step, ident) in enumerate(self.snapshots):
            
            loaded_snapshots.append(step.get_snapshot(ident, run_id='current'))
            
        compos = set()
        for recs in zip(*loaded_snapshots):
            for rec in recs:
                compos.update(rec.components)
            
        ordered_compos = []
        for c in 'wesndu':
            if c in compos:
                ordered_compos.append(c)
                
        nrecs = len(zip(*loaded_snapshots))
        
        plural = { 'seismogram': 'seismograms',
                   'spectrum': 'spectra' }
                   
        for typ in 'seismogram', 'spectrum':
        
            if config.show_progress:
                widgets = ['Plotting %s' % plural[typ], ' ',
                        progressbar.Bar(marker='-',left='[',right=']'), ' ',
                        progressbar.Percentage(), ' ',]
                
                pbar = progressbar.ProgressBar(widgets=widgets, maxval=nrecs).start()
            filez = []
            dummy = (num.arange(1), num.arange(1))
            for irec, recs in enumerate(zip(*loaded_snapshots)):
                data_by_compo = []
                data_range = [None,None]
                x_range = [None,None]
                for c in ordered_compos:
                    data = []
                    for r in recs:
                        icomp = r.components.find(c)
                        if icomp >= 0:
                            if typ == 'seismogram':
                                dsyn = r.syn_seismograms[icomp]
                                dref = r.ref_seismograms[icomp]
                            elif typ == 'spectrum':
                                dsyn = r.syn_spectra[icomp]
                                dref = r.ref_spectra[icomp]
                            data.append(dsyn)
                            data.append(dref)
                            grow( data_range, ra(dsyn[1]), ra(dref[1]) )
                            grow( x_range,  nonzero_range( dsyn ), nonzero_range( dref ))
                        else:
                            data.append(dummy)
                            data.append(dummy)
                        
                    data_by_compo.append((c, data))
                
                conf = {}
                proto = recs[0]
                conf['title'] = 'Receiver %i' % (irec+1)
                if all([r.name == proto.name for r in recs]):
                    conf['title'] += ': %s' % proto.name
                    
                conf['yrange'] = data_range
                conf['xrange'] = x_range
                    
                filename = pjoin(dirname, '%s_%i.pdf' % (typ,irec+1))
                
                plotting.seismogram_plot(data_by_compo, filename, conf_overrides=conf)
                if config.show_progress: pbar.update(irec+1)
                filez.append(filename)
            
            filename = pjoin(dirname, '%s_all.pdf' % plural[typ])
            plotting.pdfjoin(filez, filename)
            
            if config.show_progress: pbar.finish()
        
        
class MisfitGridStats:
    
    def str_best_and_confidence(self):
        lw = ''
        uw = ''
        if self.percentile16_warn: lw = ' (?)'
        if self.percentile84_warn: uw = '(?) '
           
        return '%s = %g   (confidence interval 68%%) = [ %g%s, %g %s]' % \
                (self.paramname, self.best, self.percentile16, lw, self.percentile84, uw)

    def str_mean_and_stddev(self):
        return '%(paramname)s = %(mean)g +- %(std)g' % self.__dict__

class MisfitGrid:
    '''Brute force grid search minimizer with builtin bootstrapping.'''
    
    def __init__( self, base_source,
                        param_ranges,
                        source_constraints=None):
        
        self.base_source = copy.copy(base_source)
        self.param_ranges = param_ranges
        self.sources = self.base_source.grid( self.param_ranges, 
                                              source_constraints=source_constraints)

        self.sourceparams = [ x[0] for x in self.param_ranges ]

        # will be set by compute()
        self.misfits_by_src = None
        self.norms_by_src = None
        self.receivers = None
        self.nreceivers = None
        
        # will be set by postprocess()
        self.best_source = None
        self.misfits_by_s = None
        self.bootstrap_sources = None
        self.stats = None
        
    def compute(self, seis):
        '''Let seismosizer do the calculations'''
        
        nsources = len(self.sources)
        nreceivers = len(seis.receivers)
        ncomponents = max([ len(r.components) for r in seis.receivers ])
        
        
        # results gathered by (source,receiver,component)
        misfits_by_src = num.zeros( (nsources, nreceivers, ncomponents), dtype=num.float)
        norms_by_src = num.zeros(  (nsources, nreceivers, ncomponents), dtype=num.float)
        
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
                    for icomp, comp in enumerate(receiver.components):
                        misfits_by_src[isource, ireceiver, icomp] = receiver.misfits[icomp]
                        norms_by_src[isource, ireceiver, icomp] = receiver.misfit_norm_factors[icomp]
                        
            except seismosizer.SeismosizersReturnedErrors:
                failings.append(isource)
            
            if config.show_progress: pbar.update(isource+1)
            
        if config.show_progress: pbar.finish()
        
        self.misfits_by_src = misfits_by_src
        self.norms_by_src = norms_by_src
        self.nreceivers = nreceivers
        self.receivers = seis.receivers
        
        self.best_source = None
        self.misfits_by_s = None
        self.bootstrap_sources = None
        self.stats = None
    
    def postprocess(self, **outer_misfit_config):
        self.best_source, self.misfits_by_s = self._best_source(**outer_misfit_config)
        self.bootstrap_sources = self._bootstrap(**outer_misfit_config)
        stepsizes = [ x[3] for x in self.param_ranges ]
        self.stats = self._stats(self.param_ranges, self.best_source, self.bootstrap_sources)
    
    def get_mean_misfits_by_r(self):
        mean_misfits_by_r = num.zeros( self.nreceivers, dtype=num.float )
        for irec in range(self.nreceivers):
            mean_misfits_by_r[irec] = num.sum(self.misfits_by_src[:,irec,:])/(len(self.receivers[irec].components)*len(self.sources))
            
        return mean_misfits_by_r
        
    def _best_source(self, **outer_misfit_config):
        misfits_by_s = self._make_global_misfits(**outer_misfit_config)
        return self.sources[num.argmin(misfits_by_s)], misfits_by_s
        
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
        
    def _make_global_misfits(self, receiver_weights=1., outer_norm='l2norm', anarchy=False, bootstrap=False, **kwargs):
        
        if isinstance(receiver_weights, float):
            rweights = receiver_weights
        else:
            rweights = receiver_weights[num.newaxis,:].copy()
        
        if bootstrap:
            bweights = num.zeros(self.nreceivers, dtype=num.float)
            bweights_x = num.bincount(num.random.randint(0,self.nreceivers,self.nreceivers))
            bweights[:len(bweights_x)] = bweights_x[:]
        
        if outer_norm == 'l1norm':
            misfits_by_sr = num.sum(self.misfits_by_src,2)
            norms_by_sr   = num.sum(self.norms_by_src,2)
            
            if anarchy:
                rweights /= num.where( norms_by_sr != 0., norms_by_sr, -1.)
                rweights = num.maximum(rweights, 0.)
            
            if bootstrap:
                rweights *= bweights
            
            ms = num.sum( rweights*misfits_by_sr, 1 )
            ns = num.sum( rweights*norms_by_sr, 1 )
            
            misfits_by_s  = num.where(ns > 0., ms/ns, -1.)
            maxm = num.amax(misfits_by_s)
            misfits_by_s = num.where(misfits_by_s<0, maxm, misfits_by_s)
            
        elif outer_norm == 'l2norm':
            misfits_by_sr = num.sqrt(num.sum(self.misfits_by_src**2,2))
            norms_by_sr   = num.sqrt(num.sum(self.norms_by_src**2,2))
            
            if anarchy:
                rweights /= num.where( norms_by_sr != 0., norms_by_sr, -1.)
                rweights = num.maximum(rweights, 0.)
                
            if bootstrap:
                rweights *= num.sqrt(bweights)
            
            ms = num.sum( (rweights*misfits_by_sr)**2, 1 )
            ns = num.sum( (rweights*norms_by_sr)**2, 1 )
            
            misfits_by_s  = num.where(ns > 0., num.sqrt(ms/ns), -1.)
            maxm = num.amax(misfits_by_s)
            misfits_by_s = num.where(misfits_by_s<0, maxm, misfits_by_s)
        
        else:
            raise Exception('unknown norm method: %s' % outer_norm)
        
        return misfits_by_s
        
    def plot(self, dirname):
        
        if os.path.exists( dirname ): 
            shutil.rmtree( dirname )
        os.makedirs( dirname )
        
        best_source = self.best_source
        bootstrap_sources = self.bootstrap_sources
        
        for iparam, param in enumerate(self.sourceparams):
            
            #
            # 1D misfit cross section
            #
            xdata = num.array([ s[param] for s in self.sources ], dtype=num.float)
            conf = dict( xlabel = param.title(),
                         xunit = self.base_source.sourceinfo(param).unit )
            
            plotting.km_hack(conf)
            plotting.misfit_plot_1d( [(xdata, self.misfits_by_s)],
                                     pjoin(dirname, 'misfit-%s.pdf' % param),
                                     conf )
                                     
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
            plotting.histogram_plot_1d( [(xdata, hist)],
                                        pjoin(dirname, 'histogram-%s.pdf' % param),
                                        conf )
            
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
                plotting.misfit_plot_2d( [(ax, ay, az)],
                                         pjoin(dirname, 'misfit-%s-%s.pdf' % (xparam, yparam)),
                                         conf )
                
                conf['xrange'] = (num.amin(ax), num.amax(ax))
                conf['yrange'] = (num.amin(ay), num.amax(ay))
                
                plotting.histogram_plot_2d( [(axc, ayc, ac)],
                                         pjoin(dirname, 'histogram-%s-%s.pdf' % (xparam, yparam)),
                                         conf )
                
                conf['zrange'] = (num.amin(az), num.amax(az))
                plotting.misfogram_plot_2d( [(ax, ay, az), best_loc, (axc, ayc, ac)],
                                        pjoin(dirname, 'misfogram-%s-%s.pdf' % (xparam, yparam)),
                                        conf )
                                        
                                        
                                        