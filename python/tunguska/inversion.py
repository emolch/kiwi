import receiver
import gfdb
import seismosizer
import config
import source
import util
import plotting

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

from os.path import join as pjoin    

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


class Step:
    def __init__(self, workdir, name):
        self.baseworkdir = workdir
        self.stepname = name
        self.in_config = None
        self.out_config = None
        self.stepdir = pjoin(self.baseworkdir, self.stepname)
        
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
    
    def pre_work(self):
        assert(self.in_config is not None)
        rundir = self.make_rundir_path('incomplete')
        if os.path.exists( rundir ): 
            shutil.rmtree( rundir )
        os.makedirs( rundir )
        self.in_config.dump( pjoin(rundir,'config-in.pickle') )
        self.out_config = config.Config()
                
    def post_work(self):
        rundir = self.make_rundir_path('incomplete')
        self.out_config.dump(pjoin(rundir, 'config-out.pickle'))
        if os.path.exists( self.make_rundir_path('current') ):
            shutil.move(self.make_rundir_path('current'), self.next_available_rundir())
        shutil.move(rundir, self.make_rundir_path('current'))
        
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

class AdjustPointsource(Step):
    
    def __init__(self, workdir, name='adjustment'):
        Step.__init__(self, workdir, name)
        
        self.timefinder_method = dict( anarchy=True, norm='l1norm', niterations=1000 )
        self.momentfinder_method = dict( anarchy=True, norm='l1norm', niterations=1000 )
        self.sdrfinder_method = dict( anarchy=True, norm='l1norm', niterations=1000 )
        
    def work(self):
        self.pre_work()
        conf = self.in_config.get_config()
        seis = standard_setup( **conf )
        
        seis.set_taper(conf['taper'])
        seis.set_filter(conf['filter'])
        
        sourcetype = 'eikonal'
        base_source = source.Source( sourcetype, 
                               {"strike":float(conf['strike']),
                                "dip":float(conf['dip']),
                                "slip-rake":float(conf['slip_rake']),
                                "depth": float(conf['depth']),
                                "bord-radius": 0.,
                                "moment": float(conf['moment'])/100.,
                                "rise-time": float(conf['duration']) } )
        
        
        #
        # Time
        # 
        
        seis.set_misfit_method('l1norm')
        timefinder = MisfitGrid( base_source, [( 'time', )+conf['time_range']])
        timefinder.compute(seis)
        self.dump(timefinder, 'time')
        
        r = timefinder.bootstrap(**self.timefinder_method)[0]
        time = r['time'].mean
        self.out_config.time = time
        str_result = 'time = %(mean)g +- %(std)g\n' % r['time'].__dict__
        self.result( str_result, 'time' )
        print str_result
        base_source['time'] = time
        
        
        #
        # Rough Moment
        #
        seis.set_misfit_method('ampspec_l1norm')
        
        base_source['moment'] =float(conf['moment'])
        moment = base_source['moment']
        moment_grid    = ('moment', moment*0.1, moment*2., moment*0.1 )
        momentfinder = MisfitGrid( base_source, [moment_grid])
        momentfinder.compute(seis)
        self.dump(momentfinder, 'rough_moment')
        
        r = momentfinder.bootstrap(**self.momentfinder_method)[0]
        str_result = 'rough moment = %(mean)g +- %(std)g\n' % r['moment'].__dict__
        print str_result
        base_source['moment'] = r['moment'].mean
        
                
        #
        # Strike, Dip, Rake, Moment
        #
        seis.set_misfit_method('ampspec_l1norm')
        
        strike, dip, slip_rake, moment = [ base_source[param] for param in ('strike', 'dip', 'slip-rake', 'moment') ]
        
        strike_grid    = ('strike', strike - 20., strike + 20., 5.)
        dip_grid       = ('dip', dip - 30., dip + 30., 5.)
        slip_rake_grid = ('slip-rake', slip_rake-40., slip_rake+40., 10.)
        moment_grid    = ('moment', moment*0.5, moment*2.0, moment*0.1 )
        #sdrfinder = MisfitGrid( base_source, [strike_grid, dip_grid] )
        sdrfinder = MisfitGrid( base_source, [strike_grid, dip_grid, slip_rake_grid, moment_grid] )
        sdrfinder.compute(seis)
        self.dump(sdrfinder, 'sdr')
        
        seis.close()
        self.post_work()
        
    def plot( self, run_id ):
        rundir = self.make_rundir_path(run_id)
        timefinder = self.load('time', run_id=run_id)
        timefinder.plot( pjoin(rundir,'timefinder'), **self.timefinder_method )
        
        momentfinder = self.load('rough_moment', run_id=run_id)
        momentfinder.plot( pjoin(rundir, 'momentfinder'), **self.momentfinder_method )
        
        sdrfinder = self.load('sdr', run_id= run_id)
        sdrfinder.plot( pjoin(rundir,'sdrfinder'), **self.sdrfinder_method )
        
        
class RetrieveExtension(Step):
    
    def __init__(self, workdir, name='extension'):
        Step.__init__(self, workdir, name)
        
    def work(self):
        self.pre_work()
        self.post_work()
        
    def plot( self, run_id ):
        pass
    
    

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

        self.misfits_by_src = None
        self.norms_by_src = None
        self.nreceivers = None
        
    def compute(self, seismosizer):
        '''Let seismosizer do the calculations'''
        
        nsources = len(self.sources)
        nreceivers = len(seismosizer.receivers)
        ncomponents = max([ len(r.components) for r in seismosizer.receivers ])
        
        
        # results gathered by (source,receiver,component)
        misfits_by_src = num.zeros( (nsources, nreceivers, ncomponents), dtype=num.float)
        norms_by_src = num.zeros(  (nsources, nreceivers, ncomponents), dtype=num.float)
        
        if config.show_progress:
            if len(self.sourceparams) == 1:
                title = 'Grid search param: ' + ', '.join(self.sourceparams) + ' '
            else:
                title = 'Grid search params: ' + ', '.join(self.sourceparams) + ' '
            widgets = [title,
                    progressbar.Bar(marker='-',left='[',right=']'), ' ',
                    progressbar.Percentage(), ' ',]
            
            pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(self.sources)).start()
            
        for isource, source in enumerate(self.sources):
            #
            # ??? loop over different taper-filter-misfit setups ???
            #
            seismosizer.make_misfits_for_source( source )
            for ireceiver, receiver in enumerate(seismosizer.receivers):
                for icomp, comp in enumerate(receiver.components):
                    misfits_by_src[isource, ireceiver, icomp] = receiver.misfits[icomp]
                    norms_by_src[isource, ireceiver, icomp] = receiver.misfit_norm_factors[icomp]
                    
            if config.show_progress: pbar.update(isource+1)
            
        if config.show_progress: pbar.finish()
        
        self.misfits_by_src = misfits_by_src
        self.norms_by_src = norms_by_src
        self.nreceivers = nreceivers
         
    
        
    def best_source(self, **outer_misfit_config):
        misfits_by_s = self.make_global_misfits(**outer_misfit_config)
        return self.sources[num.argmin(misfits_by_s)]
        
    def bootstrap(self, niterations=1000, **outer_misfit_config):
        best_sources = []
        for i in xrange(niterations):
            best_sources.append( self.best_source(bootstrap=True, **outer_misfit_config) )
        
        return self._bootstrap_stats(best_sources, self.sourceparams), best_sources
        
    def _bootstrap_stats(self, sources, sourceparams):
        
        class Stats:
            pass
        
        results = {}
        for param in self.sourceparams:
            param_results = num.zeros(len(sources), dtype=num.float)
            for i, source in enumerate(sources):
                param_results[i] = source[param]
                
            result = Stats()
            result.mean = num.mean(param_results)
            result.std = num.std(param_results)
            result.median = num.median(param_results)
            result.percentile16 = scipy.stats.scoreatpercentile(param_results, 16.)
            result.percentile84 = scipy.stats.scoreatpercentile(param_results, 84.)
            result.distribution = param_results
            results[param] = result
            
        return results
        
    def make_global_misfits(self, weights_by_r=1., norm='l2norm', anarchy=False, bootstrap=False, **kwargs):
        
        if isinstance(weights_by_r, float):
            rweights = weights_by_r
        else:
            rweights = weights_by_r[num.newaxis,:].copy()
        
        if bootstrap:
            bweights = num.zeros(self.nreceivers, dtype=num.float)
            bweights_x = num.bincount(num.random.randint(0,self.nreceivers,self.nreceivers))
            bweights[:len(bweights_x)] = bweights_x[:]
        
        if norm == 'l1norm':
            misfits_by_sr = num.sum(self.misfits_by_src,2)
            norms_by_sr   = num.sum(self.norms_by_src,2)
            
            if anarchy:
                rweights /= num.where( norms_by_sr != 0., norms_by_sr, -1.)
                rweights = num.maximum(rweights, 0.)
            
            if bootstrap:
                rweights *= bweights
            
            misfits_by_s  = (num.sum( rweights*misfits_by_sr, 1 ) /
                             num.sum( rweights*norms_by_sr, 1 ))
        
        elif norm == 'l2norm':
            misfits_by_sr = num.sqrt(num.sum(self.misfits_by_src**2,2))
            norms_by_sr   = num.sqrt(num.sum(self.norms_by_src**2,2))
            
            if anarchy:
                rweights /= num.where( norms_by_sr != 0., norms_by_sr, -1.)
                rweights = num.maximum(rweights, 0.)
                
            if bootstrap:
                rweights *= num.sqrt(bweights)
            
            misfits_by_s  = num.sqrt(num.sum( (rweights*misfits_by_sr)**2, 1 ) /
                                     num.sum( (rweights*norms_by_sr)**2, 1 ))
        
        else:
            raise Exception('unknown norm method: %s' % norm)
        
        return misfits_by_s
        
    def plot(self, dirname, **misfit_setup_method):
        
        if os.path.exists( dirname ): 
            shutil.rmtree( dirname )
        os.makedirs( dirname )
        misfits_by_s = self.make_global_misfits(**misfit_setup_method)
        best_source = self.sources[num.argmin(misfits_by_s)]
        
        r,best_sources = self.bootstrap( **misfit_setup_method)
        
        for iparam, param in enumerate(self.sourceparams):
            
            #
            # 1D misfit cross section
            #
            xdata = num.array([ s[param] for s in self.sources ], dtype=num.float)
            conf = dict( xlabel = param.title(),
                         xunit = self.base_source.sourceinfo(param).unit )
            
            plotting.km_hack(conf)
            plotting.misfit_plot_1d( [(xdata, misfits_by_s)],
                                     pjoin(dirname, 'misfit-%s.pdf' % param),
                                     conf )
                                     
            #
            # 1D histogram
            #
            # mi = num.min(r[param].distribution)
            # ma = num.max(r[param].distribution)
            mi, ma, step = self.param_ranges[iparam][1:]
            n = int(round((ma-mi)/step))+1
            hist, edges = num.histogram(r[param].distribution,
                                        bins=n,
                                        range=(mi-step/2., ma+step/2.),
                                        new=True)
            
            hist = hist/float(misfit_setup_method['niterations'])
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
                maxmisfit = num.amax(misfits_by_s)
                for (source, misfit) in zip(self.sources, misfits_by_s):
                    vx = source[xparam]
                    vy = source[yparam]
                    vmisfits[(vx,vy)] = min(vmisfits.get((vx,vy), maxmisfit), misfit)
                    
                vcounts = {}
                for source in best_sources:
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
                ac = ac/len(best_sources)
                
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
                                        
                                        
                                        