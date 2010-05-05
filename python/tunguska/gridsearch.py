import config
import plotting
import seismosizer
import util
import logging

import re, copy
import progressbar 
import numpy as num
import scipy.stats
from os.path import join as pjoin

def mimainc_to_gvals(mi,ma,inc):
    vmin, vmax, vinc = float(mi), float(ma), float(inc)
    n = int(round((vmax-vmin)/vinc))+1
    vinc = (vmax-vmin)/(n-1)
    return num.array([ vmin+i*vinc for i in xrange(n) ], dtype=num.float)

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

def mimainc_to_gvals(mi,ma,inc):
    vmin, vmax, vinc = float(mi), float(ma), float(inc)
    n = int(round((vmax-vmin)/vinc))+1
    vinc = (vmax-vmin)/(n-1)
    return num.array([ vmin+i*vinc for i in xrange(n) ], dtype=num.float)


class MisfitGridStats:
    def __init__(self, paramname, best, distribution, tested_values=None):
        self.paramname = paramname
        self.best = best
        self.distribution = distribution
        self.tested_values = tested_values
        self.mean = num.mean(distribution)
        self.std = num.std(distribution)
        self.median = num.median(distribution)
        self.percentile16 = scipy.stats.scoreatpercentile(distribution, 16.)
        self.percentile84 = scipy.stats.scoreatpercentile(distribution, 84.)
        
        if tested_values is not None:
            self.percentile16 -= step_at( tested_values, self.percentile16 )/2.
            self.percentile84 += step_at( tested_values, self.percentile84 )/2.
            self.percentile16_warn = self.percentile16 < num.min(tested_values)
            self.percentile84_warn = self.percentile84 > num.max(tested_values)
        else:
            self.percentile16_warn = False
            self.percentile84_warn = False
        
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
    
    def as_xml(self):
        tmpl = util.unindent('''
        <parameter>
            <name>%s</name>
            <value>%e</value>
            <confidenceinterval>
                <interval>68</interval>
                <low>%e</low>
                <high>%e</high>
                <low_unclear>%i</low_unclear>
                <high_unclear>%i</high_unclear>
            </confidenceinterval>
        </parameter>
        ''')
        return tmpl % (self.paramname.title(), self.best, 
                       self.percentile16, self.percentile84,
                       self.percentile16_warn, self.percentile84_warn)
                       
    def converted(self, paramname, function):
        best = function(self.best)
        distribution = function(self.distribution)
        if self.tested_values is not None:
            tested_values = function(self.tested_values)
        else:
            tested_values = None
            
        return MisfitGridStats(paramname, best, distribution, tested_values=tested_values)


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
            self.param_values = []
            for param, gvalues in param_values:
                self.param_values.append((param, num.asarray(gvalues)))
            
        else:
            self.param_values = []
            for param, mi, ma, inc in param_ranges:
                self.param_values.append( (param, mimainc_to_gvals(mi,ma,inc)) )
            
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
        self.nreceivers_enabled = None
        
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
        nreceivers_enabled = len( [ rec for rec in seis.receivers if rec.enabled ] )
        receiver_mask = num.array( [ rec.enabled for rec in seis.receivers ], dtype=num.bool )
        
        # results, gathered by (source,receiver,component)
        misfits_by_src, norms_by_src, failings = seis.make_misfits_for_sources( 
                                                    self.sources,
                                                    show_progress=config.show_progress,
                                                    progress_title=progress_title)
              
        # results for reference source, gathered by (source,receiver,component)
        ref_misfits_by_src, ref_norms_by_src, failings = seis.make_misfits_for_sources([self.ref_source])
        
        self.misfits_by_src = misfits_by_src
        self.norms_by_src = norms_by_src
        self.ref_misfits_by_src = ref_misfits_by_src
        self.ref_norms_by_src = ref_norms_by_src
        
        self.nreceivers = nreceivers
        self.nreceivers_enabled = nreceivers_enabled
        self.receiver_mask = receiver_mask
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
            rec = self.receivers[irec]
            ncomps = len(rec.components)
            if ncomps != 0:
                x = num.sum(self.misfits_by_src[:,irec,:ncomps]) /( ncomps*len(self.sources))
            else:
                x = -1.0
            
            mean_misfits_by_r[irec] = x
        return mean_misfits_by_r
        
    def get_median_of_misfits_by_r(self):
        '''Get the median of the misfits per station of best source.'''
        misfits = []
        for irec in range(self.nreceivers):
            if self.receivers[irec].enabled:
                misfits.append(self.misfits_by_r[irec])
        
        if len(misfits) == 0:
            raise Exception('No receivers are enabled')
        
        return num.median(misfits)
        
    def get_best_misfit(self):
        return num.nanmin(self.misfits_by_s)
        
    def _best_source(self, return_misfits_by_r=False, **outer_misfit_config):
        misfits_by_s, misfits_by_sr = seismosizer.make_global_misfits( 
            self.misfits_by_src, self.norms_by_src, 
            receiver_mask=self.receiver_mask,
            **outer_misfit_config )
        
        ibest = num.nanargmin(misfits_by_s)
        if not return_misfits_by_r:
            return self.sources[ibest], misfits_by_s

        else:
            # misfit variability by receiver
            misfits_varia_by_r = num.std(misfits_by_sr,0)
            return self.sources[ibest], misfits_by_s, misfits_by_sr[ibest,:], misfits_varia_by_r
        
    def _ref_misfits_by_r(self, **outer_misfit_config):
        misfits_by_s, misfits_by_sr = seismosizer.make_global_misfits(
            self.ref_misfits_by_src, self.ref_norms_by_src, 
            receiver_mask=self.receiver_mask, **outer_misfit_config)
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

            distribution = num.array([ source[param] for source in bootstrap_sources ], dtype=num.float)
            result = MisfitGridStats( param, best_source[param], distribution, tested_values=gvalues )
            results[param] = result
            
        return results
        

        
    def plot(self, dirname, nsets=1, source_model_infos=None, conf_overrides=None):
        
        best_source = self.best_source
        bootstrap_sources = self.bootstrap_sources
        
        
        plot_files = []
        
        param_min_misfits = {} 
        
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
            if conf_overrides:
                conf.update(conf_overrides)
            
            plotting.km_hack(conf)
            fn = 'misfit-%s.pdf' % param
            plotting.misfit_plot_1d( [(xdata,ydata)],
                                     pjoin(dirname, fn),
                                     conf, apply_moment_to_magnitude_hack=True )
            plot_files.append(fn)
            
            
            
            #
            # 1D histogram
            #
            gvalues = self.param_values[iparam][1]
            gedges = values_to_bin_edges(gvalues)
            
            kwargs = {}
            if util.cmp_version(num.version.short_version, '1.3.0') < 0:
                kwargs = {'new': True}
                
            hist, edges = num.histogram(self.stats[param].distribution,
                                        bins=gedges, **kwargs)
            
            hist = hist/float(len(bootstrap_sources))
            
            conf = dict( xlabel = param.title(),
                         xunit = self.base_source.sourceinfo(param).unit)
            if conf_overrides:
                conf.update(conf_overrides)
            
            plotting.km_hack(conf)
            fn = 'histogram-%s.pdf' % param
            plotting.histogram_plot_1d( edges, hist, 
                                        pjoin(dirname, fn),
                                        conf, apply_moment_to_magnitude_hack=True)
            
            plot_files.append( fn )
            
            #
            #
            # gather data for 1D projected misfit cross sections
            indi = num.digitize(xdata, gedges) - 1
            for_minmisfits = []
            for ival in range(gvalues.size):
                vydata = num.extract(indi == ival, ydata)
                if vydata.size != 0:
                    for_minmisfits.append((gvalues[ival], vydata.min()))
                        
            param_min_misfits[param] = num.array(for_minmisfits,dtype=num.float).transpose()
        
        # data range for projected misfit cross section
        mi, ma = num.inf, -num.inf
        for xdata,ydata in param_min_misfits.values():
            mi = min(mi, ydata.min())
            ma = max(ma, ydata.max())
            
        param_min_misfits_range = (mi,ma)
        
        # 1D projected misfit cross sections
        for iparam, param in enumerate(self.sourceparams):
            mini, maxi = param_min_misfits_range
            
            if mini < 0.95:
                maxi = 1.
                       
            else:
                if mini < 1.95:
                    maxi = 2.
                       
            conf = dict( xlabel = param.title(),
                         xunit = self.base_source.sourceinfo(param).unit,
                         ylimits = (mini, maxi),
                         ymasking=False)
            
            if conf_overrides:
                conf.update(conf_overrides)

            plotting.km_hack(conf)
            fn = 'min-misfit-%s.pdf' % param
                        
            plotting.misfit_plot_1d( [param_min_misfits[param]],
                                     pjoin(dirname, fn),
                                     conf,apply_moment_to_magnitude_hack=True )
            plot_files.append(fn)
            
            fn_table = 'min-misfit-%s.table' % param
            num.savetxt(pjoin(dirname, fn_table),  param_min_misfits[param].transpose())
                        
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
                if len(bootstrap_sources) > 1:
                    bootstrap_data = (axc, ayc, ac)
                else:
                    bootstrap_data = ([],[],[])
                
                if az.min() < 0.95:
                    zlimits = (az.min(),min(1.,az.max()))
                else:
                    zlimits = (az.min(),az.max())
                
                if 'misfit_limits' in dir(config):
                    zlimits = config.misfit_limits
                
                conf = dict( xlabel = xparam.title(),
                             xunit = self.base_source.sourceinfo(xparam).unit,
                             ylabel = yparam.title(),
                             yunit = self.base_source.sourceinfo(yparam).unit,
                             zlimits = zlimits,
                             zmasking=False,
                             zsnap = True
                         )
                
                if conf_overrides:
                    for k in conf_overrides:
                        if not k.endswith('snap'):
                            conf[k] = conf_overrides[k]
                    
                
                best_loc = (num.array([ best_source[xparam]]), num.array([best_source[yparam]]))
                plotting.km_hack(conf)
                plotting.nukl_hack(conf)
                
                fn = 'misfogram-%s-%s.pdf' % (xparam, yparam)
                plotting.misfogram_plot_2d_gmtpy( [(ax, ay, az), best_loc, bootstrap_data],
                                        pjoin(dirname, fn),
                                        conf, apply_moment_to_magnitude_hack=True )
                plot_files.append(fn)
        
        #
        # Station plot
        #
        lats = num.array( [ r.lat for r in self.receivers ], dtype='float' )
        lons = num.array( [ r.lon for r in self.receivers ], dtype='float' )
        dists = num.array( [ r.distance_deg for r in self.receivers ], dtype='float' )
        rnames = [ ' '.join(r.name.split('.')) for r in self.receivers ]
        slat, slon = self.source_location[:2]
        station_misfits = self.misfits_by_r-self.ref_misfits_by_r
        #station_varia = self.variability_by_r / num.sum(self.variability_by_r) * len(self.receivers)
        station_varia = self.misfits_by_r / num.sum(self.misfits_by_r) * len(self.receivers)
        plotting.station_plot( slat, slon, lats, lons, rnames, station_misfits, station_varia, 
                              best_source, num.amax(dists)*1.05, pjoin(dirname, 'stations.pdf'), {}, nsets=nsets)
        plot_files.append('stations.pdf')
        
        indicate_plane = 0
        if source_model_infos:
            source_size = 0.
            if 'bord-radius' in best_source.keys():
                source_size = best_source['bord-radius']
                
            delta_lat = 1.5*max(source_size*2.,50000)/(20000.*1000.)*180
            
            si = None
            if 'bord-radius' in best_source.keys() and best_source['bord-radius'] > 0.0:
                si = source_model_infos
                indicate_plane = 1
            
            plotting.location_map(pjoin(dirname, 'location.pdf'), slat, slon, delta_lat, {}, source=best_source,
                                  source_model_infos=si, receivers=(lats,lons,rnames))
                                  
            plot_files.append('location.pdf')
        
        fns = plotting.beachball(best_source, pjoin(dirname, 'beachball.pdf'), indicate_plane=indicate_plane)
        if fns:
            plot_files.append('beachball.pdf')
        
        
        
        return plot_files
    