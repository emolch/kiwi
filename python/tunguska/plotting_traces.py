

import gmtpy
from gmtpy import cm
import sys
import config
from os.path import join as pjoin
import plotting
import scipy.stats
import numpy as num
import math
from pyrocko import orthodrome, trace

import tracy
        
class PileTracy(tracy.Tracy):
    
    def set_pile(self, pile):
        traces = pile.all()
        tmin = min([trace.tmin for trace in traces ])
        for trace in traces:
            trace.shift(-tmin)
        
        self.set_traces(traces)
    
    def xminmax(self, traces):
        return trace.minmaxtime(traces, key=self.map_xscaling)
    
    def yminmax(self, traces):
        return trace.minmax(traces, key=self.map_yscaling)
        
    def data(self, trace):
        return trace.get_xdata(), trace.get_ydata()

class MyTracy(tracy.Tracy):
    
    def __init__(self, source_location=(0.,0.,0.), plot_type='seismograms', *args, **kwargs):
        tracy.Tracy.__init__(self, *args, **kwargs)
        self.mode = 3.
        self.percentile = 10.
        self.nxauxwidgets = 1
        self.slat, self.slon, self.stime = source_location
        self.plot_type = plot_type
        
    def map_xgroup(self, trace):
        return trace.channel
    
    def map_ygroup(self, trace):
        return trace.station, trace.network, trace.azimuth, trace.distance_deg
    
    def map_yscaling(self, trace):
        return None #self.map_xgroup(trace)
    
    def order_xgroup(self, a,b):
        return cmp(b,a)
    
    def order_ygroup(self, a,b):
        return cmp((math.floor(a[2]/30.),a[3]),(math.floor(b[2]/30.),b[3]))
    
    def yminmax1(self, trace):
        mean = trace.ydata.mean()
        std = trace.ydata.std()
        mi, ma = mean-std*self.mode, mean+std*self.mode
        return mi, ma
    
    def yminmax(self, traces):
        return tracy.minmax(traces, self.yminmax1, key=self.map_yscaling, 
                mincombine = lambda mins: scipy.stats.scoreatpercentile(mins,self.percentile),
                maxcombine = lambda maxs: scipy.stats.scoreatpercentile(maxs, 100.-self.percentile) )
    
    def map_color(self, trace):
        if trace.set == 'references':
            return -1
        else:
            return trace.snapshot

    def colors(self, color_key):
        ind = self.color_index[color_key]
        if ind == 0:
            return '0.5p,%s,6_4:2p' % gmtpy.color('black')
        else:
            return '0.5p,%s' % gmtpy.color(ind-1)

    def label_xgroup(self, a):
        return a
    
    def label_ygroup(self, a):
        return '.'.join(a[:2]).rstrip('.')
    
    def draw_xstuff(self, gmt, widget, xscaler, trace, ivisit):
        if ivisit == 0:
            p = xscaler.get_params()
            ma, mi = p['xmax'], p['xmin']
            span = ma - mi
            if self.plot_type == 'seismograms':
                beg = mi, 0., self.font_size, 0., 1, 'LM', '%g s' % mi
                end = ma, 0., self.font_size, 0., 1, 'RM', '+ %g s' % span
            elif self.plot_type == 'spectra':
                beg = mi, 0., self.font_size, 0., 1, 'LT', '%.3g Hz' % mi
                end = ma, 0., self.font_size, 0., 1, 'RT', '%.3g Hz' % ma
                
            eps = span/10000.
            if self.plot_type == 'seismograms':
                gmt.pstext(in_rows=[beg,end], N=True, *(widget.JXY() + xscaler.R()))
                gmt.psxy(W='1p', in_rows=[(mi+eps,0.3),(mi+eps,0.7)], *(widget.JXY() + xscaler.R()))
                gmt.psxy(W='1p', in_rows=[(ma-eps,0.3),(ma-eps,0.7)], *(widget.JXY() + xscaler.R()))
            elif self.plot_type == 'spectra':
                gmt.pstext(in_rows=[beg,end], N=True, D='%gp/%gp' % (0., -self.font_size/4.), *(widget.JXY() + xscaler.R()))
                gmt.psxy(W='1p', in_rows=[(mi+span/20.,0.0),(mi+eps,0.0),(mi+eps,0.3)], *(widget.JXY() + xscaler.R()))
                gmt.psxy(W='1p', in_rows=[(ma-span/20.,0.0),(ma-eps,0.0),(ma-eps,0.3)], *(widget.JXY() + xscaler.R()))
            
            
    def draw_trace(self, gmt, widget, scaler, trace, ivisit):
        tracy.Tracy.draw_trace(self, gmt, widget, scaler, trace, ivisit)
        
    def draw_xaux(self, gmt, widget, trace):
        p = self.azidist_scaler.get_params()
        widget['J'] = '-JE%g/%g/%g' % (self.slon, self.slat, min(p['ymax'],180.)) + '/%(width)gp'

        phi = num.arange(361, dtype=num.float)
        circle = phi, num.ones(361)*trace.distance_deg
        direction = num.array([[trace.azimuth, trace.azimuth], [0., p['ymax']]], dtype=num.float)
        
        circle = orthodrome.azidist_to_latlon(self.slat, self.slon, *circle)
        direction = orthodrome.azidist_to_latlon(self.slat, self.slon, *direction)
        
        gmt.pscoast( R='g', B=True, D='c', A=10000, G=(200,200,200), *widget.JXY())
        gmt.psxy('-:', in_columns=circle, R='g', W='1p', *widget.JXY())
        gmt.psxy('-:', in_columns=direction, R='g', W='1p', *widget.JXY())
        
    def set_traces(self, traces):
        tracy.Tracy.set_traces(self, traces)
        dists = [ tr.distance_deg for tr in traces ]
        azimuths = [ tr.azimuth for tr in traces ]
        conf = dict(xmode='off', xlimits=(0.,360.), xinc=45., masking=False, ymode='0-max', yspace=0.1)
        axes = [ gmtpy.simpleconf_to_ax(conf,x) for x in 'xy' ]
        self.azidist_scaler = gmtpy.ScaleGuru([(azimuths, dists)], axes)
        
class UTrace:
    def __init__(self, **kwargs):
        for k,v in kwargs.iteritems():
            self.__dict__[k] = v

def multi_seismogram_plot2(snapshots, source_locations, plotdir):
    
    plural = { 'seismogram': 'seismograms',
                'spectrum': 'spectra' }
                
    datasource = { ('synthetics', 'spectrum'): 'syn_spectra',
                    ('references', 'spectrum'): 'ref_spectra',
                    ('synthetics', 'seismogram'): 'syn_seismograms',
                    ('references', 'seismogram'): 'ref_seismograms' }
    
    source_location = source_locations[0]
    
    fns_all = []
    for typ in 'seismogram', 'spectrum':
        traces = []
        fns = []
        for isnap, receivers in enumerate(snapshots):
            for rec in receivers:
                for icomp, comp in enumerate(rec.components):
                    
                    for set in 'synthetics', 'references':
                        data = rec.__dict__[datasource[set,typ]][icomp] 
                                
                        if data is None: continue
                        
                        xdata, ydata = data
                        if typ == 'spectrum':
                            indices = num.where(ydata != 0)[0]
                            if indices.size == 0:
                                continue
                            
                            ibeg = indices[0]
                            iend = indices[-1]+1
                            xdata = xdata[ibeg:iend]
                            ydata = ydata[ibeg:iend]
                        
                        trace = UTrace(
                            station = rec.get_station(),
                            network = rec.get_network(),
                            set = set,
                            snapshot = isnap,
                            channel = config.component_names[comp],
                            distance_deg = rec.distance_deg,
                            azimuth = rec.azimuth,
                            lat = rec.lat,
                            lon = rec.lon,
                            misfit = rec.misfits[icomp],
                            misfit_norm_factor = rec.misfit_norm_factors[icomp],
                            weight = rec.weight,
                            xdata = xdata, 
                            ydata = ydata*rec.weight)
                        
                        traces.append(trace)
                    
        axconfig = {}
        if typ == 'spectrum':
            axconfig['ymode'] = '0-max'
        
        plotter = MyTracy(
            width = 14*gmtpy.cm, 
            height = 18*gmtpy.cm,
            margins = (2.*gmtpy.cm,0.5*gmtpy.cm, 1*gmtpy.cm, 1*gmtpy.cm), 
            source_location = source_location,
            font_size = 8,
            plot_type = plural[typ],
            axconfig = axconfig)
        
        plotter.set_traces(traces)
        
        fns.extend(plotter.save(pjoin(plotdir,'%s_%s.pdf' % (plural[typ], '%i'))))
        fn_all = pjoin(plotdir, '%s_%s.pdf' % (plural[typ], 'all'))
        plotting.pdfjoin_gs(fns, fn_all)
        fns.append(fn_all)
        
        fns_all.extend(fns)
        
    return fns_all
    