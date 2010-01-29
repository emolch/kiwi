

import gmtpy
from gmtpy import cm
import pymseed
import sys
import config
from os.path import join as pjoin

def minmax(traces, minmaxtrace, key=lambda tr: None):
    
    ranges = {}
    for trace in traces:
        mi, ma = minmaxtrace(trace)
        k = key(trace)
        if k not in ranges:
            ranges[k] = mi, ma
        else:
            tmi, tma = ranges[k]
            ranges[k] = min(tmi,mi), max(tma,ma)
    
    return ranges



class Tracy:
    
    def __init__(self,
            width=10*cm,
            height=10*cm*gmtpy.golden_ratio,
            axconfig=None,
            gmtconfig=None):
        
        self.gmtconfig = {'TICK_PEN': '1.25p',
            'TICK_LENGTH': '0.2c',
            'ANNOT_FONT_PRIMARY': '1',
            'ANNOT_FONT_SIZE_PRIMARY': '11p',
            'LABEL_FONT': '1',
            'LABEL_FONT_SIZE': '11p',
            'CHAR_ENCODING': 'ISOLatin1+'}
        
        if gmtconfig is not None:
            self.gmtconfig.update(gmtconfig)
            
        self.axconfig = dict(xmode='min-max', ymode='symmetric')
            
        if axconfig is not None:
            self.axconfig.update(axconfig)
            
        self.width = width
        self.height = height
        
            
    def gather(self, traces, mapping, ordering):
        keys = set()
        for tr in traces:
            keys.add(mapping(tr))
        return sorted(keys, cmp=ordering)

    def map_xgroup(self, trace):
        return None

    def map_ygroup(self, trace):
        return None
    
    def map_zgroup(self, trace):
        return None
    
    def order_xgroup(self, a,b):
        return cmp(a,b)
    
    def order_ygroup(self, a,b):
        return cmp(a,b)
    
    def order_zgroup(self, a,b):
        return cmp(a,b)
    
    def map_xscaling(self, trace):
        '''Default: each subplot has its own x-scaling.'''
        return self.map_xgroup(trace), self.map_ygroup(trace), self.map_zgroup(trace)
        
    def map_yscaling(self, trace):
        '''Default: each subplot has its own y-scaling.'''
        return self.map_xgroup(trace), self.map_ygroup(trace), self.map_zgroup(trace)
    
    def map_color(self, trace):
        return None
    
    def order_color(self, a,b):
        return cmp(a,b)
    
    def colors(self, color_key):
        return gmtpy.color(self.color_index[color_key])
    
    def data(self, trace):
        return trace.xdata, trace.ydata
    
    def xminmax1(self, trace):
        return trace.xdata.min(), trace.xdata.max()
    
    def yminmax1(self, trace):
        return trace.ydata.min(), trace.ydata.max()
    
    def xminmax(self, traces):
        return minmax(traces, self.xminmax1, key=self.map_xscaling)
    
    def yminmax(self, traces):
        return minmax(traces, self.yminmax1, key=self.map_yscaling)
    
    def nwidgets(self):
        return self.nxgroups, self.nygroups
    
    def npages(self):
        return self.nzgroups
    
    def group_to_widget_and_page(self, xgroup, ygroup, zgroup):
        return self.xgroup_index[xgroup], self.ygroup_index[ygroup], self.zgroup_index[zgroup]
    
    def set_traces(self, traces):
        self.traces = traces
        self._update_groups()
        self._update_scalers()
        self._update_widgets()

    def _update_groups(self):
        self.xgroup_keys = self.gather(self.traces, self.map_xgroup, self.order_xgroup)
        self.ygroup_keys = self.gather(self.traces, self.map_ygroup, self.order_ygroup)
        self.zgroup_keys = self.gather(self.traces, self.map_zgroup, self.order_zgroup)
        self.color_keys = self.gather(self.traces, self.map_color, self.order_color)
        self.xgroup_index = dict([ (key, i) for (i,key) in enumerate(self.xgroup_keys) ])
        self.ygroup_index = dict([ (key, i) for (i,key) in enumerate(self.ygroup_keys) ])
        self.zgroup_index = dict([ (key, i) for (i,key) in enumerate(self.zgroup_keys) ])
        self.color_index = dict([ (key, i) for (i,key) in enumerate(self.color_keys) ])
        self.nxgroups = len(self.xgroup_keys)
        self.nygroups = len(self.ygroup_keys)
        self.nzgroups = len(self.zgroup_keys)
        self.ncolors = len(self.color_keys)
    
    def _update_scalers(self):
        xranges = self.xminmax(self.traces)
        yranges = self.yminmax(self.traces)
        
        scalekeys = []
        for trace in self.traces:
            scalekeys.append( (self.map_xscaling(trace), self.map_yscaling(trace)) )
            
        self.scalers = {}
        for xscale_key, yscale_key in scalekeys:
            xr, yr = xranges[xscale_key], yranges[yscale_key]
            axes = [ gmtpy.simpleconf_to_ax(self.axconfig,x) for x in 'xy' ]
            self.scalers[xscale_key, yscale_key] = gmtpy.ScaleGuru([(xr,yr)], axes=axes)
        
    def _update_widgets(self):
        
        nxwidgets, nywidgets = self.nwidgets()
        grid = gmtpy.GridLayout(nxwidgets, nywidgets)
        mw=0.1*cm
        widgets = {}
        for iywidget in range(nywidgets):
            for ixwidget in range(nxwidgets):
                frame = gmtpy.FrameLayout()
                frame.set_fixed_margins(mw*2.,mw*2.,mw/2.,mw/2.)
                frame.get_widget('center').set_vertical(0., 1.)
                grid.set_widget(ixwidget,iywidget, frame)
                widgets[(ixwidget, iywidget)] = frame.get_widget('center')
                
        self.inner_layout = grid
        self.widgets = widgets
        
    def _plot_traces(self, gmt, ipage):
        
        for trace in self.traces:
            ixwidget, iywidget, ipage_ = self.group_to_widget_and_page(
                self.map_xgroup(trace),
                self.map_ygroup(trace),
                self.map_zgroup(trace))
            
            if ipage != ipage_: continue
            
            widget = self.widgets[ixwidget,iywidget]
            scaler = self.scalers[self.map_xscaling(trace),self.map_yscaling(trace)]
            xdata, ydata = self.data(trace)
            gmt.psxy( 
                in_columns=(xdata, ydata),
                N=True, 
                W=self.colors(self.map_color(trace)),
                *(widget.JXY() + scaler.R()) )
            
    def save(self, filename_tmpl):
        self.gmtconfig['PAPER_MEDIA'] = 'Custom_%ix%i' % (self.width,self.height)
        fns = []
        for ipage in range(self.npages()):
            gmt = gmtpy.GMT(config=self.gmtconfig)
            layout = gmt.default_layout()
            layout.set_widget('center', self.inner_layout)
            mw=0.2*cm
            layout.set_fixed_margins(mw,mw,mw,mw)
            self._plot_traces(gmt, ipage)
            fn = filename_tmpl % ipage
            gmt.save( fn )
            fns.append(fn)
            
        return fns
        
class PileTracy(Tracy):
    
    def set_pile(self, pile):
        traces = pile.all()
        self.set_traces(traces)
        
    def xminmax(self, traces):
        return pymseed.minmaxtime(traces, key=self.map_xscaling)
    
    def yminmax(self, traces):
        return pymseed.minmax(traces, key=self.map_yscaling)
        
    def data(self, trace):
        return trace.make_xdata(), trace.ydata()

class MyTracy(Tracy):
    
    def map_ygroup(self, trace):
        return trace.station, trace.network, trace.azimuth, trace.distance_deg
    
    def map_xgroup(self, trace):
        return trace.channel
    
    def order_xgroup(self, a,b):
        return cmp(b,a)
    
    def order_ygroup(self, a,b):
        return cmp((round(a[3]/30.),a[2]),(round(a[3]/30.),b[2]))
    
   # def map_yscaling(self, trace):
   #     return round(trace.distance_deg/10.)
    
    def map_color(self, trace):
        return trace.location


class UTrace:
    def __init__(self, **kwargs):
        for k,v in kwargs.iteritems():
            self.__dict__[k] = v

def multi_seismogram_plot2(snapshots, plotdir):
    
    plural = { 'seismogram': 'seismograms',
                'spectrum': 'spectra' }
                
    datasource = { ('synthetics', 'spectrum'): 'syn_spectra',
                    ('references', 'spectrum'): 'ref_spectra',
                    ('synthetics', 'seismogram'): 'syn_seismograms',
                    ('references', 'seismogram'): 'ref_seismograms' }
                
    fns = []
    for typ in 'seismogram', 'spectrum':
        traces = []
        for receivers in snapshots:
            for rec in receivers:
                for icomp, comp in enumerate(rec.components):
                    
                    for set in 'synthetics', 'references':
                        data = rec.__dict__[datasource[set,typ]][icomp] 
                                
                        if data is None: continue
                        
                        trace = UTrace(
                            station = rec.get_station(),
                            network = rec.get_network(),
                            location = set,
                            channel = config.component_names[comp],
                            distance_deg = rec.distance_deg,
                            azimuth = rec.azimuth,
                            misfit = rec.misfits[icomp],
                            misfit_norm_factor = rec.misfit_norm_factors[icomp],
                            xdata = data[0], ydata = data[1])
                        
                        traces.append(trace)
                    
        plotter = MyTracy(height=40*cm)
        plotter.set_traces(traces)
        
        fns.extend(plotter.save(pjoin(plotdir,'%s_%s.pdf' % (plural[typ], '%i'))))
        
    return fns
    