import gmtpy
from gmtpy import cm
import copy
import numpy as num

def afloat(l):
    return num.array(l, dtype=num.float)

def minmax(traces, minmaxtrace, key=lambda tr: None, 
                   mincombine=num.amin, maxcombine=num.amax):
    limits = {}
    for trace in traces:
        mi, ma = minmaxtrace(trace)
        k = key(trace)
        if k not in limits:
            limits[k] = [], []
        limits[k][0].append(mi)
        limits[k][1].append(ma)
    
    ranges = {}
    for k,(mins,maxs) in limits.iteritems():
        mins, maxs = afloat(mins), afloat(maxs)
        mi = mincombine(mins)
        ma = maxcombine(maxs)
        ranges[k] = mi, ma
    
    return ranges
    
class Tracy:
    
    def __init__(self,
            width=10*cm,
            height=10*cm*gmtpy.golden_ratio,
            margins=(0.5*cm,0.5*cm,0.5*cm,0.5*cm),
            subplot_margins=(0.2*cm,0.2*cm,0.1*cm,0.1*cm),
            axconfig=None,
            font_size = 10,
            nypaginate=10,
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
        self.margins = margins
        self.nypaginate = nypaginate
        self.nxauxwidgets = 0
        self.font_size = font_size
        self.subplot_margins =subplot_margins
            
    def gather(self, traces, mapping, ordering):
        keys = set()
        for tr in traces:
            keys.add(mapping(tr))
        return sorted(keys, cmp=ordering)

    def label_xgroup(self, a):
        return a
    
    def label_ygroup(self, a):
        return a
    
    def label_zgroup(self, a):
        return a

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
        return self.map_xgroup(trace), self.map_ygroup(trace)
        
    def map_yscaling(self, trace):
        '''Default: each subplot has its own y-scaling.'''
        return self.map_xgroup(trace), self.map_ygroup(trace)
    
    def map_color(self, trace):
        return self.map_zgroup(trace)
    
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
        return self.nxgroups, self.nypaginate
    
    def npages(self):
        return (len(self.ygroup_index)-1) / self.nypaginate + 1
    
    def group_to_widget_and_page(self, xgroup, ygroup, zgroup):
        return (self.xgroup_index[xgroup], 
                self.ygroup_index[ygroup] % self.nypaginate, 
                self.ygroup_index[ygroup] / self.nypaginate)
    
    def set_traces(self, traces):
        self.traces = traces
        self._update_groups()
        self._update_scalers()
        self._update_widgets()
        
    def draw_trace(self, gmt, widget, scaler, trace, ivisit):
        xdata, ydata = self.data(trace)
        gmt.psxy( 
            in_columns=(xdata, ydata),
            N=True, 
            W=self.colors(self.map_color(trace)),
            *(widget.JXY() + scaler.R()) )
            
    def draw_xstuff(self, gmt, widget, xscaler, trace, ivisit):
        pass
    
    def draw_ystuff(self, gmt, widget, yscaler, trace, ivisit):
        pass

    def draw_xaux(self, gmt, widget, trace):
        pass

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
    
    def _scaling_extra(self, scaler):
        
        scaler_x = copy.deepcopy(scaler)
        scaler_x.data_ranges[1] = (0.,1.)
        scaler_x.axes[1].mode = 'off'
        
        scaler_y = copy.deepcopy(scaler)
        scaler_y.data_ranges[0] = (0.,1.)
        scaler_y.axes[0].mode = 'off'
        
        return scaler_x, scaler_y

    def _update_scalers(self):
        xranges = self.xminmax(self.traces)
        yranges = self.yminmax(self.traces)
        
        scalekeys = set()
        for trace in self.traces:
            scalekeys.add( (self.map_xscaling(trace), self.map_yscaling(trace)) )
            
        self.scalers = {}
        for xscale_key, yscale_key in scalekeys:
            xr, yr = xranges[xscale_key], yranges[yscale_key]
            axes = [ gmtpy.simpleconf_to_ax(self.axconfig,x) for x in 'xy' ]
            scaler = gmtpy.ScaleGuru([(xr,yr)], axes=axes)
            xscaler, yscaler = self._scaling_extra(scaler)
            self.scalers[xscale_key, yscale_key] = scaler, xscaler, yscaler
            
    def _update_widgets(self):
        
        nxwidgets, nywidgets = self.nwidgets()
        grid = gmtpy.GridLayout(nxwidgets+self.nxauxwidgets, nywidgets)
        mw=0.1*cm
        widgets = {}
        xauxwigets = {}
        for iywidget in range(nywidgets):
            for ixwidget in range(nxwidgets):
                frame = gmtpy.FrameLayout()
                frame.set_fixed_margins(*self.subplot_margins)
                frame.get_widget('center').set_vertical(0., 1.)
                grid.set_widget(ixwidget,iywidget, frame)
                widgets[(ixwidget, iywidget)] = frame.get_widget('center')
                
            for ixauxwidget in range(self.nxauxwidgets):
                frame = gmtpy.FrameLayout()
                frame.set_fixed_margins(*self.subplot_margins)
                frame.set_horizontal(0.,0.25)
                frame.get_widget('center').set_vertical(0., 1.)
                grid.set_widget(nxwidgets+ixauxwidget, iywidget, frame)
                widget = frame.get_widget('center')
                widget.set_aspect(1.)
                xauxwigets[(ixauxwidget, iywidget)] = widget

        self.inner_layout = grid
        self.widgets = widgets
        self.xauxwidgets = xauxwigets
        
    def _plot_traces(self, gmt, ipage):
        
        ivisits = {}
        for trace in self.traces:
            ixwidget, iywidget, ipage_ = self.group_to_widget_and_page(
                self.map_xgroup(trace),
                self.map_ygroup(trace),
                self.map_zgroup(trace))
            
            if ipage != ipage_: continue
            
            widget = self.widgets[ixwidget,iywidget]
            
            if widget not in ivisits:
                ivisits[widget] = 0
            
            scaler, xscaler, yscaler = self.scalers[self.map_xscaling(trace), self.map_yscaling(trace)]
            
            self.draw_xstuff(gmt, widget, xscaler, trace, ivisits[widget])
            self.draw_ystuff(gmt, widget, yscaler, trace, ivisits[widget])
            self.draw_trace(gmt, widget, scaler, trace, ivisits[widget])
            
            ivisits[widget] += 1
            
    def _plot_labels(self, gmt, ipage):
        
        xhave = {}
        yhave = {}
        
        nxwidgets, nywidgets = self.nwidgets()

        for trace in self.traces:
            ixwidget, iywidget, ipage_ = self.group_to_widget_and_page(
                self.map_xgroup(trace),
                self.map_ygroup(trace),
                self.map_zgroup(trace))
            
            if ipage != ipage_: continue
            
            if iywidget not in yhave:
                left_widget = self.widgets[0,iywidget]
                left_scaler = gmtpy.ScaleGuru([([0, left_widget.width()],[-1,1])])
                
                mleft = self.margins[0]
                smleft = self.subplot_margins[0]
                x,y,size,angle,fontno,justify,text = (
                    -0.5*mleft-smleft, 0.0, self.font_size, 0., 1, 'MC',
                    self.label_ygroup(self.map_ygroup(trace)))
                    
                gmt.pstext(
                    in_rows=[(x,y,size,angle,fontno,justify,text)],
                    N=True,
                    *(left_widget.JXY() + left_scaler.R()) )
                
                yhave[iywidget] = True
                
                for ixauxwidget in range(self.nxauxwidgets):
                    widget = self.xauxwidgets[ixauxwidget, iywidget]
                    self.draw_xaux(gmt, widget, trace)
            
            if ixwidget not in xhave:
                top_widget = self.widgets[ixwidget,0]
                top_scaler = gmtpy.ScaleGuru([([-1,1],[0,top_widget.height()])])
                
                mtop = self.margins[3]
                smtop = self.subplot_margins[3]
                x,y,size,angle,fontno,justify,text = (
                    0.0, top_widget.height()+0.5*mtop+smtop, 10., 0., 1, 'MC',
                    self.label_xgroup(self.map_xgroup(trace)))
                    
                gmt.pstext(
                    in_rows=[(x,y,size,angle,fontno,justify,text)],
                    N=True,
                    *(top_widget.JXY() + top_scaler.R()) )
            
                xhave[ixwidget] = True
            
    def save(self, filename_tmpl):
        self.gmtconfig['PAPER_MEDIA'] = 'Custom_%ix%i' % (self.width,self.height)
        self.gmtconfig['FRAME_PEN'] = '1p/150/150/150'
        self.gmtconfig['FRAME_WIDTH'] = '0.5p'
        fns = []
        for ipage in range(self.npages()):
            gmt = gmtpy.GMT(config=self.gmtconfig)
            layout = gmt.default_layout()
            layout.set_widget('center', self.inner_layout)
            layout.set_fixed_margins(*self.margins)
            self._plot_traces(gmt, ipage)
            self._plot_labels(gmt, ipage)
            fn = filename_tmpl % ipage
            gmt.save( fn )
            fns.append(fn)
            
        return fns
