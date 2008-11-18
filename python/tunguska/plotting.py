import util
import config
import copy
import subprocess
import gmt
import numpy as num
import progressbar
from os.path import join as pjoin

def ra(a):
    return a.min(), a.max()

def all(x):
    for e in x:
        if not e: return False
    return True


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

def km_hack(conf):
    
    if 'xunit' in conf and conf['xunit'] == 'm':
        conf.pop('xunit')
        conf['xfunit'] = 'km'
        conf['xexp'] = 3
        
    if 'yunit' in conf and conf['yunit'] == 'm':
        conf.pop('yunit')
        conf['yfunit'] = 'km'
        conf['yexp'] = 3
        
        
def pdfjoin(files, outfile):
    cmd = ['pdfjoin', '--outfile', outfile ]
    cmd.extend(files)
    subprocess.call(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
def misfit_plot_1d( data, filename, conf_overrides ):
    
    conf = dict(**config.misfit_plot_1d_config)
    symbols = conf.pop('symbols_SGW')
    
    conf.update( conf_overrides )
    conf['argopts'] = [ '%s symbol' % symbol for symbol in symbols[:len(data)] ]
    
    args = tuple(data) + (filename,)
    
    util.autoplot( *args, **conf )
                    
                    
                    
def histogram_plot_1d( data, filename, conf_overrides ):
    
    conf = dict(**config.histogram_plot_1d_config)
    conf.update( conf_overrides )
    
    symbols = conf.pop('symbols_GW')
    
    barwidth = conf['width']/len(data[0][0])*0.8
    
    conf['argopts'] = [ '-Sb%gi %s symbol' % (barwidth, symbol) for symbol in symbols[:len(data)] ]
    args = tuple(data) + (filename,)
    
    util.autoplot( *args, **conf )
                    
                    
                    
def misfit_plot_2d( data, filename, conf_overrides ):
    
    conf = dict(**config.misfit_plot_2d_config)
    conf.update( conf_overrides )
    
    conf['argopts'] = ['density']
    
    args = tuple(data) + (filename,)
    util.autoplot( *args, **conf )
    
def histogram_plot_2d( data, filename, conf_overrides ):
    
    conf = dict(**config.histogram_plot_2d_config)
    conf.update( conf_overrides )
    
    symbols = conf.pop('symbols_SGW')
    
    conf['argopts'] = [ '%s symbol' % symbol for symbol in symbols[:len(data)] ]
    
    args = tuple(data) + (filename,)
    util.autoplot( *args, **conf )
    
    
def misfogram_plot_2d( data, filename, conf_overrides ):
    
    conf = dict(**config.misfogram_plot_2d_config)
    conf.update( conf_overrides )

    symbols = conf.pop('symbols_SGW')
    conf['argopts'] = [ 'density' ]
    conf['argopts'].extend( [ '%s symbol' % symbol for symbol in symbols[:len(data)] ] )
    
    args = tuple(data) + (filename,)
    util.autoplot( *args, **conf )
    
def seismogram_plot( data_by_component, filename, conf_overrides, are_spectra=False ):
    
    if are_spectra:
        conf = dict(**config.spectrum_plot_config)
    else:
        conf = dict(**config.seismogram_plot_config)
        
    conf.update( conf_overrides )
    
    conf_master = conf
    
    ncomps = len(data_by_component)
    for icomp, (comp, data) in enumerate(data_by_component):
        conf = copy.copy(conf_master)
        
        conf["ynlayout"] = ncomps
        conf["yilayout"] = icomp+1
        
        conf['ylabel'] = config.component_names[comp]
        
        axes = 'eW'
        if icomp == 0:
            axes = axes + 'S'
        if icomp == ncomps-1:
            axes = axes + 'n'
        else:
            if 'title' in conf:
                del conf['title']
    
        if icomp > 0:
            conf['O'] = True
        if icomp < ncomps-1:
            conf['K'] = True
            
        symbols = conf.pop('symbols_SGW')
        conf['argopts'] = [ '%s symbol' % symbol for symbol in symbols[:len(data)] ]
        
        conf['axes'] = axes
        
        args = tuple(data) + (filename,)
        util.autoplot( *args, **conf )
        
def multi_seismogram_plot( snapshots, plotdir ):
    
    nrecs = len(snapshots[0])
    for snap in snapshots:
        assert(len(snap) == nrecs)
    
    compos = set()
    for recs in zip(*snapshots):
        for rec in recs:
            compos.update(rec.components)
        
    ordered_compos = []
    for c in 'wesnducalr':
        if c in compos:
            ordered_compos.append(c)
    
    plural = { 'seismogram': 'seismograms',
                'spectrum': 'spectra' }
    allfilez = []
    for typ in 'seismogram', 'spectrum':
    
        if config.show_progress:
            widgets = ['Plotting %s' % plural[typ], ' ',
                    progressbar.Bar(marker='-',left='[',right=']'), ' ',
                    progressbar.Percentage(), ' ',]
            
            pbar = progressbar.ProgressBar(widgets=widgets, maxval=nrecs).start()
        filez = []
        dummy = (num.arange(1), num.arange(1))
        for irec, recs in enumerate(zip(*snapshots)):
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
            
            if None in data_range or None in x_range:
                continue
            
            conf = {}
            proto = recs[0]
            conf['title'] = 'Receiver %i' % (irec+1)
            if all([r.name == proto.name for r in recs]):
                conf['title'] += ': %s' % proto.name
            else:
                conf['title'] += ': %s' % 'comparing different stations'
                
            conf['yrange'] = data_range
            conf['xrange'] = x_range
                
            filename = pjoin(plotdir, '%s_%i.pdf' % (typ,irec+1))
            
            seismogram_plot(data_by_compo, filename, conf_overrides=conf, are_spectra = typ == 'spectrum')
            if config.show_progress: pbar.update(irec+1)
            filez.append(filename)
        
        filename = pjoin(plotdir, '%s_all.pdf' % plural[typ])
        pdfjoin(filez, filename)
        allfilez.extend( filez )
        
        if config.show_progress: pbar.finish()
        
    return allfilez
        
def station_plot( slat, slon, lat, lon, rnames, station_color, station_size, source, maxdist, filename, conf_overrides, zexpand=1.0, nsets=1):
    conf = dict(**config.station_plot_config)
    conf.update( conf_overrides )
    
    width = conf.pop('width')
    height = conf.pop('height')
    margins = conf.pop('margins')
    
    plot = gmt.GMT(width, width, margins)
    plot.pscoast( R='g', J='E%g/%g/%g/%gi' % (slon, slat, maxdist, plot.width), B='', 
                  D='c', A=10000, S=(114,159,207), G=(233,185,110), W='thinnest')
    
    # shift source location slightly, so that projection does not crash...
    plot.psmeca( rows=[[float(slon+0.01), float(slat+0.01), 1., source['strike'], source['dip'], source['slip-rake'], 6., 0.,0., '' ]],
                 S='a0.3' )
                 
    mi = num.nanmin(station_color)
    ma = num.nanmax(station_color)
    ma = max(abs(mi),abs(ma)) * zexpand 
    mi = -ma
    inc = (ma-mi)/128.
    
    f_cpt, fn_cpt = plot.tempfile()
    plot.makecpt( C='polar', T=(mi,ma,inc), Z=True, output=f_cpt )
    
    plot.psxy( columns=(lon[::nsets],lat[::nsets], station_color[::nsets], num.sqrt(station_size[::nsets])/1.5 ),
               C=fn_cpt, S='klflag', W='1p/black', G='white', N=True )
    if nsets == 2:
        plot.psxy( columns=(lon[1::nsets],lat[1::nsets], station_color[1::nsets], num.sqrt(station_size[1::nsets])/1.5 ), 
                   C=fn_cpt, S='krflag', W='1p/black', G='white', N=True )
    
    nr = len(lat)
    size = [9]*nr
    angle = [0]*nr
    fontno = [1]*nr
    justify = ['MC']*nr
    plot.pstext( columns=(lon,lat,size,angle,fontno,justify,rnames), D=(0.,-0.1), N=True )
    
    
    plot.save( filename )
    
    