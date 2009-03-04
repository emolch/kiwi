import util
import config
import copy
import subprocess
import gmt              # <- must get rid of this old version of gmtpy
import gmtpy
import numpy as num
import progressbar
import math
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
    
    
def point_in_region(p,r):
    p = [ num.mod(x,360.) for x in p ]
    r = [ num.mod(x,360.) for x in r ]
    if r[0] <= r[1]:
        blon = r[0] <= p[0] <= r[1] 
    else:
        blon = not (r[1] < p[0] < r[0])
    if r[2] <= r[3]:
        blat = r[2] <= p[1] <= r[3] 
    else:
        blat = not (r[3] < p[1] < r[2])
        
    return blon and blat

    
def location_map( filename, lat, lon, lat_delta, conf_overrides, source=None, source_model_info=None ):

    conf = dict(**config.location_map_config)
    conf.update( conf_overrides )
    
    w = conf.pop('width')
    h = conf.pop('height')
    margins = conf.pop('margins')
    
    topo_img_file_1m = conf.pop('topo_img_file_1m')
    topo_grd_file_5m = conf.pop('topo_grd_file_5m')
    cptfile = conf.pop('topocpt_sealand')
    cptfile1 = conf.pop('topocpt_sea')
    cptfile2 = conf.pop('topocpt_land')    
    
    gmt = gmtpy.GMT( config={ 'ANNOT_FONT_PRIMARY': 'Helvetica',
                              'PLOT_DEGREE_FORMAT':'dddF',
                              'PAPER_MEDIA':'Custom_%ix%i' % (w,h),
                              'GRID_PEN_PRIMARY': 'thinnest/0/0/0' } )
    
    d2r = math.pi/180.0
    
    lat_delta = 5.
    
    south = lat - 0.5*lat_delta
    north = lat + 0.5*lat_delta
    if lat_delta > 20. or south < -80. or north > 80.:
        resolution = 5
        coastline_resolution = 'i'
        rivers=[]
    else:
        resolution = 1
        coastline_resolution = 'f'
        rivers = ['-Ir',]
        
    lon_delta = lat_delta/math.cos(lat*d2r)
    
    delta = lat_delta/360.*config.earthradius*2.*math.pi
    scale_km = gmtpy.nice_value(delta/10.)/1000.
    
    west = lon - 0.5*lon_delta
    east = lon + 0.5*lon_delta
    
    x,y,z = (west, east), (south,north), (-6000.,4500.)
    xax = gmtpy.Ax(mode='min-max', approx_ticks=4.)
    yax = gmtpy.Ax(mode='min-max', approx_ticks=4.)
    zax = gmtpy.Ax(mode='min-max', inc=1000., label='Height', scaled_unit='km', scaled_unit_factor=0.001)
    scaler = gmtpy.ScaleGuru( data_tuples=[(x,y,z)], axes=(xax,yax,zax))    
    
    layout = gmt.default_layout(with_palette=True)
    layout.set_fixed_margins(*margins)
    widget = layout.get_widget().get_widget(0,0)
    palette_layout = gmtpy.FrameLayout()
    palette_layout.set_policy( *layout.get_widget().get_widget(2,0).get_policy() )
    layout.get_widget().set_widget(2,0,palette_layout)
    inset = h-margins[2]-margins[3]
    palette_layout.set_fixed_margins(0.,0.,inset/6., inset/6.)
    palette_widget = palette_layout.get_widget()

    widget['J'] =  ('-JT%g/%g'  % (lon, lat)) + '/%(width)gp'
    aspect = gmtpy.aspect_for_projection( *(widget.J() + scaler.R()) )
    widget.set_aspect(aspect)
    if resolution == 1:
        grdfile = gmt.tempfilename()
        gmt.img2grd( topo_img_file_1m, T=1, S=1, m=1, D=True, G=grdfile, out_discard=True, suppress_defaults=True, *scaler.R())
    else:
        grdfile = topo_grd_file_5m
    
    # work around GMT bug... detect if region contains coastlines
    checkfile = gmt.tempfilename()
    gmt.pscoast( M=True, D=coastline_resolution, W='thinnest/black', A=10., out_filename=checkfile, *(widget.JXY()+scaler.R()))
    f = open(checkfile, 'r')
    has_coastlines = False
    for line in f:
        ls = line.strip()
        if ls.startswith('#') or ls.startswith('>') or ls == '': continue
        clon, clat = [ float(x) for x in ls.split() ]
        
        if point_in_region( (clon,clat), (west,east,south,north) ):
            has_coastlines = True
            break
        else:
            continue
        
    f.close()
    
    if has_coastlines:
        gmt.pscoast( D=coastline_resolution, S='c', A=10., *(widget.JXY()+scaler.R()) )
        gmt.grdimage( grdfile, C=cptfile1, *(widget.JXY()+scaler.R()) )
        gmt.pscoast( Q=True, *(widget.JXY()+scaler.R()))
    
        gmt.pscoast(  D=coastline_resolution, G='c', A=10., *(widget.JXY()+scaler.R()) )
        gmt.grdimage( grdfile, C=cptfile2, *(widget.JXY()+scaler.R()) )
        gmt.pscoast( Q=True, *(widget.JXY()+scaler.R()))
    else:
        gmt.grdimage( grdfile, C=cptfile, *(widget.JXY()+scaler.R()) )
    
    gmt.pscoast( D=coastline_resolution, W='thinnest/black', A=10., *(rivers+widget.JXY()+scaler.R()))
    
    if source:
        gmt.psmeca( in_rows=[[lon, lat, 1., source['strike'], source['dip'], source['slip-rake'], 6., 0.,0., '' ]],
                 S='a0.3', *(widget.JXY()+scaler.R()) )

    if lat > 0:
        axes_layout = 'WSen'
    else:
        axes_layout = 'WseN'
        
    print lon, lat
    gmt.psbasemap( B=('%(xinc)gg%(xinc)g:%(xlabel)s:/%(yinc)gg%(yinc)g:%(ylabel)s:' % scaler.get_params())+axes_layout,
                   L=('x%gp/%gp/%g/%g/%gk' % (widget.width()/2., widget.height()/7.,lon,lat,scale_km) ),
                   *(widget.JXY()+scaler.R()) )

    gmtpy.nice_palette(gmt, palette_widget, scaler, cptfile, innerticks=False, zlabeloffset=1.5*gmtpy.cm)
    
    gmt.save(filename)
    
    
    