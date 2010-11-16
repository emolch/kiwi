import util
import config
import copy
import subprocess
import gmtpy
import numpy as num
import progressbar
import math
from os.path import join as pjoin
import orthodrome
import moment_tensor
from moment_tensor import moment_to_magnitude
from plotting_traces import multi_seismogram_plot2


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

def subst(conf, param, oldvalue, newvalue):
    if param in conf and conf[param].lower() == oldvalue:
            conf[param] = newvalue
            
def nukl_hack(conf):
    for x in 'xyz':
        param = x+'label'
        subst(conf,param, 'nukl-shift-x', 'Pos. of Nucleation along Strike')
        subst(conf,param, 'nukl-shift-y', 'Pos. of Nucleation down Dip')
        subst(conf,param, 'bord-radius', 'Border Radius')
        subst(conf,param, 'rel-rupture-velocity', 'Rel. Rupture Velocity')
        
def pdfjoin(files, outfile):
    cmd = ['pdfjoin', '--outfile', outfile ]
    cmd.extend(files)
    subprocess.call(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
def pdfjoin_gs(files, outfile):
    cmd = ['gs', '-dBATCH', '-dNOPAUSE', '-q', '-sDEVICE=pdfwrite', '-sOutputFile=%s' % outfile ]
    cmd.extend(files)
    subprocess.call(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
def misfit_plot_1d( data, filename, conf_overrides, apply_moment_to_magnitude_hack=False ):
    conf = dict(**config.misfit_plot_1d_config)
    conf.update( conf_overrides )
    data = list(data)
    if apply_moment_to_magnitude_hack:
        moment_to_magnitude_hack( conf, data, 'x')

    autoplot_gmtpy( data, filename, **conf )
                    
def histogram_plot_1d( edges, hist, filename, conf_overrides, apply_moment_to_magnitude_hack=False ):
    
    conf = dict(**config.histogram_plot_1d_config)
    conf.update( conf_overrides )
    
    xdata = num.zeros(len(hist)*4,dtype=num.float)
    ydata = num.zeros(len(hist)*4,dtype=num.float)
    xdata[0::4] = edges[:-1]
    xdata[1::4] = edges[:-1]
    xdata[2::4] = edges[1:]
    xdata[3::4] = edges[1:]
    ydata[0::4] = 0.
    ydata[1::4] = hist
    ydata[2::4] = hist
    ydata[3::4] = 0.
   
    data = [(xdata,ydata)]
    if apply_moment_to_magnitude_hack:
        moment_to_magnitude_hack( conf, data, 'x')

    autoplot_gmtpy( data, filename, **conf )
    
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
        
def multi_seismogram_plot( snapshots, plotdir, summary=False):
    
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
    
        if config.show_progress and not summary:
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
            if config.show_progress and not summary: pbar.update(irec+1)
            filez.append(filename)
        
        filename = pjoin(plotdir, '%s_all.pdf' % plural[typ])
        pdfjoin(filez, filename)
        allfilez.extend( filez )
        
        if config.show_progress and not summary: pbar.finish()
        
    return allfilez
    
    
        
def station_plot( slat, slon, lat, lon, rnames, station_color, station_size, source,
                 maxdist, filename, conf_overrides, zexpand=1.0, nsets=1, symbols=('klflag','krflag')):
    conf = dict(**config.station_plot_config)
    conf.update( conf_overrides )
    
    width = conf.pop('width')
    height = conf.pop('height')
    margins = conf.pop('margins')
    
    gmt = gmtpy.GMT( config={ 'LABEL_FONT_SIZE': '12p',
                               'PLOT_DEGREE_FORMAT':'dddF',
                               'PAPER_MEDIA':'Custom_%ix%i' % (width,height),
                               'GRID_PEN_PRIMARY': 'thinnest/0/0/0' } )
    
    layout = gmt.default_layout()
    layout.set_fixed_margins(*margins)
    widget = layout.get_widget()
    widget['J'] = '-JE%g/%g/%g' % (slon, slat, min(maxdist,180.)) + '/%(width)gp'
    aspect = gmtpy.aspect_for_projection( R='g', *widget.J() )
    widget.set_aspect(aspect)
    
    if maxdist > 5.:
        draw_topo(gmt, JXY=widget.JXY(), region='g', resolution=2, coastline_resolution='l', rivers=[], conf=conf, no_clipping=True)

    #gmt.pscoast( R='g', B=True, D='c', A=10000, S=(114,159,207), G=(233,185,110), W='thinnest', *widget.JXY())
    
    # shift source location slightly, so that projection does not crash...
    sdelta = maxdist*0.0001
    
    gmt.psbasemap( R='g', B=True, *widget.JXY() )
    if source and 'strike' in source.keys() and 'dip' in source.keys() and 'slip-rake' in source.keys():
        gmt.psmeca( in_rows=[[float(slon+sdelta), float(slat+sdelta), 1., source['strike'], source['dip'], source['slip-rake'], 6., 0.,0., '' ]],
                        R='g', S='a0.3', *widget.JXY())
    else:
        gmt.psxy( in_rows=[(float(slon+sdelta), float(slat+sdelta))], S='c20p', W='2p/0',  R='g', *widget.JXY())
    
    mi = num.nanmin(station_color)
    ma = num.nanmax(station_color)
    
    colorscale = gmtpy.AutoScaler()
    colorscale.mode = 'symmetric'

    fn_cpt = gmt.tempfilename()
    gmt.makecpt( C='polar', T=colorscale.make_scale((mi,ma)), Z=True, out_filename=fn_cpt )
    
    gmt.psxy( in_columns=(lon[::nsets],lat[::nsets], station_color[::nsets], num.sqrt(station_size[::nsets])/1.5 ),
               R='g', C=fn_cpt, S=symbols[0], W='1p/black', G='white', N=True, *widget.JXY())
    if nsets == 2:
        gmt.psxy( in_columns=(lon[1::nsets],lat[1::nsets], station_color[1::nsets], num.sqrt(station_size[1::nsets])/1.5 ), 
                   R='g', C=fn_cpt, S=symbols[1], W='1p/black', G='white', N=True, *widget.JXY())
    
    nr = len(lat)
    size = [7]*nr
    angle = [0]*nr
    fontno = [1]*nr
    justify = ['MC']*nr
    gmt.pstext( in_columns=(lon,lat,size,angle,fontno,justify,rnames), R='g', D=(0.,-0.1), N=True, *widget.JXY())    
    gmt.save( filename )
    
    
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

def draw_topo(gmt, JXY, region, resolution, coastline_resolution, rivers, conf, no_clipping=False, iluminate=True):
    
    if region != 'g':
        R = [ '-R%g/%g/%g/%g' % region ]
    else:
        R = [ '-Rg' ]
    
    conf = dict(conf)
    
    topo_img_file_1m = conf.pop('topo_img_file_1m')
    topo_grd_file_5m = conf.pop('topo_grd_file_5m')
    cptfile = conf.pop('topocpt_sealand')
    cptfile1 = conf.pop('topocpt_sea')
    cptfile2 = conf.pop('topocpt_land')    

    if 'custom_grd_file' in dir(config) and region != 'g':
        grdfile = config.custom_grd_file
    else:
        if resolution == 1:
            grdfile = gmt.tempfilename()
            gmt.img2grd( topo_img_file_1m, T=1, S=1, m=1, D=True, G=grdfile, out_discard=True, suppress_defaults=True, *R)
        else:
            grdfile = topo_grd_file_5m
    
    # work around GMT bug... detect if region contains coastlines
    if region != 'g':
        checkfile = gmt.tempfilename()
        
        # newer versions of gmt expect -m instead of -M
        if gmtpy.cmp_version(gmt.installation['version'], '4.5.0') < 0:
            mtrue = { 'M':True }
        else:
            mtrue = { 'm':True }
        
        gmt.pscoast( D=coastline_resolution, W='thinnest,black', A=10., out_filename=checkfile, *(JXY+R), **mtrue)
        f = open(checkfile, 'r')
        has_coastlines = False
        for line in f:
            ls = line.strip()
            if ls.startswith('#') or ls.startswith('>') or ls == '': continue
            clon, clat = [ float(x) for x in ls.split() ]
            
            if point_in_region( (clon,clat), region ):
                has_coastlines = True
                break
            else:
                continue
            
        f.close()
    else:
        has_coastlines = True
    
    if iluminate:
        ilumfn = gmt.tempfilename()
        gmt.grdgradient( grdfile, N='t1', A=-45, G=ilumfn, out_discard=True)
        ilumarg = [ '-I%s' % ilumfn ]
    else:
        ilumarg = []
    
    if has_coastlines and not no_clipping:
        gmt.pscoast( D=coastline_resolution, S='c', A=10., *(JXY+R) )
        gmt.grdimage( grdfile, C=cptfile1, *(ilumarg+JXY+R) )
        gmt.pscoast( Q=True, *(JXY+R))
    
        gmt.pscoast(  D=coastline_resolution, G='c', A=10., *(JXY+R) )
        gmt.grdimage( grdfile, C=cptfile2, *(ilumarg+JXY+R) )
        gmt.pscoast( Q=True, *(JXY+R))
    else:
        gmt.grdimage( grdfile, C=cptfile, *(ilumarg+JXY+R) )
    
    gmt.pscoast( D=coastline_resolution, W='thinnest/black', A=10., *(rivers+JXY+R))
    return cptfile
    
def draw_coastlines(gmt, JXY, region, coastline_resolution, rivers):
    
    if region != 'g':
        R = [ '-R%g/%g/%g/%g' % region ]
    else:
        R = [ '-Rg' ]
    
    gmt.pscoast( D=coastline_resolution, W='thinnest/black', A=10., *(rivers+JXY+R))

def draw_shakemap( gmt, widget, scaler, axes, shakemap_range, shakemap_cpt, lat, lon, *datasets):
    
    zax = gmtpy.Ax(mode='0-max', limits=shakemap_range, scaled_unit_factor=100., scaled_unit='cm/s^2', label='Peak Ground Acceleration', masking=False )
    
    zscaler = gmtpy.ScaleGuru([ (lat,lon,dataset) for dataset in datasets ],
        axes=(axes[0],axes[1],zax), 
        aspect=widget.height()/widget.width() )

    grdfile =  gmt.tempfilename()

    R = scaler.R()
    par = scaler.get_params()
    inc_interpol = (0.1*(par['xmax']-par['xmin'])/math.sqrt(len(datasets[0])),
                    0.1*(par['ymax']-par['ymin'])/math.sqrt(len(datasets[0])))
    rxyj = R + widget.XYJ()
        
    colors = gmtpy.color(0), gmtpy.color(1)

    zpar = zscaler.get_params()
    
    clip = (zpar['zinc']/2., zpar['zmax'])
        
    fn_cpt = gmt.tempfilename()
    gmt.makecpt( C=shakemap_cpt, out_filename=fn_cpt,  suppress_defaults=True, *zscaler.T())
    
    if len(datasets[0]) > 3:
        for dataset, color in zip(datasets, colors):
        
            gmt.surface( 
                in_columns=(lon,lat,dataset), 
                T=1,
                G=grdfile, 
                I=inc_interpol, 
                out_discard=True, 
                *R )
        
            gmt.grdimage( 
                grdfile,
                C=fn_cpt,
                #W='2p,%s' % color, 
                #G='d5c', 
                #L=clip,
                *rxyj )
                
    return fn_cpt, zscaler
        
def location_map( filename, lat, lon, lat_delta, conf_overrides, source=None, 
                source_model_infos=None, receivers=None, 
                with_palette=False, shakemap_data=None, shakemap_range=None, shakemap_cpt=None, show_topo=True):

    conf = dict(**config.location_map_config)
    conf.update( conf_overrides )
    
    w = conf.pop('width')
    h = conf.pop('height')
    margins = conf.pop('margins')
        
    d2r = math.pi/180.0
    if source:
        slat, slon = lat, lon
        lat, lon = orthodrome.ne_to_latlon(lat,lon, source['north-shift'], source['east-shift'])
    
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
    
    curcfg = {'TICK_PEN': '1.25p',
                                'TICK_LENGTH': '0.2c',
                                'ANNOT_FONT_PRIMARY': '1',
                                'ANNOT_FONT_SIZE_PRIMARY': '14p',
                                'LABEL_FONT': '1',
                                'LABEL_FONT_SIZE': '14p',
                                'CHAR_ENCODING': 'ISOLatin1+',
                                'D_FORMAT': '%.1f',
                                'PLOT_DEGREE_FORMAT':'DF',
                                'PAPER_MEDIA':'Custom_%ix%i' % (w,h),
                                'GRID_PEN_PRIMARY': 'thinnest/0/0/0' }
    
    #degree_format = 'dddF'
    #if scaler.get_params()['xinc'] < 1. or scaler.get_params()['yinc'] < 1.:
    #    degree_format = 'ddd:mmF'
        
    gmt = gmtpy.GMT( config=curcfg )
    
    if with_palette:
        layout = gmt.default_layout(with_palette=True)
        layout.set_fixed_margins(*margins)
        widget = layout.get_widget().get_widget(0,0)
        palette_layout = gmtpy.FrameLayout()
        palette_layout.set_policy( *layout.get_widget().get_widget(2,0).get_policy() )
        layout.get_widget().set_widget(2,0,palette_layout)
        inset = h-margins[2]-margins[3]
        palette_layout.set_fixed_margins(0.,0.,inset/6., inset/6.)
        palette_widget = palette_layout.get_widget()
    else:
        layout = gmt.default_layout()
        layout.set_fixed_margins(*margins)
        widget = layout.get_widget()
        
    widget['J'] =  ('-JT%g/%g'  % (lon, lat)) + '/%(width)gp'
    aspect = gmtpy.aspect_for_projection( *(widget.J() + scaler.R()) )
    widget.set_aspect(aspect)
    
    if show_topo:
        cptfile = draw_topo(gmt, widget.JXY(), (west,east,south,north), resolution, coastline_resolution, rivers, conf)
    
    zscaler = scaler
    if shakemap_data is not None:
        cptfile, zscaler = draw_shakemap( gmt, widget, scaler, (xax,yax,zax), 
                shakemap_range, shakemap_cpt, *shakemap_data)
        draw_coastlines(gmt, widget.JXY(), (west,east,south,north),  coastline_resolution, rivers)
    
    if source_model_infos and source:
        try:
            draw_rupture(gmt, widget, source_model_infos, axes=(xax,yax,zax), view='from-top', source_loc=(slat,slon), outer_scaler=scaler, with_centroids=False)
        except NoOutlineFound:
            gmt.psxy( in_rows=[[lon,lat]], S='c20p', W='2p,0/0/0',  *(widget.JXY()+scaler.R()))
    
    else:
        if source and 'strike' in source.keys() and 'dip' in source.keys() and 'slip-rake' in source.keys():
            gmt.psmeca( in_rows=[[lon, lat, 1., source['strike'], source['dip'], source['slip-rake'], 6., 0.,0., '' ]],
                     S='a0.3', *(widget.JXY()+scaler.R()) )
        else:
            gmt.psxy( in_rows=[[lon,lat]], S='c20p', W='2p,0/0/0',  *(widget.JXY()+scaler.R()))
    
    if receivers:
        rlats, rlons, rnames = receivers
        gmt.psxy( in_columns=(rlons, rlats), S='t12p', W='1p,black', G='red', *(widget.JXY()+scaler.R()))
        
        nr = len(rnames)
        size = [7]*nr
        angle = [0]*nr
        fontno = [1]*nr
        justify = ['MC']*nr
        gmt.pstext( in_columns=(rlons,rlats,size,angle,fontno,justify,rnames), D=(0.,-0.1), *(widget.JXY()+scaler.R()))
        
    if lat > 0:
        axes_layout = 'WSen'
    else:
        axes_layout = 'WseN'
        
    gmt.psbasemap( B=('%(xinc)gg%(xinc)g:%(xlabel)s:/%(yinc)gg%(yinc)g:%(ylabel)s:' % scaler.get_params())+axes_layout,
                   L=('x%gp/%gp/%g/%g/%gk' % (widget.width()/2., widget.height()/7.,lon,lat,scale_km) ),
                   *(widget.JXY()+scaler.R()) )
                   
    if with_palette and (show_topo or shakemap_data is not None):
        gmtpy.nice_palette(gmt, palette_widget, zscaler, cptfile, innerticks=False, zlabeloffset=1.5*gmtpy.cm)
    gmt.save(filename)
    
class NoOutlineFound(Exception):
    pass
    
def draw_rupture(gmt, widget, source_infos, axes, view='rupture-plane', source_loc=None, outer_scaler=None, with_centroids=True):
    
    conf = dict(**config.rupture_vis_config)
    
    symbol_nucleation_point = conf.pop('symbol_nucleation_point')
    rupture_cpt = conf.pop('rupture_cpt')

    if 'outline' not in source_infos: raise NoOutlineFound()
     
    sec_outline = source_infos['outline']
    
    if (num.all(sec_outline[1][:,0] == sec_outline[1][0,0]) and
        num.all(sec_outline[1][:,1] == sec_outline[1][0,1]) and
        num.all(sec_outline[1][:,2] == sec_outline[1][0,2])):
        raise NoOutlineFound()
    
    if view == 'rupture-plane':
        outline = sec_outline[1][:,3], sec_outline[1][:,4]
    elif view == 'from-top':
        lats, lons = orthodrome.ne_to_latlon(source_loc[0],source_loc[1], sec_outline[1][:,0], sec_outline[1][:,1])
        outline = (lons,lats)
        
    eigrid = source_infos['eikonal-grid'][1].transpose()
    eigrid_valid = num.compress(eigrid[5,:] >= 0., eigrid, axis=1)
    if view == 'rupture-plane':
        ruptime = eigrid_valid[3:6,:]
    elif view == 'from-top':
        lats, lons = orthodrome.ne_to_latlon(source_loc[0],source_loc[1], eigrid_valid[0], eigrid_valid[1])
        ruptime = (lons, lats, eigrid_valid[5])

    eigrid_valid_d = num.compress(eigrid[2,:] > 0., eigrid, axis=1)
    depth = eigrid_valid_d[2,:]
    depth1 = (-(depth-num.min(depth))/(num.max(depth)-num.min(depth))*2.+1.)*0.8
    if view == 'rupture-plane':
        rupdepth = (eigrid_valid_d[3,:],eigrid_valid_d[4,:],depth1)
    elif view == 'from-top':
        lats, lons = orthodrome.ne_to_latlon(source_loc[0],source_loc[1], eigrid_valid_d[0],eigrid_valid_d[1])
        rupdepth = (lons,lats,depth1)
    
    nucleation_point = None
    if 'nucleation-point' in source_infos:
        nucl = source_infos['nucleation-point'][1][0]
        if view == 'rupture-plane':
            nucleation_point = (nucl[3],nucl[4])
        elif view == 'from-top':
            lat, lon = orthodrome.ne_to_latlon(source_loc[0],source_loc[1], nucl[0],nucl[1])
            nucleation_point = (lon,lat)
        
    zax = gmtpy.Ax(snap=True, label='Time', unit='s')

    scaler = gmtpy.ScaleGuru( [ outline, ruptime ], axes=(axes[0],axes[1],zax), aspect=widget.height()/widget.width() )
    grdfile = gmt.tempfilename()
    grdfile_extra = gmt.tempfilename()
    igrdfile = gmt.tempfilename()
    cptfile = gmt.tempfilename()
    
    if outer_scaler is not None:
        R = outer_scaler.R()
    else:
        R = scaler.R()
    
    if outer_scaler is None: outer_scaler = scaler
    par = outer_scaler.get_params()
    inc_interpol = ((par['xmax']-par['xmin'])/(widget.width()/gmtpy.inch*150.),
                    (par['ymax']-par['ymin'])/(widget.height()/gmtpy.inch*150.))
    inc_interpol_extra = ((par['xmax']-par['xmin'])/(widget.width()/gmtpy.inch*20.),
                          (par['ymax']-par['ymin'])/(widget.height()/gmtpy.inch*20.))
    
    rxyj = R + widget.XYJ()
        
    dark_color = '%g/%g/%g' % gmtpy.tango_colors['aluminium6']
    gmt.makecpt( I=True, C=rupture_cpt, Z=True, out_filename=cptfile, *scaler.T())

    if len(ruptime[0]) > 3:
        gmt.triangulate(in_columns=ruptime, G=grdfile, I=inc_interpol, out_discard=True, *R )
        gmt.surface(in_columns=ruptime, G=grdfile_extra, I=inc_interpol_extra, T=0.3, out_discard=True, *R )
        #gmt.triangulate(in_columns=ruptime, G=grdfile_extra, I=inc_interpol_extra, out_discard=True, *R )
        gmt.surface(in_columns=rupdepth, G=igrdfile, I=inc_interpol_extra, T=0.3, out_discard=True, *R )
        gmt.psclip( in_columns=outline, *rxyj )
        gmt.grdimage( grdfile_extra, I=igrdfile, C=cptfile, *rxyj)
        if with_centroids:
            gmt.psxy( in_columns=(ruptime[0],ruptime[1]), S='c3p', G=dark_color, *rxyj)
        gmt.grdcontour( grdfile, W='1p,%s' % dark_color, G='d5c', A='%g+g%s+kwhite+us+o' % (scaler.get_params()['zinc'], dark_color), *rxyj )
        gmt.psclip( C=True, *widget.XY() )
        gmt.psxy( in_columns=outline,  L=True, W='1p,%s' % dark_color, *rxyj )
        
    if nucleation_point:
        gmt.psxy( in_rows=[nucleation_point], *(symbol_nucleation_point+rxyj))
    return scaler, cptfile
    
def rupture_plot(filename, source_infos, conf_overrides=None, gmtconfig=None):
    conf = dict(**config.rupture_plot_config)
    if conf_overrides is not None:
        conf.update( conf_overrides )
    w = conf.pop('width')
    h = conf.pop('height')
    margins = conf.pop('margins')
    
    if gmtconfig is None:
        gmtconfig = {'LABEL_FONT_SIZE': '12p'}
    
    gmtconfig['PAPER_MEDIA'] = 'Custom_%ix%i' % (w,h)
    
    gmt = gmtpy.GMT( config=gmtconfig ) 
    
    layout = gmt.default_layout(with_palette=True)
    layout.set_fixed_margins(*margins)

    widget = layout.get_widget().get_widget(0,0)
    palette_widget = layout.get_widget().get_widget(2,0)
    
    xax, yax, zax = [ gmtpy_ax_from_autoplot_conf(conf,x) for x in ('x','y','z') ]    
    widget['J'] = '-JX%(width)gp/-%(height)gp'
    
    try:
        scaler, cptfile = draw_rupture(gmt, widget, source_infos, axes=(xax,yax,zax))
        gmt.psbasemap( *(widget.JXY() + scaler.RB(ax_projection=True)) )
        gmtpy.nice_palette( gmt, palette_widget, scaler, cptfile, zlabeloffset=1.5*gmtpy.cm )
    except NoOutlineFound:
        scaler = gmtpy.ScaleGuru([([0],[0])], axes=(xax,yax))
        gmt.pstext( in_rows=[(0,0,12,0,0,'MC','No rupture data!')],  *(scaler.R() + widget.XYJ()) )
    
    gmt.save(filename)
    
def gmtpy_ax_from_autoplot_conf(conf, axname):
    c = {}
    x = axname
    for x in ('', axname):
        if x+'limits' in conf:      c['limits'] = conf[x+'limits']
        if x+'label' in conf:       c['label'] = conf[x+'label']
        if x+'unit' in conf:        c['unit'] = conf[x+'unit']
        if x+'funit' in conf:       c['scaled_unit'] = conf[x+'funit']
        if x+'exp' in conf:         c['scaled_unit_factor'] = 10.**(-conf[x+'exp'])
        if x+'expand' in conf:      c['space'] = conf[x+'expand']
        if x+'autoscale' in conf:   c['mode'] = conf[x+'autoscale']
        if x+'approxticks' in conf: c['approx_ticks'] = conf[x+'approxticks']
        if x+'snap' in conf:        c['snap'] = conf[x+'snap']
            
    return gmtpy.Ax( **c )

def to_01(x):
    xo = num.min(x)
    xs = num.max(x)-xo
    if xs == 0.: xs = 1
    return (x-xo)/xs, xo, xs
    
def moment_to_magnitude_hack( conf, data, xy ):
    icomp = {'x':0,'y':1}[xy]

    if conf[xy+'label'].lower() == 'moment':
        for i,dset in enumerate(data):
            data[i] = list(dset)
            data[i][icomp] = moment_to_magnitude(num.asarray(dset[icomp]))
            
        conf[xy+'label'] = 'M@-W@-'
        if xy+'unit' in conf: del conf[xy+'unit']
        if xy+'inc' in conf:  del conf[xy+'inc']
        if xy+'limits' in conf:
            conf[xy+'limits'] = [ moment_to_magnitude(x) for x in conf[xy+'limits'] ]

def misfogram_plot_2d_gmtpy(data, filename, conf_overrides, apply_moment_to_magnitude_hack=False):
    
    conf = dict(**config.misfogram_plot_2d_gmtpy_config)
    if conf_overrides is not None:
        conf.update( conf_overrides )
        
    if apply_moment_to_magnitude_hack:
        data = list(data)
        moment_to_magnitude_hack( conf, data, 'x')
        moment_to_magnitude_hack( conf, data, 'y')            
                
    w = conf.pop('width')
    h = conf.pop('height')
    margins = conf.pop('margins')
    symbols_SGW = conf.pop('symbols_SGW')
    misfit_cpt = conf.pop('misfit_cpt')
    zlabeloffset = conf.pop('zlabeloffset')
    
    gmtconfig = {}
    for k in conf.keys():
        if k.upper() == k:
            gmtconfig[k] = conf[k]
    
    gmtconfig['PAPER_MEDIA'] = 'Custom_%ix%i' % (w,h)
    
    gmt = gmtpy.GMT( config=gmtconfig ) 
    
    layout = gmt.default_layout(with_palette=True)
    layout.set_fixed_margins(*margins)

    widget = layout.get_widget().get_widget(0,0)
    palette_widget = layout.get_widget().get_widget(2,0)
    
    axes = [ gmtpy_ax_from_autoplot_conf(conf,x) for x in ('x','y','z') ]
    scaler = gmtpy.ScaleGuru( data[0:2], axes=axes )
    
    grdfile = gmt.tempfilename()
    cptfile = gmt.tempfilename()
    
    par = scaler.get_params()
    inc_interpol = ((par['xmax']-par['xmin'])/(widget.width()/gmtpy.inch*150.),
                    (par['ymax']-par['ymin'])/(widget.height()/gmtpy.inch*150.))
    inc_interpol = (1./(widget.width()/gmtpy.inch*150.),
                    1./(widget.height()/gmtpy.inch*150.))
    
    r = scaler.R()
    rxyj = scaler.R() + widget.JXY()
    x,y,z = data[0]
    x, xo, xs = to_01(x)
    y, yo, ys = to_01(y)
    xyz_sorted = sorted([ xyz_ for xyz_ in zip(x,y,z) ])
    xyz_sorted.sort()
        
    gmt.triangulate(in_rows=xyz_sorted, G=grdfile, I=inc_interpol, out_discard=True, R=(0,1,0,1))#*r)
    gmt.grdedit( grdfile, R=(xo,xo+xs,yo,yo+ys), out_discard=True)
    gmt.makecpt( I=True, C=misfit_cpt, Z=True, out_filename=cptfile, *scaler.T())
    gmt.grdimage( grdfile, C=cptfile, *rxyj)
    gmt.grdcontour( grdfile, W='1p', C=cptfile, *rxyj )
    gmtpy.nice_palette( gmt, palette_widget, scaler, cptfile, zlabeloffset=zlabeloffset)
    gmt.psxy( in_columns=data[2], *(symbols_SGW[1].split()+rxyj) )
    gmt.psxy( in_columns=data[1], *(symbols_SGW[0].split()+rxyj) )
    gmt.psbasemap( *(widget.JXY() + scaler.RB(ax_projection=True)) )

    gmt.save(filename)
    
def autoplot_gmtpy(data, filename, **conf ):
    
    w = conf.pop('width')
    h = conf.pop('height')
    margins = conf.pop('margins')
    symbols_SGW = conf.pop('symbols_SGW')
    
    gmtconfig = {}
    for k in conf.keys():
        if k.upper() == k:
            gmtconfig[k] = conf[k]
    
    gmtconfig['PAPER_MEDIA'] = 'Custom_%ix%i' % (w,h)
    
    gmt = gmtpy.GMT( config=gmtconfig ) 
    
    layout = gmt.default_layout()
    layout.set_fixed_margins(*margins)

    widget = layout.get_widget()
    
    axes = [ gmtpy_ax_from_autoplot_conf(conf,x) for x in ('x','y') ]
    scaler = gmtpy.ScaleGuru( data[0:2], axes=axes )
    
    gmt.psbasemap( *(widget.JXY() + scaler.RB(ax_projection=True)) )
    rxyj = scaler.R() + widget.JXY()
    
    for dat, sym in zip(data,symbols_SGW):
        bdata = num.ascontiguousarray(num.asarray( dat ).transpose(), dtype=num.float)
        gmt.psxy( in_string=bdata.tostring(), b='i%i' % len(dat), *(sym.split()+rxyj) )
    
    gmt.save(filename)


def beachball(source, filename, **conf_overrides):
    
    if 'strike' in source.keys() and 'dip' in source.keys() and 'slip-rake' in source.keys():
       
        r2d = 180./math.pi
        d2r = 1./r2d

        strike, dip, rake = source['strike'], source['dip'], source['slip-rake']        
        dip, strike, mrake = [r2d*xx for xx in moment_tensor.unique_euler(dip*d2r, strike*d2r, -rake*d2r)]
        rake = -mrake
                
        return beachball_mt(None, filename,  sdr=(strike,dip,rake), **conf_overrides)
    
    else:
        return []
        
def beachball_mt(mt, filename, sdr=None, **conf_overrides):
    
    conf = dict(**config.beachball_config)
    if conf_overrides is not None:
        conf.update( conf_overrides )
    
    w = conf.pop('width')
    h = conf.pop('height')
    margins = conf.pop('margins')
    
    indicate_plane = 0
    if 'indicate_plane' in conf:
        indicate_plane = conf.pop('indicate_plane')
    
    gmtconfig = {}
    for k in conf.keys():
        if k.upper() == k:
            gmtconfig[k] = conf[k]
    
    gmtconfig['PAPER_MEDIA'] = 'Custom_%ix%i' % (w,h)
    gmtconfig['PS_MITER_LIMIT'] = '180'
    
    gmt = gmtpy.GMT( config=gmtconfig ) 
    
    layout = gmt.default_layout()
    layout.set_fixed_margins(*margins)

    widget = layout.get_widget()
    
    if mt is not None:
        strike, dip, rake = mt.both_strike_dip_rake()[0]
    else:
        strike, dip, rake = sdr
    
    kwargs = dict( R=(-1.,1.,-1.,1.), S='a%gc' % (((w-margins[0]-margins[1])/gmtpy.cm)-2./28.34), in_rows=[[
                        0., 0., 1.,
                        strike, dip, rake, 
                        5., 0.,0., '' ]]) 
    
    #gmt.psmeca( *widget.JXY(), **kwargs )
    if indicate_plane == 0:
        gmt.psmeca( L='2p,black', G=conf['fillcolor'], *widget.JXY(), **kwargs )
    else:
        gmt.psmeca( L='1p,black,.', G=conf['fillcolor'], *widget.JXY(), **kwargs )        
        gmt.psmeca( T='%i/2p,%s' % (indicate_plane, gmtpy.color('black')), *widget.JXY(), **kwargs )
                
    gmt.save(filename)
    return filename
    
    