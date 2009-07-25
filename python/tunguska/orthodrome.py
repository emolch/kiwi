import math
import numpy as num
import config


d2r = math.pi/180.
r2d = 1./d2r

def clip(x, mi, ma):
    return num.minimum(num.maximum(mi,x),ma)
        
def wrap(x, mi, ma):
    return x - num.floor((x-mi)/(ma-mi)) * (ma-mi)

def ne_to_latlon( lat0, lon0, north_m, east_m ):
    
    '''Transform local carthesian coordinates to latitude and longitude.
    
    lat0, lon0:      Origin of the carthesian coordinate system.
    north_m, east_m: 1D numpy arrays with distances from origin in meters.
    
    Returns: lat, lon: 1D numpy arrays with latitudes and longitudes
    
    The projection used preserves the azimuths of the input points.
    '''

    b = math.pi/2.-lat0*d2r
    a = num.sqrt(north_m**2+east_m**2)/config.earthradius

    gamma = num.arctan2(east_m,north_m)
    alphasign = 1.
    alphasign = num.where(gamma < 0, -1., 1.)
    gamma = num.abs(gamma)
    
    c = num.arccos( clip(num.cos(a)*math.cos(b)+num.sin(a)*math.sin(b)*num.cos(gamma),-1.,1.) )
    alpha = num.arcsin( clip(num.sin(a)*num.sin(gamma)/num.sin(c),-1.,1.) )
    
    alpha = num.where(num.cos(a)-num.cos(b)*num.cos(c) < 0, 
        num.where(alpha > 0,  math.pi-alpha, -math.pi-alpha), alpha)
    
    lat = r2d * (math.pi/2. - c)
    lon = wrap(lon0 + r2d*alpha*alphasign,-180.,180.)
    
    return lat, lon


def ne_to_latlon_alternative_method( lat0, lon0, north_m, east_m ):

    '''Like ne_to_latlon(), but this method, although it should be numerically
    more stable, suffers problems at points which are 'across the pole' as seen
    from the carthesian origin.'''

    b = math.pi/2.-lat0*d2r
    a = num.sqrt(north_m**2+east_m**2)/config.earthradius

    
    gamma = num.arctan2(east_m,north_m)
    alphasign = 1.
    alphasign = num.where(gamma < 0., -1., 1.)
    gamma = num.abs(gamma)
    
    z1 = num.cos((a-b)/2.)*num.cos(gamma/2.)
    n1 = num.cos((a+b)/2.)*num.sin(gamma/2.)
    z2 = num.sin((a-b)/2.)*num.cos(gamma/2.)
    n2 = num.sin((a+b)/2.)*num.sin(gamma/2.)
    t1 = num.arctan2( z1,n1 )
    t2 = num.arctan2( z2,n2 )
    
    alpha = t1 + t2
    beta  = t1 - t2
    
    sin_t1 = num.sin(t1)
    sin_t2 = num.sin(t2)           
    c = num.where( num.abs(sin_t1)>num.abs(sin_t2), 
                num.arccos(z1/sin_t1)*2.,
                num.arcsin(z2/sin_t2)*2. )
            
    lat = r2d * (math.pi/2. - c)
    lon = wrap(lon0 + r2d*alpha*alphasign,-180.,180.)
    return lat, lon


if __name__ == '__main__':
    import sys
    import gmtpy
    import random
    import subprocess
    import time
    
    while True:
        w,h = 20,15
    
        km = 1000.
        gsize = random.uniform(0.,1.)*4.*10.**random.uniform(4.,7.)
        north_grid, east_grid = num.meshgrid(num.linspace(-gsize/2.,gsize/2.,11) ,
                                             num.linspace(-gsize/2.,gsize/2.,11) )
        
        north_grid = north_grid.flatten()
        east_grid = east_grid.flatten()
        
        lat_delta = gsize/config.earthradius*r2d*2.
        lon = random.uniform(-180.,180.)
        lat = random.uniform(-90.,90.)
        
        print gsize/1000.
        
        lat_grid, lon_grid = ne_to_latlon(lat, lon, north_grid, east_grid)
        lat_grid_alt, lon_grid_alt = ne_to_latlon_alternative_method(lat, lon, north_grid, east_grid)
    
    
        maxerrlat = num.max(num.abs(lat_grid-lat_grid_alt))
        maxerrlon = num.max(num.abs(lon_grid-lon_grid_alt))
        eps = 1.0e-8
        if maxerrlon > eps or maxerrlat > eps:
            print lat, lon, maxerrlat, maxerrlon
        
            gmt = gmtpy.GMT( config={ 'PLOT_DEGREE_FORMAT':'ddd.xxxF',
                                    'PAPER_MEDIA':'Custom_%ix%i' % (w*gmtpy.cm,h*gmtpy.cm),
                                    'GRID_PEN_PRIMARY': 'thinnest/0/50/0' } )
        
            south = max(-85., lat - 0.5*lat_delta)
            north = min(85., lat + 0.5*lat_delta)
                
            lon_delta = lat_delta/math.cos(lat*d2r)
            
            delta = lat_delta/360.*config.earthradius*2.*math.pi
            scale_km = gmtpy.nice_value(delta/10.)/1000.
            
            west = lon - 0.5*lon_delta
            east = lon + 0.5*lon_delta
            
            x,y = (west, east), (south,north)
            xax = gmtpy.Ax(mode='min-max', approx_ticks=4.)
            yax = gmtpy.Ax(mode='min-max', approx_ticks=4.)
            scaler = gmtpy.ScaleGuru( data_tuples=[(x,y)], axes=(xax,yax))    
            scaler['R'] = '-Rg'
            layout = gmt.default_layout()
            mw = 2.5*gmtpy.cm
            layout.set_fixed_margins(mw,mw,mw/gmtpy.golden_ratio,mw/gmtpy.golden_ratio)
            widget = layout.get_widget()
            # widget['J'] =  ('-JT%g/%g'  % (lon, lat)) + '/%(width)gp'
            widget['J'] =  ('-JE%g/%g/%g'  % (lon, lat, min(lat_delta/2.,180.))) + '/%(width)gp'
            aspect = gmtpy.aspect_for_projection( *(widget.J() + scaler.R()) )
            widget.set_aspect(aspect)
            
            if lat > 0:
                axes_layout = 'WSen'
            else:
                axes_layout = 'WseN'
        
            gmt.psbasemap( #B=('%(xinc)gg%(xinc)g:%(xlabel)s:/%(yinc)gg%(yinc)g:%(ylabel)s:' % scaler.get_params())+axes_layout,
                           B='5g5',
                        L=('x%gp/%gp/%g/%g/%gk' % (widget.width()/2., widget.height()/7.,lon,lat,scale_km) ),
                        *(widget.JXY()+scaler.R()) )
            
            gmt.psxy( in_columns=(lon_grid,lat_grid), S='x10p', W='1p/200/0/0', *(widget.JXY()+scaler.R()) )
            gmt.psxy( in_columns=(lon_grid_alt,lat_grid_alt), S='c10p', W='1p/0/0/200', *(widget.JXY()+scaler.R()) )
            
            gmt.save('orthodrome.pdf')
            subprocess.call( [ 'xpdf', '-remote', 'ortho', '-reload' ] )
            time.sleep(2)
        else:
            print 'ok', gsize, lat, lon
            