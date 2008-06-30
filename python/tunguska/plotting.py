import util
import config

def km_hack(conf):
    
    if 'xunit' in conf and conf['xunit'] == 'm':
        conf.pop('xunit')
        conf['xfunit'] = 'km'
        conf['xexp'] = 3
        
    if 'yunit' in conf and conf['yunit'] == 'm':
        conf.pop('yunit')
        conf['yfunit'] = 'km'
        conf['yexp'] = 3
        
        
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
    