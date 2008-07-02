import util
import config
import copy
import subprocess


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
    
def seismogram_plot( data_by_component, filename, conf_overrides ):
    
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
        
        