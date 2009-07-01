import pickle
import sys
import logging
from util import gmt_color
import gmtpy
from os.path import join as pjoin
from util import kiwi_aux_file

earthradius = 6371.*1000.

show_progress = True
default_loglevel = logging.INFO

source_info_prog = 'source_info'
gfdb_info_prog = 'gfdb_info'
gfdb_extract_prog = 'gfdb_extract'
sys.stderr.write('warning: config.py:  /bonus/src/kiwi/minimizer')
seismosizer_prog = '/bonus/src/kiwi/minimizer'

topo_dir = '/bonus/topo'
topo_img_file_1m = pjoin(topo_dir,'topo_11.1.img')
topo_grd_file_5m = pjoin(topo_dir,'ETOPO1_Ice_g_gmt4_downsampled_5m.grd')
topocpt_sealand = pjoin(topo_dir, 'light_topo_continuous.cpt')
topocpt_sea = pjoin(topo_dir, 'light_topo_continuous_sea_only.cpt')
topocpt_land = pjoin(topo_dir, 'light_topo_continuous_land_only.cpt')


component_names = {     'a':'R@-+@-',
                        'c':'R@--@-',
                        'r':'T@-+@-',
                        'l':'T@--@-',
                        'd':'Z@-+@-',
                        'u':'Z@--@-',
                        'n':'N',
                        'e':'E',
                        's':'S',
                        'w':'W' }
                        
tango_colors = {
'butter1': (252, 233,  79),
'butter2': (237, 212,   0),
'butter3': (196, 160,   0),
'chameleon1': (138, 226,  52),
'chameleon2': (115, 210,  22),
'chameleon3': ( 78, 154,   6),
'orange1': (252, 175,  62),
'orange2': (245, 121,   0),
'orange3': (206,  92,   0),
'skyblue1': (114, 159, 207),
'skyblue2': ( 52, 101, 164),
'skyblue3': ( 32,  74, 135),
'plum1': (173, 127, 168),
'plum2': (117,  80, 123),
'plum3': ( 92,  53, 102),
'chocolate1': (233, 185, 110),
'chocolate2': (193, 125,  17),
'chocolate3': (143,  89,   2),
'scarletred1': (239,  41,  41),
'scarletred2': (204,   0,   0),
'scarletred3': (164,   0,   0),
'aluminium1': (238, 238, 236),
'aluminium2': (211, 215, 207),
'aluminium3': (186, 189, 182),
'aluminium4': (136, 138, 133),
'aluminium5': ( 85,  87,  83),
'aluminium6': ( 46,  52,  54)
}

graph_colors = [ tango_colors[x] for x in  [ 'skyblue2', 'scarletred2', 'chameleon3', 'chocolate2', 'orange2', 'butter3', 'plum2', 'aluminium4' ]]
plot_background_color = (255,255,255)

symbol_colors = graph_colors
symbol_markers = [ 'c5p', 'd5p', 't5p', 'i5p' ]

def taco(s):
    return gmt_color(tango_colors[s])

symbol_best_result = '-Sa20p  -G%s -W1p/%s' % (taco('scarletred1'), taco('scarletred3'))

#
# Plot configurations using autoplot
#


zebra = [ b or c for c in graph_colors for b in None, (0,0,0) ]
seismogram_plot_config = dict(
    fit = True,
    width    = 4,
    yspacing = 0.1,
    symbols_SGW = ['-W1.5p/%s' % gmt_color(c) for c in zebra ],
    yannotevery = 0,
    yfunit = 'off',
    ylabelofs = 0.15,
    xlabel = 'Time',
    xunit = 's',
    xautoscale = 'min-max',
)

spectrum_plot_config = dict(
    fit = True,
    width    = 4,
    yspacing = 0.1,
    symbols_SGW = ['-W1.5p/%s' % gmt_color(c) for c in zebra ],
    yannotevery = 0,
    yfunit = 'off',
    ylabelofs = 0.15,
    xlabel = 'Frequency',
    xunit = 'Hz',
    xautoscale = 'min-max',
)

# 
# plot configurations using the gmtpy module
#

pw = 15.*gmtpy.cm
mw = 2.5*gmtpy.cm
misfogram_plot_2d_gmtpy_config = dict(
    width = pw,
    height = pw - 1.*gmtpy.cm - 2.*(mw-mw/gmtpy.golden_ratio),
    margins = (mw,mw,mw/gmtpy.golden_ratio,mw/gmtpy.golden_ratio),    
    zapproxticks = 7,
    zlimits = (0.,1.),
    autoscale = 'min-max',
    symbols_SGW = [ symbol_best_result, '-Sc -W1p/black', ],
    zlabel = 'Misfit',
    misfit_cpt = 'ocean',
)

histogram_plot_1d_config = dict(
    width=pw,
    height = pw/gmtpy.golden_ratio,
    margins = (mw,mw,mw/gmtpy.golden_ratio,mw/gmtpy.golden_ratio),
    yexpand = 0.05,
    symbols_SGW = [ '-G%s -W1p/black' % gmt_color(symbol_colors[0]) ],
    xautoscale = 'min-max',
    ylabel = 'Probability of Result'
)

misfit_plot_1d_config = dict(
    yexpand = 0.05,
    width = pw,
    height = pw/gmtpy.golden_ratio,
    margins = (mw,mw,mw/gmtpy.golden_ratio,mw/gmtpy.golden_ratio),
    
    symbols_SGW = ['-S%s -G%s -W1p/black' % (m, gmt_color(c)) for (m,c) in zip(symbol_markers, symbol_colors) ],
    approxticks = 5,
    xautoscale = 'min-max',
    expand = 0.05,
    ylabel =  'Misfit',
)

mw = 0.5*gmtpy.cm
station_plot_config = dict(
    width = pw,
    height = pw,
    margins = (mw,mw,mw,mw),
    topo_img_file_1m = topo_img_file_1m,
    topo_grd_file_5m = topo_grd_file_5m,
    topocpt_sealand = topocpt_sealand,
    topocpt_sea = topocpt_sea,
    topocpt_land = topocpt_land
)

mw = 2.5*gmtpy.cm
location_map_config = dict(
    width = pw,
    height = pw,
    margins = (mw,mw,mw,mw),
    topo_img_file_1m = topo_img_file_1m,
    topo_grd_file_5m = topo_grd_file_5m,
    topocpt_sealand = topocpt_sealand,
    topocpt_sea = topocpt_sea,
    topocpt_land = topocpt_land
)

rupture_plot_config = dict(
    width = pw,
    height = pw/gmtpy.golden_ratio,
    margins = (mw,mw,mw/gmtpy.golden_ratio,mw/gmtpy.golden_ratio),
)

rupture_vis_config = dict(
    symbol_nucleation_point = symbol_best_result.split(),
    rupture_cpt = kiwi_aux_file('cpt', 'rupture.cpt'),    
)

class Config:
    def __init__(self, *configs):
        '''Create configuration opject, based on config files, other config
           objects and dicts. Later entries override earlier.'''
        
        self.configs = configs
    
    def update(self, config):
        self.configs.append( config )
    
    def get_config(self, keys=None):
        
        configdict = {}
        for config in self.configs:
            # either load from file
            if isinstance(config, str):
                try:
                    cd = self.load(config)
                    configdict.update(cd)
                except:
                    sys.exit('failed to load config from file: %s' % config)
                    
            # from another config object / any object which has o.get_config()
            elif hasattr(config, 'get_config'):
                configdict.update(config.get_config())
                
            # or from a mapping type
            else:
                configdict.update(config)
        
        for k in self.__dict__:
            if k not in 'configs':
                configdict[k] = self.__dict__[k]
        if keys == None:
            return configdict
        else:
            sub_configdict = {}
            for k in keys:
                if k in configdict:
                    sub_configdict[k] = configdict[k]
                
            return sub_configdict
    
    def dump(self, filename):
        f = open(filename, 'w')
        pickle.dump(self.get_config(), f)
        f.close()
        
    def load(self, filename):
        f = open(filename, 'r')
        d = pickle.load(f)
        f.close()
        return d
    