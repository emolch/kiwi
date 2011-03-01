
from distutils.core import setup

setup( name    = 'kiwi',
       version = '0.1',
       description = 'Python interface to the Kiwi codes.',
       py_modules = ['ugly_minimizer', 'QVTKRenderWindowInteractor', 'pysacio'],
       packages = ['tunguska'],
       scripts = ['kinherd_sourceview',
                  'kinherd_sourceview_pres', 'gfdb_phaser', 'gfdb_scale',
                   'autokiwi'],
       author = 'Sebastian Heimann',
       author_email = 'sebastian.heimann@zmaw.de',
       url = 'http://www.kinherd.org/',
       )

