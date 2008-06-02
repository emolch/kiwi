
from distutils.core import setup

setup( name    = 'ugly_minimizer',
       version = '0.1',
       description = 'Python interface to the invearthquake codes.',
       py_modules = ['ugly_minimizer', 'moment_tensor', 'QVTKRenderWindowInteractor'],
       scripts = ['kinherd_prepare', 'kinherd_sourceview', 'gfdb_phaser'],
       author = 'Sebastian Heimann',
       author_email = 'sebastian.heimann@zmaw.de',
       url = 'http://www.kinherd.org/',
       )

