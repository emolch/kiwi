
from distutils.core import setup

setup( name    = 'ugly_minimizer',
       version = '0.1',
       description = 'Python interface to the kiwi codes.',
       py_modules = ['ugly_minimizer', 'moment_tensor', 'QVTKRenderWindowInteractor', 'pysacio'],
       packages = ['tunguska'],
       scripts = ['kiwi_prepare', 'kinherd_sourceview', 'kinherd_sourceview_pres', 'gfdb_phaser', 'gfdb_scale', 'kiwi_autopilot'],
       author = 'Sebastian Heimann',
       author_email = 'sebastian.heimann@zmaw.de',
       url = 'http://www.kinherd.org/',
       )

