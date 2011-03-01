
from distutils.core import setup

from os.path import join as pjoin

setup(
    name    = 'kiwi',
    version = '0.1',
    description = 'Python interface to the Kiwi codes.',
    
    packages = ['tunguska'],
    
    scripts = [ pjoin('scripts', x) for x in [
        'kinherd_sourceview',
        'kinherd_sourceview_pres', 
        'gfdb_phaser',
        'gfdb_downsample',
        'autokiwi', 
        'snufflek',
    ]],
    author = 'Sebastian Heimann',
    author_email = 'sebastian.heimann@zmaw.de',
    url = 'http://www.kinherd.org/',
)

