#!/usr/bin/env python
from tunguska import gfdb, receiver, seismosizer, source
from pyrocko import orthodrome

from subprocess import call, Popen, PIPE
import numpy as num
import sys

import logging

logformat =  '[%(asctime)s] %(levelname)-8s %(message)s'
dateformat = '%Y-%m-%d %H:%M:%S'

logging.basicConfig( level   = logging.DEBUG,
                         format  = logformat,
                         datefmt = dateformat)

def feed(args, data=''):
    p = Popen(args, stdin=PIPE)
    p.communicate(data)

def dump(fn, s):
    f = open(fn,'w')
    f.write(s)
    f.close()
    

command = sys.argv[1]

if command == 'makedb':

    dump('material.table', 
    '''2300.0  3200.0  1600.0
    ''')

    dump('stf.table', 
    '''0 0
    0.1 0
    0.2 0
    0.3 0
    0.4 0
    0.5 0
    0.6 0.1
    0.7 0.2
    0.8 0.3
    0.9 0.4
    1 0.5
    1.1 0.6
    1.2 0.7
    1.3 0.8
    1.4 0.9
    1.5 1
    1.6 1
    1.7 1
    1.8 1
    1.9 1
    ''')

    feed([str(x) for x in '../gfdb_build', 'benchdb', 1, 
            200, 200, 10, 0.1, 50., 50., 50., 0. ])

    db = gfdb.Gfdb('benchdb')

    p = Popen(['../gfdb_build_ahfull', 'benchdb', 'material.table',
                 'stf.table'], stdin=PIPE)

    for ix in xrange(db.nx):
        print 'distance',  db.firstx + ix*db.dx
        for iz in xrange(db.nz):
            x = db.firstx + ix*db.dx
            z = db.firstz + iz*db.dz
            p.stdin.write("%g %g T T\n" % (x,z))
            p.stdin.flush()
            
    p.stdin.close()
    p.wait()


elif command == 'syntheseis':

    olat, olon = 30., 70.
    receivers = []
    distances = num.linspace(3000, 4000, 2)
    
    for dist in distances:
        lat, lon = orthodrome.ne_to_latlon(olat, olon, dist, 0.)
        r = receiver.Receiver(lat,lon, components='ned')
        receivers.append(r)
        print r
    
    db = gfdb.Gfdb('benchdb')
    
    seis = seismosizer.Seismosizer(hosts=['localhost'])
    seis.set_database(db)
    seis.set_effective_dt(0.1)
    seis.set_local_interpolation('bilinear')
    seis.set_receivers(receivers)
    seis.set_source_location( olat, olon, 0.0)
    
    s = source.Source('bilateral', 
        sourceparams_str='0 0 0 5000 1.0  91 87 164 0 900 700 1000 2500 0.2')
        
    seis.set_source(s)
    seis.output_seismograms('seis', 'table', 'synthetics', 'plain')
    seis.close()
    

