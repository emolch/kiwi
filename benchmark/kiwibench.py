#!/usr/bin/env python
from tunguska import gfdb, receiver, seismosizer, source, config
from pyrocko import orthodrome

from subprocess import call, Popen, PIPE
import numpy as num
import sys, time, signal

import logging
from os.path import join as pj

logformat =  '[%(asctime)s] %(levelname)-8s %(message)s'
dateformat = '%Y-%m-%d %H:%M:%S'

logging.basicConfig( level   = logging.WARN,
                         format  = logformat,
                         datefmt = dateformat)

def feed(args, data=''):
    p = Popen(args, stdin=PIPE)
    p.communicate(data)

def dump(fn, s):
    f = open(fn,'w')
    f.write(s)
    f.close()
    
    
    
bindir = sys.argv[1]
command = sys.argv[2]



def pbin(bin):
    return pj(bindir, bin)


config.source_info_prog = pbin('source_info')
config.gfdb_info_prog = pbin('gfdb_info')
config.gfdb_extract_prog = pbin('gfdb_extract')
config.gfdb_build_prog = pbin('gfdb_build')
config.seismosizer_prog =pbin('minimizer')

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

    feed([str(x) for x in pbin('gfdb_build'), 'benchdb', 1, 
            200, 200, 10, 0.1, 50., 50., 50., 0. ])

    db = gfdb.Gfdb('benchdb')

    p = Popen([ pbin('gfdb_build_ahfull'), 'benchdb', 'material.table',
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
    distances = num.linspace(3000, 4000, 100)
    
    for dist in distances:
        lat, lon = orthodrome.ne_to_latlon(olat, olon, dist, 0.)
        r = receiver.Receiver(lat,lon, components='ned')
        receivers.append(r)
    
    db = gfdb.Gfdb('benchdb')
    
    seis = seismosizer.Seismosizer(hosts=['localhost']*4)
    seis.set_database(db)
    seis.set_effective_dt(0.1)
    seis.set_local_interpolation('bilinear')
    seis.set_receivers(receivers)
    seis.set_source_location( olat, olon, 0.0)
    
   # s = source.Source('bilateral', 
   #     sourceparams_str='0 0 0 5000 1e12  91 87 164 0 0 0 0 2500 0.2')
    
    s = source.Source('bilateral', 
        sourceparams_str='0 0 0 5000 1e12  91 87 164 0 600 600 600 2500 0.2')
    
    seis.set_source(s)
    seis.set_synthetic_reference()
    seis.set_ignore_sigint('T')
    
    seis.output_seismograms( 'seis', 'mseed', 'synthetics', 'plain')
    
    done = []
    def got_sigint(*args):
        done.append(1)
        
    signal.signal(signal.SIGINT, got_sigint)
    
    start = time.time()
    times = [ start ]
    strikes = num.linspace(0.,360.,361)
    for istrike, strike in enumerate(strikes):
        s['strike'] = strike
        seis.set_source(s)
        seis.do_get_misfits()
        now = time.time()
        total_mps = float(istrike+1) / (now - start)
        last_mps = 1. / (now - times[-1])
        last_10_mps = min(10, len(times)) / (now-times[-min(10,len(times))])
        
        print 'MPS   total avg: %.3g   last 10 avg: %.3g   last: %.3g' % (total_mps, last_10_mps, last_mps)
        
        times.append(now)
        if done:
            break
        
    seis.close()
    

