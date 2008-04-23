#!/usr/bin/env python

#
#     Copyright 2007 Sebastian Heimann
#  
#     Licensed under the Apache License, Version 2.0 (the "License");
#     you may not use this file except in compliance with the License.
#     You may obtain a copy of the License at
#  
#         http://www.apache.org/licenses/LICENSE-2.0
#  
#     Unless required by applicable law or agreed to in writing, software
#     distributed under the License is distributed on an "AS IS" BASIS,
#     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#     See the License for the specific language governing permissions and
#     limitations under the License.
# 

# A tool wrapping gfdb_redeploy to redeploy a complete database, or to 
# extract particular phases from a database.

import os, sys
from subprocess import Popen, PIPE
import ugly_minimizer as minimizer
pjoin = os.path.join


def get_arrival(phases, distance):
    for phase in phases:
        t = phase(distance)
        if t is not None:
            return t

# argument processing

if len(sys.argv) < 5:
    sys.exit('usage: gfdb_phaser in_db_path out_db_path [ time_begin time_end phase_names ... ]')

in_db_path = sys.argv[1]
out_db_path = sys.argv[2]

cut_mode = False
if len(sys.argv) >= 3:
    cut_mode = True
    trange = [float(x) for x in sys.argv[3:5]]
    phase_names = sys.argv[5:]


# create database if it does not exist already

in_db = minimizer.get_gfdb_infos( in_db_path )
if not os.path.isfile( out_db_path + '.index' ):
    cmd = [str(x) for x in ['gfdb_build',   out_db_path, 
                                            in_db.nchunks, 
                                            in_db.nx, 
                                            in_db.nz, 
                                            in_db.ng, 
                                            in_db.dt, 
                                            in_db.dx, 
                                            in_db.dz, 
                                            in_db.firstx, 
                                            in_db.firstz ]]
    p = Popen( cmd, stdin=PIPE )
    p.communicate()
    

if cut_mode:
    phases = [ minimizer.Phase(phase_name) for phase_name in phase_names]
    beg = minimizer.Phase('begin')
    end = minimizer.Phase('end')

cmd = [str(x) for x in ['gfdb_redeploy', in_db_path, out_db_path ]]
p = Popen( cmd, stdin=PIPE )

for ix in xrange(in_db.nx):
    x = in_db.firstx + ix * in_db.dx
    if cut_mode:
        arrival = get_arrival( phases+[beg], x )
        twin_begin = max(arrival+trange[0], get_arrival((beg,),x))
        twin_end   = min(arrival+trange[1], get_arrival((end,),x))
        print 'distance: %g  time window: %g %g' % (x, twin_begin, twin_end)
    else:
        print 'distance: %g' % x
        
    for iz in xrange(in_db.nz):
        z = in_db.firstz + iz * in_db.dz
        if cut_mode:
            p.stdin.write( '%g %g %g %g\n' % ( x, z, twin_begin, twin_end ) )
        else:
            p.stdin.write( '%g %g\n' % ( x, z ) )
        
        p.stdin.flush()
        
