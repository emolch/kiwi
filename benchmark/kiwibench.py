#!/usr/bin/env python
from tunguska import gfdb

from subprocess import call, Popen, PIPE

def feed(args, data=''):
    p = Popen(args, stdin=PIPE)
    p.communicate(data)

def dump(fn, s):
    f = open(fn,'w')
    f.write(s)
    f.close()
    
    
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

feed([str(x) for x in '../gfdb_build', '/tmp/testdb', 1, 200, 200, 10, 0.1, 50., 50., 50., 0. ])

db = gfdb.Gfdb('testdb')

p = Popen(['../gfdb_build_ahfull', '/tmp/testdb', 'material.table', 'stf.table'], stdin=PIPE)

for ix in xrange(db.nx):
    print 'distance',  db.firstx + ix*db.dx
    for iz in xrange(db.nz):
        x = db.firstx + ix*db.dx
        z = db.firstz + iz*db.dz
        p.stdin.write("%g %g T T\n" % (x,z))
        p.stdin.flush()
        
p.stdin.close()
p.wait()


