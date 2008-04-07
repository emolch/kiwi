import os, sys
from subprocess import Popen, PIPE
import ugly_minimizer as minimizer
pjoin = os.path.join


def get_arrival(phases, distance):
    for phase in phases:
        t = phase(distance)
        if t is not None:
            return t

in_db_path = sys.argv[1]
out_db_path = sys.argv[2]
trange = [float(x) for x in sys.argv[3:5]]
phase_names = sys.argv[5:]

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
    

phases = [ minimizer.Phase(phase_name) for phase_name in phase_names]
beg = minimizer.Phase('begin')
end = minimizer.Phase('end')

cmd = [str(x) for x in ['gfdb_redeploy', in_db_path, out_db_path ]]
p = Popen( cmd, stdin=PIPE )
delta = 100000.
for ix in xrange(in_db.nx):
    x = in_db.firstx + ix * in_db.dx
    arrival = get_arrival( phases+[beg], x )
    twin_begin = max(arrival+trange[0], get_arrival((beg,),x))
    twin_end   = min(arrival+trange[1], get_arrival((end,),x))
    print x, twin_begin, twin_end
    for iz in xrange(in_db.nz):
        z = in_db.firstz + iz * in_db.dz
        p.stdin.write( '%g %g %g %g\n' % ( x, z, twin_begin, twin_end ) )
        p.stdin.flush()
        
