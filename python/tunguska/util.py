
import os
import subprocess
import tempfile

import pickle
import shutil
import sys
import numpy as num
pjoin = os.path.join

def kiwi_aux_dir():
    ieq_home = os.getenv('KIWI_HOME')
    if ieq_home is None:
         sys.exit('KIWI_HOME environment variable not set')
    d = pjoin(ieq_home, 'aux')
    if not os.path.isdir(d):
        sys.exit('directory not found: "%s"' % d)
    return d

def kiwi_aux_file(*p):
    return pjoin(kiwi_aux_dir(), *p)

def unindent(s):
    lines = s.splitlines()
    if lines and not lines[0].strip(): lines.pop(0)
    if lines and not lines[-1].strip(): lines.pop()
    minindent = min([ len(line) - len(line.lstrip()) for line in lines if line.strip() ])
    eat = ' '*minindent
    outlines = []
    for line in lines:
        if line.startswith(eat):
            outlines.append(line[minindent:])
        else:
            outlines.append(line)
    return '\n'.join(outlines)


def gform( number, significant_digits=3 ):
    '''Pretty print floating point numbers.
    
Align floating point numbers at the decimal dot.
       
      |  -d.dde+xxx|
      |  -d.dde+xx |
      |-ddd.       |
      | -dd.d      |
      |  -d.dd     |
      |  -0.ddd    |
      |  -0.0ddd   |
      |  -0.00ddd  |
      |  -d.dde-xx |
      |  -d.dde-xxx|'''
    
    no_exp_range = (pow(10.,-1), 
                    pow(10.,significant_digits))
    width = significant_digits+significant_digits-1+1+1+5
    
    if (no_exp_range[0] <= abs(number) < no_exp_range[1]) or number == 0.:
        s = ('%#.*g' % (significant_digits, number)).rstrip('0')
    else:
        s = '%.*E' % (significant_digits-1, number)
    s = (' '*(-s.find('.')+(significant_digits+1))+s).ljust(width)
    return s
    
    
def gmt_color( rgb ):
    return '%i/%i/%i' % rgb

def autoplot( *args, **kwargs ):

    outfile = args[-1]
    things = args[:-1]
    infiles = []
    
    if 'O' in kwargs:
        mode = 'a'
    else:
        mode = 'w'
        
    f = open('%s.autoplot' % outfile, mode)
    pickle.dump((args,kwargs), f)
    f.close()
    
    argopts = None
    if 'argopts' in kwargs:
        argopts = kwargs.pop('argopts')
        assert(len(argopts) == len(things)) 
    
    tempdir = tempfile.mkdtemp("","autplot-")
    
    for ithing, thing in enumerate(things):
                
        if type(thing) is str:
            infiles.append(str)
            
        elif isinstance(thing, num.ndarray):
            tfn = pjoin(tempdir, 'plotdata-%i.table' % ithing)
            num.savetxt(tfn, thing)
            infiles.append(tfn)
        
        elif type(thing) is tuple or type(thing) is list:
            
            tfn = pjoin(tempdir, 'plotdata-%i.table' % ithing)
                        
            if thing[-1].ndim == 2 and len(thing) == 3:     # (x[:],y[:],z[:,:])
                x,y,z = thing
                nx = x.size
                ny = y.size
                assert( ny == z.shape[0] and nx == z.shape[1] )
                nz = nx*ny
                tab = num.zeros((nx*ny,3),dtype=num.float)
                tab[:,2] = z.reshape((nz,))
                tab[:,1] = y.repeat(nx)
                for i in range(ny): tab[i*nx:i*nx+nx,0] = x
                num.savetxt(tfn,tab)
            else:
                nrows = thing[0].size
                ncols = len(thing)
                tab = num.zeros((nrows,ncols), dtype=num.float)
                for icol,col in enumerate(thing):
                    assert(col.size == nrows)
                    tab[:,icol] = col
                num.savetxt(tfn,tab)
            infiles.append(tfn)
            
        

    options = []
    for k,v in kwargs.items():
        if type(v) is bool:
            if v:
                options.append( '--%s' % k )
        elif type(v) is tuple or type(v) is list:
            options.append( '--%s=' % k + '/'.join([ str(x) for x in v]) )
        else:
            options.append( '--%s=%s' % (k,str(v)) )
    
    args = [ 'autoplot' ]
    
    infiles_with_argopts = []
    if argopts:
        for infile, argopt in zip(infiles, argopts):
            if argopt:
                infiles_with_argopts.append( infile+'['+argopt+']' )
            else:
                infiles_with_argopts.append( infile )
    else:
        infiles_with_argopts.extend(infiles)
    
    args.extend( infiles_with_argopts )
    args.append( outfile )
    args.extend( options )
    subprocess.call( args )
    
    shutil.rmtree(tempdir)


if __name__ == '__main__':
    import random
    for i in range(0,15):
        n = 10**(i-7)
        print '|%s|' % gform(n)
    for i in range(0,15):
        n = 10**(i-7)*(random.random()*2.-1.)
        print '|%s|' % gform(n)
    print '|%s|' % gform(0.)
    print '|%s|' % gform(0.1)


