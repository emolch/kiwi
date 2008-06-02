
import os
pjoin = os.path.join

def invearthquake_aux_dir():
    ieq_home = os.getenv('INVEARTHQUAKE_HOME')
    if ieq_home is None:
         sys.exit('INVEARTHQUAKE_HOME environment variable not set')
    d = pjoin(ieq_home, 'aux')
    if not os.path.isdir(d):
        sys.exit('directory not found: "%s"' % d)
    return d


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
    