import util

import sys, os

class Phase:
    def __init__(self,name,filename=None):
    
        self.name = name
        if filename is None:
            filename = os.path.join(util.invearthquake_aux_dir(), 'phases', name)
        
        f = open(filename,'r')
        self.ref_points = []
        for line in f:
            distance, time = [float(x) for x in line.split()]
            self.ref_points.append( (distance, time) )
        f.close()
        
    def __call__(self, distance):
        for (low,high) in zip( self.ref_points[0:-1],self.ref_points[1:len(self.ref_points)]):
            if low[0] <= distance <= high[0]:
                return low[1] + (distance-low[0])/(high[0]-low[0])*(high[1]-low[1])
        return None

if __name__ == '__main__':
    p = Phase('P')
    pn = Phase('Pn')
    for i in range(15):
        print p(i*1000000.), pn(i*1000000.)
