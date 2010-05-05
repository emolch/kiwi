
class Filter:
    '''Simple filter definition to be used with sesimosizer'''
    
    def __init__(self, frequencies):
        '''Setup cosine frequency tapering filter with given 4 frequencies.'''
        assert(len(frequencies)==4)
        self.frequencies = frequencies
        
    def set(self, i, f):
        self.frequencies = list(self.frequencies)
        self.frequencies[i] = f
        
    def __call__(self):
        '''Returns representation which can be used for Seismosizer.do_set_misfit_filter().'''
        return ( self.frequencies[0], 0.,
                 self.frequencies[1], 1.,
                 self.frequencies[2], 1.,
                 self.frequencies[3], 0. )


    def __repr__(self):
        return 'Filter( frequencies=(%g, %g, %g, %g) )' % tuple(self.frequencies)
     