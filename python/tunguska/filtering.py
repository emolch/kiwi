
class Filter:
    '''Simple filter definition to be used with sesimosizer'''
    
    def __init__(self, frequencies):
        '''Setup cosine frequency tapering filter with given 4 frequencies.'''
        assert(len(frequencies)==4)
        self.frequencies = frequencies
        
    def __call__(self):
        '''Returns representation which can be used for Seismosizer.do_set_misfit_filter().'''
        return ( self.frequencies[0], 0.,
                 self.frequencies[1], 1.,
                 self.frequencies[2], 1.,
                 self.frequencies[3], 0. )
