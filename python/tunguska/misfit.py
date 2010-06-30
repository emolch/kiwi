class InnerMisfitSetup:
    def __init__(self, inner_norm, tapers_by_set, filter):
        self._inner_norm = inner_norm
        self._tapers_by_set = tapers_by_set
        self._filter = filter
        
    def setup(self, seis, depth):
        tapers = [ self._tapers_by_set[i%len(self._tapers_by_set)] for i in range(len(seis.receivers)) ]
        seis.set_taper(tapers, depth)
        seis.set_filter(self._filter)
        seis.set_misfit_method(self._inner_norm)
        
class OuterMisfitSetup:
    def __init__(self, outer_norm='l1norm', bootstrap_iterations=1000, anarchy=False, receiver_weights=None):
        self._outer_norm = outer_norm
        self._bootstrap_iterations = bootstrap_iterations
        self._anarchy = anarchy
        self._receiver_weights = receiver_weights
        
    def set_receiver_weights(self, receiver_weights):
        self._receiver_weights = receiver_weights
        
    def get_params(self):
        return dict(outer_norm=self._outer_norm,
                    bootstrap_iterations=self._bootstrap_iterations,
                    anarchy=self._anarchy,
                    receiver_weights=self._receiver_weights)
