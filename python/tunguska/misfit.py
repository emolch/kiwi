class InnerMisfitSetup:
    def __init__(self, inner_norm, tapers_by_set=None, filters_by_set=None, taper=None, filter=None):
        self._inner_norm = inner_norm
        self._tapers_by_set = tapers_by_set
        self._filters_by_set = filters_by_set
        self._filter = filter
        self._taper = taper
        
    def setup(self, seis, depth):
        tapers, filters = [], []
        for i in range(len(seis.receivers)):
            taper = self._taper
            if self._tapers_by_set is not None:
                taper = self._tapers_by_set[i%len(self._tapers_by_set)]

            
            filter = self._filter
            if self._filters_by_set is not None:
                filter = self._filters_by_set[i%len(self._filters_by_set)]

            tapers.append(taper)
            filters.append(filter)
            
        seis.set_taper(tapers, depth)
        seis.set_filters(filters)
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
