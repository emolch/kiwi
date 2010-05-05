from pyrocko import pile, trace, eventdata, util, model, pz
import glue, receiver


class KiwiDataAccess(eventdata.EventDataAccess):
    
    def __init__(self, receivers_fn, event_fn, reference_fn_base):
        eventdata.EventDataAccess.__init__(self)
        self._receivers_fn = receivers_fn
        self._event_fn = event_fn
        self._reference_fn_base = reference_fn_base
        self._receivers = None
        self._pile = None
        
    def iter_displacement_traces(self):
        return self.iter_traces()
        
    def get_pile(self):
        if not self._pile:
            fns = []
            for irec, rec in enumerate(self._get_receivers()):
                irec_fortran = irec+1
                for comp in rec.components:
                    fn = '%s-%i-%s.%s' % (self._reference_fn_base, irec_fortran, comp, 'mseed')
            
            self._pile = pile.Pile()
            self._pile.add_files(fns)
            
        return self._pile
        
    def _get_receivers(self):
        if not self._receivers:
            self._receivers = receiver.load_table(self._receivers_fn)
            
        return self._receivers

    def _get_stations_from_file(self):
       
        stations = []
        for rec in self._get_receivers():
            sta = glue.receiver_to_station(rec)
            stations.append(sta)

        return stations
        
    def _get_events_from_file(self):
        ev = model.Event(load=self._event_fn)
        return [ev]
        