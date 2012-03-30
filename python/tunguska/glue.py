
from pyrocko import eventdata
from pyrocko import model
from pyrocko import io

from tunguska import source
from tunguska import receiver
from tunguska import gfdb
from tunguska import seismosizer
from tunguska import gridsearch

import numpy as num

def get_nsl(x):
    return x.network, x.station, x.location

kiwi_channels = {
    'u': model.Channel('u', azimuth=0., dip=-90.),
    'd': model.Channel('d', azimuth=0., dip=90.),
    'e': model.Channel('e', azimuth=90., dip=0.),
    'w': model.Channel('w', azimuth=-90., dip=0.),
    'n': model.Channel('n', azimuth=0., dip=0.),
    's': model.Channel('s', azimuth=180., dip=0.),
    'a': model.Channel('a'),
    'c': model.Channel('c'),
    'r': model.Channel('r'),
    'l': model.Channel('l'),
}

kiwi_component_map_default = {
    'N': 'n',
    'S': 's',
    'E': 'e',
    'W': 'w',
    'Z': 'u',
    'R': 'a',
    'T': 'r'
}

def station_to_receiver(station, wanted_components=None, kiwi_component_map=None):
    '''Convert pyrocko-style station into kiwi-style receiver.'''
    
    if kiwi_component_map == 'default':
        kiwi_component_map = kiwi_component_map_default

    components = ''
    for channel in station.get_channels():
        cname = channel.name
        if kiwi_component_map is not None:
            if cname not in kiwi_component_map:
                continue

            cname = kiwi_component_map[cname]
            
        assert len(cname) == 1
        if wanted_components is None or cname in wanted_components:
            components += cname
    
    rname = '%s.%s.%s' % (station.network, station.station, station.location)
    depth = station.depth
    if station.depth is None:
        depth = 0.0
    
    r = receiver.Receiver(station.lat, station.lon, depth=depth, components=components, name=rname)
    return r
    
def receiver_to_station(rec):
    channels = []
    for comp in rec.components:
        channels.append(kiwi_channels[comp])
        
    sta = model.Station(rec.get_network(), rec.get_station(), rec.get_location(), rec.lat, rec.lon, 0.0, depth=rec.depth, channels=channels)
    return sta
    
def start_seismosizer( gfdb_or_path, event, stations=None, receivers=None,
                    local_interpolation     = 'bilinear',
                    spacial_undersampling   = [ 1, 1 ],
                    effective_dt            = 1,
                    crustal_thickness_limit = None,
                    constraining_planes     = None,
                    hosts                   = ['localhost'],
                    balance_method          = '123321',
                    verbose                 = False,
                    gfdb_path               = None, # for backward compatibility
                    ):
        
        assert stations==None or receivers==None
       
        if gfdb_path is not None:
            gfdb_or_path = gfdb_path

        if isinstance(gfdb_or_path, str):
            database = gfdb.Gfdb(gfdb_or_path)
        else:
            database = gfdb_or_path
        
        seis = seismosizer.Seismosizer(hosts, balance_method)
        if verbose: seis.set_verbose('T')
        seis.set_database(database)
        seis.set_effective_dt(effective_dt)
        seis.set_local_interpolation(local_interpolation)
        seis.set_spacial_undersampling(*spacial_undersampling)
        seis.set_source_location(event.lat, event.lon, event.time)
        
        if crustal_thickness_limit is not None:
            seis.set_source_crustal_thickness_limit( crustal_thickness_limit )
        
        if constraining_planes is not None:
            values = []
            for plane in constraining_planes:
                for vect in plane:
                    values.extend(vect)
                    
            seis.set_source_constraints( *values )
        
        if stations is not None:
            receivers = [ station_to_receiver(station) for station in stations ]
        
        seis.set_receivers(receivers)
     
        return seis

    
class EventDataToKiwi:
    '''Interface to convert from pyrocko-style eventdata accessor to kiwi structures'''
    
    def __init__(self, accessor, 
                 station_order      = lambda a,b: cmp(a.dist_deg, b.dist_deg),
                 station_filter     = lambda sta: True,
                 station_splitting  = ['nsewudacrl'],
                 kiwi_component_map = kiwi_component_map_default,
                 trace_factor       = 1.0,
                 trace_time_zero    = ('system','event')[0]
                ):
        
        self._acc = accessor
        self._station_splitting = station_splitting
        self._station_order = station_order
        self._station_filter = station_filter
        self._kiwi_component_map = kiwi_component_map
        self._trace_factor = trace_factor
        self._dataset = None
        self._have_observations = None
        self._update()
        if trace_time_zero == 'system':
            self._zero_time = 0.
        else:
            self._zero_time = self.get_event().time
        
    def make_seismosizer(self, 
                    gfdb_or_path,
                    local_interpolation     = 'bilinear',
                    spacial_undersampling   = [ 1, 1 ],
                    effective_dt            = 1,
                    crustal_thickness_limit = None,
                    constraining_planes     = None,
                    shifts                  = None,
                    blacklist               = None,
                    xblacklist              = None,
                    hosts                   = ['localhost'],
                    balance_method          = '123321',
                    verbose                 = False,
                    gfdb_path               = None, # for backward compatibility
                    ):
       
        if gfdb_path is not None:
            gfdb_or_path = gfdb_path

        if isinstance(gfdb_or_path, str):
            database = gfdb.Gfdb(gfdb_or_path)
        else:
            database = gfdb_or_path

        seis = seismosizer.Seismosizer(hosts, balance_method)
        if verbose: seis.set_verbose('T')
        seis.set_database(database)
        seis.set_effective_dt(effective_dt)
        seis.set_local_interpolation(local_interpolation)
        seis.set_spacial_undersampling(*spacial_undersampling)
        seis.set_source_location(*self.get_source_location())
        
        if crustal_thickness_limit is not None:
            seis.set_source_crustal_thickness_limit( crustal_thickness_limit )
        
        if constraining_planes is not None:
            values = []
            for plane in constraining_planes:
                for vect in plane:
                    values.extend(vect)
                    
            seis.set_source_constraints( *values )
        
        seis.set_receivers(self.get_receivers())
       
        if self._have_observations:
            self._ref_seismogram_stem = seis.tempdir + '/reference-from-eventdata'
            self.put_ref_seismograms()
            seis.set_ref_seismograms( self._ref_seismogram_stem, 'mseed')
        
        if blacklist:
            seis.blacklist_receivers( blacklist )
        if xblacklist:
            seis.xblacklist_receivers( xblacklist )
            
        # apply reference seismograms shifts
        if self._have_observations:
            if shifts is not None:
                seis.shift_ref_seismograms( shifts )
        
        return seis
        
    def make_receiver_weights(self, seis, base_source):
        
        base_source = source.Source( 'circular', 
                                time=base_source['time'],
                                depth=base_source['depth'],
                                radius= 0.,
                                moment=7.0e18,
                                rise_time=base_source['rise-time'])
                                
        zero_source = base_source.clone()
        zero_source['moment'] = 0.0
        seis.set_source(zero_source)
        seis.set_synthetic_reference()
        seis.set_source(base_source)
        
        strike_grid    = ('strike',  -180., 150., 30.)
        dip_grid       = ('dip',  0., 90., 30.)
        slip_rake_grid = ('slip-rake', -180., 150., 30.)
        sdr_grid = gridsearch.MisfitGrid( base_source, [strike_grid, dip_grid, slip_rake_grid])
        
        sdr_grid.compute( seis )
        
        means = sdr_grid.get_mean_misfits_by_r()
        means /= num.mean(means[means>0.])
        dweights = num.where(means>0., 1./means, 0.)
        
        # reset reference seismograms
        seis.set_ref_seismograms( self._ref_seismogram_stem, 'mseed' )
        
        return dweights
        
    def get_event(self):
        return self._acc.get_events()[0]
        
    def get_receivers(self):
        receivers = []
        for (station, receiver, traces) in self._dataset:
            receivers.append(receiver)
        return receivers
    
    def get_source_location(self):
        event = self.get_event()
        return event.lat, event.lon, event.time-self._zero_time
    
    def put_ref_seismograms(self):
        traces_path = self._ref_seismogram_stem+ '-%(ireceiver)i-%(component)s.mseed'
        for irec, (station, receiver, traces) in enumerate(self._dataset):
            for tr in traces:
                fn = traces_path % {
                        'ireceiver': irec+1, 
                        'component': self._kiwi_component_map[tr.channel]}
                        
                tr = tr.copy()
                if self._zero_time != 0.0:
                    tr.shift(-self._zero_time)
                ydata = tr.get_ydata()
                if self._trace_factor != 1.0:
                    ydata *= self._trace_factor
                io.save([tr], fn)

    def get_preprocessed_traces(self):
        all = []
        for trs in self._acc.iter_displacement_traces():
            all.extend(trs)
        return all

    def _update(self):
        
        event = self.get_event()
        
        stations = self._acc.get_stations(relative_event=event).values()
        stations.sort(self._station_order)
        
        kcm = self._kiwi_component_map
       
        traces = self.get_preprocessed_traces() 
        
        # gather traces by station, apply station splitting, and  select desired 
        # components in each set
        
        dataset = []
        for station in stations:
            if not self._station_filter(station):
                continue

            dataset_station = []
            anydata = False
            for set_components in self._station_splitting: 
                
                if traces:
                    station_traces = []
                    trace_ids = {}
                    for tr in traces:
                        if get_nsl(tr) == get_nsl(station) and \
                           tr.channel in kcm and \
                           kcm[tr.channel] in set_components and \
                           tr.nslc_id not in trace_ids:
                           
                            station_traces.append(tr)
                            trace_ids[tr.nslc_id] = True  # cannot handle gappy traces
            
                    station_traces.sort(lambda a,b: cmp(kcm[a.channel], kcm[b.channel]))
                    
                    kiwi_components = ''
                    for tr in station_traces:
                        kiwi_components += kcm[tr.channel]
                        anydata = True
                
                else:
                    kiwi_components = set_components
                    station_traces = []

                receiver = station_to_receiver(station, 
                            wanted_components=kiwi_components, 
                            kiwi_component_map=kcm)
                
                dataset_station.append( (station, receiver, station_traces))
                    
            if anydata or not traces:
                # keep empty datasets if another set from this station has data or if only forward modelling is wanted (no data at all)
                dataset.extend(dataset_station)
                
        self._have_observations = any( [ x[2] for x in dataset ] )
        self._dataset = dataset
        
