
from pyrocko import eventdata

from tunguska import source
from tunguska import receiver
from tunguska import gfdb
from tunguska import seismosizer

def station_to_receiver(station, kiwi_component_map=None):
    '''Convert pyrocko-style station into kiwi-style receiver.'''
    
    components = ''
    for channel in station.channels:
        cname = channel.name
        if kiwi_component_map is not None:
            cname = kiwi_component_map[cname]
        assert len(cname) == 1
        components += cname
        
    rname = '%s.%s.%s' % (station.network, station.station, station.location)
    r = receiver.Receiver(station.lat, station.lon, components, name=rname)
    
def start_seismosizer(  acc,
                        gfdb_path, 
                        effective_dt, 
                        hosts=['localhost'], 
                        local_interpolation='bilinear'
                        kiwi_component_map=None):

    database = gfdb.Gfdb(gfdb_path)

    events = acc.get_events()
    assert len(events) == 1
    event = events[0]
    
    stations = acc.get_stations(relative_event=event)
    receivers = []
    for station in stations.values():
        receivers.append(station_to_receiver(station))
    
    seis = seismosizer.Seismosizer(hosts)
    seis.set_database(database)
    seis.set_effective_dt(effective_dt)
    seis.set_local_interpolation(local_interpolation)
    seis.set_source_location(event.lat, event.lon, event.time)
    seis.set_receivers(receivers)
    return seis
    
    
    

