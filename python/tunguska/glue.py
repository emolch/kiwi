
from pyrocko import eventdata
from pyrocko import model

from tunguska import source
from tunguska import receiver
from tunguska import gfdb
from tunguska import seismosizer

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
    return r
    
def receiver_to_station(rec):
    channels = []
    for comp in rec.components:
        channels.append(kiwi_channels[comp])
        
    sta = model.Station(rec.get_network(), rec.get_station(), rec.get_location(), rec.lat, rec.lon, 0.0, channels=channels)
    return sta
    
def start_seismosizer(  gfdb_path, 
                        effective_dt, 
                        event,
                        stations,
                        hosts=['localhost'], 
                        local_interpolation='bilinear',
                        kiwi_component_map=None ):

    database = gfdb.Gfdb(gfdb_path)
    receivers = []
    for station in stations:
        receivers.append(station_to_receiver(station))
    
    seis = seismosizer.Seismosizer(hosts)
    seis.set_database(database)
    seis.set_effective_dt(effective_dt)
    seis.set_local_interpolation(local_interpolation)
    seis.set_source_location(event.lat, event.lon, event.time)
    seis.set_receivers(receivers)
    return seis
    

