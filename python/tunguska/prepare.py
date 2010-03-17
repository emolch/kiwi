from pyrocko import util, io, trace, moment_tensor, rdseed

from tunguska import phase, gfdb, edump_access

import sys, os, logging, shutil, time, copy
from os.path import join as pjoin

logger = logging.getLogger('tunguska.prepare')

def get_nsl(x):
    return x.network, x.station, x.location

def get_ns(x):
    return x.network, x.station

def save_rapid_station_table(stations_path, stations):
    '''Save station table in format for rapidinv'''
    util.ensuredirs(stations_path)
    f = open(stations_path, 'w')
    for station in stations:
        nsl = '.'.join((station.network,station.station,station.location))
        f.write('%-10s %15.8e %15.8e\n' % (nsl, station.lat, station.lon))
    f.close()
    
def save_rapid_input_file(rapid_input_path, event, gfdb):
    mag = event.magnitude
    if mag is None:
        mag = 6.0
    
    depth = event.depth
    if depth is None:
        depth = 20000.
    
    if gfdb is not None:
        depth_range = (gfdb.firstz, gfdb.firstz+gfdb.dz*(gfdb.nz-1))
    else:
        depth_range = (1000.,600000.)
    
    x = dict(
        lat=event.lat,
        lon=event.lon,
        eventname=event.name,
        eventmom=moment_tensor.magnitude_to_moment(mag),
        depth=depth,
        depth_upperlim=depth_range[0]/1000.+1,
        depth_bottomlim=depth_range[1]/1000.-1)
        
    f = open(rapid_input_path, 'w')
    f.write('''#FileSystem
INVERSION_DIR    results
DATA_DIR         data
DATA_FILE        %(eventname)s

#SourceStartingValues
LATITUDE_NORTH   %(lat)e
LONGITUDE_EAST   %(lon)e
DEPTH_1          %(depth)e
DEPTH_2          %(depth)e
DEPTH_UPPERLIM   %(depth_upperlim)e
DEPTH_BOTTOMLIM  %(depth_bottomlim)e
SCAL_MOM_1       %(eventmom)e
SCAL_MOM_2       %(eventmom)e
SCAL_MOM_STEP    1e16

#SpecificFocalMecInversion
BP_F1_STEP1      0.01
BP_F2_STEP1      0.01
BP_F3_STEP1      0.03
BP_F4_STEP1      0.03

#SpecificCentroidInversion
REL_TIME_1       0
REL_TIME_2       16
REL_TIME_STEP    2
BP_F1_STEP2      0.01
BP_F2_STEP2      0.01
BP_F3_STEP2      0.03
BP_F4_STEP2      0.03

#SpecificKinematicInversion
EFFECTIVE_DT_ST3 2.0
KIN_RISETIME     2
BP_F1_STEP3      0.01
BP_F2_STEP3      0.01
BP_F3_STEP3      0.10
BP_F4_STEP3      0.10
''' % x )
    f.close()
    
    
def save_kiwi_dataset( stations, traces, event, config):
    if hasattr(config, 'skeleton_dir'): 
        copy_files(config.path('skeleton_dir'), config.path('main_dir'))

    # gather traces by station
    dataset = []
    for station in stations:
        station_traces = []
        for tr in traces:
            if get_nsl(tr) == get_nsl(station):
                if tr.channel in config.wanted_components:
                    station_traces.append(tr)
                        
        station_traces.sort(lambda a,b: cmp(a.channel, b.channel))
        kiwi_components = ''
        for tr in station_traces:
            kiwi_components += config.kiwi_component_map[tr.channel]
        if station_traces:
            dataset.append( (station, kiwi_components, station_traces) )
    
    fpath = config.path('receivers_path')
    util.ensuredirs(fpath)
    f = open(fpath, 'w')
    iref = 1
    for station, components, traces in dataset:
        nsl = '.'.join((get_nsl(station)))
        for i in range(config.nsets):
            f.write('%15.8e %15.8e %3s %-15s\n' % (station.lat, station.lon, components, nsl) )
            for tr in traces:
                tr = tr.copy()
                if config.trace_time_zero == 'event':
                    tr.shift(-event.time)
                ydata = tr.get_ydata()
                ydata *= config.trace_factor
                fn = config.path('trace_path', {
                    'ireceiver': iref, 
                    'component': config.kiwi_component_map[tr.channel]})
                io.save([tr], fn)
                
            iref += 1
    
    f.close()
    
    fpath = config.path('reference_time_path')
    f = open(fpath, 'w')
    f.write('%i %s\n' % (event.time, 
                    time.strftime('%Y/%m/%d %H:%M:%S', 
                                    time.gmtime(event.time))))
    f.close()
    
    fpath = config.path('source_origin_path')
    f = open(fpath, 'w')
    f.write('%e %e 0\n' % (event.lat, event.lon))
    f.close()
    
    
def save_rapid_dataset(stations, traces, event, gfdb, config):
    
    if hasattr(config, 'skeleton_dir'): 
        copy_files(config.path('skeleton_dir'), config.path('main_dir'))
    
    save_rapid_station_table(config.path('stations_path'), stations)
    save_rapid_input_file(config.path('rapid_input_path'), event, gfdb)
    
    used_traces = []
    for station in stations:
        for tr in traces:
            if get_nsl(tr) == get_nsl(station):
                tr = tr.copy()
                if config.trace_time_zero == 'event':
                    tr.shift(-event.time)
                ydata = tr.get_ydata() 
                ydata *= config.trace_factor
                used_traces.append(tr)
    
    io.save(used_traces, config.path('trace_path'))

    
def copy_files(source_dir, dest_dir):
    if os.path.exists(dest_dir):
        shutil.rmtree(dest_dir)
    shutil.copytree(source_dir, dest_dir)

def prepare(config, kiwi_config, rapid_config, event_names):

    if hasattr(config, 'gfdb_path'):
        db = gfdb.Gfdb(config.path('gfdb_path'))
        deltat = db.dt
        min_dist = db.firstx + config.gfdb_margin
        max_dist = db.firstx + (db.nx-1)*db.dx - config.gfdb_margin
    else:
        db = None
        deltat = config.deltat
        min_dist = None
        max_dist = None
    
    
    
    for event_name in event_names:
        
        logger.info('Preparing event %s' % event_name)
        
        config.event_name = event_name
        
        sw = util.Stopwatch()
        if hasattr(config, 'seed_volume'):
            try:
                acc = rdseed.SeedVolumeAccess(config.path('seed_volume'))
            except rdseed.SeedVolumeNotFound:
                logger.error('SEED volume not found for event %s' % event_name)
                continue
               
        elif hasattr(config, 'edump_data_dir'):
            acc = edump_access.EventDumpAccess(config.path('edump_data_dir'))
        else:
            sys.exit('config has neither entry "seed_volume" nor "edump_data_dir"')
            
        events = acc.get_events()
        if not events:
            logger.error('No event metainformation found for %s\n' % event_name)
            continue
        
        event = events[0]
        event.name = event_name
        
        stations = acc.get_stations(relative_event=event)
        get_angle = lambda tr: stations[tr.network, tr.station, tr.location].backazimuth + 180.
        
        
        chan_count = {}
        
        processed_traces = []
        
        if hasattr(config, 'rotation_table'):
            rotate = (get_angle, config.rotation_table)
        else:
            rotate = None
        
        displacement_limit = None
        if hasattr(config, 'displacement_limit'):
            displacement_limit = config.displacement_limi
        
        for traces in acc.iter_displacement_traces(
                config.restitution_fade_time, 
                config.restitution_frequencyband, 
                deltat=deltat,
                rotate=rotate,
                maxdisplacement=displacement_limit,
                allowed_methods=config.restitution_methods,
                crop=config.restitution_crop):
            
            for tr in traces:
                station = stations[get_nsl(tr)]
                
                if min_dist is not None and station.dist_m < min_dist:
                    continue
                if max_dist is not None and station.dist_m > max_dist:
                    continue
                
                span_complete = True
                
                timings = []
                if hasattr(config, 'check_span'):
                    timings = config.check_span
                    
                if hasattr(config, 'cut_span'):
                    timings.extend(config.cut_span)
                
                for timing in timings:
                    tt = timing(station.dist_m, event.depth)
                    if tt is None:
                        logger.warn('Trace does not contain all required arrivals: %s.%s.%s.%s (timing not present)' % tr.nslc_id)
                        span_complete = False
                        break
                    
                    arrival_time = event.time + tt
                    if not (tr.tmin <= arrival_time and arrival_time <= tr.tmax):
                        what = arrival_time
                        logger.warn('Trace does not contain all required arrivals: %s.%s.%s.%s (timing not in trace)' % tr.nslc_id)
                        span_complete = False
                        break
                                
                if not span_complete:
                    continue
                
                if hasattr(config, 'cut_span'):
                    cs = config.cut_span
                    tmin, tmax = (event.time+cs[0](station.dist_m, event.depth),
                                event.time+cs[1](station.dist_m, event.depth))
                                
                    tr.chop(tmin, tmax, inplace=True)
                
                processed_traces.append(tr)
                
                if tr.channel not in chan_count:
                    chan_count[tr.channel] = 0
                    
                chan_count[tr.channel] += 1
                    
        # use only one sensor at each station; use lexically first
        sstations = stations.values()
        sstations.sort( lambda a,b: cmp(get_nsl(a),  get_nsl(b)) )
        xstations = {}
        for station in sstations:
            if get_ns(station) not in xstations:
                xstations[get_ns(station)] = station
        
        # stations by distance
        dstations = xstations.values()
        dstations.sort( lambda a,b: cmp(a.dist_m, b.dist_m) )
        
        if kiwi_config is not None:
            save_kiwi_dataset(dstations, processed_traces, event, kiwi_config)
        
        if rapid_config is not None:
            save_rapid_dataset(dstations, processed_traces, event, db, rapid_config)
        
        for k,v in chan_count.iteritems():
            logger.info( 'Number of displacement traces for channel %s: %i\n' % (k,v) )
            
        if hasattr(config, 'raw_trace_path'):
            io.save(acc.get_pile().all(), config.path('raw_trace_path'))
            
        logger.info( 'Stopwatch: %5.1f s' % sw() )
    