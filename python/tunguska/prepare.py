from pyrocko import util, io, trace, moment_tensor, rdseed, model

from tunguska import phase, gfdb, edump_access

import sys, os, logging, shutil, time, copy
from os.path import join as pjoin
import numpy as numpy

logger = logging.getLogger('tunguska.prepare')

def get_nsl(x):
    return x.network, x.station, x.location

def get_ns(x):
    return x.network, x.station

def get_badness(streams_badness_dir, time):
    fns = os.listdir(streams_badness_dir)
    fns.sort()
    bef = {}
    for fn in fns:
        toks = fn.split('_')
        if len(toks) != 5 or toks[0] != 'badness': continue
        beg_d, beg_t, end_d, end_t = toks[1:]
        begin = util.ctimegm(beg_d+' '+beg_t.replace('-',':'))
        end = util.ctimegm(end_d+' '+end_t.replace('-',':'))
        dist = 0.
        if time < begin:
            dist = begin-time
        if time > end:
            dist = time-end
        
        bef[dist] = fn
        
    fn = bef[sorted(bef.keys())[0]]
    
    f = open(pjoin(streams_badness_dir, fn), 'r')
    badness = {}
    for line in f:
        toks = line.split()
        if len(toks) != 2: continue
        nslc = toks[0].split('.')
        if len(nslc) != 4: continue
        val = float(toks[1])
        badness[tuple(nslc)] = val
    f.close()
    return badness

def save_rapid_station_table(stations_path, stations):
    '''Save station table in format for rapidinv'''
    util.ensuredirs(stations_path)
    f = open(stations_path, 'w')
    for station in stations:
        nsl = '.'.join((station.network,station.station,station.location))
        f.write('%-10s %15.8e %15.8e\n' % (nsl, station.lat, station.lon))
    f.close()
    
def save_event_info_file(event_info_path, event):
    event.dump(event_info_path)
    
def save_kiwi_dataset(acc, stations, traces, event, config):
    
    if config.has('data_dir'):
        data_dir = config.path('data_dir')
        if os.path.exists(data_dir):
            shutil.rmtree(data_dir)
    
    if config.has('skeleton_dir'): 
        copy_files(config.path('skeleton_dir'), config.path('main_dir'))

    if config.has('raw_trace_path'):
        trace_selector = None
        if config.has('station_filter'):
            trace_selector = lambda tr: get_nsl(tr) in stations and config.station_filter(stations[get_nsl(tr)])
        
        for raw_traces in acc.iter_traces(trace_selector=trace_selector):
            io.save(raw_traces, config.path('raw_trace_path'))
    
    eventfn = config.path('event_info_path')
    util.ensuredirs(eventfn)
    save_event_info_file(eventfn, event)
    
    dstations = stations.values()
    dstations.sort( lambda a,b: cmp(a.dist_m, b.dist_m) )
    
    # gather traces by station
    dataset = []
    used_channels = {}
    for station in dstations:
        station_traces = []
        for tr in traces:
            if get_nsl(tr) == get_nsl(station):
                if tr.channel in config.wanted_channels:
                    station_traces.append(tr)
                    if get_nsl(tr) not in used_channels:
                        used_channels[get_nsl(tr)] = []
                    
                    used_channels[get_nsl(tr)].append(station.get_channel(tr.channel))
                        
        station_traces.sort(lambda a,b: cmp(a.channel, b.channel))
        kiwi_components = ''
        for tr in station_traces:
            kiwi_components += config.kiwi_component_map[tr.channel]
        if station_traces:
            dataset.append( (station, kiwi_components, station_traces) )
    
    if config.has('stations_path'):
        fpath = config.path('stations_path')
        util.ensuredirs(fpath)
        ddstations = copy.deepcopy(dstations)
        for sta in ddstations:
            if get_nsl(sta) in used_channels:
               sta.set_channels(used_channels[get_nsl(sta)])
        model.dump_stations(ddstations, fpath)
    
    if config.has('receivers_path'):
        fpath = config.path('receivers_path')
        util.ensuredirs(fpath)
        f = open(fpath, 'w')
        
    iref = 1
    nsets = 1
    if config.has('nsets'):
        nsets = config.nsets
    for station, components, traces in dataset:
        nsl = '.'.join((get_nsl(station)))
        for i in range(nsets):
            depth = 0.0
            if station.depth is not None:
                depth = station.depth
            if config.has('receivers_path'):
                f.write('%15.8e %15.8e %15.8e %3s %-15s\n' % (station.lat, station.lon, depth, components, nsl) )
            for tr in traces:
                tr = tr.copy()
                if config.trace_time_zero == 'event':
                    tr.shift(-event.time)
                ydata = tr.get_ydata()
                ydata *= config.trace_factor
                fn = config.path('displacement_trace_path', {
                    'ireceiver': iref, 
                    'component': config.kiwi_component_map[tr.channel],
                    'network': tr.network,
                    'station': tr.station,
                    'location': tr.location,
                    'channel': tr.channel})
                io.save([tr], fn)
                
            iref += 1
    if config.has('receivers_path'):
        f.close()
    
    if config.has('reference_time_path'):
        fpath = config.path('reference_time_path')
        f = open(fpath, 'w')
        f.write('%i %s\n' % (event.time, 
                        time.strftime('%Y/%m/%d %H:%M:%S', 
                                        time.gmtime(event.time))))
        f.close()
    
    if config.has('source_origin_path'):
        fpath = config.path('source_origin_path')
        f = open(fpath, 'w')
        f.write('%e %e 0\n' % (event.lat, event.lon))
        f.close()
            

    
def save_rapid_dataset(acc, stations, traces, event, config):
    
    if config.has('data_dir'):
        data_dir = config.path('data_dir')
        if os.path.exists(data_dir):
            shutil.rmtree(data_dir)

    if config.has('skeleton_dir'): 
        copy_files(config.path('skeleton_dir'), config.path('main_dir'))
    
    if config.has('raw_trace_path'):
        trace_selector = None
        if config.has('station_filter'):
            trace_selector = lambda tr: get_nsl(tr) in stations and config.station_filter(stations[get_nsl(tr)])
        
        for raw_traces in acc.iter_traces(trace_selector=trace_selector):
            io.save(raw_traces, config.path('raw_trace_path'))

    dstations = stations.values()
    dstations.sort( lambda a,b: cmp(a.dist_m, b.dist_m) )
    
    save_rapid_station_table(config.path('stations_path'), dstations)
    save_event_info_file(config.path('event_info_path'), event)
    
    used_traces = []
    for station in dstations:
        for tr in traces:
            if get_nsl(tr) == get_nsl(station):
                tr = tr.copy()
                if config.trace_time_zero == 'event':
                    tr.shift(-event.time)
                ydata = tr.get_ydata() 
                ydata *= config.trace_factor
                used_traces.append(tr)
    
    io.save(used_traces, config.path('displacement_trace_path'))
    
def copy_files(source_dir, dest_dir):
    if os.path.exists(dest_dir):
        shutil.rmtree(dest_dir)
    shutil.copytree(source_dir, dest_dir)

def prepare(config, kiwi_config, rapid_config, event_names):

    if config.has('gfdb_path'):
        db = gfdb.Gfdb(config.path('gfdb_path'))
        deltat = db.dt
        min_dist = db.firstx + config.gfdb_margin
        max_dist = db.firstx + (db.nx-1)*db.dx - config.gfdb_margin
    else:
        if config.has('deltat'):
            deltat = config.deltat
        else:
            deltat = None
        db = None
        min_dist = None
        max_dist = None
    
    
    for event_name in event_names:
        
        logger.info('Preparing event %s' % event_name)
        
        config.event_name = event_name
        
        sw = util.Stopwatch()
        if config.has('seed_volume'):
            try:
                acc = rdseed.SeedVolumeAccess(config.path('seed_volume'))
            except rdseed.SeedVolumeNotFound:
                logger.error('SEED volume not found for event %s' % event_name)
                continue
            
            fne = config.path('seed_volume') + '.event'
            if os.path.exists(fne):
                override_event = model.Event(load=fne)
               
        elif config.has('edump_data_dir'):
            acc = edump_access.EventDumpAccess(config.path('edump_data_dir'))
            
            
        elif config.has('custom_accessor'):
            gargs = []
            for arg in config.custom_accessor_args:
                gargs.append( config.mkpath(arg) )
            
            if config.has('plugins_dir'):
                pd = config.path('plugins_dir')
                if pd not in sys.path: sys.path[0:0] = [ pd ]
                
            module_name, class_name = config.custom_accessor
            module = __import__(module_name)
            acc_class = getattr(module, class_name)
            acc = acc_class(*gargs)
            
        else:
            sys.exit('config has neither entry "seed_volume" nor "edump_data_dir" nor "custom_data_accessor"')
            
       
        if override_event:
            event = override_event
        else:
            events = acc.get_events()
            if not events:
                logger.error('No event metainformation found for %s\n' % event_name)
                continue
       
            event = events[0]
        
        event.name = event_name
        
        stations = acc.get_stations(relative_event=event)
        
        chan_count = {}
        
        processed_traces = []
        
        displacement_limit = None
        if config.has('displacement_limit'):
            displacement_limit = config.displacement_limit
                
        extend = None
        if config.has('restitution_pre_extend'):
            extend = config.restitution_pre_extend
        
        crop = True
        if config.has('restitution_crop'):
            crop = config.restitution_crop
        
        projections = None
        if config.has('projection_functions'):
            projections = config.projection_functions
            
        rotations = None
        if config.has('rotation_functions'):
            rotations = config.rotation_functions
        
        whitelist = lambda tr: True
        if config.has('streams_badness_dir') and config.has('streams_badness_limit'):
            badness_dir = config.path('streams_badness_dir')
            badness_limit = config.streams_badness_limit
            badness = get_badness(badness_dir, event.time)
            whitelist = lambda tr: tr.nslc_id in badness and badness[tr.nslc_id] <= badness_limit
            
        station_filter = lambda tr: True
        if config.has('station_filter'):
            station_filter = lambda tr: config.station_filter(stations[get_nsl(tr)])
            
        restitution_off_hack = False
        if config.has('restitution_off_hack'):
            restitution_off_hack = config.restitution_off_hack
            
        trace_selector = lambda tr: station_filter(tr) and whitelist(tr)
        out_stations = {}
        for traces in acc.iter_displacement_traces(
                config.restitution_fade_time, 
                config.restitution_frequencyband,
                deltat=deltat,
                rotations=rotations,
                projections=projections,
                relative_event=event,
                maxdisplacement=displacement_limit,
                allowed_methods=config.restitution_methods,
                trace_selector=trace_selector,
                extend=extend,
                crop=crop,
                out_stations=out_stations,
                restitution_off_hack=restitution_off_hack):
                redundant_channel_priorities=config.get_or_none('redundant_channel_priorities')):
                        
            for tr in traces:
                
                station = stations[get_nsl(tr)]
                
                if min_dist is not None and station.dist_m < min_dist:
                    logger.warn('Station %s is too close to the source (distance = %g m, limit = %g m' % (station.nsl_string(), station.dist_m, min_dist) )
                    continue
                if max_dist is not None and station.dist_m > max_dist:
                    logger.warn('Station %s is too far from the source (distance = %g m, limit = %g m' % (station.nsl_string(), station.dist_m, max_dist) )
                    continue
                span_complete = True
                
                timings = []
                if config.has('check_span'):
                    timings = config.check_span
                    
                if config.has('cut_span'):
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
                        acc.problems().add('gappy', tr.full_id)
                        span_complete = False
                        break
                
                if not span_complete:
                    continue
                
                if config.has('cut_span'):
                    cs = config.cut_span
                    tmin, tmax = (event.time+cs[0](station.dist_m, event.depth),
                                event.time+cs[1](station.dist_m, event.depth))
                                
                    tr.chop(tmin, tmax, inplace=True)
                processed_traces.append(tr)
                
                if tr.channel not in chan_count:
                    chan_count[tr.channel] = 0
                    
                chan_count[tr.channel] += 1
                
        # use only one sensor at each station; use lexically first
        sstations = out_stations.values()
        sstations.sort( lambda a,b: cmp(get_nsl(a),  get_nsl(b)) )
        xstations = {}
        have = {}
        for station in sstations:
            if get_ns(station) not in have:
                have[get_ns(station)] = 1
                xstations[get_nsl(station)] = station
        
        if kiwi_config is not None:
            save_kiwi_dataset(acc, xstations, processed_traces, event, kiwi_config)
        
        if rapid_config is not None:
            save_rapid_dataset(acc, xstations, processed_traces, event, rapid_config)
        
        for k,v in chan_count.iteritems():
            logger.info( 'Number of displacement traces for channel %s: %i\n' % (k,v) )
            
        if config.has('raw_trace_path'):
            io.save(acc.get_pile().all(trace_selector=trace_selector), config.path('raw_trace_path'))
            
        if config.has('problems_file'):
            acc.problems().dump(config.path('problems_file'))
            
        logger.info( 'Stopwatch: %5.1f s' % sw() )
    