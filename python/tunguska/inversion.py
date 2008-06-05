import sys
import receiver
import gfdb



def standard_setup( datadir,
                    gfdb_path,
                    effective_dt,
                    spacial_undersampling = [ 1, 1 ],
                    hosts = ['localhost'],
                    crustal_thickness_limit = None,
                    # not so important:
                    local_interpolation = 'bilinear',
                    source_origin_file = 'source-origin.table',
                    receivers_file = 'receivers.table' 
                    ref_seismogram_stem = 'reference',
                    ref_seismogram_format = 'sac',
                    )
                    
    '''Start seismosizers and setup source location, Greens functions, 
       receivers, and reference seismograms.'''

    source_origin_file      = pjoin(datadir, source_origin_file)
    ref_seismogram_stem     = pjoin(datadir, ref_seismogram_stem)
    receivers_file          = pjoin(datadir, receivers_file)

    # setup database
    database = gfdb.Gfdb(gfdb_path)
    seis = Seismosizer(hosts)
    seis.set_database(database)
    seis.set_effective_dt(effective_dt)
    seis.set_local_interpolation(local_interpolation)
    seis.set_spacial_undersampling(*spacial_undersampling)
    
    # setup source origin
    f = open( source_origin_file, 'r' )
    (slat, slon, stime) = [ float(x) for x in f.read().split() ]
    f.close()
    seis.set_source_location(slat, slon, stime)
    seis.set_source_crustal_thickness_limit( crustal_thickness_limit )
    
    # setup receivers
    receivers = receiver.load_table(receivers_file)
    seis.set_receivers(receivers)
    seis.set_ref_seismograms( ref_seismogram_stem, ref_seismogram_format )
    seis.fill_distazi(receivers)
    distances = [ r.distance_m for r in receivers ]
    dist_range = [min(distances), max(distances) ]
    dist_center = (dist_range[0]+dist_range[1])/2.
    dist_delta = (dist_range[1]-dist_range[0])/len(seis)/2.
    
    for irec, dist in enumerate( distances ):
        if dist < dist_center:
            iproc = max(min(int((dist-dist_range[0])/dist_delta),len(seis)-1),0)
        else:
            iproc = max(min(int((dist_range[2]-dist)/dist_delta),len(seis)-1),0)
        
        seis.set_receiver_process(irec, iproc)
        
    return seis
    

