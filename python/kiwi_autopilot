from tunguska import wilber
from os.path import join as pjoin
from subprocess import call
import calendar, time, os, logging, shutil, sys

class GlobalDataRequest(wilber.IrisWilber):
    
    def station_filter(self, station):
        return station.dist > 3.0
        
    def event_filter(self, event):
        return  event.depth < 50. and event.mag > 6.

def event_dir_name(event):
    region = '-'.join(event.region.replace('.','').lower().split())
    date_time = time.strftime('%Y-%m-%d_%H-%M-%S', time.gmtime(event.timestamp))
    return '%s_%s' % (region, date_time)

def get_next_event(events_dir):
    request = GlobalDataRequest(username='sebastian', email='sebastian.heimann@zmaw.de')
    now = time.time()
    recently = now - 60*60*3
    tr = ( calendar.timegm((2009,1,1,0,0,0)), recently )
    
    events = request.get_events(time_range=tr)
    events.sort(lambda a,b: cmp(a.timestamp, b.timestamp))
    events.reverse()
    
    for event in events:
        dirname = pjoin(events_dir, event_dir_name(event))
        if not os.path.exists(dirname):
            os.mkdir(dirname)
            request.get_data(event, outfilename=pjoin(dirname, 'data.seed'), vnetcodes=['_GSN-BROADBAND'])
            return dirname
        
    logging.info('Already have all data for given time range and event selection.')

def prepare_event(skeleton_dir, event_dir):
    oldcwd = os.getcwd()
    files_to_copy = [ 'kiwi_prepare.config', 'kiwi_ampspec', 'report.html' ]
    for fn in files_to_copy:
        shutil.copy(pjoin(skeleton_dir,fn),
                    pjoin(event_dir,fn))

    os.chdir(event_dir)        
    call('kiwi_prepare', 'kiwi_prepare.config')
    os.chdir(oldcwd)

def process_event(event_dir):
    oldcwd = os.getcwd()
    os.chdir(event_dir)
    call('./kiwi_ampspec', 'work')
    call('./kiwi_ampspec', 'report')
    os.chdir(oldcwd)

def post_event(event_dir):
    report_dir = pjoin(event_dir, 'report')
    event_d = os.path.split(event_dir.rstrip('/'))[-1]
    remote_dir = 'www/kinherd/reports'
    host = 'emolchrsync'
    if os.path.isdir(report_dir):
        subprocess.call(['rsync', '-d', event_dir, host+':'+remote_dir ])
        subprocess.call(['rsync', '-a', report_dir+'/', host+':'+pjoin(remote_dir,event_d) ])

logging.basicConfig(level=logging.INFO,)

if len(sys.argv) >= 2:
    event_dir = sys.argv[1]
else:
    event_dir = get_next_event('.')

prepare_event(event_dir)
process_event(event_dir)
post_event(event_dir)


