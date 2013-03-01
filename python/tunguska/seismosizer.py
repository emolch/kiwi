
import os
import sys
import time
import threading
from Queue import Queue, Empty
import tempfile
import subprocess
import signal
from os.path import join as pjoin
import shutil
import logging
import numpy as num
import copy
import progressbar

import config
import phase

logger = logging.getLogger('kiwi.seismosizer')


def load_table(fn):
    f = open(fn, 'r')
    values = [ float(x) for x in f.read().split() ]
    f.close()
    v = num.array(values, dtype=num.float)
    x = v[::2].copy()
    y = v[1::2].copy()
    return (x,y)

def getsigdict():
    r = {}
    for name in dir(signal):
        if name.startswith("SIG"):
           r[getattr(signal, name)] = name
    return r

class Fatal(Exception):
    
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)
        if config.exit_on_fatal:
            sys.exit(args[0])

class SeismosizerBase:
    '''Controls a group of seismosizer processes'''
    
    commands = ['set_database',
                'set_local_interpolation',
                'set_spacial_undersampling',
                'set_receivers',
                'switch_receiver',
                'set_ref_seismograms',
                'set_source_location',
                'set_source_crustal_thickness_limit',
                'get_source_crustal_thickness',
                'set_source_params',
                'set_source_params_mask',
                'set_source_subparams',
                'set_source_constraints',
                'set_effective_dt',
                'set_misfit_method',
                'set_misfit_filter',
                'set_misfit_filter_1',
                'set_misfit_taper',
                'set_synthetics_factor',
                'minimize_lm',
                'output_seismograms',
                'output_seismogram_spectra',
                'output_source_model',
                'get_source_subparams',
                'get_global_misfit',
                'get_misfits',
                'get_peak_amplitudes',
                'get_arias_intensities',
                'get_principal_axes',
                'output_distances',
                'output_cross_correlations',
                'shift_ref_seismogram',
                'autoshift_ref_seismogram',
                'set_floating_shiftrange',
                'get_cached_traces_memory',
                'set_cached_traces_memory_limit',
                'set_verbose',
                'set_ignore_sigint',]
  
    def __init__(self, hosts):
        
        # start processes
        processes = []
        for host in hosts:
            p = SeismosizerProcess( host )
            processes.append(p)

        self.tempdir = tempfile.mkdtemp("","seismosizer-")
            
        self.processes = processes
        signal.signal(signal.SIGINT, self.sighandler)
        signal.signal(signal.SIGTERM, self.sighandler)
        signal.signal(signal.SIGQUIT, self.sighandler)

    def __del__(self):
        self.close()
        shutil.rmtree(self.tempdir)
        
    def sighandler(self, *args):
        logger.warn('Caught signal %s, license to kill' % getsigdict()[args[0]])
        self.close()
        time.sleep(1.0)
        self.terminate()
        
    def close(self):
        for p in self.processes:
            p.stop()
            
    def terminate(self):
        for p in self.processes:
            p.terminate()
            
    def __len__(self):
        return len(self.processes)
                
    def do(self, cmd, where=None):
        '''Parallel execution of same command on selected processes.
        
           everywhere:            s.do( ['command','arg' ] )
           on single process:     s.do( ['command','arg' ], where=2 )
           on group of processes: s.do( ['command','arg' ], where=range(1,4) )
        '''
        
        if where is None:
            processes = self.processes
        else:
            if isinstance(where, int):
                processes = [ self.processes[where] ]
            else:
                processes = [ self.processes[i] for i in where ]
        
        strcommand  =  ' '.join( [str(arg) for arg in cmd ] )
        
        # distribute command to each process
        runners = []
        ready = Queue()
        for p in processes:
            logger.debug('Do (%i): %s' % (p.tid, strcommand))
            p.push(strcommand, ready)
            runners.append(p)
        
        # gather answers
        answers = {}
        errors = {}
        fatal = False
        try:
            while runners:
                p = ready.get()
                if not p:
                    raise ThreadIsDead()

                runners.remove(p)
                try:
                    answer = p.poll()
                    answers[p.tid] = answer
                    logger.debug('Answer (%i): %s' % (p.tid, answer))
                
                except SeismosizerReturnedError, error:
                    errors[p.tid] = error.args[1]
                
        except ThreadIsDead:
            self.close()
            time.sleep(1.0)
            self.terminate()
            raise Fatal('Shit happens!')
       
        if errors:
            raise SeismosizersReturnedErrors(errors)
        
        return [ answers[i] for i in sorted(answers.keys()) ]

# add commands as methods
def gen_do_method(command):
    def func(self, *args, **kwargs):
        return self.do( (func.command,)+args, **kwargs )
        
    func.command = command
    return func

for command in SeismosizerBase.commands:
    method = gen_do_method(command)
    setattr( SeismosizerBase, 'do_'+command, method )

class SeismosizersReturnedErrors(Exception):
    pass

class SeismosizerReturnedError(Exception):
    pass

class SeismosizerDied(Exception):
    pass
    
class ThreadIsDead(Exception):
    pass

class SeismosizersDied(Exception):
    pass

class SeismosizerInsane(Exception):
    pass

class SeismosizerProcess(threading.Thread):
    '''Controls a single seismosizer process in a thread.'''
    
    tid = 0
    
    def __init__(self, host):
        threading.Thread.__init__(self)
        
        self.tid = SeismosizerProcess.tid
        SeismosizerProcess.tid += 1
        self.host = host
        
        cmd = [config.source_info_prog]
        
        self.tempdir = tempfile.mkdtemp("","seismosizer-process-")
        
        if host == 'localhost':
            cmd = [config.seismosizer_prog]
        else:
            cmd = [ 'ssh', host, config.seismosizer_prog ]
        
        try:
            p = subprocess.Popen( cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True )
            self.p = p
            self.to_p = p.stdin
            self.from_p = p.stdout
        except:
            raise Fatal("cannot start %s on %s" % (config.seismosizer_prog, self.host))
            
        self.commands = Queue()
        self.answers = Queue()
        self.the_end_has_come = False
        self.have_sent_term_signal = False
        self.start()
        
  
    def __del__(self):
        self.stop()
        self.join()
        shutil.rmtree(self.tempdir)
            
    def stop(self):
       self.the_end_has_come = True
       self.commands.put(('close',None))
        
    def terminate(self):
        if not self.have_sent_term_signal:
            if self.p.poll() is None:
                try:
                    os.kill(self.p.pid, signal.SIGTERM)
                    logger.warn( 'Sent SIGTERM to seismosizer %i, pid %i' % (self.tid, self.p.pid) )
                except OSError, error:
                    logger.warn( error )

            self.have_sent_term_signal = True
        
    def run(self):
        logger.debug('Starting seismosizer %i.' % self.tid)
        ready = None
        try:
            while True:
                retcode = self.p.poll()
                if retcode is not None: 
                    raise SeismosizerDied('Seismosizer %i died or has been killed. Exit status: %i' % (self.tid, retcode))
                
                self.check_end()
                cmd, ready = self.commands.get()
                
                if cmd == 'close':
                    break

                else:
                    answer = self._do(cmd)
                    self.answers.put(answer)
                    ready.put(self)
                    ready = None
            
                
                
        except (IOError, SeismosizerDied, SeismosizerInsane), e:
            logger.warn( e )
            self.the_end_has_come = True
            if ready is not None:
                ready.put(None)
        
        self.to_p.close()
        self.from_p.close()
        retcode = self.p.wait()
        if retcode != 0:
            logger.error('Seismosizer %i exited with nonzero exit status: %i' % (self.tid, retcode))
         
        logger.debug('Seismosizer %i finished.' % self.tid)
               
    def _do(self, command):
        '''Put command to minimizer and return the results'''
               
        answer = None
        line = command.strip()
        
        self.to_p.write(line+"\n")
        self.to_p.flush()
        
        answer = ''
        retval = self.from_p.readline().rstrip()
        self.check_end()
        
        if retval.endswith('nok'):
            return SeismosizerReturnedError("%s on %s (tid=%i) failed doing command: %s" %
                        (config.seismosizer_prog, self.host, self.tid, line), '' )
                        
        elif retval.endswith('nok >'):
            error_str = self.from_p.readline().rstrip()
            self.check_end()
            return SeismosizerReturnedError("%s on %s (tid=%i) failed doing command: %s" %
                        (config.seismosizer_prog, self.host, self.tid, line), error_str )
            
        elif retval.endswith('ok >'):
            answer = self.from_p.readline().rstrip()
            self.check_end()
            return answer
        
        elif retval.endswith('ok'):
            return ''
            
        else:
            raise SeismosizerInsane('Seismosizer %i did not answer correctly.' % self.tid)
    
    def check_end(self):
        if self.the_end_has_come: 
            mess = '"You gotta go when you gotta go!", says process %i' % self.tid
                    
            raise SeismosizerInsane(mess)
        
    def push(self, cmd, ready):
        '''Enqueue command for seismosizer.'''
        self.commands.put((cmd, ready))
        
    def poll(self):
        '''Check for results of commands.'''
        if not self.isAlive(): raise ThreadIsDead()
        answer = self.answers.get()
        if isinstance(answer, SeismosizerReturnedError):
            raise answer
        return answer

class NoValidSources(Exception):
    pass

class Seismosizer(SeismosizerBase):

    # these commands are passed directly to all running processes
    plain_commands = ['set_local_interpolation',
                      'set_spacial_undersampling',
                      'set_source_crustal_thickness_limit',
                      'set_ref_seismograms',
                      'set_source_params',
                      'set_source_params_mask',
                      'set_source_subparams',
                      'set_source_constraints',
                      'set_effective_dt',
                      'set_synthetics_factor',
                      'output_seismograms',
                      'output_seismogram_spectra',
                      'output_source_model',
                      'output_distances',
                      'output_cross_correlations',
                      'set_floating_shiftrange',
                      'get_cached_traces_memory',
                      'set_cached_traces_memory_limit',
                      'set_verbose',
                      'set_ignore_sigint']
                      
    def __init__(self, hosts, balance_method='123321'):
        SeismosizerBase.__init__(self, hosts)
        self.database = None
        self.receivers = None
        self.source = None
        self.source_location = None
        self.shifts = None
        self.local_interpolation = None
        self.balance_method = balance_method
        
    def set_database(self, gfdb, **kwargs):
        self.database = gfdb
        self.do_set_database( gfdb.path, **kwargs )
      
    def set_source_location(self, *args, **kwargs ):
        self.do_set_source_location(*args,**kwargs)
        self.source_location = args
        self._locations_changed()
        
    def get_source_location(self):
        return self.source_location
        
    def set_receivers(self, receivers, **kwargs ):
        self.receivers = receivers
        
        receiverfn = pjoin(self.tempdir, "receivers")
        file = open(receiverfn, "w")
        for r in receivers:
            file.write( "%s\n" % str(r) )
        file.close()
        self.do_set_receivers(receiverfn, 'has_depth')
        os.remove(receiverfn)
        self._locations_changed()
    
    def blacklist_receivers(self, blacklist):
        for irec, r in enumerate(self.receivers):
            sid = r.get_station()
            if sid in blacklist:
                self.switch_receiver(irec+1, 'off')
    
    def xblacklist_receivers(self, xblacklist):
        '''This blacklist contains the numbers of the receivers to be disabled.
           Yes, I know that this needs to be cleaned up. But I don't have time to
           do that.'''
        for irec in xblacklist:
            self.switch_receiver(irec+1, 'off')
    
    def switch_receiver(self, irec, onoff, **kwargs):
        # irec counted from 1 !!!
        assert onoff in ('on', 'off')
        assert 1 <= irec <= len(self.receivers) 
        if 'where' in kwargs: raise Exception("Can't use kwarg 'where' in switch_receiver")
        kwargs['where'] = self.receivers[irec-1].proc_id
        args = (onoff,)
        self.do_switch_receiver( irec, *args, **kwargs)
        self.receivers[irec-1].enabled = onoff == 'on'
        
    def set_source( self, source, **kwargs ):
        self.source = source
        self.do_set_source_params( source, **kwargs )
        
    def set_taper( self, taper, sourcedepth=None ):
        '''Set a general taper.
           If the taper is built on phases, which are not available at 
           a receivers distance, this will switch off the receiver in question.'''
           
        if isinstance(taper, phase.Taper):
            taper = [ taper ] * len(self.receivers)
            
        for irec, rec in enumerate(self.receivers):
            if taper[irec] is not None:
                #if rec.depth != 0.0:
                #    dist = num.sqrt(rec.distance_m**2 + (sourcedepth - rec.depth)**2)
                #    values = taper[irec](dist, 0.0)
                #else:
                values = taper[irec](rec.distance_m, sourcedepth)
            else:
                values = []
                
            if None in values:
                # Phase(s) not existant at this distance
                self.switch_receiver(irec+1, 'off')
            else:
                self.do_set_misfit_taper( irec+1, *values )
    
    def set_filter( self, filter):
        self.filter = filter
        self.do_set_misfit_filter( *self.filter() )

    def set_filters(self, filters):
        for irec, rec in enumerate(self.receivers):
            filter = filters[irec]
            if filter is not None:
                self.do_set_misfit_filter_1( irec+1, *filter() )
            else:
                self.do_set_misfit_filter_1( irec+1 )
        
    def set_misfit_method( self, method ):
        self.inner_misfit_method = method
        self.do_set_misfit_method( method )
    
    def autoshift_ref_seismograms(self, shift_range, irec_range=None, **kwargs):
        if 'where' in kwargs: raise Exception("Can't use kwarg 'where' in autoshift_ref_seismogram")
        
        # irecs are counted here from 1
        
        if irec_range is None:
            irec_range = range(1,len(self.receivers)+1)
        shifts = []
        for irec in irec_range:
            iproc = self.receivers[irec-1].proc_id
            shift = float(self.do_autoshift_ref_seismogram(irec, shift_range[0], shift_range[1], where=iproc, **kwargs)[0])
            shifts.append(shift)
            if len(self)>1:
                # shift reference seismograms in other seismosizer processes
                others = range(len(self))
                others.pop(iproc)
                self.do_shift_ref_seismogram(irec, shift, where=others, **kwargs)
                
            self.receivers[irec-1].cumulative_shift += shift
            
            
        return shifts
        
    def set_synthetic_reference(self):
        """Calculate seismograms and use these as reference"""
        tempfnbase = self.tempdir + "/syntref"
        self.do_output_seismograms(tempfnbase, "mseed", "synthetics", "plain")
        self.do_set_ref_seismograms(tempfnbase, "mseed")
        
    def shift_ref_seismograms( self, shifts, irec_range=None, **kwargs):
        if 'where' in kwargs: raise Exception("Can't use kwarg 'where' in shift_ref_seismogram")
        
        # irecs are counted here from 1
        
        if irec_range is None:
            irec_range = range(1,len(self.receivers)+1)
        
        for irec, shift in zip(irec_range, shifts):
            self.do_shift_ref_seismogram(irec, shift, **kwargs)
            self.receivers[irec-1].cumulative_shift += shift
        
    def get_receivers_snapshot( self,
                                which_seismograms = ('syn','ref'), 
                                which_spectra     = ('syn','ref'),
                                which_processing  = 'filtered'):
        
        tdir = pjoin(self.tempdir, 'get_receivers_copy')
        if os.path.isdir(tdir):
            logger.warn('Found stale traces output dir')
            shutil.rmtree(tdir)
        os.mkdir(tdir)
        
        tempfnbase1 = pjoin(tdir, "ref_seismogram")
        tempfnbase2 = pjoin(tdir, "syn_seismogram")
        tempfnbase3 = pjoin(tdir, "ref_spectrum")
        tempfnbase4 = pjoin(tdir, "syn_spectrum")
        
        extension = 'table'
        if 'ref' in which_seismograms: self.do_output_seismograms(tempfnbase1, extension, "references", which_processing)
        if 'syn' in which_seismograms: self.do_output_seismograms(tempfnbase2, extension, "synthetics", which_processing)
        if 'ref' in which_spectra: self.do_output_seismogram_spectra(tempfnbase3, "references", which_processing)
        if 'syn' in which_spectra: self.do_output_seismogram_spectra(tempfnbase4, "synthetics", which_processing)
        
        receivers = copy.deepcopy(self.receivers)
        for irec, rec in enumerate(receivers):
            irec_fortran = irec+1
            for icomp, comp in enumerate(rec.components):
                fn1 = '%s-%i-%s.%s' % (tempfnbase1, irec_fortran, comp, extension)
                fn2 = '%s-%i-%s.%s' % (tempfnbase2, irec_fortran, comp, extension)
                fn3 = '%s-%i-%s.%s' % (tempfnbase3, irec_fortran, comp, extension)
                fn4 = '%s-%i-%s.%s' % (tempfnbase4, irec_fortran, comp, extension)
                
                if os.path.isfile(fn1):
                    rec.ref_seismograms[icomp] = load_table( fn1 )
                else:
                    rec.ref_seismograms[icomp] = None
                
                if os.path.isfile(fn2):
                    rec.syn_seismograms[icomp] = load_table( fn2 )
                else:
                    rec.syn_seismograms[icomp] = None
                
                if os.path.isfile(fn3):
                    rec.ref_spectra[icomp] = load_table( fn3 )
                else:
                    rec.ref_spectra[icomp] = None
                    
                if os.path.isfile(fn4):
                    rec.syn_spectra[icomp] = load_table( fn4 )
                else:
                    rec.syn_spectra[icomp] = None
                        
        shutil.rmtree(tdir)
        return receivers
                
    def get_dsm_infos( self ):
        si_base = pjoin(self.tempdir, "source-info")
        self.do_output_source_model( si_base, where=0 )
        si_tdsm = si_base+'-tdsm.info'
        f = open(si_tdsm,'r')
        sect = ''
        infos = {}
        for line in f:
            sline = line.strip()
            if sline in ['ncentroids', '']: 
                sect = sline
                continue
            
            if sect == 'ncentroids':
                infos['ncentroids'] = int(sline)
        return infos
            
    def get_psm_infos( self ):
        si_base = pjoin(self.tempdir, "source-info")
        self.do_output_source_model( si_base, where=0 )
        si_psm = si_base+'-psm.info'
        sections = set(["center","outline","rupture","slip","eikonal-grid", "nucleation-point"])
        f = open(si_psm)
        atsec = ''
        points = []
        data = {}
        def fill(atsec, points):
            if atsec == 'eikonal-grid': nhead = 1
            else: nhead = 0
            data[atsec] = (num.array(points[:nhead],dtype=num.float), num.array(points[nhead:], dtype=num.float))
            
        for line in f:
            sline = line.strip()
            
            if sline == '':    # at a section end
                if atsec != '':
                    fill(atsec, points)
                atsec = ''
                points = []
                continue
            if sline in sections:
                atsec = sline
                continue
            if atsec != '':
                points.append( sline.split() )
        if atsec != '':
            fill(atsec, points)
        
        f.close()
        return data
        
    def _gather_misfits_into_receivers(self, results):
        
        values = [[ float(x) for x in result.split() ] for result in results ]
        ipos = [ 0 ] * len(results)
        for irec, rec in enumerate(self.receivers):
            if rec.enabled:
                iproc = rec.proc_id
                for icomp in xrange(len(rec.components)):
                    rec.misfits[icomp] = values[iproc][ipos[iproc]]
                    ipos[iproc] += 1
                    rec.misfit_norm_factors[icomp] = values[iproc][ipos[iproc]]
                    ipos[iproc] += 1
        
        for iproc in range(len(results)):
            assert ipos[iproc] == len(values[iproc]), "if this is printed, there is a bug in make_misifits_for_source()"
        
    def make_misfits_for_source( self, source ):
        """Calculate misfits for given source and fill these into the receivers datastructure."""
        
        self.set_source(source)
        results = self.do_get_misfits()
        self._gather_misfits_into_receivers(results)
    
    def make_misfits_for_sources(self, sources, show_progress=False, progress_title='grid search'):
        '''Get misfits for many sources gathered by (source,receiver,component).'''
        
        nsources = len(sources)
        nreceivers = len(self.receivers)
        ncomponents = max([ len(r.components) for r in self.receivers ])
        
        # results gathered by (source,receiver,component)
        misfits_by_src = num.zeros( (nsources, nreceivers, ncomponents), dtype=num.float)
        norms_by_src = num.zeros(  (nsources, nreceivers, ncomponents), dtype=num.float)
        
        if nsources == 0: show_progress=False
        
        if show_progress:
            widgets = [progress_title, ' ',
                    progressbar.Bar(marker='-',left='[',right=']'), ' ',
                    progressbar.Percentage(), ' ',]
            
            pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(sources)).start()
        
        failings = []
        for isource, source in enumerate(sources):
            try:
                self.set_source(source)
                results = self.do_get_misfits()
                
                self._gather_misfits_into_receivers(results)
                for ireceiver, receiver in enumerate(self.receivers):
                    if receiver.enabled:
                        for icomp, comp in enumerate(receiver.components):
                            misfits_by_src[isource, ireceiver, icomp] = receiver.misfits[icomp]
                            norms_by_src[isource, ireceiver, icomp] = receiver.misfit_norm_factors[icomp]
                        
            except SeismosizersReturnedErrors:
                failings.append(isource)
                
            if show_progress: pbar.update(isource+1)
                        
        if show_progress: pbar.finish()
        
        return misfits_by_src, norms_by_src, failings
    
    def make_global_misfits_for_sources(self, sources, **outer_misfit_config):
        receiver_mask = num.array([ rec.enabled for rec in self.receivers ], dtype=num.bool)

        misfits_by_src, norms_by_src, failings = self.make_misfits_for_sources(sources)
        misfits_by_s, misfits_by_sr = make_global_misfits( misfits_by_src, norms_by_src, 
            reveiver_mask=receiver_mask, **outer_misfit_config)
            
        misfits_by_s = num.where( misfits_by_s > 0, misfits_by_s, num.NaN)
        
        return misfits_by_s, failings
    
    def best_source(self, sources, return_failings=False, **outer_misfit_config):
        misfits_by_s, failings = self.make_global_misfits_for_sources( sources, **outer_misfit_config)
        ibest = num.nanargmin(misfits_by_s)
        if num.isnan(ibest) or num.isnan(misfits_by_s[ibest]): raise NoValidSources()
        if return_failings:
            return sources[ibest], misfits_by_s[ibest], failings
        return sources[ibest], misfits_by_s[ibest]
    
    def get_peak_amplitudes_for_source( self, ndiff, source ):
        """Calculate peak amplitudes at receivers for given source."""
        
        self.set_source(source)
        results = self.do_get_peak_amplitudes(ndiff)
        values = [[ float(x) for x in result.split() ] for result in results ]
        
        maxabs = [ 0.0 ] * len(self.receivers)
        
        ipos = [ 0 ] * len(results)
        for irec, rec in enumerate(self.receivers):
            iproc = rec.proc_id
            maxabs[irec] = values[iproc][ipos[iproc]]
            ipos[iproc] += 1
    
        return maxabs
        
    def get_arias_intensities_for_source( self, source ):
        """Calculate peak amplitudes at receivers for given source."""
        
        self.set_source(source)
        results = self.do_get_arias_intensities()
        values = [[ float(x) for x in result.split() ] for result in results ]
        
        intensities = [ 0.0 ] * len(self.receivers)
        
        ipos = [ 0 ] * len(results)
        for irec, rec in enumerate(self.receivers):
            iproc = rec.proc_id
            
            intensities[irec] = values[iproc][ipos[iproc]]
            ipos[iproc] += 1
    
        return intensities
    
    # "private" methods:
    
    def _locations_changed(self):
        if self.source_location is None or self.receivers is None: return
        self._fill_distazi()
        self.balance(self.balance_method)
        
    def balance(self, method='123321'):
        
        if method == 'serial':
            for irec in range(len(self.receivers)):
                self._set_receiver_process(irec+1, None)
                
        else:


            #assert len(self.receivers) >= len(self)

            if len(self) > 1:
                distances = [ r.distance_m for r in self.receivers ]
                distances_active = [ r.distance_m for r in self.receivers if r.enabled ]
                dist_range = [ min(distances_active), max(distances_active) ]
                dist_center = (dist_range[0]+dist_range[1])/2.
                dist_delta = (dist_range[1]-dist_range[0])/len(self)
                if method == '123123':
                    for idist, irec in enumerate(num.argsort(distances)):
                        self._set_receiver_process(irec+1, idist % len(self))
                        
                else:
                    for irec, dist in enumerate( distances ):
                        if method == '123321':
                            if dist < dist_center:
                                iproc = max(min(int((dist-dist_range[0])/(dist_delta/2.)),len(self)-1),0)
                            else:
                                iproc = max(min(int((dist_range[1]-dist)/(dist_delta/2.)),len(self)-1),0)
                        elif method == '112233':
                            iproc = max(min(int((dist-dist_range[0])/dist_delta),len(self)-1),0)
                        
                        self._set_receiver_process(irec+1, iproc)
                
            else:
                for irec in range(len(self.receivers)):
                    self._set_receiver_process(irec+1, 0)
            
    def _set_receiver_process(self, irec, iproc):
        self.receivers[irec-1].proc_id = iproc
        if len(self) > 1:
            self.do_switch_receiver(irec, 'off')
            if self.receivers[irec-1].enabled:
                self.do_switch_receiver(irec, 'on', where=iproc)
        
    def _fill_distazi( self ):
        fn = self.tempdir + '/distances'
        self.do_output_distances( fn, where=0 )
        f = open(fn,'r')
        distances_deg = []
        distances_m = []
        azimuths = []
        for irec,line in enumerate(f):
            toks = [ float(x) for x in line.split() ]
            self.receivers[irec].set_distazi( *toks )
            
        f.close()
        os.remove(fn)

def make_global_misfits(misfits_by_src, norms_by_src, receiver_mask=None, receiver_weights=1., outer_norm='l2norm', anarchy=False, bootstrap=False, **kwargs):
    
    nreceivers = misfits_by_src.shape[1]
    
    if isinstance(receiver_weights, float):
        mask2 = None
        rweights = receiver_weights
    else:
        mask2 = receiver_weights != 0
        rweights = receiver_weights[num.newaxis,:].copy()
        
    if bootstrap:
        
        mask = None
        if receiver_mask is not None:
            mask = num.asarray(receiver_mask, dtype=num.bool)
        
        if mask2 is not None:
            if mask is not None:
                mask = num.logical_and(mask2, receiver_mask)
            else:
                mask = mask2
                
        if mask is not None:
            nenabled = sum(mask)
            enabled = num.arange(nreceivers, dtype=num.int)[mask]
            bweights = num.zeros(nreceivers, dtype=num.float)
            _bweights = num.bincount(num.take(enabled, num.random.randint(0,nenabled,nenabled)))
            bweights[:_bweights.size] = _bweights
        else:
            bweights = num.zeros(nreceivers, dtype=num.float)
            _bweights = num.bincount(num.random.randint(0,nreceivers,nreceivers))
            bweights[:_bweights.size] = bweights_x[:]
                    
    if outer_norm == 'l1norm':
        misfits_by_sr = num.sum(misfits_by_src,2)
        norms_by_sr   = num.sum(norms_by_src,2)
        
        if anarchy:
            xrweights = num.zeros(norms_by_sr.shape, dtype=float)
            xrweights[:,:] = rweights
            xrweights /= num.where( norms_by_sr != 0., norms_by_sr, -1.)
            rweights = num.maximum(xrweights, 0.)
        
        if bootstrap:
            rweights *= bweights
        
        misfits_by_sr *= rweights
        norms_by_sr *= rweights
        
        ms = num.sum( misfits_by_sr, 1 )
        ns = num.sum( norms_by_sr, 1 )
        
        misfits_by_s  = num.where(ns > 0., ms/ns, -1.)
        misfits_by_s = num.where(misfits_by_s<0., num.NaN, misfits_by_s)
        
    elif outer_norm == 'l2norm':
        misfits_by_sr = num.sqrt(num.sum(misfits_by_src**2,2))
        norms_by_sr   = num.sqrt(num.sum(norms_by_src**2,2))
        
        if anarchy:
            rweights /= num.where( norms_by_sr != 0., norms_by_sr, -1.)
            rweights = num.maximum(rweights, 0.)
            
        if bootstrap:
            rweights *= num.sqrt(bweights)
        
        misfits_by_sr *= rweights
        norms_by_sr *= rweights
        
        ms = num.sum( (misfits_by_sr)**2, 1 )
        ns = num.sum( (norms_by_sr)**2, 1 )
        
        misfits_by_s  = num.where(ns > 0., num.sqrt(ms/ns), -1.)
        misfits_by_s = num.where(misfits_by_s<0, num.NaN, misfits_by_s)
    
    else:
        raise Exception('unknown norm method: %s' % outer_norm)
    
    return misfits_by_s, misfits_by_sr

for command in Seismosizer.plain_commands:
    method = gen_do_method(command)
    setattr( Seismosizer, command, method )

if __name__ == '__main__':
    s = Seismosizer(['localhost', 'localhost', 'localhost'])
    print s.do_set_effective_dt(0.5)
    print s.do_output_distances() # should fail nicely
    s.close()
    
        
