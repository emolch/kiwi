
import os
import sys
import time
import threading, Queue
import tempfile
import subprocess
import signal
from os.path import join as pjoin
import shutil
import logging

import config

runners_sleep = 0.001
pollers_sleep = 0.001

        
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
                'set_effective_dt',
                'set_misfit_method',
                'set_misfit_filter',
                'set_misfit_taper',
                'set_synthetics_factor',
                'minimize_lm',
                'output_seismograms',
                'output_seismogram_spectra',
                'output_source_model',
                'get_source_subparams',
                'get_global_misfit',
                'get_misfits',
                'get_principal_axes',
                'output_distances',
                'output_cross_correlations',
                'shift_ref_seismogram',
                'autoshift_ref_seismogram',
                'set_verbose']
  
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
        self.close()
    
    def close(self):
        for p in self.processes:
            p.stop()
            
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
        
        strcommand  = ' '.join( [str(arg) for arg in cmd ] )
        
        # distribute command to each process
        runners = []
        for p in processes:
            logging.debug('do (%i): %s' % (p.tid, strcommand))
            p.push( cmd )
            runners.append(p)
            
        # gather answers
        seismosizer_went_mad = False
        answers = {}
        errors = {}
        while runners:
            p = runners.pop(0)
            try:
                
                answer = p.poll()
                answers[p.tid] = answer
                logging.debug('answer (%i): %s' % (p.tid, answer))
                
            except NonePending:
                runners.append(p)
            
            except SeismosizerReturnedError, error:
                errors[p.tid] = error.args[1]
            
            except ThreadIsDead:
                sys.exit('shit happens!')
            
            time.sleep(pollers_sleep)
            
        if errors:
            raise SeismosizersReturnedErrors(errors);
            #self.close()
            #sys.exit('\n'.join([ error.args[0] for error in errors]) + '\n')
            
            
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
    
class NonePending(Exception):
    pass

class ThreadIsDead(Exception):
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
            p = subprocess.Popen( cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE )
            self.p = p
            self.to_p = p.stdin
            self.from_p = p.stdout
        except:
            sys.exit("cannot start %s on %s" % (config.seismosizer_prog, self.host))
            
        self.commands = Queue.Queue()
        self.answers = Queue.Queue()
        self.the_end_has_come = False
        self.start()
  
    def __del__(self):
        self.stop()
        self.join()
        shutil.rmtree(self.tempdir)
            
    def stop(self):
        self.the_end_has_come = True
            
    def run(self):
        while True:
            retcode = self.p.poll()
            if retcode is not None: 
                raise SeismosizerDied('Seismosizer died with an exit status of %i' % retcode)
            try:
                if self.the_end_has_come: break
                cmd = self.commands.get_nowait()
                if type(cmd) is str:
                    cmd = [cmd]
                answer = self._do(cmd)
                self.answers.put(answer)
            except SeismosizerReturnedError, inst:
                self.answers.put(inst)
            except Queue.Empty:
                time.sleep(runners_sleep)
            
        self.to_p.close()
        self.from_p.close()
        self.p.wait()
       
    def _do(self, cmd):
        '''Put command to minimizer and return the results'''
        
        command  = ' '.join( [str(arg) for arg in cmd ] )
       
        answer = None
            
        self.to_p.write(command+"\n")
        self.to_p.flush()
        retval = self.from_p.readline().rstrip()
        
        if retval.endswith('nok'):
            raise SeismosizerReturnedError("%s on %s (tid=%i) failed doing command: %s" %
                        (config.seismosizer_prog, self.host, self.tid, command), '' )
                        
        if retval.endswith('nok >'):
            error_str = self.from_p.readline().rstrip()
            raise SeismosizerReturnedError("%s on %s (tid=%i) failed doing command: %s" %
                        (config.seismosizer_prog, self.host, self.tid, command), error_str )
            
        if retval.endswith('ok >'):
            answer = self.from_p.readline().rstrip()
                
        return answer

    def push(self, cmd):
        '''Enqueue command for seismosizer.'''
        self.commands.put(cmd)
        
    def poll(self):
        '''Check for results of commands.'''
        if not self.isAlive(): raise ThreadIsDead()
        try:
            answer = self.answers.get_nowait()
            if isinstance(answer, SeismosizerReturnedError):
                raise answer
            return answer
            
        except Queue.Empty:
            raise NonePending()

class Seismosizer(SeismosizerBase):

    # these commands are passed directly to all running processes
    plain_commands = ['set_local_interpolation',
                      'set_spacial_undersampling',
                      'set_source_crustal_thickness_limit',
                      'set_ref_seismograms',
                      'set_source_params',
                      'set_source_params_mask',
                      'set_source_subparams',
                      'set_effective_dt',
                      'set_synthetics_factor',
                      'output_seismograms',
                      'output_seismogram_spectra',
                      'output_source_model',
                      'output_distances',
                      'output_cross_correlations',
                      'set_verbose']
                      
    def __init__(self, hosts):
        SeismosizerBase.__init__(self, hosts)
        self.database = None
        self.receivers = None
        self.source = None
        self.source_location = None
        self.taper = None
        self.shifts = None
        self.local_interpolation = None
        
    def set_database(self, gfdb, **kwargs):
        self.database = gfdb
        self.do_set_database( gfdb.path, **kwargs )
      
    def set_source_location(self, *args, **kwargs ):
        self.do_set_source_location(*args,**kwargs)
        self.source_location = args
        self._locations_changed()
        
    def set_receivers(self, receivers, **kwargs ):
        self.receivers = receivers
        
        receiverfn = pjoin(self.tempdir, "receivers")
        file = open(receiverfn, "w")
        for r in receivers:
            file.write( ' '.join( (str(r.lat), str(r.lon), r.components) ) + "\n" )
        file.close()
        self.do_set_receivers(receiverfn)
        os.remove(receiverfn)
        self._locations_changed()
    
    def switch_receiver(self, irec, *args, **kwargs):
        if 'where' in kwargs: raise Exception("Can't use kwarg 'where' in switch_receiver")
        kwargs['where'] = self.receivers[irec].proc_id
        self.do_switch_receiver( irec, *args, **kwargs)
        
    def set_source( self, source, **kwargs ):
        self.source = source
        self.do_set_source_params( source, **kwargs )
        
    def set_taper( self, taper ):
        '''Set a general taper.
           If the taper is built on phases, which are not available at 
           a receivers distance, this will switch off the receiver in question.'''
        self.taper = taper
        for irec, rec in enumerate(self.receivers):
            values = self.taper(rec.distance_m)
            if None in values:
                # Phase(s) not existant at this distance
                self.switch_receiver(irec+1, 'off')
            else:
                self.do_set_misfit_taper( irec+1, *values )
    
    def set_filter( self, filter):
        self.filter = filter
        self.do_set_misfit_filter( *self.filter() )
        
    def set_misfit_method( self, method ):
        self.inner_misfit_method = method
        self.do_set_misfit_method( method )
    
    def autoshift_ref_seismograms(self, irec_range=None, **kwargs):
        if 'where' in kwargs: raise Exception("Can't use kwarg 'where' in autoshift_ref_seismogram")
        
        # irecs are counted here from 1
        
        if irec_range is None:
            irec_range = range(1,len(self.receivers)+1)
        
        for irec in irec_range:
            iproc = self.receivers[irec-1].proc_id
            shift = float(self.do_autoshift_ref_seismogram(irec, where=ipro, **kwargsc))[0]
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
        
        for irec in zip(irec_range, shifts):
            self.do_shift_ref_seismogram(irec, shift, **kwargs)
            self.receivers[irec-1].cumulative_shift += shift
        
            
    def make_misfits_for_source( self, source ):
        """Calculate misfits for given source and fill these into the receivers datastructure."""
        
        self.set_source(source)
        results = self.do_get_misfits()
        values = [[ float(x) for x in result.split() ] for result in results ]
        
        ipos = 0
        for irec, rec in enumerate(self.receivers):
            iproc = rec.proc_id
            for icomp in xrange(len(rec.components)):
                rec.misfits[icomp] = values[iproc][ipos]
                ipos += 1
                rec.misfit_norm_factors[icomp] = values[iproc][ipos]
                ipos += 1
    
    # "private" methods:
    
    def _locations_changed(self):
        if self.source_location is None or self.receivers is None: return
        self._fill_distazi()
        
        if len(self) > 1:
            distances = [ r.distance_m for r in self.receivers ]
            dist_range = [min(distances), max(distances) ]
            dist_center = (dist_range[0]+dist_range[1])/2.
            dist_delta = (dist_range[1]-dist_range[0])/len(self)/2.
            
            for irec, dist in enumerate( distances ):
                if dist < dist_center:
                    iproc = max(min(int((dist-dist_range[0])/dist_delta),len(self)-1),0)
                else:
                    iproc = max(min(int((dist_range[1]-dist)/dist_delta),len(self)-1),0)
                
                self._set_receiver_process(irec+1, iproc)
        
        else:
            for irec in range(len(self.receivers)):
                self._set_receiver_process(irec+1, 0)
            
    def _set_receiver_process(self, irec, iproc):
        self.receivers[irec-1].proc_id = iproc
        if len(self) > 1:
            self.do_switch_receiver(irec, 'off')
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




for command in Seismosizer.plain_commands:
    method = gen_do_method(command)
    setattr( Seismosizer, command, method )

if __name__ == '__main__':
    s = Seismosizer(['localhost', 'localhost', 'localhost'])
    print s.do_set_effective_dt(0.5)
    print s.do_output_distances() # should fail nicely
    s.close()
    
        