
import sys
import time
import threading
import tempfile
import subprocess

import config

runners_sleep = 0.1
pollers_sleep = 0.1

class Seismosizer(SeismosizerBase):

    # these commands are passed directly to all running processes
    plain_commands = ['set_local_interpolation',
                      'set_spacial_undersampling',
                      'set_source_crustal_thickness_limit',
                      'set_ref_seismograms',
                      'set_source_location',
                      'set_source_params',
                      'set_source_params_mask',
                      'set_source_subparams',
                      'set_effective_dt',
                      'set_misfit_method',
                      'set_misfit_filter',
                      'set_misfit_taper',
                      'set_synthetics_factor',
                      'output_seismograms',
                      'output_seismogram_spectra',
                      'output_source_model',
                      'output_distances',
                      'output_cross_correlations',
                      'shift_ref_seismogram',
                      'set_verbose']
                      
    def __init__(self, hosts):
        SeismosizerBase.__init__(self, hosts)
        self.database = None
        self.receivers = None
        self.source = None
        
    def set_database(self, gfdb, **kwargs)
        self.database = gfdb
        self.do_set_database( gfdb.path, **kwargs )
        
    def set_receivers(self, receivers, **kwargs )
        self.receivers = receivers
        
        receiverfn = pjoin(self.tempdir, "receivers")
        file = open(receiverfn, "w")
        for r in receivers:
            file.write( ' '.join( (str(r.lat), str(r.lon), r.components) ) + "\n" )
            
        file.close()

        self.do_set_receivers(receiversfn)
        os.remove(receiversfn)
   
    def set_receiver_process(self, irec, iproc)
        self.receivers[irec].proc_id = iproc
        self.switch_receiver(irec, 'off')
        self.switch_receiver(irec, 'on', where=iproc)
   
    def switch_receiver(self, *args, **kwargs)
        if 'where' in kwargs: raise Exception("Can't use kwarg 'where' in switch_receiver")
        self.do_switch_receiver(*args, where=self.receivers.proc_id, **kwargs)
        
    def set_source( self, source, **kwargs )
        self.source = source
        self.do_set_source_params( source, **kwargs )
        
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
                processes = self.processes[where]
            else:
                processes = [ self.processes[i] for i in where ]
        
        # distribute command to each process
        runners = []
        for p in processes:
            p.push( cmd )
            runners.append(p)
            
        # gather answers
        seismosizer_went_mad = False
        answers = {}
        errors = []
        while runners:
            p = runners.pop(0)
            try:
                answer = p.poll()
                answers[p.tid] = answer
                
            except NonePending:
                runners.append(p)
            
            except SeismosizerReturnedError, error:
                errors.append(error)
            
            time.sleep(pollers_sleep)
            
        if errors:
            self.close()
            sys.exit('\n'.join([ error.args[0] for error in errors]) + '\n')
            
        return [ answers[i] for i in sorted(answers.keys()) ]
    
    def close(self):
        self.do('close')

# add commands as methods
def gen_do_method(command):
    def func(self, *args, **kwargs):
        return self.do( (func.command,)+args, **kwargs )
        
    func.command = command
    return func

for command in Seismosizer.commands:
    method = gen_do_method(command)
    setattr( Seismosizer, 'do_'+command, method )


class SeismosizerReturnedError(Exception):
    pass
    
class NonePending(Exception):
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
            
        self.commands = []
        self.answers = []
        self.start()
            
    def run(self):
        while True:
            if self.commands:
                cmd = self.commands.pop(0)
                if type(cmd) is str:
                    cmd = [cmd]
                if cmd[0] == 'close':
                    self.answers.append('closed') 
                    break
                try:
                    answer = self._do(cmd)
                    self.answers.append(answer)
                except SeismosizerReturnedError, inst:
                    self.answers.append(inst)
            time.sleep(runners_sleep)
        self.close()

    def close(self):
        '''Shut down connected seismosizer process.'''
        self.to_p.close()
        self.from_p.close()
        self.p.wait()
                
    def _do(self, cmd):
        '''Put command to minimizer and return the results'''
        
        command  = ' '.join( [str(arg) for arg in cmd ] )
       
        answer = None
        if config.verbose >= 2:
            sys.stderr.write( 'tid=%i: %s\n' % (self.tid, command) )
            
        self.to_p.write(command+"\n")
        self.to_p.flush()
        retval = self.from_p.readline().rstrip()
        
        if config.verbose >= 3:
            sys.stderr.write( retval+"\n" )
        
        if retval.endswith('nok'):
            raise SeismosizerReturnedError("%s on %s (tid=%i) failed doing command: %s" %
                        (config.seismosizer_prog, self.host, self.tid, command) )
            
        if retval.endswith('ok >'):
            answer = self.from_p.readline().rstrip()
            if config.verbose >= 3:
                sys.stderr.write( answer+"\n" )
                
        return answer

    def push(self, cmd):
        '''Enqueue command for seismosizer.'''
        self.commands.append(cmd)
        
    def poll(self):
        '''Check for results of commands.'''
        if self.answers:
            answer = self.answers.pop(0)
            if isinstance(answer, SeismosizerReturnedError):
                raise answer
        else:
            raise NonePending()


if __name__ == '__main__':
    s = Seismosizer(['localhost', 'localhost', 'localhost'])
    print s.do_set_effective_dt(0.5)
    print s.do_output_distances() # should fail nicely
    s.close()
    
        