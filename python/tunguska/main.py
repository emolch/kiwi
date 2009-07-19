from optparse import OptionParser
import logging, os, sys, shutil, re

from Cheetah.Template import Template

pjoin = os.path.join

def progress_off(option, opt_str, value, parser):
    config.show_progress = False

def install(src, dst):
    dirs = os.path.split(dst)
    d,x = os.path.split(dst)
    dirs = []
    while d and not os.path.isdir(d):
        dirs.append(d)
        d,x = os.path.split(d)
        
    dirs.reverse()
    
    for d in dirs:
        if not os.path.exists(d):
            os.mkdir(d)
    
    shutil.copy(src, dst)
    
def kiwi_main(steps):
    parser = OptionParser()
    parser.add_option('--loglevel', action='store', dest='loglevel', type='choice', choices=('info', 'debug'), default='info')
    parser.add_option('--no-progress', action='callback', callback=progress_off)
    parser.add_option('--no-search', action='store_false', dest='do_search', default=True)
    parser.add_option('--no-forward', action='store_false', dest='do_forward', default=True)
    parser.add_option('--run-id', action='store', dest='run_id', type='string', default='current')
    
    (options, args) = parser.parse_args()
    
    logformat =  '[%(asctime)s] %(levelname)-8s %(message)s'
    dateformat = '%Y-%m-%d %H:%M:%S'
    levels = { 'info' : logging.INFO, 'debug' : logging.DEBUG }
    logging.basicConfig( level   = levels[options.loglevel],
                         format  = logformat,
                         datefmt = dateformat,
                         filename='kiwi.log' )
    console = logging.StreamHandler()
    formatter = logging.Formatter(logformat, datefmt=dateformat)
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)
    
    if options.loglevel == 'debug': progress_off
   
    commands = 'work', 'replot', 'show-steps', 'show-in-config', 'show-out-config', 'show-active-config', 'report'
    try:
        command = args.pop(0)
        
        assert(command in  commands)
    except:
        parser.error('available commands are: '+' '.join(commands))
    
    if command == 'show-steps':
        for step in steps:
            print step.stepname
        sys.exit()
    
    stepnames_all = [ step.stepname for step in steps ]
    
    if command == 'report':
        data = {}
        for step in steps: 
            data[step.stepname] = step
            
        report_dir = 'report'
        templates_dir = 'report_templates'
        template_filenames = []
        for entry in os.listdir(templates_dir):
            if entry.startswith('.') or entry.startswith('_'): continue
            template_filename = pjoin(templates_dir, entry)
            t = Template(file=template_filename, searchList=[ data ])
            
            page = str(t)
            files = [ x[1] or x[3] for x in re.findall(r'("([^"]+\.(png|pdf))"|\'([^\']+\.(png|pdf))\')', page) ]
        
            for file in files:
                install(file, pjoin(report_dir, file))
            
            f = open(pjoin(report_dir,entry),'w')
            f.write( page )
            f.close()
            
        sys.exit()
    
    if args:
        for arg in args:
            if arg != '-' and arg not in stepnames_all:
                parser.error('unknown stepname: %s\n' % arg + 'available stepnames are: '+' '.join(stepnames_all))
                
        stepnumber = dict([ (stepname,i) for i,stepname in enumerate(stepnames_all) ])
        xargs = []
        iarg = 0
        if args[0] == '-': args.insert(0,stepnames_all[0])
        if args[-1] == '-': args.append(stepnames_all[-1])
        while iarg < len(args)-2:
            if args[iarg+1] == '-':
                xargs.extend(stepnames_all[stepnumber[args[iarg]]:stepnumber[args[iarg+2]]+1])
                iarg += 3
            else:
                xargs.append(args[iarg])
                iarg += 1
        while iarg < len(args):
            xargs.append(args[iarg])
            iarg += 1
        stepnames_to_do = xargs
        
    else:
        stepnames_to_do = list(stepnames_all)
    
    for stepname in stepnames_to_do:
        if stepname not in stepnames_all:
            parser.error('unknown stepname: %s\n' % stepname + 'available stepnames are: '+' '.join(stepnames_all))
    
    for step in steps:
        if step.stepname in stepnames_to_do:
            if command == 'work':
                step.work(search=options.do_search, forward=options.do_forward, run_id=options.run_id)
                step.plot(run_id=options.run_id)
                if step.stepname != stepnames_to_do[-1]: logging.info('---')
            if command == 'replot':
                step.plot(run_id=options.run_id)
                if step.stepname != stepnames_to_do[-1]: logging.info('---')
            if command == 'show-in-config':
                step.show_in_config(run_id=options.run_id)
            if command == 'show-out-config':
                step.show_out_config(run_id=options.run_id)
            if command == 'show-active-config':
                step.show_active_config(run_id=options.run_id)
    
