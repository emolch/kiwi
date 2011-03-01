import re
import sys
import os
import time
import tempfile
from struct import pack
import copy
import math
import random
from subprocess import Popen, PIPE, call
import numpy as num
#import matplotlib
#matplotlib.use('PDF')
#import matplotlib.pylab as lab
import pickle
import scipy.optimize

pjoin = os.path.join

def kiwi_aux_dir():
    ieq_home = os.getenv('KIWI_HOME')
    if ieq_home is None:
         sys.exit('KIWI_HOME environment variable not set')
    d = pjoin(ieq_home, 'aux')
    if not os.path.isdir(d):
        sys.exit('directory not found: "%s"' % d)
    return d
    

def dump(data, filename):
    f = open(filename, 'w')
    pickle.dump( data, f )
    f.close()

    
class Config:
    seismosizer_prog  = "minimizer"
    source_info_prog  = "source_info"
    gfdb_info_prog = "gfdb_info"
    plotstyle = [ "1p/0/0/200", "1p/200/0/0", "1p/0/150/50", "1p/226/113/0", "1p/226/0/113" ]
    verbose = 0
    component_names = { 'a':'R@-+@-',
                        'c':'R@--@-',
                        'r':'T@-+@-',
                        'l':'T@--@-',
                        'd':'Z@-+@-',
                        'u':'Z@--@-',
                        'n':'N',
                        'e':'E',
                        's':'S',
                        'w':'W' }
    earthradius = 6371.*1000.
    
class Minimizer:
    """This is a wrapper to the seismogram calculation part of a minimizer process (see minimizer.f90),
       allowing on-the-fly conversion of the seismosizer output files to vtk files.
       Communication with the child seismosizer process is done via pipes,
       connected to stdin and stdout of the child process.
       Furthermore, at startup, it queries source_info (source_info.f90) for
       information about the possible source model parameterizations, which it
       stores in self.params"""
       
    
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
                
            
    def __init__(self):
        """start minimizer program"""
        
        cmd = [Config.source_info_prog]
                
        self.tempdir = tempfile.mkdtemp("","minimizer-")
        self.tempfilebase = self.tempdir + "/s"
        
        cmd = [Config.seismosizer_prog]
                
        try:
            p = Popen( cmd, stdin=PIPE, stdout=PIPE )
            self.p = p
            self.to_p = p.stdin
            self.from_p = p.stdout
        except:
            sys.exit("cannot start minimizer")
            
    def close(self):
        self.to_p.close()
        self.from_p.close()
        
            
    def do(self, command_name, *args):
        """Put commands to minimizer and return the results"""
        
        command  = command_name+' '+' '.join( [str(arg) for arg in args ] )
       
        answer = None
        if Config.verbose >= 2:
            sys.stderr.write( command+"\n" )
            
        self.to_p.write(command+"\n")
        self.to_p.flush()
        retval = self.from_p.readline().rstrip()
        
        if Config.verbose >= 3:
            sys.stderr.write( retval+"\n" )
        
        if retval.endswith('nok'):
            sys.exit("minimizer failed doing command: "+ command)
            
        if retval.endswith('ok >'):
            answer = self.from_p.readline().rstrip()
            if Config.verbose >= 3:
                sys.stderr.write( answer+"\n" )
                
        return answer
    
    def set_receivers( self, lat=None, lon=None, components=None, filename=None ):
        """Set minimizers receiver list to given lists or from file"""
        names = []
        if filename != None:
            file = open(filename, "r")
            lat = []
            lon = []
            components = []
            for line in file:
                toks = line.split()
                if len(toks) >= 3:
                    lat.append(float(toks[0]))
                    lon.append(float(toks[1]))
                    components.append(toks[2])
                if len(toks) == 4:
                    names.append(toks[3])
            file.close()
        
        receiverfn = self.tempdir+"/receivers"
        
        self.lat = lat 
        self.lon = lon
        self.components = components
        self.receiver_names = names
        self.nreceivers = len(lat)
        
        file = open(receiverfn, "w")
        for i in xrange(len(lat)):
            file.write( ' '.join( (str(lat[i]), str(lon[i]), components[i]) ) + "\n" )
            
        file.close()
        
        self.do("set_receivers", receiverfn)
        os.remove( receiverfn )
        
    def set_source( self, source ):
        self.do("set_source_params", source )
        
    def set_source_receiver_circle( self, nreceivers, epidist, comps='ard' ):
        """set source to pole and layout n receivers around it at distance epidist"""
        
        self.do("set_source_location", 90, 0, 0)
        lat = []
        lon = []
        compo = []
        
        epidist = 5. # degrees
        for i in xrange(nreceivers):
            lat.append( 90. - epidist )
            lon.append( i*360./nreceivers )
            compo.append( comps )
        
        self.set_receivers( lat, lon, compo )

    def set_source_receiver_zigzag( self, nreceivers, epidist, comps='ard' ):
        """set source to pole and layout n receivers around it at distance epidist"""
        
        self.do("set_source_location", 90, 0, 0)
        lat = []
        lon = []
        compo = []
        
        epidist = 5. # degrees
        for i in xrange(nreceivers):
            if i % 2:
                lat.append( 90. - epidist )
            else:
                lat.append( 90. - epidist*1.5 )
            lon.append( i*360./nreceivers )
            compo.append( comps )
        
        self.set_receivers( lat, lon, compo )

    def set_source_receiver_random( self, nreceivers, mindist, maxdist, comps='ard' ):
        """set source to pole and layout n receivers around it at distance epidist"""
        
        self.do("set_source_location", 90, 0, 0)
        lat = []
        lon = []
        compo = []
        
        epidist = 5. # degrees
        for i in xrange(nreceivers):
            epidist = ( (maxdist-mindist)*random.random()+mindist ) * 360. / (2.0 * math.pi * Config.earthradius)
            lat.append( 90. - epidist )
            lon.append( random.random()*360. )
            compo.append( comps )
        
        self.set_receivers( lat, lon, compo )

        
    def set_synthetic_reference( self, source ):
        """Calculate seismograms for given source and use these as reference"""
        
        tempfnbase = self.tempdir + "/reference"
        self.set_source( source )
        self.do("output_seismograms", tempfnbase, "mseed", "synthetics", "plain")
        self.do("set_ref_seismograms", tempfnbase, "mseed")
    
    def get_misfits_for_source( self, source, source_has_been_set=False ):
        """Calculate misfit for given source"""
        
        if not source_has_been_set:
            self.set_source( source )
        values = [ float(x) for x in self.do_get_misfits().split() ]
        misfits_by_r = []
        norms_by_r = []
        for irec in xrange(len(self.components)):
            misfits_by_c = []
            norms_by_c = []
            for icomp in xrange(len(self.components[irec])):
                misfits_by_c.append( values.pop(0) )
                norms_by_c.append( values.pop(0) )
                
            misfits_by_r.append(misfits_by_c)
            norms_by_r.append(norms_by_c)
            
        return source, misfits_by_r, norms_by_r
        
    def get_misfits_for_sources( self, sources,  prefunc=None ):
        results = []
        for source in sources:
            if prefunc is not None: 
                source = prefunc(self,source)
            results.append( self.get_misfits_for_source( source ) )
            
        return results
    
    def global_misfit( self, misfits_by_r, norms_by_r ):
        smisfit = 0.
        for misfit_by_c, norms_by_c in zip(misfits_by_r, norms_by_r):
            sm = 0.
            sn = 0.
            for misfit, norm  in zip(misfit_by_c, norms_by_c):
                if norm > 0.:
                    sm += misfit**2
                    sn += norm**2
        
            if sn > 0.0:
                smisfit += sm / sn
                
        return math.sqrt(smisfit)
        
    
        
    def get_best_factor_for_source( self, source ):
        
        def func(factor):
            self.do_set_synthetics_factor( factor )
            misfits_by_r, norms_by_r = self.get_misfits_for_source(source, source_has_been_set=True )[1:]
            return self.global_misfit(misfits_by_r, norms_by_r)
        
        self.set_source(source)
        result = scipy.optimize.brent(func, brack=(0.,1.), tol=1.0e-4)
        
        self.do_set_synthetics_factor(1.)
        return result
        
         
    def get_grid_minimum( self, sm_grid ):
        smg_grid = []
        for (source, misfits, norms) in sm_grid:
            glob_misfit = self.global_misfit( misfits, norms )
            smg_grid.append( (source, glob_misfit) )
            
        misfit, s =  min( [(xm,xs) for (xs,xm) in smg_grid ] )
        return s, misfit
        
    def get_grid_maximum( self, sm_grid ):
        smg_grid = []
        for (source, misfits, norms) in sm_grid:
            glob_misfit = self.global_misfit( misfits, norms )
            smg_grid.append( (source, glob_misfit) )
            
        misfit, s =  max( [(xm,xs) for (xs,xm) in smg_grid ] )
        return s, misfit
        
        
    def get_grid_minimum_bootstrap( self, sm_grid ):
        
        # randomization
        nrecs = len(sm_grid[0][1])
        sample = []
        for i in xrange(nrecs):
            sample.append( random.randint(0,nrecs-1) )
        
        # construct misfit and norm vectors for this randomization
        smb_grid = []
        for (source, misfits, norms) in sm_grid:
            bmisfits = []
            bnorms = []
            for irec in xrange(nrecs):
                bmisfits.append( misfits[sample[irec]] )
                bnorms.append( norms[sample[irec]] )
            smb_grid.append((source, bmisfits, bnorms))
        
        # get minimum in the same way as usual    
        return self.get_grid_minimum( smb_grid )
        
    def grid_bootstrap( self, sm_grid, niter,  sourceparams=None ):
        sources = []
        for i in xrange(niter):
            source, misfit = self.get_grid_minimum_bootstrap( sm_grid )
            sources.append( source )
       
        if sourceparams==None:
            sourceparams = sources[0].params.keys()
        
        meansource = Source(sources[0].sourcetype)
        meansourceparams = {}
        stddevsourceparams = {}
        bresults = {}
        for param in sourceparams:
            param_results = num.zeros(niter, dtype=num.float)
            param_histogram = {}
            for i, source in enumerate(sources):
                value = source.params[param]
                param_results[i] = value
                param_histogram[value] = param_histogram.get(value, 0) + 1
            
            meansourceparams[param] = param_results.mean()
            stddevsourceparams[param] = param_results.std()
            bresults[param] = param_results
            
        return meansourceparams, stddevsourceparams, bresults
        
    def get_misfit_for_source( self, source ):
        """Calculate misfit for given source"""
        
        self.set_source( source )
        misfit = float(self.do("get_global_misfit"))
        return misfit
        
    def get_misfit_for_sources( self, sources, misfit_setup_functions=None ):
        """Calculate misfit for many sources"""
        
        results = []
        for source in sources:
            self.set_source( source )
            if misfit_setup_functions:
                for setup in misfit_setup_functions:
                    setup()
                    misfit = float(self.do_get_global_misfit())
                    results.append( ((source,setup), misfit) )
            else:
                misfit = float(self.do_get_global_misfit())
                results.append( (source, misfit) )
            
        return results
        
    def brute_force_minimize(self, 
                    base_source,
                    param_ranges,
                    source_constraints=None, 
                    prefunc=None,
                    misfit_method='l2norm', 
                    bootstrap_iterations=1000,
                    histogram_filename='histogram-%s.pdf', 
                    result_filename='results-%s.txt',
                    misfit_filename='misfit-%s.txt'):
        
        sourcetype = base_source.sourcetype
        params = []
        pmin = {}
        pmax = {}
        pinc = {}
        plen = {}
        for param, vmin, vmax, vinc in param_ranges:
            params.append( param )
            pmin[param] = vmin
            pmax[param] = vmax
            pinc[param] = vinc
            plen[param] = int(round((vmax-vmin)/vinc))+1
        
        self.do_set_misfit_method( misfit_method )
        
        if prefunc is not None: base_source=prefunc(self,base_source)
        base_misfit = self.global_misfit(
                          *self.get_misfits_for_source(base_source)[1:] )
        
        sources = base_source.make_source_grid( 
            param_ranges,
            source_constraints=source_constraints)
            
        results = self.get_misfits_for_sources( sources, prefunc=prefunc )
       
        opt_source, misfit = self.get_grid_minimum( results )
        worst_source, worst_misfit = self.get_grid_maximum( results )
        
        mean, stddev, bresults = self.grid_bootstrap( results, 
                                                        bootstrap_iterations,
                                                        params )
        
        
        for param in params:
            
            # plot histogram            
            count = {}
            for r in bresults[param]:
                count[r] = count.get(r,0) + 1
                
            counted = ([],[])
            for r in sorted(count.keys()):
                counted[0].append(r)
                counted[1].append(float(count[r]) / bootstrap_iterations)
            
            unit = source_infos(sourcetype)[param].unit
            lab.clf()
            lab.bar(counted[0], counted[1], 
                facecolor=(0.74,0.81,0.74),
                align='center',
                width=pinc[param]*0.7 )
            lab.xlabel('%s [%s]' % (param,unit) )
            lab.ylabel('probability of result')
            lab.xlim(pmin[param]-pinc[param]/2.,pmax[param]+pinc[param]/2.)
            fn = histogram_filename % param
            lab.savefig(fn)
            dump( [counted[0], counted[1], param, unit, pmin[param], pmax[param], pinc[param] ], fn+'.data' )

            fout = open(result_filename % param, 'w')
            fout.write( '--- Grid search and bootstrap results for parameter: %s ---\n' % param )
            fout.write( 'Best value (using all stations):         %g\n' % opt_source.params[param] )
            fout.write( 'Mean value of bootstrap results:         %g\n' % mean[param] )
            fout.write( 'Standard deviation of bootstrap results: %g\n' % stddev[param] )
            fout.write( 'Misfit method used:                      %s\n' % misfit_method )
            fout.write( 'Base misfit value:                       %g\n' % base_misfit )
            fout.write( 'Worst misfit value:                      %g\n' % worst_misfit )
            fout.write( 'Best misfit value:                       %g\n' % misfit )
            fout.write( 'Base source:  %s\n' % str(base_source) )
            fout.write( 'Worst source: %s\n' % str(worst_source) )
            fout.write( 'Best source:  %s\n' % str(opt_source) )
            fout.close()
            
            # plot misfit cross-section
            aparams = num.zeros(len(results), dtype=num.float)
            amisfits = num.zeros(len(results), dtype=num.float)
            for i, (source, misfits, norms) in enumerate(results):
                aparams[i] = source.params[param]
                amisfits[i] = self.global_misfit( misfits, norms )

            lab.clf()
            xunit = source_infos(sourcetype)[param].unit
            if xunit == 'm':
                xunit = 'km'
                aparams /= 1000.
                    
            lab.scatter( aparams, amisfits, facecolor=(0.74,0.81,0.74) )
            lab.xlabel('%s [%s]' % (param,xunit) )
            lab.ylabel('norm. misfit')
            fn = misfit_filename % param
            lab.savefig(fn)
            
            dump( [aparams, amisfits, param, xunit], fn+'.data' )
            
            
        for ixparam in xrange(len(params)):
            for iyparam in xrange(ixparam+1,len(params)):
            
                # plot 2D histograms
                xparam = params[ixparam]
                yparam = params[iyparam]
                
                xunit = source_infos(sourcetype)[xparam].unit
                yunit = source_infos(sourcetype)[yparam].unit
                xunitfactor = 1.
                yunitfactor = 1.
                if xunit == 'm':
                    xunit = 'km'
                    xunitfactor = 1000.
                    
                if yunit == 'm':
                    yunit = 'km'
                    yunitfactor = 1000.
                    
                count = {}
                    
                for r,s in zip(bresults[xparam],bresults[yparam]):
                    count[(r,s)] = count.get((r,s),0) + 1

                counted = ([],[],[])
                for (r,s) in sorted(count.keys()):
                    counted[0].append(r)
                    counted[1].append(s)
                    counted[2].append(count[(r,s)])
                
                lab.clf()
                lab.scatter(num.array(counted[0])/xunitfactor, num.array(counted[1])/yunitfactor, counted[2], facecolor=(0.74,0.81,0.74))
                lab.xlabel('%s [%s]' % (xparam,xunit) )
                lab.ylabel('%s [%s]' % (yparam,yunit) )
                fn = histogram_filename % '-'.join((xparam,yparam))
                lab.savefig(fn)
                dump( [counted, xunitfactor, yunitfactor, xparam, xunit, yparam, yunit ], fn+'.data' )
                
                # plot 2D misfit cross-section
                axparams = num.linspace(pmin[xparam], pmax[xparam], plen[xparam])
                ayparams = num.linspace(pmin[yparam], pmax[yparam], plen[yparam])
                
                amisfits = num.zeros((plen[yparam],plen[xparam]),dtype=num.float) - 1.0
                
                for i, (source, misfits, norms) in enumerate(results):
                    ix = round((source.params[xparam] - pmin[xparam])/pinc[xparam])
                    iy = round((source.params[yparam] - pmin[yparam])/pinc[yparam])
                    if amisfits[iy,ix] < 0.:
                        amisfits[iy,ix] = self.global_misfit( misfits, norms )
                    else:
                        amisfits[iy,ix] = min(amisfits[iy,ix], self.global_misfit( misfits, norms ))

                amisfits_masked = num.ma.masked_less( amisfits, 0.0 )
                    
                lab.clf()
                lab.contourf( axparams/xunitfactor, ayparams/yunitfactor, amisfits_masked )
                lab.xlabel('%s [%s]' % (xparam,xunit) )
                lab.ylabel('%s [%s]' % (yparam,yunit) )
                fn = misfit_filename % '-'.join((xparam,yparam))
                lab.savefig(fn)
                
                dump( [axparams, ayparams, amisfits_masked, xunitfactor, yunitfactor, xparam, xunit, yparam, yunit ], fn+'.data' )
                
        return opt_source, misfit
        
    def get_distazi( self ):
        fn = self.tempdir + '/distances'
        self.do_output_distances( fn )
        f = open(fn,'r')
        distances_deg = []
        distances_m = []
        azimuths = []
        for line in f:
            toks = [ float(x) for x in line.split() ]
            distances_deg.append( toks[0] )
            distances_m.append( toks[1] )
            azimuths.append( toks[2] )
        f.close()
        return distances_deg, distances_m, azimuths
        
    def make_seismograms( self, sources, fnbase ):
        """Calculate seismograms for many sources."""
        
        for i,source in enumerate(sources):
            self.set_source( source )
            self.do("output_seismograms", fnbase+"-"+str(i), "table", "synthetics", "plain")
        
    def show_seismograms_gmt( self, which_set="synthetics", which_processing="plain" ):
        tempfnbase = self.tempdir + "/seismogram"
        outfile = tempfnbase+".pdf"
        self.do("output_seismograms ",tempfnbase, "table", which_set, which_processing)
        fns = self.file_names( tempfnbase, "table" )
        cmd = ["autoplot", " ".join(fns), outfile, "--ydistrib", "--fit", "--width=10" ]
        os.system( " ".join(cmd) )
        return outfile
    
    def show_seismograms_gmt_new_dgg( self, comp='d', which_processing='plain', 
                                    scaling='per-component', addopts=[], plotstyle=None,
                                    timeranges=None, ncols=2, phases=None, filenamepattern=None):
        
        outfile     = filenamepattern
        tempfnbase1 = self.tempdir + "/reference"
        tempfnbase2 = self.tempdir + "/synthetic"
        dummyfile   = self.tempdir + "/dummy"
        
        extension = 'table'

        f = open(dummyfile, 'w')
        f.write( '0 0\n')
        f.close()
        
        if plotstyle is None:
            plotstyle = Config.plotstyle
        
        nrecs = len(self.components)
        
        self.do("output_seismograms", tempfnbase1, "table", "references", which_processing)
        self.do("output_seismograms", tempfnbase2, "table", "synthetics", which_processing)
        
        # layout
        ncomps = len(comp)
        outfiles = []
        nrows = int(math.ceil(nrecs/float(ncols)))
        frames = []
        for irow in xrange(nrows):
            for icol in xrange(ncols):
                for icomp in xrange(ncomps):
                    frames.append( ( "--xnlayout=" + str(1),
                                     "--xilayout=" + str(1),
                                     "--ynlayout=" + str(ncomps),
                                     "--yilayout=" + str(icomp + 1),
                                     "--yspacing=0.1"  ) )
            
        # gather ranges acording to scaling method
        if scaling == 'per-component':
            yranges = []
            for icomp in xrange(ncomps):
                infiles = []
                infiles.extend(  self.file_names( tempfnbase1, extension, components=comp[icomp] )  )
                infiles.extend(  self.file_names( tempfnbase2, extension, components=comp[icomp] )  )
                ymin, ymax = self.autoplot_get_yrange( infiles, addopts )                
                yranges.append( (ymin,ymax) )
        
        if scaling == 'all-the-same':
            yranges = []
            infiles = []
            infiles.extend(  self.file_names( tempfnbase1, extension, components=comp )  )
            infiles.extend(  self.file_names( tempfnbase2, extension, components=comp )  )
            ymin, ymax = self.autoplot_get_yrange( infiles, addopts )                
            yranges.append( (ymin,ymax) )
        
        if scaling == 'per-receiver':
            yranges = []
            infiles = []
            for irec in xrange(nrecs):
                infiles = []
                infiles.extend(  self.file_names( tempfnbase1, extension, components=comp, receivers=[ irec ] )  )
                infiles.extend(  self.file_names( tempfnbase2, extension, components=comp, receivers=[ irec ] )  )
                ymin, ymax = self.autoplot_get_yrange( infiles, addopts )                
                yranges.append( (ymin,ymax) )
            
        if phases:
            distances_deg, distances_m, azi = self.get_distazi()
            
            
        # plot it
        for irec in xrange(nrecs):
            currentoutfile = outfile % (irec+1)
            for icomp in xrange(ncomps):
                if comp[icomp] in self.components[irec]:
                    f1 = tempfnbase1 + "-" + str(irec+1) + "-" + comp[icomp] + "." + extension
                    f2 = tempfnbase2 + "-" + str(irec+1) + "-" + comp[icomp] + "." + extension
                else:
                    f1 = dummyfile
                    f2 = dummyfile
                
                cmd = [ "autoplot", 
                    f1+"["+plotstyle[0]+"]", f2+"["+plotstyle[1]+"]", 
                    currentoutfile, "--fit",
                     "--ylabel="+Config.component_names[comp[icomp]], 
                    "--yannotevery=0", "--yfunit=off", "--ylabelofs=0.15" ]
                cmd.extend( frames[irec*ncomps+icomp] )
                cmd.extend( addopts )
                yrange = None
                if scaling == 'per-component':
                    yrange = yranges[icomp]
                if scaling == 'all-the-same':
                    yrange = yranges[0]
                if scaling == 'per-receiver':
                    yrange = yranges[irec]
                if yrange:
                    cmd.append("--yrange=%g/%g" % yrange )
                axes = 'eW'
                if icomp == 0:
                    axes = axes + 'S'
                if icomp == ncomps-1:
                    axes = axes + 'n'
                    title = "Receiver %i" % (irec+1)
                    if len(self.receiver_names) == nrecs:
                        title += ': ' + self.receiver_names[irec]
                    cmd.append( "--title='%s'" % title )
                    if phases:
                        anotations = open(currentoutfile+'.anot','w')
                        for phase in phases:
                            arrivaltime = phase(distances_m[irec])
                            if arrivaltime is not None:
                                anotations.write( "%f %f %s\n" % (arrivaltime, yrange[0]/2, phase.name) )
                        anotations.close()
                        
                    
                cmd.append( "--axes="+axes )
                if icomp > 0: 
                    cmd.append( "-O" )
                if icomp < ncomps-1:
                    cmd.append( "-K" )
                
                if icomp == 0:
                    cmd.extend( ["--xlabel=time", "--xunit=s"] )
                    
                if not timeranges is None:
                    cmd.append( "--xrange=%f/%f" % timeranges[irec] )
                    
                call(cmd)
                
                if (os.path.isfile(currentoutfile+'.anot')):
                    os.remove(currentoutfile+'.anot')
                
            outfiles.append(currentoutfile)
        
        cmd = ['pdfjoin', '--outfile', outfile % 'all' ]
        cmd.extend(outfiles)
        call(cmd)
        
        return outfile
        
        
    def autoplot_get_yrange(self, infiles, addopts=[] ):
        cmd = [ "autoplot" ]
        cmd.extend( infiles )
        cmd.extend( [ "nofile.pdf", "--fit", "--printenv" ] )
        cmd.extend( addopts )
        plotenv = os.popen( " ".join(cmd), "r")
        penv = {}
        for line in  plotenv:
            k,v = line.strip().split("=")
            if (re.match( r'[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?$', v ) ):
                v = float(v)
                penv[k] = v
        
        return penv['ymin'], penv['ymax']
                        
    def file_names( self, basefn, extension, components=None, receivers=None ):
        fns = []
        
        if receivers is None:
            receivers = range(len(self.components))
            
        for irec in receivers:
            if components == None:
                components = self.components[irec]
            for comp in components:
                if comp in self.components[irec]:
                    fns.append( basefn + "-" + str(irec+1) + "-" + comp + "." + extension)
            
        return fns
        

    def show_grid_results_1d_gmt( self, sourcetype, paramname, traces, plotstyle=None, setupname='rel.misfit' ):
    
        if plotstyle is None:
            plotstyle = Config.plotstyle
            
        xyfns = []
        for i, trace in enumerate(traces):
            tempfnbase = self.tempdir + "/xy" + str(i)
            xyfn = tempfnbase+".table"
            file = open(xyfn, "w")
            
            for point in trace:
                (x, y) = point
                file.write( str(x)+" "+str(y)+"\n" )
                
            file.close()
            xyfns.append( xyfn + "[" + plotstyle[i] + "]" )
            
        outfilename = self.tempdir+"/"+paramname+".pdf"
        
        unit = source_infos(sourcetype)[paramname].unit
        
        cmd = [ "autoplot", ' '.join(xyfns), outfilename, "--fit",
                "--width=4", "--xlabel="+paramname, "--ylabel="+setupname,
                "--yapproxticks=3" ]
        
        if unit == 'm':
            cmd.extend( [ "--xfunit=km", "--xexp=3" ] )
        elif unit == 'm/s':
            cmd.extend( [ "--xfunit=km", "--xexp=3" ] )
        else:
            cmd.extend( [ "--xunit="+unit ] )
        
        os.system( " ".join(cmd) )
        return outfilename
        

    def make_misfit_cross_sections_1d( self, gfdb, s, steps, params, misfit_setup_functions=None ):
        
        results = []
        for param in params:
            if Config.verbose >= 1:
                sys.stderr.write( "-- making misfit cross-section for: " + param + "\n" )
                
            vmin = source_infos(s.sourcetype)[param].soft_min
            vmax = source_infos(s.sourcetype)[param].soft_max
            
            # exceptions:
            if param == 'depth':
                depth_ext = math.sqrt( 0.25*(s.params['length-a']+s.params['length-b'])**2 + 
                                        0.25*(s.params['width'])**2 )
                                        
                if vmin < gfdb.firstz + depth_ext:
                    vmin = gfdb.firstz + depth_ext
                    
                if vmax > gfdb.firstz + (gfdb.nz-1) * gfdb.dz - depth_ext:
                    vmax = gfdb.firstz + (gfdb.nz-1) * gfdb.dz - depth_ext
                    
            if param == 'radius':    
                vmax_limit = s.params['depth']
                if vmax > vmax_limit:
                    vmax = vmax_limit
        
            if param == 'length-a':    
                vmax_limit = 2.0*math.sqrt(s.params['depth']**2 - 0.25*s.params['width']**2 ) - s.params['length-b']
                if vmax > vmax_limit:
                    vmax = vmax_limit
            
            if param == 'length-b':    
                vmax_limit = 2.0*math.sqrt(s.params['depth']**2 - 0.25*s.params['width']**2 ) - s.params['length-a']
                if vmax > vmax_limit:
                    vmax = vmax_limit
            
            if param == 'width':    
                vmax_limit = 2.0*math.sqrt(s.params['depth']**2 - 0.25*(s.params['length-a'] + s.params['length-b'])**2 )
                if vmax > vmax_limit:
                    vmax = vmax_limit
            
            if param == 'moment':
                vmax = s.params['moment']*2.
                vmin = s.params['moment']/2.
            
            vdelta = (vmax - vmin) / float(steps)
            grid = s.make_source_grid( [ ( param, vmin, vmax, vdelta )] )
            results_this_param = self.get_misfit_for_sources( grid, misfit_setup_functions )
            
            for result in results_this_param:
                if misfit_setup_functions:
                    ((source,setup), misfit) = result
                    results.append( ((param,source,setup), misfit) )
                else:
                    (source, misfit) = result
                    results.append( ((param,source), misfit) )
            
        return results

# add commands as methods
def gen_do_method(command):
    def func(self,*args):
        return self.do( func.command, *args )
        
    func.command = command
    return func

for command in Minimizer.commands:
    method = gen_do_method(command)
    setattr( Minimizer, 'do_'+command, method )

class Source:
    def __init__(self, sourcetype='bilateral', sourceparams=None, sourceparams_str=None):
    
        self.sourcetype = sourcetype
        self.params = {}
        
        if not sourceparams_str is None:
            sourceparams_float = [ float(s) for s in sourceparams_str.split() ]
            for i,sparam in enumerate(param_names(self.sourcetype)):
                self.params[sparam] = sourceparams_float[i]
        
        if sourceparams is None:
            sourceparams = {}
            
        for sparam in param_names(sourcetype):
            self.params[sparam] = source_infos(sourcetype)[sparam].default
        
        self.set_params( sourceparams )
        
    def set_params(self, sourceparams):
        for sparam in sourceparams.keys():
            if sparam in param_names(self.sourcetype):
                self.params[sparam] = sourceparams[sparam]
            else:
                raise Exception('invalid source parameter: "%s" for source type: "%s"'% (sparam, self.sourcetype) )
                
    def __str__(self):
    
        values = []
        for sparam in param_names(self.sourcetype):
            values.append( self.params[sparam] )
            
        return self.sourcetype + ' ' + ' '.join( [ str(f) for f in values ] )

    def get_params_as_list(self):
        values = []
        for sparam in param_names(self.sourcetype):
            values.append( self.params[sparam] )
        return values
        
    def set_params_from_list(self, values):
        for sparam, value in zip( param_names(self.sourcetype), values):
            self.params[sparam] = value

    def make_source_grid( self, sourceparams, irecurs=0, paramlist=None, sourceslist=None, source_constraints=None ):
        """create a grid of sources by """
        
        if paramlist is None:
            paramlist = [0.] * len(sourceparams)
            
        if sourceslist is None:
            sourceslist = []
            
        key, vmin, vmax, vstep = sourceparams[ irecurs ]
        n = int(round((vmax-vmin)/vstep))+1
        for i in xrange(n):
            v = vmin + i * vstep
            paramlist[irecurs] = v
            if irecurs < len(sourceparams)-1:
                self.make_source_grid( sourceparams, irecurs+1, paramlist, sourceslist, source_constraints=source_constraints )
            else:
                v = {}
                for i,elem in enumerate(sourceparams):
                    v[elem[0]] = paramlist[i]
                
                s = copy.deepcopy(self)
                s.set_params(v)
                if source_constraints == None or source_constraints(s):
                    sourceslist.append(s)
        
        return sourceslist
        
        
    def make_source_randomize( self, sourceparams, nsources ):
        """create a grid of sources by """
       
        sourceslist = []
        for isource in xrange(nsources):
            v = {}
            for iparam in xrange(len(sourceparams)):
                key, vmin, vmax = sourceparams[ iparam ]
                v[key] = random.uniform(vmin, vmax)
            
            s = copy.deepcopy(self)
            s.set_params(v)
            sourceslist.append(s)
            
        return sourceslist

class Phase:
    def __init__(self,name,filename=None):
    
        self.name = name
        if filename is None:
            filename = os.path.join(kiwi_aux_dir(), 'phases', name)
            
        f = open(filename,'r')
        self.ref_points = []
        for line in f:
            distance, time = [float(x) for x in line.split()]
            self.ref_points.append( (distance, time) )
            
    def __call__(self, distance):
        for (low,high) in zip( self.ref_points[0:-1],self.ref_points[1:len(self.ref_points)]):
            if low[0] <= distance <= high[0]:
                return low[1] + (distance-low[0])/(high[0]-low[0])*(high[1]-low[1])
        return None
        #raise Exception("distance %f out of range [%f,%f] for phase %s" % 
        #                 (distance, self.ref_points[0][0], self.ref_points[-1][0], self.name) )
        
        
class SourceInfo:
    info = None
    info_flat = None
    
def source_types():
    """Get available source types."""
    
    if SourceInfo.info is None:
        SourceInfo.info, SourceInfo.info_flat = get_source_infos()
     
    return SourceInfo.info.keys()
     

def source_infos( sourcetype ):
        """get some information about possible sources"""
    
        if SourceInfo.info is None:
            SourceInfo.info, SourceInfo.info_flat = get_source_infos()
        
        return SourceInfo.info[sourcetype]

def source_infos_flat( sourcetype ):
        """get some information about possible sources"""
    
        if SourceInfo.info is None:
            SourceInfo.info, SourceInfo.info_flat = get_source_infos()
        
        return SourceInfo.info_flat[sourcetype]

def param_names( sourcetype ):
        """returns param names in correct order for given source type"""
    
        if SourceInfo.info is None:
            SourceInfo.info, SourceInfo.info_flat = get_source_infos()
            
            
        return SourceInfo.info_flat[sourcetype]["name"]
            
def get_source_infos():
        """get some information about possible sources in minimizer by asking source_info"""
        info = {}
        info_flat = {}
        
        class Entry:
            pass
        
        # get avail. source types
        cmd = [ Config.source_info_prog ]
        
        source_info = os.popen( ' '.join(cmd), 'r' )
        
        for line in source_info:
            if re.match(r'\s*source types: ', line):
                sourcetypes = re.sub(r'\s*source types: ','',line).split()
        
        
        # get parameter names for source types
        for sourcetype in sourcetypes:
            cmd = [ Config.source_info_prog, sourcetype ]
           
            source_info = os.popen( ' '.join(cmd), 'r' )
            
            key_translate = { 'parameter names: ': 'name',
                              'parameter units: ': 'unit',
                              'parameter hard min: ': 'hard_min',
                              'parameter hard max: ': 'hard_max',
                              'parameter soft min: ': 'soft_min',
                              'parameter soft max: ': 'soft_max',
                              'parameter defaults: ': 'default' }
                              
            params_flat = {}
            try:
              for line in source_info:
                  # string fields
                  for key in ['parameter names: ',
                              'parameter units: ']:
                              
                      if re.match(r'\s*'+key, line):
                          pars = re.sub(r'\s*'+key, '', line).split()
                          pkey = key_translate[key]
                          params_flat[pkey] = pars
                          
                  for key in ['parameter hard min: ',
                              'parameter hard max: ',
                              'parameter soft min: ',
                              'parameter soft max: ',
                              'parameter defaults: ' ]:
                              
                      if re.match(r'\s*'+key, line):
                          pars = [ float(s) for s in re.sub(r'\s*'+key, '', line).split() ]
                          pkey = key_translate[key]
                          params_flat[pkey] = pars
            except:
                print 'uuuups'
            
            source_info.close()
            info_flat[sourcetype] = params_flat
            
            params = {}
            for i,pname in enumerate(params_flat['name']):
                entry = Entry()
                for key in params_flat.keys():
                    setattr(entry, key, params_flat[key][i])
                params[pname] = entry
                        
            info[sourcetype] = params
            
        return info, info_flat
        
def get_gfdb_infos( gfdbpath ):
    
    class Gfdb:
        pass
        
    gfdb_infos_str = {}
    cmd = [ Config.gfdb_info_prog, gfdbpath ]
    gfdb_info = os.popen( ' '.join(cmd), 'r' )
    for line in gfdb_info:
        k,v = line.split('=')
        gfdb_infos_str[k] = v.strip()
    
    gfdb_infos = Gfdb()
    for k in [ 'dt', 'dx', 'dz', 'firstx', 'firstz' ]:
        setattr(gfdb_infos, k, float( gfdb_infos_str[k] ))
        
    for k in [ 'nchunks', 'nx', 'nz', 'ng' ]:
        setattr(gfdb_infos, k, int( gfdb_infos_str[k] ))

    return gfdb_infos

    
def table_to_bin(ifn, ofn):
    """make binary file which vtk can understand from tabular output from seismosizer"""
    
    i = open(ifn)
    o = open(ofn,"w")
    for line in i:
        vals = line.split()
        nvals = len(vals)
        val = vals[-1]
        for ival in range(nvals,4):
            vals.append(val)
            
        data = pack("ffff", *([float(x) for x in vals[:4]]))
        o.write(data)
    
def psm_info_to_vtk(infofilename, outfilenamebase):
    """convert PSM info file to VTK format files"""
    sections = set(["center","outline","rupture","slip","eikonal-grid"])
    i = open(infofilename)
    atsec = ''
    points = []
    for line in i:
        sline = line.strip()
        if sline == '':    # at a section end
            if atsec != '':
                ofn = outfilenamebase+ "-" + atsec +".vtk"
                psm_info_to_vtk_section( atsec, points, ofn )
            atsec = ''
            points = []
            continue
        if sline in sections:
            atsec = sline
            continue
        if atsec != '':
            points.append( sline.split() )
    if atsec != '':
        psm_info_to_vtk_section( atsec, points )
    
    i.close()
    
def psm_info_to_vtk_section(atsec, points, filename ):
    """called by PsmInfoToVtk() for every section <atsec> in the PSM infofile
        with point data in <points>, this then writes a vtk file for this section."""
    npoints = len(points)
    vtk_head = """# vtk DataFile Version 3.1 
generated by minimizer.py
ASCII
DATASET POLYDATA
"""
        
    vtk_head_ugr = """# vtk DataFile Version 3.1 
generated by minimizer.py
ASCII
DATASET UNSTRUCTURED_GRID
"""

    vtk_head_sgr = """# vtk DataFile Version 3.1 
generated by minimizer.py
ASCII
DATASET STRUCTURED_GRID
"""
    
    o = open(filename,"wb")
    
    def vecstr(v):
        return ' '.join([str(e) for e in v])
    
    if atsec == "outline":
        o.write(vtk_head)
        o.write("POINTS %i FLOAT\n" % (npoints*2))
        for p in points:
            o.write(vecstr(p[0:3]) + "\n")
        for p in points:
            o.write(vecstr(p[0:2]) + " 0\n")
        o.write("\nPOLYGONS 2 %i\n" % ((npoints+1)*2))
        o.write(str(npoints)+" ")
        o.write(vecstr(range(0,npoints)) + "\n")
        o.write(str(npoints) + " ")
        o.write(vecstr(range(npoints,npoints*2)) + "\n")
        o.write("\nLINES 2 %i\n" % ((npoints+1)*2))
        o.write(str(npoints)+" ")
        o.write(vecstr(range(0,npoints)) + "\n")
        o.write(str(npoints) + " ")
        o.write(vecstr(range(npoints,npoints*2)) + "\n")

    if atsec == "center":
        o.write(vtk_head)
        o.write("POINTS 2 FLOAT\n")
        o.write(vecstr(points[0][0:2]) + " 0\n" )
        o.write(vecstr(points[0]) + "\n")
        o.write("\nLINES 1 3\n")
        o.write("2 0 1\n")
        
    if atsec == "rupture" or atsec == "slip":
        o.write(vtk_head_ugr)
        o.write("POINTS %i FLOAT\n" % (npoints/2))
        for i in range(0,npoints,2):
            o.write(vecstr(points[i])+"\n")
        o.write("\nPOINT_DATA %i\n" % (npoints/2))
        o.write("VECTORS "+atsec+"vector FLOAT\n")
        for i in xrange(1,npoints,2):
            o.write(vecstr(points[i])+"\n")
            
    if atsec == "eikonal-grid":
        o.write(vtk_head_sgr)
        gridsize = (int(points[0][0]), int(points[0][1]), 1)
        o.write("DIMENSIONS %i %i %i\n" % gridsize)
        o.write("POINTS %i FLOAT\n" % (gridsize[0]*gridsize[1]))
        for i in xrange(1,npoints):
            o.write(vecstr(points[i][0:3])+"\n")
        o.write("\nPOINT_DATA %i\n" % (gridsize[0]*gridsize[1]))
        o.write("SCALARS rupturetime FLOAT 1\n")
        o.write("LOOKUP_TABLE default\n")
        
        itimecol = 3
        if len(points[1]) > 4:
            itimecol = 5
        for i in xrange(1,npoints):
            o.write('%s\n' % points[i][itimecol])
            
    o.close()

            