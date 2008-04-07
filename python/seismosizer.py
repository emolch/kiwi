#
#     Copyright 2007 Sebastian Heimann
#  
#     Licensed under the Apache License, Version 2.0 (the "License");
#     you may not use this file except in compliance with the License.
#     You may obtain a copy of the License at
#  
#         http://www.apache.org/licenses/LICENSE-2.0
#  
#     Unless required by applicable law or agreed to in writing, software
#     distributed under the License is distributed on an "AS IS" BASIS,
#     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#     See the License for the specific language governing permissions and
#     limitations under the License.
# 

import re
import sys
import os
import time
from tempfile import *
from struct import pack

# only available in python 2.4:
#from subprocess import *

# only needed in 2.3:
import sets
import popen2


def rrmdir(dirname):
    """recursively delete directory"""
    for root, dirs, files in os.walk(dirname, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    os.rmdir(dirname)

class Seismosizer:
    """This is a wrapper to the seismogram calculation part of a minimizer process (see minimizer.f90),
       allowing on-the-fly conversion of the seismosizer output files to vtk files.
       Communication with the child seismosizer process is done via pipes,
       connected to stdin and stdout of the child process.
       Furthermore, at startup, it queries source_info (source_info.f90) for
       information about the possible source model parameterizations, which it
       stores in self.params"""
    
    seismosizer_prog  = "./minimizer"
    source_info_prog  = "./source_info"
    source_params_dir = "source_params"
    
    def __init__(self, gfdb, effective_dt, source_loc, receiverfile):
        """start seismosizer program, see minimizer.f90 for help on the parameters"""
        
        cmd = [Seismosizer.source_info_prog]
        
        self.params = self.GetSourceInfos()
        
        self.tempdir = mkdtemp("","seismosizer-")
        self.tempfilebase = self.tempdir + "/s"
        
        cmd = [Seismosizer.seismosizer_prog]
        
        
        initialization = "\n".join([
            "set_database " + gfdb,
            "set_effective_dt " + str(effective_dt),
            "set_source_location " + str(source_loc[0]) + " " +  str(source_loc[1]) + " 0",
            "set_receivers " + receiverfile,
	    "set_local_interpolation bilinear"])+"\n"
            
        self.calc_and_save = "\n".join([
            "output_seismograms "+self.tempfilebase+" table", 
            "output_source_model "+self.tempfilebase ])
        
        try:
            # the following should be used in python >= 2.4
            # p = Popen( cmd, stdin=PIPE, stdout=PIPE )
            # self.p = p
            # self.to_p = p.stdin
            # self.from_p = p.stdout
            (self.to_p,self.from_p) = os.popen2( cmd )
        except:
            sys.exit("cannot start seismosizer")
    
            
        self.Do(initialization)
            
    def GetSourceInfos(self):
        """get some information about possible sources in seismosizer by asking source_info"""
    
        # get avail. source types
        cmd = [Seismosizer.source_info_prog]
        #try:
            # python 2.4:
            # p = Popen( cmd, stdout=PIPE )
            # pout = p.stdout
        (pout,pin) = popen2.popen2( cmd )
        #except:
        #    sys.exit("cannot start source_info")
            
        for line in pout:
            if re.match(r'\s*source types: ', line):
                types = re.sub(r'\s*source types: ','',line).split()
        
        pout.close()
        # python 2.4:
        # p.wait()
        
        # get parameter names for source types
        params = {}
        for t in types:
            cmd = [Seismosizer.source_info_prog,t]
            try:
                # python 2.4:
                # p = Popen( cmd, stdout=PIPE )
                # pout = p.stdout
                (pout,pin) = popen2.popen2( cmd )
            except:
                sys.exit("cannot start source_info")
                
            for line in pout:
                if re.match(r'\s*parameter names: ', line):
                    pars = re.sub(r'\s*parameter names: ','',line).split()
        
            params[t] = {'names' : pars}
            pout.close()
            # python 2.4:
            # p.wait()
             
        
        # get default, min and max parameters
        for t in types:
            fn = Seismosizer.source_params_dir+"/"+t
            f = open(fn)
            for line in f:
                for k in [ "min", "max", "default" ]:
                    if re.match(r'\s*'+k+': ', line):
                        params[t][k] = [ float(s) for s in 
                                        re.sub(r'\s*'+k+': ','',line).split() ]
                        nparams = len(params[t]['names'])
                        if len(params[t][k]) != nparams:
                            sys.exit("need "+str(nparams)+" parameters in line "+k+
                                     ":... in file "+fn)
        
        return params
        
        
            
    def Calculate(self, sourcetype, params):
        """let seismosizer do it's calculations for a specific source,
           and then convert the output into something VTK can understand"""
        
        if sourcetype not in self.params:
            sys.stderr.write("unknown source type: "+sourcetype+"\n")
            return
            
        nparams =len(self.params[sourcetype]['names'])
        if len(params) != nparams:
            sys.stderr.write(sourcetype+" source needs "+str(nparams)+" parameters\n")
            return
        
        set_source_params = "set_source_params "+sourcetype+" "+' '.join([str(n) for n in params])+"\n"
        self.Do(set_source_params)
        self.Do(self.calc_and_save)
        self.PsmInfoToVtk()
        self.SeismogramsAndDsmToBin()
    
    def Do(self,args):
        """Put commands to minimizer and return the results"""
        
        commands = args.strip().splitlines()
        answers = []
        for command in commands:
            self.to_p.write(command+"\n")
            self.to_p.flush()
            retval = self.from_p.readline().rstrip()
            #print command
            if retval.endswith('nok'):
                sys.exit("minimizer failed doing command: "+ command)
                
            if retval.endswith('ok >'):
                answers.append( self.from_p.readline().rstrip() )
                
        return answers
    
    def Shutdown(self):
        """stop seismosizer and remove temporary files"""
        self.to_p.close()
        self.from_p.close()
        # python 2.4:
        #self.p.wait()
        rrmdir( self.tempdir )
        
        
    def PsmInfoToVtk(self):
        """convert PSM info file to VTK format files"""
        infofilename = self.tempfilebase+"-psm.info"
        # python 2.4:
        #sections = set(["center","outline","rupture","slip"])
        sections = sets.Set(["center","outline","rupture","slip"])
        i = open(infofilename)
        atsec = ''
        points = []
        for line in i:
            sline = line.strip()
            if sline == '':    # at a section end
                if atsec != '':
                    self.PsmInfoToVtk_( atsec, points )
                atsec = ''
                points = []
                continue
            if sline in sections:
                atsec = sline
                continue
            if atsec != '':
                points.append( sline.split() )
        if atsec != '':
            self.PsmInfoToVtk_( atsec, points )
        
        i.close()
        
    def PsmInfoToVtk_(self, atsec, points ):
        """called by PsmInfoToVtk() for every section <atsec> in the PSM infofile
           with point data in <points>, this then writes a vtk file for this section."""
        npoints = len(points)
        vtk_head = """# vtk DataFile Version 3.1 
generated by seismosizer.py
ASCII
DATASET POLYDATA
"""
        
        vtk_head_ugr = """# vtk DataFile Version 3.1 
generated by seismosizer.py
ASCII
DATASET UNSTRUCTURED_GRID
"""
        ofn = self.tempfilebase + "-psm-" + atsec +".vtk"
        o = open(ofn,"w")
        
        def vecstr(v):
            return ' '.join([str(e) for e in v])
        
        if atsec == "outline":
            o.write(vtk_head)
            o.write("POINTS %i FLOAT\n" % (npoints*2))
            for p in points:
                o.write(vecstr(p) + "\n")
            for p in points:
                o.write(vecstr(p[0:2]) + " 0\n")
            o.write("\nPOLYGONS 2 %i\n" % ((npoints+1)*2))
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
            for i in range(1,npoints,2):
                o.write(vecstr(points[i])+"\n")
    
    def SeismogramsAndDsmToBin(self):
        """make binary file which vtk can understand from tabular output from seismosizer"""
        
        for compo in [ "-1-n", "-1-e", "-1-d", "-dsm" ]:
            ifn = self.tempfilebase+compo+".table"
            ofn = self.tempfilebase+compo+".bin"
        
            try:
                i = open(ifn)
            except:
                sys.exit("failed to open file "+ifn)
            try:
                o = open(ofn,"w")
            except:
                sys.exit("failed to open file "+ofn)
            
            for line in i:
                vals = line.split()
                nvals = len(vals)
                val = vals[-1]
                for ival in range(nvals,4):
                    vals.append(val)
                data = pack("ffff", *([float(x) for x in vals]))
                o.write(data)

