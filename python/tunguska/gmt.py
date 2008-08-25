# A thin python wrapper to GMT
#
# GMT commands are translated into method calls by the GMT class.
# GMT options are given as keyword arguments.
# 
# * The following GMT options should not be touched and are handled
#   automatically by the GMT class:
#
#    -P -K -O
#
# * -R and -J default to be appended, use R=False or J=False to turn off.
#
# * Output of a series of GMT commands is accumulated in memory and can be 
#   saved as PS or PDF file with the save() method.
# 
# * Any .gmtdefaults files are ignored. The GMT class uses a fixed set of
#   defaults, which may be altered via an argument to the constructor.
#
# * The figure is always centered on the output page.
# 
# * A corrected bounding box is put at specified margins around the figure.
#
#
# An example usage is at the bottom of this script.

import subprocess
import cStringIO
import re
import os
import shutil
from itertools import izip
from os.path import join as pjoin
import tempfile
import random

golden_ratio   = 1.61803
point          = 1.0/72.0

alphabet = 'abcdefghijklmnopqrstuvwxyz'

gmt_defaults_by_version = {}
gmt_defaults_by_version['4.2.1'] = r'''
#
#       GMT-SYSTEM 4.2.1 Defaults file
#
#-------- Plot Media Parameters -------------
PAGE_COLOR              = 255/255/255
PAGE_ORIENTATION        = portrait
PAPER_MEDIA             = a4+
#-------- Basemap Annotation Parameters ------
ANNOT_MIN_ANGLE         = 20
ANNOT_MIN_SPACING       = 0
ANNOT_FONT_PRIMARY      = Helvetica
ANNOT_FONT_SIZE         = 12p
ANNOT_OFFSET_PRIMARY    = 0.075i
ANNOT_FONT_SECONDARY    = Helvetica
ANNOT_FONT_SIZE_SECONDARY       = 16p
ANNOT_OFFSET_SECONDARY  = 0.075i
DEGREE_SYMBOL           = ring
HEADER_FONT             = Helvetica
HEADER_FONT_SIZE        = 36p
HEADER_OFFSET           = 0.1875i
LABEL_FONT              = Helvetica
LABEL_FONT_SIZE         = 14p
LABEL_OFFSET            = 0.1125i
OBLIQUE_ANNOTATION      = 1
PLOT_CLOCK_FORMAT       = hh:mm:ss
PLOT_DATE_FORMAT        = yyyy-mm-dd
PLOT_DEGREE_FORMAT      = +ddd:mm:ss
Y_AXIS_TYPE             = hor_text
#-------- Basemap Layout Parameters ---------
BASEMAP_AXES            = WESN
BASEMAP_FRAME_RGB       = 0/0/0
BASEMAP_TYPE            = plain
FRAME_PEN               = 1.25p
FRAME_WIDTH             = 0.075i
GRID_CROSS_SIZE_PRIMARY = 0i
GRID_CROSS_SIZE_SECONDARY       = 0i
GRID_PEN_PRIMARY        = 0.25p
GRID_PEN_SECONDARY      = 0.5p
MAP_SCALE_HEIGHT        = 0.075i
TICK_LENGTH             = 0.075i
POLAR_CAP               = 85/90
TICK_PEN                = 0.5p
X_AXIS_LENGTH           = 9i
Y_AXIS_LENGTH           = 6i
X_ORIGIN                = 1i
Y_ORIGIN                = 1i
UNIX_TIME               = FALSE
UNIX_TIME_POS           = -0.75i/-0.75i
#-------- Color System Parameters -----------
COLOR_BACKGROUND        = 0/0/0
COLOR_FOREGROUND        = 255/255/255
COLOR_NAN               = 128/128/128
COLOR_IMAGE             = adobe
COLOR_MODEL             = rgb
HSV_MIN_SATURATION      = 1
HSV_MAX_SATURATION      = 0.1
HSV_MIN_VALUE           = 0.3
HSV_MAX_VALUE           = 1
#-------- PostScript Parameters -------------
CHAR_ENCODING           = ISOLatin1+
DOTS_PR_INCH            = 300
N_COPIES                = 1
PS_COLOR                = rgb
PS_IMAGE_COMPRESS       = none
PS_IMAGE_FORMAT         = ascii
PS_LINE_CAP             = butt
PS_LINE_JOIN            = miter
PS_MITER_LIMIT          = 0
PS_VERBOSE                      = FALSE
GLOBAL_X_SCALE          = 1
GLOBAL_Y_SCALE          = 1
#-------- I/O Format Parameters -------------
D_FORMAT                = %lg
FIELD_DELIMITER         = tab
GRIDFILE_SHORTHAND      = FALSE
GRID_FORMAT             = nf
INPUT_CLOCK_FORMAT      = hh:mm:ss
INPUT_DATE_FORMAT       = yyyy-mm-dd
IO_HEADER               = FALSE
N_HEADER_RECS           = 1
OUTPUT_CLOCK_FORMAT     = hh:mm:ss
OUTPUT_DATE_FORMAT      = yyyy-mm-dd
OUTPUT_DEGREE_FORMAT    = +D
XY_TOGGLE               = FALSE
#-------- Projection Parameters -------------
ELLIPSOID               = WGS-84
MAP_SCALE_FACTOR        = default
MEASURE_UNIT            = inch
#-------- Calendar/Time Parameters ----------
TIME_FORMAT_PRIMARY     = full
TIME_FORMAT_SECONDARY   = full
TIME_EPOCH              = 2000-01-01T00:00:00
TIME_IS_INTERVAL        = OFF
TIME_INTERVAL_FRACTION  = 0.5
TIME_LANGUAGE           = us
TIME_SYSTEM             = other
TIME_UNIT               = d
TIME_WEEK_START         = Sunday
Y2K_OFFSET_YEAR         = 1950
#-------- Miscellaneous Parameters ----------
HISTORY                 = TRUE
INTERPOLANT             = akima
LINE_STEP               = 0.01i
VECTOR_SHAPE            = 0
VERBOSE                 = FALSE'''

presets = {}
presets['test'] = dict(
    BASEMAP_FRAME_RGB       = '255/255/255',
    PAGE_COLOR              = '0/0/0',
    COLOR_BACKGROUND        = '255/255/255',
    COLOR_FOREGROUND        = '0/0/0',
)

gmt_defaults = gmt_defaults_by_version['4.2.1']

gmt_config_defaults = {}
for line in gmt_defaults.splitlines():
    sline = line.strip()
    if not sline or sline.startswith('#'): continue
    k,v = sline.split('=',1)
    gmt_config_defaults[k.strip()] = v.strip()


paper_sizes_a = '''A0 2380 3368
                      A1 1684 2380
                      A2 1190 1684
                      A3 842 1190
                      A4 595 842
                      A5 421 595
                      A6 297 421
                      A7 210 297
                      A8 148 210
                      A9 105 148
                      A10 74 105
                      B0 2836 4008
                      B1 2004 2836
                      B2 1418 2004
                      B3 1002 1418
                      B4 709 1002
                      B5 501 709
                      archA 648 864
                      archB 864 1296
                      archC 1296 1728
                      archD 1728 2592
                      archE 2592 3456
                      flsa 612 936
                      halfletter 396 612
                      note 540 720
                      letter 612 792
                      legal 612 1008
                      11x17 792 1224
                      ledger 1224 792'''
                      

paper_sizes = {}
for line in paper_sizes_a.splitlines():
    k, w, h = line.split()
    paper_sizes[k.lower()] = float(w)*point, float(h)*point


def make_bbox( width, height, gmt_config, margins=(0.8,0.8,0.8,0.8)):
    
    leftmargin, topmargin, rightmargin, bottommargin = margins
    portrait = gmt_config['PAGE_ORIENTATION'].lower() == 'portrait'
    
    paper_size = paper_sizes[gmt_config['PAPER_MEDIA'].lower().rstrip('+')]
    if not portrait: paper_size = paper_size[1], paper_size[0]
            
    xoffset = (paper_size[0] - (width + leftmargin + rightmargin)) / 2.0 + leftmargin;
    yoffset = (paper_size[1] - (height + topmargin + bottommargin)) / 2.0 + bottommargin;
    
    if portrait:
        bb1 = int((xoffset - leftmargin) / point);
        bb2 = int((yoffset - bottommargin) / point);
        bb3 = bb1 + int((width+leftmargin+rightmargin) / point);
        bb4 = bb2 + int((height+topmargin+bottommargin) / point);
    else:
        bb1 = int((yoffset - topmargin) / point);
        bb2 = int((xoffset - leftmargin) / point);
        bb3 = bb1 + int((height+topmargin+bottommargin) / point);
        bb4 = bb2 + int((width+leftmargin+rightmargin) / point);
    
    return xoffset, yoffset, (bb1,bb2,bb3,bb4)

class GMT:

    def __init__(self, width=4., height=None, margins=(0.8,0.8,0.8,0.8), gmt_config=None, preset=None ):
        
        self.gmt_config = dict(gmt_config_defaults)
        
        if preset:
            self.gmt_config.update(presets[preset])
        
        if gmt_config:
            self.gmt_config.update(gmt_config)
        
        self.tempdir = tempfile.mkdtemp("","pygmt-")
        self.gmt_config_filename = pjoin(self.tempdir, 'gmtdefaults')
        f = open(self.gmt_config_filename,'w')
        for k,v in self.gmt_config.iteritems():
            f.write( '%s = %s\n' % (k,v) )
        f.close()
        
        if height is None:
            height = width/golden_ratio
        
        self.output = cStringIO.StringIO()
        self.needstart = True
        self.finished = False
        
        self.width = width
        self.height = height
        self.xoffset, self.yoffset, self.bbox = make_bbox( width, height, self.gmt_config, margins )
       
    
    def __del__(self):
        import shutil
        shutil.rmtree(self.tempdir)
        
        
    def _gmtcommand(self, command, data=None, columns=None, rows=None, finish=False, R=True, J=True, output=None, **kwargs):
        
        assert(not self.finished)
        
        options = []
        
        if not output:
            kwargs['R'] = R
            kwargs['J'] = J
        
        for k,v in kwargs.items():
            if type(v) is bool:
                if v:
                    options.append( '-%s' % k )
            elif type(v) is tuple or type(v) is list:
                options.append( '-%s' % k + '/'.join([ str(x) for x in v]) )
            else:
                options.append( '-%s%s' % (k,str(v)) )
                
        if not output:
            if not finish:
                options.append('-K')
            else:
                self.finished = True
            
            if not self.needstart:
                options.append('-O')
            else:
                options.append('-X%gi' % self.xoffset )
                options.append('-Y%gi' % self.yoffset )
                self.needstart = False
    
        args = [ command ]
        args.extend( options )
        args.append( '+'+self.gmt_config_filename )
                
        p = subprocess.Popen( args, stdin=subprocess.PIPE, stdout=subprocess.PIPE )
        out = p.stdin
        if data is not None:
            out.write( data )
        
        if columns:
            for row in izip(*columns):
                out.write(' '.join([str(x) for x in row]))
                out.write('\n')
            
        if rows:
            for row in rows:
                out.write(' '.join([str(x) for x in row]))
                out.write('\n')
            
        out.close()
        if not output:
            self.output.write( p.stdout.read() )
        else:
            output.write( p.stdout.read() )
            output.close()
 
        p.wait()
    
    def tempfile(self, name=None):
        '''Create and open file in the plots tempdir.'''
        
        if not name: 
            name = ''.join( [ random.choice(alphabet) for i in range(10) ])
        
        fn = pjoin(self.tempdir, name)
        f = open(fn, 'w')
        return f, fn
    
    def save(self, filename=None):
        '''Finish and save figure as PDF or PS file.
           
           If filename ends in '.pdf' a PDF file is created by piping the 
           GMT output through epstopdf.
           
           The bounding box is set according to the margins and size specified
           in the constructor.'''
        
        if not self.finished:
            self.psxy(finish=True)
        
        if self.bbox:
            oldbb = re.compile(r'%%BoundingBox:((\s+\d+){4})')
            newbb = '%%%%BoundingBox: %s' % ' '.join([ str(int(x)) for x in self.bbox ])
        
        if filename:
            tempfn = pjoin(self.tempdir, 'incomplete')
            out = open(tempfn, 'w')
        else:
            out = sys.stdout
            
        if self.bbox:
            out.write(oldbb.sub(newbb,self.output.getvalue()))
        else:
            out.write(self.output.getvalue())
        
        out.close()
        
        if filename.endswith('.pdf'):
            subprocess.call([ 'epstopdf', '--outfile='+filename, tempfn])
        else:
            shutil.move(tempfn, filename)
                
    def __getattr__(self, command):
        def f(*args, **kwargs):
            return self._gmtcommand(command, *args, **kwargs)
        return f
         
         

if __name__ == '__main__':
    plot = GMT(width=4, height=4, margins=[0.5,0.5,0.5,0.5])
    plot.pscoast( R='g', J='E32/30/%gi' % plot.width, B='10g10', D='c', A=10000, S=(114,159,207), G=(233,185,110), W='thinnest')
    plot.save('GMT_example1.pdf')
    plot.save('GMT_example1.ps')


    plot = GMT(width=8., height=8., margins=[0.5,0.5,0.5,0.5])
    plot.pscoast( R='g', J='E%g/%g/%g/%gi' % (0., 0., 180., plot.width), B='0g0', 
                  D='c', A=10000, S=(114,159,207), G=(233,185,110), W='thinnest')
    for i in range(500):
        strike = random.random()*360.
        dip = random.random()*90.
        rake =random.random()*360.-180.
        print strike, dip, rake
        lat = random.random() * 180.-90.
        lon = random.random() * 360.-180.
        plot.psmeca( 
                    rows=[[ lon, lat, 0., strike, dip, rake, 4., 0.,0., '%.3g, %.3g, %.3g' % (strike, dip, rake)]], S='a0.5' )
    
    
    plot.save('Beach.ps')
    plot.save('Beach.pdf')
    
    

