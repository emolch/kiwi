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

import os
import sys
import vtk
import time
import signal
from seismosizer import *
from qt import *


class MySlider(QSlider):

    def __init__(self, *args):
        apply(QSlider.__init__, (self,) + args)
        self.connect( self, SIGNAL("valueChanged(int)"), self.myValueChange )
    
    def myValueChange(self):
        self.emit(PYSIGNAL("slideto"), (self.value(), ))

class MyLineEdit(QLineEdit):

    def __init__(self, *args):
        apply(QLineEdit.__init__, (self,) + args)
        self.connect( self, SIGNAL("returnPressed()"), self.myReturnPressed )
        self.connect( self, SIGNAL("lostFocus()"), self.myReturnPressed )
        self.setValidator( QDoubleValidator(self) )
        
    def myReturnPressed(self):
        try:
            x = float(str(self.text()))
            self.emit(PYSIGNAL("edited"), (x, ))
        except:
            pass

class ValControl(QFrame):

    def __init__(self, *args):
        apply(QFrame.__init__, (self,) + args)
        self.layout = QHBoxLayout( self )
        self.layout.setSpacing(5)
        self.lname = QLabel( "name", self )
        self.lname.setFixedWidth(100)
        self.lvalue = MyLineEdit( "value", self )
        self.lvalue.setFixedWidth(100)
        self.slider = MySlider(Qt.Horizontal, self)
        self.slider.setMaxValue( 1000000 )
        self.slider.setLineStep( 10000 )
        self.slider.setPageStep( 100000 )
        self.layout.addWidget( self.lname )
        self.layout.addWidget( self.lvalue )
        self.layout.addWidget( self.slider )
        #self.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        self.connect( self.slider, PYSIGNAL("slideto"),
                      self.slided )
        self.connect( self.lvalue, PYSIGNAL("edited"),
                      self.edited )
                          
    def setup(self, name, mi, ma, cur, ind):
        self.lname.setText( name )
        self.mi = mi
        self.ma = ma
        self.cur = cur
        self.cursl = float(cur-mi)/(ma-mi) * 1000000.
        self.ind = ind
        self.slider.setValue( self.cursl )
        self.adjtext()
        
    def slided(self,val):
        if self.cursl != val:
            self.cursl = val
            self.cur = self.mi+(self.ma-self.mi)*self.cursl/1000000.
            self.adjtext()
            self.emit(PYSIGNAL("valchange"), (self.cur, self.ind, ))

    def edited(self,val):
        if self.cur != val:
            self.cur = val
            cursl = (self.cur-self.mi)/(self.ma-self.mi) * 1000000.
            if (cursl != self.cursl):
                self.slider.setValue( cursl )
            
            self.emit(PYSIGNAL("valchange"), (self.cur, self.ind, ))
        
    def adjtext(self):
        self.lvalue.setText( "%8.5g" % self.cur )
        
def printChildren(obj, indent=""):
    children=obj.children()
    if children==None:
        return
    for child in children:
        print indent, child.name(), child.__class__
        printChildren(child, indent + "  ")
        
class ValControlPan(QScrollView):
    
    def __init__(self, *args ):
        apply(QScrollView.__init__, (self,) + args  )
        self.frame = False
   
    def setup(self, names, mins, maxs, curs):
        n = len(names)
        self.setResizePolicy(QScrollView.AutoOneFit)
        if (self.frame):
            self.frame.hide()
            self.removeChild(self.frame)
            del self.frame
            
            
        self.frame = QFrame(self)
        self.addChild( self.frame )
        grid = QGridLayout(self.frame,n,1)
        grid.setMargin(5)
        grid.setSpacing(5)
        self.data = curs
        for i in range(n):
            valc = ValControl(self.frame)
            grid.addWidget(valc,i,1)
            valc.setup( names[i], mins[i], maxs[i], curs[i], i)
            self.connect( valc, PYSIGNAL("valchange"), self.valchange )
        self.frame.show()
            
    def valchange(self, value, index):
        self.data[index] = value
        self.emit( PYSIGNAL("valchange"), (self.data, ) )

class MyWindow(QMainWindow):

    def __init__(self, clargs, *args):
        apply(QMainWindow.__init__, (self,) + args )
        self.valpan = ValControlPan( self )

        self.setCaption( "Seismosizer" )
        
        self.filemenu = QPopupMenu( self )
        
        self.menuBar().insertItem( "&File", self.filemenu )
        self.filemenu.insertItem( "&Quit", MyApp.app.myquit )
        
        self.sourcemenu = QPopupMenu( self )
        self.sourcemenu.setCheckable( True )
        self.menuBar().insertItem( "Source&type", self.sourcemenu )
        
        if (len(clargs) != 6): 
            sys.exit("usage: show_seismograms.py database effective-dt origin-lat origin-lon receiverfile")
        
        (pn, gfdb, effective_dt, origin_lat, origin_lon, receiverfile) = clargs
        
        #(fid,self.receiverfile) = mkstemp()
        #f = open(self.receiverfile,'w')
        #f.write("0 10.8\n")
        #f.close()
        self.receiverfile = receiverfile
        
        self.seismosizer = Seismosizer( gfdb, effective_dt, [origin_lat,origin_lon],
                                        self.receiverfile)
        self.vis = False
        defaultsource = "bilateral"
        self.SetSourceType( defaultsource )
        
        self.setCentralWidget( self.valpan )
        self.connect( self.valpan, PYSIGNAL("valchange"),
                      self.Calculate )
        self.vis = Visualization(self.seismosizer.tempfilebase)
        
        self.mdict = {}
        
        for name in self.seismosizer.params:
            mid =self.sourcemenu.insertItem( name, self.SetSourceTypeById )
            self.mdict[mid] = name
            self.sourcemenu.setItemChecked( mid, name == defaultsource )
               
    def SetSourceTypeById(self, ident):
        name = self.mdict[ident]
        self.SetSourceType( name )
        for mid,n in self.mdict.iteritems():
            self.sourcemenu.setItemChecked( mid, name == n )
        
    def SetSourceType(self, name):
        
        self.sourcetype = name
        para = self.seismosizer.params[name]
        self.valpan.setup( para['names'], para['min'], para['max'], para['default'] )
        self.Calculate( para['default'] )
        
        
    def Calculate(self,data):
        
        self.seismosizer.Calculate( self.sourcetype, data )
        if self.vis:
            self.vis.modified()
    
class MyApp(QApplication):

    app = 0

    def __init__(self, *args):
        
        apply(QApplication.__init__, (self,) + args)
        self.timer = QTimer( self )
        self.connect( self.timer, SIGNAL("timeout()"), self.periodical ) 
        self.timer.start( 1000, False )
        self.goaway = False
        
    def myquit(self):
        self.win.seismosizer.Shutdown()
        self.quit()
    
    def periodical(self):
        self.vis.rerender()
        if self.goaway:
            self.myquit()
        
    def timetogo(self,*args):
        self.goaway = True
    
    def setwin( self, win ):
        self.win = win
                                
    def setvis( self, vis ):
        self.vis = vis
        
class Visualization:
    def __init__(self,fnbase):
    
        fpreader = vtk.vtkPolyDataReader()
        fpreader.SetFileName( fnbase+"-psm-outline.vtk" )
        self.fpreader = fpreader
        
        fpmapper = vtk.vtkPolyDataMapper()
        fpmapper.SetInput( fpreader.GetOutput() )
        
        fpmapper2 = vtk.vtkPolyDataMapper()
        fpmapper2.SetInput( fpreader.GetOutput() )
        
        fpactor = vtk.vtkActor()
        fpactor.SetMapper( fpmapper );
        fpactor.GetProperty().SetOpacity(0.3)
        fpactor.GetProperty().SetColor( 1.,1.,1. )
        
        fpactor2 = vtk.vtkActor()
        fpactor2.SetMapper( fpmapper2 );
        fpactor2.GetProperty().SetAmbientColor( 1.,1.,1. )
        fpactor2.GetProperty().SetRepresentationToWireframe()
        fpactor2.GetProperty().BackfaceCullingOff()
        fpactor2.GetProperty().SetLineWidth(2.)
        fpactor2.GetProperty().SetOpacity(0.7)
  
        
        cereader = vtk.vtkPolyDataReader()
        cereader.SetFileName( fnbase+"-psm-center.vtk" )
        self.cereader = cereader
        
        cemapper = vtk.vtkPolyDataMapper()
        cemapper.SetInput( cereader.GetOutput() )
        
        ceactor = vtk.vtkActor()
        ceactor.SetMapper( cemapper );
        ceactor.GetProperty().SetDiffuseColor( 1.,1.,1. )
        ceactor.GetProperty().SetRepresentationToWireframe()
        
        
        rupreader = vtk.vtkDataSetReader()
        rupreader.SetFileName( fnbase+"-psm-rupture.vtk" )
        self.rupreader = rupreader
        
        arrow = vtk.vtkArrowSource()
        arrow.SetShaftRadius(0.03)
        arrow.SetTipLength(0.16)
        arrow.SetTipRadius(0.05)
        
        arrows = vtk.vtkGlyph3D()
        arrows.SetInput(rupreader.GetOutput())
        arrows.SetSource(arrow.GetOutput())
        arrows.SetVectorModeToUseVector()
        arrows.SetScaleModeToScaleByVector()
        arrows.OrientOn()
        
        rupmapper = vtk.vtkPolyDataMapper()
        rupmapper.SetInput( arrows.GetOutput() )
        
        rupactor = vtk.vtkActor()
        rupactor.SetMapper( rupmapper );
        rupactor.GetProperty().SetDiffuseColor( 1.,1.,1. )
        rupactor.GetProperty().SetOpacity(0.3)

        slreader = vtk.vtkDataSetReader()
        slreader.SetFileName( fnbase+"-psm-slip.vtk" )
        self.slreader = slreader
        
        arrow2 = vtk.vtkArrowSource()
        arrow2.SetShaftRadius(0.09)
        arrow2.SetTipLength(0.48)
        arrow2.SetTipRadius(0.15)
        
        arrows2 = vtk.vtkGlyph3D()
        arrows2.SetInput(slreader.GetOutput())
        arrows2.SetSource(arrow2.GetOutput())
        arrows2.SetVectorModeToUseVector()
        arrows2.SetScaleModeToScaleByVector()
        arrows2.OrientOn()
        
        slmapper = vtk.vtkPolyDataMapper()
        slmapper.SetInput( arrows2.GetOutput() )
        
        slactor = vtk.vtkActor()
        slactor.SetMapper( slmapper );
        slactor.GetProperty().SetDiffuseColor( 1.,1.,0. )
        
        
        reader = vtk.vtkParticleReader()
        reader.SetFileName( fnbase+"-dsm.bin" )
        reader.SetDataByteOrderToLittleEndian()
        self.reader = reader
        
        nxyplots = 3
        self.creaders = []
        self.ned = ['n','e','d']
        for i in range(nxyplots):
            newreader = vtk.vtkParticleReader()
            newreader.SetFileName( fnbase+"-1-"+self.ned[i]+".bin" )
            newreader.SetDataByteOrderToLittleEndian()
            self.creaders.append( newreader )
            
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput( reader.GetOutput() )
        mapper.SetScalarRange(-1., 1.);

        actor = vtk.vtkActor()
        actor.SetMapper( mapper );
        actor.GetProperty().SetPointSize(5)
        
        plane = vtk.vtkPlaneSource()
        plane.SetOrigin(0.,0.,0.)
        plane.SetPoint1(10000.,0.,0.)
        plane.SetPoint2(0.,10000.,0.)
        plane.SetCenter(0.,0.,0.)
        
        planeMapper = vtk.vtkPolyDataMapper()
        planeMapper.SetInput(plane.GetOutput())
        planeActor = vtk.vtkActor()
        planeActor.SetMapper(planeMapper)
        planeActor.GetProperty().SetOpacity( 0.2 )
        planeActor.GetProperty().SetDiffuseColor( 1.,1.,1. )
        
        #outlineData = vtk.vtkOutlineFilter()
        #outlineData.SetInput(fpreader.GetOutput())
        #mapOutline = vtk.vtkPolyDataMapper()
        #mapOutline.SetInput(outlineData.GetOutput())
        #outline = vtk.vtkActor()
        #outline.SetMapper(mapOutline)
        #outline.GetProperty().SetColor(1.,1.,1.)
        #outline.GetProperty().SetLineWidth(2.)

        camera = vtk.vtkCamera()
        camera.SetClippingRange(10., 500000.)
        camera.SetFocalPoint(0.,0.,0.)
        camera.SetPosition(-30000.,0.,-30000.)
        camera.SetViewUp(0.,0.,-1.)
        
        xyplots = []
        rangey = []
        rangex = []
        for i in range(nxyplots):
            xyplot = vtk.vtkXYPlotActor()
            xyplot.GetAxisLabelTextProperty().ItalicOff()
            xyplot.GetAxisLabelTextProperty().BoldOn()
            xyplot.GetAxisLabelTextProperty().SetFontFamilyToArial()
            xyplot.GetAxisLabelTextProperty().SetFontSize(5)
            xyplot.SetLabelFormat("%g")
            xyplot.AddInput( self.creaders[i].GetOutput() )
            xyplot.SetXValuesToValue()
            xyplot.SetTitle("")
            xyplot.SetXTitle("t")
            xyplot.SetYTitle("")
            xyplot.GetPositionCoordinate().SetValue(0.0, 1.-(i+1)*1./3., 0)
            xyplot.GetPosition2Coordinate().SetValue(1.0, 0.33, 0) #relative to Position
            xyplots.append( xyplot )
            ry = [0,0]
            rangey.append( ry )
            rx = [1000000,0]
            rangex.append( rx )
            
        self.xyplots = xyplots
        self.rangey = rangey
        self.rangex = rangex
            
        ren = vtk.vtkRenderer()
        ren.SetActiveCamera(camera)

        ren.AddActor(ceactor)
        ren.AddActor(fpactor)
        ren.AddActor(fpactor2)
        ren.AddActor(slactor)
        ren.AddActor(rupactor)
        ren.AddActor(actor)
        ren.AddActor(planeActor)
        #ren.AddActor(outline)
        ren.SetViewport(0, 0, .5, 1)

        ren2d = vtk.vtkRenderer()
        ren2d.SetViewport(0.5, 0.0, 1.0, 1.0)
        for i in range(nxyplots):
            ren2d.AddActor(xyplots[i])

        self.nxyplots = nxyplots
            
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)
        renWin.AddRenderer(ren2d)
        
        renWin.Render()
        self.renwin = renWin
        self.reader = reader
        
        
    def modified(self):
        self.reader.Modified()
        self.fpreader.Modified()
        self.cereader.Modified()
        self.slreader.Modified()
        self.rupreader.Modified()
        for i in range(self.nxyplots):
            self.creaders[i].Modified()
            ry = [0,0]
            rx = [1000000,0,0,0,0,0]
            self.creaders[i].GetOutput().GetScalarRange( ry )
            self.creaders[i].GetOutput().GetBounds( rx )
            yabsmax = max(max( abs( self.rangey[i][0]), abs(ry[0]) ),
                          max( abs( self.rangey[i][1]), abs(ry[1]) ) )
            self.rangey[i][0] = -yabsmax
            self.rangey[i][1] = yabsmax
            self.rangex[i][0] = min( self.rangex[i][0], rx[0] )
            self.rangex[i][1] = max( self.rangex[i][1], rx[1] )
            self.xyplots[i].SetYRange( self.rangey[i] )
            self.xyplots[i].SetXRange( self.rangex[i] )
        
        self.renwin.Render()

    def rerender(self):
        self.renwin.Render()
        
 
def main(args):
    MyApp.app=MyApp(args)
    win=MyWindow(args)
    win.show()
    MyApp.app.setwin( win )
    MyApp.app.setvis( win.vis )
    
    MyApp.app.connect(MyApp.app, SIGNAL("lastWindowClosed()"),
                MyApp.app.myquit)
    signal.signal(signal.SIGINT, MyApp.app.timetogo)
    
    MyApp.app.exec_loop()

    
    
if __name__=="__main__":
    main(sys.argv)
    