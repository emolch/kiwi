#!/usr/bin/env python

import vtk
import sys

reader = vtk.vtkParticleReader()
print sys.argv[1]
reader.SetFileName( sys.argv[1] )
reader.SetDataByteOrderToLittleEndian()

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
plane.SetXResolution(1)
plane.SetYResolution(1)

planeMapper = vtk.vtkPolyDataMapper()
planeMapper.SetInput(plane.GetOutput())
planeActor = vtk.vtkActor()
planeActor.SetMapper(planeMapper)
planeActor.GetProperty().SetOpacity( 0.3 )
planeActor.GetProperty().SetDiffuseColor( 1.,1.,1. )



camera = vtk.vtkCamera()
camera.SetClippingRange(10., 50000.)
camera.SetFocalPoint(0.,0.,0.)
camera.SetPosition(-30000.,0.,-30000.)
camera.SetViewUp(0.,0.,-1.)


# Create the usual rendering stuff.
ren = vtk.vtkRenderer()
ren.SetActiveCamera(camera)

renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(300, 150)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

#ren.SetBackground(.1, .2, .4)

ren.AddActor(actor)
ren.AddActor(planeActor)
iren.Initialize()
renWin.Render()
iren.Start()
