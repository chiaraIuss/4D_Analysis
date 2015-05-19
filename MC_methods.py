#!/usr/bin/env python

import sys,math,os
import vtk
from vmtk import vtkvmtk
from vmtk import pypes

def annotatePick(object, event):
#    global picker, textActor, textMapper
#    if picker.GetCellId() < 0:
#        textActor.VisibilityOff() 
#    else:
        selPt = picker.GetSelectionPoint() 
        pickPos = picker.GetPickPosition() 
        textMapper.SetInput("(%.6f, %.6f, %.6f)"%pickPos) 
        textActor.SetPosition(selPt[:2]) 
        textActor.VisibilityOn()

def QuitRendererCallback( obj):
    PrintLog('Quit renderer')
    Renderer.RemoveActor(TextActor)
    renderWindowInteractor.ExitCallback()

def EnterTextInputMode(renderer, interactive=1):
    CurrentTextInput = ''
    renderer.AddActor(TextInputActor)
    renderer.RemoveActor(self.TextActor)
    UpdateTextInput()
    TextInputMode = 1
    render(interactive)

def AddKeyBinding(key, text, callback=None, group='1'):
    KeyBindings = {}
    KeyBindings[key] = {'text': text, 'callback': callback, 'group': group}
    
def ReadPolyData(filename):
   reader = vtk.vtkXMLPolyDataReader()
   reader.SetFileName(filename)
   reader.Update()
   return reader.GetOutput()

def ReadPolyDataDir(directoryname):
   reader = vtk.vtkXMLPolyDataReader()
   reader.SetFileName(directoryname)
   reader.Update()
   return reader.GetOutput()
   
def ReadSTL(filename):
   reader = vtk.vtkSTLReader()
   reader.SetFileName(filename)
   reader.Update()
   return reader.GetOutput()

def WritePolyData(surface,filename):
   writer = vtk.vtkXMLPolyDataWriter()
   writer.SetInput(surface)
   writer.SetFileName(filename)
   writer.Write()
   return filename

def CreateActor(oggetto):
    mapper=vtk.vtkPolyDataMapper()
    mapper.SetInput(oggetto)
    actor=vtk.vtkActor()
    actor.SetMapper(mapper)
    return actor    

def CreateRender(actor, actor1, actor2): #, centerline, actor_sphere, actor_bar):
    renderer=vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.AddActor(actor1)    
    renderer.AddActor(actor2)
    return renderer
    
def CreateRenderMulti(actor, actor1, actor2, actor3): #, centerline, actor_sphere, actor_bar):
    renderer=vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.AddActor(actor1)    
    renderer.AddActor(actor2)
    renderer.AddActor(actor3)
    return renderer

    
def CreateRenderWindow(renderer, x, y):
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetPosition(x, y) 
    renderWindow.AddRenderer(renderer)
    #renderWindow.AddRenderer(renderer1)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.LightFollowCameraOff()
    interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(interactorStyle)
    renderWindow.SetInteractor(renderWindowInteractor)
    
    renderWindowInteractor.Initialize()
    renderWindow.Render()
    renderWindowInteractor.Start()
    
def RenderWindowMultiple(ren):
    '''One render window, multiple viewports'''
    #iren_list = []
    renderWindow = vtk.vtkRenderWindow()
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.LightFollowCameraOff()
    interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(interactorStyle)
    renderWindow.SetInteractor(renderWindowInteractor)
    #renderWindowInteractor.SetRenderWindow(renderWindow)
    # Define viewport ranges
    xmins=[0,.5,0,.5]
    xmaxs=[0.5,1,0.5,1]
    ymins=[0,0,.5,.5]
    ymaxs=[0.5,0.5,1,1]
    for i in range(4):
        #ren = vtk.vtkRenderer()
        renderWindow.AddRenderer(ren)
        ren.SetViewport(xmins[i],ymins[i],xmaxs[i],ymaxs[i])
     
        ##Create a sphere
        #sphereSource = vtk.vtkSphereSource()
        #sphereSource.SetCenter(0.0, 0.0, 0.0)
        #sphereSource.SetRadius(5)
     #
        ##Create a mapper and actor
        #mapper = vtk.vtkPolyDataMapper()
        #mapper.SetInputConnection(sphereSource.GetOutputPort())
        #actor = vtk.vtkActor()
        #actor.SetMapper(mapper)
        #ren.AddActor(actor)
        #ren.ResetCamera()
    renderWindowInteractor.Initialize() 
    renderWindow.Render()
    renderWindow.SetWindowName('RW: Multiple ViewPorts')
    renderWindowInteractor.Start()
    

def CutAortaWithPlane1(surface,origin,normal):
   plane = vtk.vtkPlane()
   plane.SetOrigin(origin)
   plane.SetNormal(normal)
 
   cutter = vtk.vtkCutter()
   cutter.SetInput(surface)
   cutter.SetCutFunction(plane)
   cutter.Update()

   stripper = vtk.vtkStripper()
   stripper.SetInput(cutter.GetOutput())
   stripper.Update()

   connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
   connectivityFilter.SetInput(stripper.GetOutput())
   connectivityFilter.SetClosestPoint(origin)
   connectivityFilter.SetExtractionModeToClosestPointRegion()
   connectivityFilter.Update()

   contour = connectivityFilter.GetOutput()
   numberOfContourPoints = contour.GetNumberOfPoints()
 
   section = vtk.vtkPolyData()
   sectionPoints = vtk.vtkPoints()
   sectionCellArray = vtk.vtkCellArray() 
   sectionCellArray.InsertNextCell(numberOfContourPoints)
 
   for j in range(numberOfContourPoints):
     point = contour.GetPoint(j)

     sectionPoints.InsertNextPoint(point)
     sectionCellArray.InsertCellPoint(j)

   sectionCellArray.InsertCellPoint(0)	

   section.SetPoints(sectionPoints)
   section.SetPolys(sectionCellArray)

   return section

def CutAortaWithPlane(surface,centerline):
   areas = []
   numberOfPoints = centerline.GetNumberOfPoints()
   tangentArray = centerline.GetPointData().GetArray('FrenetTangent')
   print tangentArray.GetRange()

   print 'Id Area'
   for i in range(numberOfPoints):
      print i, 	   
      origin = centerline.GetPoint(i)
      normal = tangentArray.GetTuple3(i)

      plane = vtk.vtkPlane()
      plane.SetOrigin(origin)
      plane.SetNormal(normal)
     
      cutter = vtk.vtkCutter()
      cutter.SetInput(surface)
      cutter.SetCutFunction(plane)
      cutter.Update()

      stripper = vtk.vtkStripper()
      stripper.SetInput(cutter.GetOutput())
      stripper.Update()

      connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
      connectivityFilter.SetInput(stripper.GetOutput())
      connectivityFilter.SetClosestPoint(origin)
      connectivityFilter.SetExtractionModeToClosestPointRegion()
      connectivityFilter.Update()

      contour = connectivityFilter.GetOutput()
      numberOfContourPoints = contour.GetNumberOfPoints()
     
      section = vtk.vtkPolyData()
      sectionPoints = vtk.vtkPoints()
      sectionCellArray = vtk.vtkCellArray() 
      sectionCellArray.InsertNextCell(numberOfContourPoints)
     
      for j in range(numberOfContourPoints):
         point = contour.GetPoint(j)

         sectionPoints.InsertNextPoint(point)
         sectionCellArray.InsertCellPoint(j)
	
      sectionCellArray.InsertCellPoint(0)	

      section.SetPoints(sectionPoints)
      section.SetPolys(sectionCellArray)

      sectionArea = vtkvmtk.vtkvmtkPolyDataBranchSections.ComputeBranchSectionArea(section)
      print sectionArea

      areas.append(sectionArea)

   return areas    

def ExtractLargestSection(surface,centerline,id_array):
   tangentArray = centerline.GetPointData().GetArray('FrenetTangent')
   for i in id_array:
      id=i
      print id 
      origin = centerline.GetPoint(id)
      normal = tangentArray.GetTuple3(id)
      
      plane = vtk.vtkPlane()
      plane.SetOrigin(origin)
      plane.SetNormal(normal)
      
      cutter = vtk.vtkCutter()
      cutter.SetInput(surface)
      cutter.SetCutFunction(plane)
      cutter.Update()
      
      stripper = vtk.vtkStripper()
      stripper.SetInput(cutter.GetOutput())
      stripper.Update()
   
      connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
      connectivityFilter.SetInput(stripper.GetOutput())
      connectivityFilter.SetClosestPoint(origin)
      connectivityFilter.SetExtractionModeToClosestPointRegion()
      connectivityFilter.Update()
   
      contour = connectivityFilter.GetOutput()
      numberOfContourPoints = contour.GetNumberOfPoints()
        
      section = vtk.vtkPolyData()
      sectionPoints = vtk.vtkPoints()
      sectionCellArray = vtk.vtkCellArray() 

   
      sectionCellArray.InsertNextCell(numberOfContourPoints)
   
      for j in range(numberOfContourPoints):
          point = contour.GetPoint(j)
          sectionPoints.InsertNextPoint(point)
          sectionCellArray.InsertCellPoint(j)

      for j in range(numberOfContourPoints):
          point = contour.GetPoint(j)
          d0=point
          print dist
   
      sectionCellArray.InsertCellPoint(0)	
      section.SetPoints(sectionPoints)
      section.SetPolys(sectionCellArray)
      sectionArea = vtkvmtk.vtkvmtkPolyDataBranchSections.ComputeBranchSectionArea(section)
      areaArray = vtk.vtkDoubleArray()
      areaArray.SetNumberOfComponents(1)
      areaArray.SetNumberOfTuples(1)
      areaArray.SetName('Area')
      areaArray.SetTuple1(0,sectionArea)
      
      meanArray = vtk.vtkDoubleArray()
      meanArray.SetNumberOfComponents(1)
      meanArray.SetNumberOfTuples(1)
      meanArray.SetName('Media')
      meanArray.SetTuple1(0,10.)
      
      section.GetCellData().AddArray(areaArray)
      section.GetCellData().AddArray(meanArray)
   return section
