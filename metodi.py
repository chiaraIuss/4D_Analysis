import sys,math,os
import vtk
import MC_methods as MC
from vmtk import vtkvmtk
from vmtk import pypes
import string
import pypescript2 as pypescript
import numpy as np  
import openFile as op
import extractlargecrosssection as ES

def ComputeCenterline(fn_root, fn_line):
    
#    str0='vmtkcenterlines -ifile '+fn_root+' -endpoints 1' 
#    str1=' -resampling 1 -resamplingstep 0.5 -ofile '+fn_line
    str0='vmtkcenterlines -ifile '+fn_root
    str1=' -costfunction 1/R^2 -endpoints 1 -resampling 1 -resamplingstep 1.0 --pipe vmtkcenterlinegeometry  -smoothing 1 -iterations 30 -factor 0.01 -outputsmoothed 1 -ofile '+fn_line
    myArguments=str0+str1
    myPype=pypes.PypeRun(myArguments)
    Surface=MC.ReadPolyData(fn_root)
    actor_surface=MC.CreateActor(Surface)
    prop_surface = actor_surface.GetProperty()
    prop_surface.SetColor(1.0, 0.0, 0.0)
    prop_surface.SetOpacity(0.3)
    
    Centerline=MC.ReadPolyData(fn_line)
    actor_line=MC.CreateActor(Centerline)
    renderer=vtk.vtkRenderer()
    renderer.AddActor(actor_surface)
    renderer.AddActor(actor_line)
    
    renderWindow=vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor=vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindow.SetInteractor(renderWindowInteractor)
    interactorStyle=vtk.vtkInteractorStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(interactorStyle)
    renderWindowInteractor.Initialize()
    renderWindow.Render()
    renderWindowInteractor.Start()  

def BranchExtractor(fn_line, fn_out):

    str0='vmtkbranchextractor -ifile '+fn_line+' -radiusarray MaximumInscribedSphereRadius -ofile '+fn_out
    myArguments=str0
    myPype=pypes.PypeRun(myArguments)

def BranchClipper(fn_root, fn_out, fn_clip):

    str0='vmtkbranchclipper -ifile '+fn_root+' -centerlinesfile '+fn_out+' -groupidsarray GroupIds -radiusarray MaximumInscribedSphereRadius -blankingarray Blanking -ofile '+fn_clip
    myArguments=str0
    myPype=pypes.PypeRun(myArguments)

def BranchClipper2(fn_root, fn_out, fn_clip, Ids):
       
        #fn_clip=base+'_branch_'+str(Ids)+'.vtp'
        #fn_cl_clip=base+'_cl_'+str(Ids)+'.vtp'
    str0='vmtkbranchclipper -ifile '+fn_root
    str1=' -centerlinesfile '+fn_out
    str2=' -groupidsarray GroupIds -blankingarray Blanking  -radiusarray MaximumInscribedSphereRadius -groupids '+str(Ids)
    str3=' -ofile '+fn_clip
    myArguments=str0+str1+str2+str3
    myPype=pypes.PypeRun(myArguments)
    
def Treshold(i, Branches,  fn_cl_treshold):
    
    thresh = vtk.vtkThreshold()
    thresh.ThresholdBetween(i, i)
    #thresh.SetAttributeModeToUseCellData()
    ###         selectCells.SetInputArrayToProcess( 
    ###            0, 
    ###            0, 
    ###            0,
    ###            1, // POINTS = 0, CELLS = 1, NONE = 2, POINTS_THEN_CELLS = 3, VERTICES = 4, EDGES = 5, ROWS = 6
    ###            0  // SCALARS = 0, VECTORS = 1, NORMALS = 2, TCOORDS = 3, TENSORS = 4, GLOBALIDS = 5, PEDIGREEIDS = 6, EDGEFLAG = 7
    ###            )
    thresh.SetInputArrayToProcess(0, 0, 0, 1, 0 )
    thresh.SetInput(Branches)
    thresh.Update()
    #print thresh.GetOutput()
    extract = vtk.vtkGeometryFilter() 
    extract.SetInput(thresh.GetOutput())
    extract.Update()
    fn_threshold=MC.WritePolyData(extract.GetOutput(), fn_cl_treshold)
    return fn_threshold
    #cl_threshold=MC.ReadPolyData(fn_threshold)
    
def CenterlineGeometry(cl_threshold, fn_cl_geometry):
    
    centerlineGeometry = vtkvmtk.vtkvmtkCenterlineGeometry()
    ##centerlineGeometry.SetInput(cl_clip)
    centerlineGeometry.SetInput(cl_threshold)
    centerlineGeometry.SetLengthArrayName('Length')
    centerlineGeometry.SetCurvatureArrayName('Curvature')
    centerlineGeometry.SetTorsionArrayName('Torsion')
    centerlineGeometry.SetTortuosityArrayName('Tortuosity')
    centerlineGeometry.SetFrenetTangentArrayName('FrenetTangent')
    centerlineGeometry.SetFrenetNormalArrayName('FrenetNormal')
    centerlineGeometry.SetFrenetBinormalArrayName('FrenetBinormal')
    centerlineGeometry.SetLineSmoothing(1)
    centerlineGeometry.SetSmoothingFactor(1.0)
    centerlineGeometry.Update()
    centerline_clip = centerlineGeometry.GetOutput()
    MC.WritePolyData(centerline_clip, fn_cl_geometry)
        
def SurfaceResampling(fn_section,  fn_section_remeshed):
    
    str0='vmtksurfaceremeshing -ifile '+fn_section+' -elementsizemode edgelength -ofile '+fn_section_remeshed
    myArguments=str0
    myPype=pypes.PypeRun(myArguments)
    
    
