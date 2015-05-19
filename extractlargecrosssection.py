#!/usr/bin/env python

import sys,math,os
import vtk
from vmtk import vtkvmtk
from vmtk import pypes

def ReadPolyData(filename):
   reader = vtk.vtkXMLPolyDataReader()
   reader.SetFileName(filename)
   reader.Update()
   return reader.GetOutput()

def WritePolyData(surface,filename):
   writer = vtk.vtkXMLPolyDataWriter()
   writer.SetInput(surface)
   writer.SetFileName(filename)
   writer.Write()

def GeneratePypeArgs(icenterlinefilename,ocenterlinefilename):
   pypeargs = 'vmtkcenterlineresampling -length 0.5 ' + \
              '-ifile %s ' % icenterlinefilename + \
	      '-ofile %s ' % ocenterlinefilename						
   return pypeargs

def CutSacWithPlane(surface,centerline):
   areas = []
   numberOfPoints = centerline.GetNumberOfPoints()
   tangentArray = centerline.GetPointData().GetArray('FrenetTangent')
   print tangentArray.GetRange()
   print numberOfPoints
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

def ExtractLargestSection(surface,centerline,id):
   tangentArray = centerline.GetPointData().GetArray('FrenetTangent')

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

   sectionCellArray.InsertCellPoint(0)	

   section.SetPoints(sectionPoints)
   section.SetPolys(sectionCellArray)

   sectionArea = vtkvmtk.vtkvmtkPolyDataBranchSections.ComputeBranchSectionArea(section)
 
   areaArray = vtk.vtkDoubleArray()
   areaArray.SetNumberOfComponents(1)
   areaArray.SetNumberOfTuples(1)
   areaArray.SetName('Area')
   areaArray.SetTuple1(0,sectionArea)

   section.GetCellData().AddArray(areaArray)

   return section

#
## -----------------------------------------------------------------------------------------
#
#
### Program:	extractlargestcrosssection.py	
### Language:	Python
### Date:	2012/02/27
### Application: Cerebral Aneurysms - Aneurysm Sac Morphology 
#
## -----------------------------------------------------------------------------------------
#
#print 'USAGE:'
#print '      ./extractlargestcrosssection inputfilesDirectory caseID'
#print ''
#
#inputfiledirectory = sys.argv[1]
#ID = sys.argv[2]
#
#print 'Inputfiles Directory   ', inputfiledirectory
#print 'case ID		      ', ID
#
## input filenames
#saccenterlinefilename       = inputfiledirectory + '/' + ID + '/' + ID + '_saccl.vtp' 
#sacsurfacefilename          = inputfiledirectory + '/' + ID + '/' + ID + '_sac.vtp'
#
## output filenames
#largestcrosssectionfilename = inputfiledirectory + '/' + ID + '/' + ID + '_largestcrosssection.vtp'
#
#pypeargs = GeneratePypeArgs(saccenterlinefilename,saccenterlinefilename)
#pipe = pypes.Pype()
#pipe.ExitOnError = 0
#pipe.Arguments = pypeargs.split()
#pipe.ParseArguments()
#pipe.Execute()
#
#sacCenterline = ReadPolyData(saccenterlinefilename)
#sacSurface = ReadPolyData(sacsurfacefilename)
#
#centerlineGeometry = vtkvmtk.vtkvmtkCenterlineGeometry()
#centerlineGeometry.SetInput(sacCenterline)
#centerlineGeometry.SetLengthArrayName('Length')
#centerlineGeometry.SetCurvatureArrayName('Curvature')
#centerlineGeometry.SetTorsionArrayName('Torsion')
#centerlineGeometry.SetTortuosityArrayName('Tortuosity')
#centerlineGeometry.SetFrenetTangentArrayName('FrenetTangent')
#centerlineGeometry.SetFrenetNormalArrayName('FrenetNormal')
#centerlineGeometry.SetFrenetBinormalArrayName('FrenetBinormal')
#centerlineGeometry.SetLineSmoothing(1)
#centerlineGeometry.SetSmoothingFactor(0.1)
#centerlineGeometry.Update()
#sacCenterline = centerlineGeometry.GetOutput()
#
#print 'Slicing sac surface'
#sectionAreas = CutSacWithPlane(sacSurface,sacCenterline)
#
#print ''
#
#maxArea = 0.0
#maxAreaId = 0
#for k in range(len(sectionAreas)):
#   if (sectionAreas[k]>maxArea):
#      maxArea = sectionAreas[k]
#      maxAreaId = k
#print 'Largest cross section id, area ', maxAreaId, maxArea
#
#largestCrossSection = ExtractLargestSection(sacSurface,sacCenterline,maxAreaId)
#WritePolyData(largestCrossSection,largestcrosssectionfilename)
#
#
