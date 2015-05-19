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
import metodi as M

file_section_arch = '/Users/chiara/Documents/Utrecht_elongation/A_2/0/0_2_cl_0_section_5_'+str(176)+'.vtp'
section = MC.ReadPolyData(file_section_arch)
file_cl_ascending = '/Users/chiara/Documents/Utrecht_elongation/A_2/0/0_2cl_arch5.vtp'
cl = MC.ReadPolyData(file_cl_ascending)
normal=cl.GetPointData().GetArray('FrenetTangent').GetTuple3(14)
center = [cl.GetPoint(176)[0],cl.GetPoint(176)[1],cl.GetPoint(176)[2]]

extrude = vtk.vtkLinearExtrusionFilter()
#extrude = vtk.vtkRotationalExtrusionFilter()
extrude.SetInput(section)
extrude.CappingOff()
extrude.SetScaleFactor(1)
extrude.SetExtrusionTypeToNormalExtrusion()
extrude.SetVector(normal[0], normal[1], normal[2])
#extrude.SetVector(-normal[0], -normal[1], -normal[2])
extrude2 = vtk.vtkLinearExtrusionFilter()
#extrude = vtk.vtkRotationalExtrusionFilter()
extrude2.SetInput(section)
extrude2.CappingOff()
extrude2.SetScaleFactor(1)
extrude2.SetExtrusionTypeToNormalExtrusion()
#extrude.SetVector(normal[0], normal[1], normal[2])
extrude2.SetVector(-normal[0], -normal[1], -normal[2])

appendFilter = vtk.vtkAppendPolyData()
appendFilter.AddInput(extrude.GetOutput());
appendFilter.AddInput(extrude2.GetOutput())
appendFilter.Update() 
cleanFilter = vtk.vtkCleanPolyData()
cleanFilter.SetInput(appendFilter.GetOutput())
cleanFilter.Update()
MC.WritePolyData(cleanFilter.GetOutput(), '/Volumes/UNTITLED/Utrecht_elongation/Healthy_Elongation_Study/Output/cylinder.vtp')

def ComputeBranchSectionShape(branchSection, center):

    branchSection.BuildCells()

    sectionPolygon = branchSection.GetCell(0)
    numberOfSectionPolygonPoints = sectionPolygon.GetNumberOfPoints()
    
    poliPoints = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    
    poliPoint_min = vtk.vtkPolyData()
    points_min = vtk.vtkPoints()
    poliPoint_max = vtk.vtkPolyData()
    points_max = vtk.vtkPoints()
    
    VTK_VMTK_LARGE_DOUBLE = 1.0E+32
    minDistance = VTK_VMTK_LARGE_DOUBLE
    maxDistance = 0.0
    minDistanceId = -1
    maxDistanceId = -1
    
    point=[0.0,0.0,0.0]
    for i in range(numberOfSectionPolygonPoints):

        sectionPolygon.GetPoints().GetPoint(i,point)
        distance = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point,center))
        if (distance > maxDistance):

            maxDistance = distance
            maxDistanceId = i


    point0=[0.0,0.0,0.0]
    point1=[0.0,0.0,0.0]
    planeNormal=[0.0,0.0,0.0]
    radialVector0=[0.0,0.0,0.0]
    radialVector1=[0.0,0.0,0.0]
    cross=[0.0,0.0,0.0]
    
    for i in range(numberOfSectionPolygonPoints):

        sectionPolygon.GetPoints().GetPoint(i,point0)
        sectionPolygon.GetPoints().GetPoint((i+numberOfSectionPolygonPoints/4)%numberOfSectionPolygonPoints,point1)
        radialVector0[0] = point0[0] - center[0]
        radialVector0[1] = point0[1] - center[1]
        radialVector0[2] = point0[2] - center[2]
        radialVector1[0] = point1[0] - center[0]
        radialVector1[1] = point1[1] - center[1]
        radialVector1[2] = point1[2] - center[2]
        vtk.vtkMath.Cross(point0,point1,cross)
        planeNormal[0] += cross[0]
        planeNormal[1] += cross[1]
        planeNormal[2] += cross[2]
        
    vtk.vtkMath.Normalize(planeNormal)
    d = []
    d1 = []
    d2 = []
    index_min = 0
    index_max = 0
    d_min = 50.0
    d_max = 0.0
    pointsIntersection = []
    for j in range(numberOfSectionPolygonPoints):
    
        minDistancePoint=[0.0,0.0,0.0]
        sectionPolygon.GetPoints().GetPoint(j,minDistancePoint)
        print 'mindistpoint', minDistancePoint
        distance1 = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(minDistancePoint,center))

        minDistanceNormal=[0.0,0.0,0.0]
        
        minDistanceNormal[0] = minDistancePoint[0] - center[0]
        minDistanceNormal[1] = minDistancePoint[1] - center[1]
        minDistanceNormal[2] = minDistancePoint[2] - center[2]
        vtk.vtkMath.Normalize(minDistanceNormal)

        minDistanceOppositePoint=[0.0,0.0,0.0]

        minDistanceOppositePoint[0] = center[0] - 2.0 * maxDistance * minDistanceNormal[0]
        minDistanceOppositePoint[1] = center[1] - 2.0 * maxDistance * minDistanceNormal[1]
        minDistanceOppositePoint[2] = center[2] - 2.0 * maxDistance * minDistanceNormal[2]
        obbTree = vtk.vtkOBBTree()
        #obbTree.SetDataSet(branchSection)
        obbTree.SetDataSet(extrude.GetOutput())
        obbTree.BuildLocator()
    
        pointsVTKintersection = vtk.vtkPoints()
        code = obbTree.IntersectWithLine(minDistancePoint, minDistanceOppositePoint, pointsVTKintersection, None)
        #code = obbTree.IntersectWithLine(center, minDistancePoint, pointsVTKintersection, None)
    
        pointsVTKIntersectionData = pointsVTKintersection.GetData()
        noPointsVTKIntersection = pointsVTKIntersectionData.GetNumberOfTuples()
        
        for idx in range(noPointsVTKIntersection):
            _tup = pointsVTKIntersectionData.GetTuple3(idx)
            pointsIntersection.append(_tup)
        print 'intersection point', pointsIntersection
        
        inter = [0.0,0.0,0.0]    
        inter[0]= pointsIntersection[-1][0]
        inter[1]= pointsIntersection[-1][1]
        inter[2]= pointsIntersection[-1][2]
        d.append(math.sqrt(vtk.vtkMath.Distance2BetweenPoints(minDistancePoint,inter)))

        if math.sqrt(vtk.vtkMath.Distance2BetweenPoints(minDistancePoint,inter)) < d_min:
            d_min = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(minDistancePoint,inter))
            index_min = j
            
        if  math.sqrt(vtk.vtkMath.Distance2BetweenPoints(minDistancePoint,inter)) > d_max:
            d_max = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(minDistancePoint,inter))
            index_max = j
            
        
        #points.InsertNextPoint(minDistancePoint)
        #points.InsertNextPoint(minDistanceOppositePoint)
        #points.InsertNextPoint(center)
        points.InsertNextPoint(inter)
        
    print 'index min',index_min
    print 'index_max',index_max
    poliPoints.SetPoints(points)
    #poliPoints.SetPoints(pointsVTKintersection)
    MC.WritePolyData(poliPoints, '/Volumes/UNTITLED/Utrecht_elongation/Healthy_Elongation_Study/Output/Pointsssss.vtp')
    
    sizeRange = [0.0,0.0]
    sizeRange[0] = min(d)
    sizeRange[1] = max(d)
    sectionShape = sizeRange[0] / sizeRange[1]
    #id_max = d.index(sizeRange[1])
    #id_min = d.index(sizeRange[0])
    
    point_temp=[0.0,0.0,0.0]
    sectionPolygon.GetPoints().GetPoint(index_min,point_temp)
    print '--------------',point_temp
    points_min.InsertNextPoint(point_temp)
    poliPoints.GetPoints().GetPoint(index_min,point_temp)
    print '--------------',point_temp
    points_min.InsertNextPoint(point_temp)
    sectionPolygon.GetPoints().GetPoint(index_max,point_temp)
    print '--------------',point_temp
    points_max.InsertNextPoint(point_temp)
    poliPoints.GetPoints().GetPoint(index_max,point_temp)
    print '--------------',point_temp
    points_max.InsertNextPoint(point_temp)
        
    poliPoint_min.SetPoints(points_min)
    poliPoint_min.Update()
    poliPoint_max.SetPoints(points_max)
    poliPoint_max.Update()
        
    MC.WritePolyData(poliPoint_min, '/Volumes/UNTITLED/Utrecht_elongation/Healthy_Elongation_Study/Output/PointMINDiam.vtp')
    MC.WritePolyData(poliPoint_max, '/Volumes/UNTITLED/Utrecht_elongation/Healthy_Elongation_Study/Output/PointMAXDiam.vtp')
    print 'max_d, min_d', max(d), min(d)
ComputeBranchSectionShape(section, center)   