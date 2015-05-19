import sys,math,os
import vtk
import MC_methods as MC
from vmtk import vtkvmtk
from vmtk import pypes
import string
import pypescript2 as pypescript
import numpy as np  
import openFile as op
#import extractlargecrosssection as ES
import metodi as M

def CutAortaWithPlane_Spline(surface,origin,normal):

	plane = vtk.vtkPlane()
	plane.SetOrigin(origin)
	plane.SetNormal(normal)

	cutter = vtk.vtkCutter()
	cutter.SetInput(surface)
	cutter.SetCutFunction(plane)
	cutter.Update()
	
	cleaner = vtk.vtkCleanPolyData()
	cleaner.SetInput(cutter.GetOutput())
	cleaner.Update()
	
	stripper = vtk.vtkStripper()
	stripper.SetInput(cleaner.GetOutput())
	stripper.Update()

	splineFilter = vtk.vtkSplineFilter()
	splineFilter.SetInput(stripper.GetOutput())
	splineFilter.Update()
        
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

	#reference_contour = vtk.vtkPolyData()
	#reference_contour.DeepCopy(splineFilter.GetOutput())
	#return reference_contour

def ComputeBranchSectionArea(branchSection):
	
	branchSection.BuildCells()
	if (branchSection.GetNumberOfCells() == 0):
		return 0.0;
	
	sectionPolygon = branchSection.GetCell(0)
	trianglePointIds = vtk.vtkIdList()
	sectionPolygon.Triangulate(trianglePointIds)
	numberOfTriangles = trianglePointIds.GetNumberOfIds() / 3
	polygonArea = 0.0
	for i in range(0,numberOfTriangles,1):
		
		pointId0 = trianglePointIds.GetId(3*i)
		pointId1 = trianglePointIds.GetId(3*i+1)
		pointId2 = trianglePointIds.GetId(3*i+2)
		point0 = [0.0,0.0,0.0]
		point1 = [0.0,0.0,0.0]
		point2 = [0.0,0.0,0.0]
		sectionPolygon.GetPoints().GetPoint(pointId0,point0)
		sectionPolygon.GetPoints().GetPoint(pointId1,point1)
		sectionPolygon.GetPoints().GetPoint(pointId2,point2)
		tri = vtk.vtkTriangle()
		triangleArea = tri.TriangleArea(point0,point1,point2)
		polygonArea += triangleArea

	#trianglePointIds.Delete()
	return polygonArea


def splineInterp(fileIn):

	cleaner = vtk.vtkCleanPolyData()
	cleaner.SetInput(fileIn)
	cleaner.Update()
	
	stripper = vtk.vtkStripper()
	stripper.SetInput(cleaner.GetOutput())
	stripper.Update()

	splineFilter = vtk.vtkSplineFilter()
	splineFilter.SetInput(stripper.GetOutput())
	splineFilter.Update()

	reference_contour = vtk.vtkPolyData()
	reference_contour.DeepCopy(splineFilter.GetOutput())
	return reference_contour


if __name__ == "__main__":
    app = op.MyFrame()
    app.mainloop()
    fn_root = app.fname

if fn_root.split('.')[-1] == 'stl':
    fn_surf = fn_root.split('.')
    base = fn_surf[0]
    surf = MC.ReadSTL(fn_root)
    fn = MC.WritePolyData(surf, base+'.vtp')

elif fn_root.split('.')[-1] == 'vtp':
    fn_surf = fn_root.split('.')
    base = fn_surf[0]
    surf = MC.ReadPolyData(fn_root)

#############
## CENTERLINE
file_surface_remesh = base+'.vtp'
print file_surface_remesh
file_line = base+'_cl.vtp'
file_branches = base+'_branches.vtp'   
input = 0

## CENTERLINE -> PRESENTE
## VISUALIZZO
if os.path.exists(file_line) == True:

    actor_surface = MC.CreateActor(MC.ReadPolyData(file_surface_remesh))
    prop_surface = actor_surface.GetProperty()
    prop_surface.SetColor(1.0, 1.0, 1.0)
    prop_surface.SetOpacity(0.5)
    
    Centerline = MC.ReadPolyData(file_line)
    actor_line = MC.CreateActor(Centerline)
    prop_line = actor_line.GetProperty()
    prop_line.SetLineWidth(1.5)
    prop_line.SetColor(1.5, 1.0, 0.0)
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor_surface)
    renderer.AddActor(actor_line)
    
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindow.SetInteractor(renderWindowInteractor)
    interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(interactorStyle)
    renderWindowInteractor.Initialize()
    renderWindow.Render()
    renderWindowInteractor.Start()
    ##################################
    ## CENTERLINE --> OK?
    input = raw_input('ti va bene? (Yes-->1,No-->0)')
    input = int(input)
    
M.ComputeCenterline(file_surface_remesh,file_line)

str0 = 'vmtkcenterlineresampling -ifile '+file_line
str1 = ' -length 1.0 -ofile '+file_line
myArguments = str0+str1
myPype=pypes.PypeRun(myArguments)

if os.path.exists(file_branches) == False:
    M.BranchExtractor(file_line, file_branches)
else:
    pass

file_surface_clip = base+'_surface_clip.vtp'
if os.path.exists(file_surface_clip) == False:
    M.BranchClipper(file_surface_remesh, file_branches, file_surface_clip)
else:
    pass

centerline_branches = MC.ReadPolyData(file_branches)
#='vmtkcenterlineresampling -ifile '+file_branches
#str1=' -length 1.0 -ofile '+file_branches
#myArguments=str0+str1
#myPype=pypes.PypeRun(myArguments)

GroupIdsArray=centerline_branches.GetCellData().GetArray('GroupIds')
GroupIds=[]

for i in range(GroupIdsArray.GetSize()):
    GroupIds.append(GroupIdsArray.GetComponent(i, 0))
    GroupIds[i]=int(GroupIds[i])

for el in GroupIds:
    if GroupIds.count(el)>=2:
        GroupIds.remove(el)   
GroupIds.sort()

for Ids in GroupIds:
        ##Scrivo i branch
        file_surface_groupids = base+'_surfaceGroupIds_'+str(Ids)+'.vtp'
        M.BranchClipper2(file_surface_remesh, file_branches, file_surface_groupids, Ids)
        
floatGroupIds=[]
for el in GroupIds:
    floatGroupIds.append(float(el))
    
centerlineCellData=centerline_branches.GetCellData()
centerlineCellData.SetActiveScalars('CenterlineIds')

file_centerline_brachio1 = base+'cl_bachio1.vtp'
M.Treshold(0, centerline_branches, file_centerline_brachio1)

file_centerline_brachio2 = base+'cl_bachio2.vtp'
M.Treshold(1, centerline_branches, file_centerline_brachio2)

file_centerline_LCCA = base+'cl_LCCA.vtp'
M.Treshold(2, centerline_branches, file_centerline_LCCA)

file_centerline_LSA = base+'cl_LSA.vtp'
M.Treshold(3, centerline_branches, file_centerline_LSA)

file_centerline_celiac = base+'cl_celiac.vtp'
M.Treshold(4, centerline_branches, file_centerline_celiac)

file_centerline_arch = base+'cl_arch.vtp'
M.Treshold(5, centerline_branches, file_centerline_arch)

M.CenterlineGeometry(MC.ReadPolyData(file_centerline_arch), file_centerline_arch)
centerline_arch = MC.ReadPolyData(file_centerline_arch)

centerlineArchPointSize = centerline_arch.GetNumberOfPoints()
fn_clip0=base+'_branch_0_'+str(int(i))+'.vtp'
surface=MC.ReadPolyData(file_surface_remesh)

area = []
file_area = []
for i in range(centerlineArchPointSize):
    
    origin=centerline_arch.GetPoint(i)
    normal=centerline_arch.GetPointData().GetArray('FrenetTangent').GetTuple3(i)

    section = CutAortaWithPlane_Spline(surface,origin,normal)
    file_section_arch = base+'_section_arch'+str(i)+'.vtp'
    
    MC.WritePolyData(section, file_section_arch)
    section_read = MC.ReadPolyData(file_section_arch)

    #triangleFilter = vtk.vtkTriangleFilter()
    #triangleFilter.SetInput(section_read)
    #triangleFilter.Update()
    
    #section_properties = vtk.vtkMassProperties()
    #section_properties.SetInputArrayToProcess(0, 0, 0, 3, 0 )
    #section_properties.SetInput(triangleFilter.GetOutput())
    #section_properties.Update()
    
    #area.append(section_properties.GetSurfaceArea()) 
    area.append(ComputeBranchSectionArea(section_read))
    
    file_area.append(str(ComputeBranchSectionArea(section_read)))
        
    #fn_section_prova=base+'_cl_0_section_'+str(int(i))+'_'+str(j)+'.vtp'
    #fn_section_remeshed=base+'_cl_0_section_'+str(int(i))+'_'+str(j)+'_remeshed.vtp'

out_file0 = open("/Volumes/UNTITLED/Utrecht_elongation/Healthy_Elongation_Study/Output/area.txt","w")
for el in file_area:
    out_file0.write(el+"\n")
out_file0.close()

diff_area=[]
for k in range(0,  30, 1):
    diff_area.append(area[k+1]-area[k])
print 'diff_area' , diff_area

STJ=[]
for j in range(len(diff_area)-1):
    if abs(diff_area[j]) < abs(diff_area[j+1]):
        STJ.append(j)
print STJ

distance_4_5=[]
index=0
for k in range(centerlineArchPointSize):
    
    celiac_centerline = MC.ReadPolyData(file_centerline_celiac)
    cp = celiac_centerline.GetPoint(k)
    ap = centerline_arch.GetPoint(k)
    
    distSquared = vtk.vtkMath.Distance2BetweenPoints(cp,ap)
    dist = math.sqrt(distSquared)
    distance_4_5.append(dist)
    
    if (distance_4_5[k] > 0.50 and distance_4_5[k] < 2.30):
        index = k;

#print distance_4_5
print index

centerlineCellData=centerline_branches.GetCellData()
centerlineCellData.SetActiveScalars('GroupIds')

GroupIdsDescending = 10
file_centerline_descending = base+'cl_descending.vtp'
M.Treshold(GroupIdsDescending, centerline_branches, file_centerline_descending)

cl_descending = MC.ReadPolyData(file_centerline_descending)
nOfPointsDescending = cl_descending.GetNumberOfPoints()
LSA = []
for i in range(nOfPointsDescending):
    
    dp = cl_descending.GetPoint(i)
    
    for j in range(centerlineArchPointSize):
        
        ap = centerline_arch.GetPoint(j)
        
        if math.sqrt(vtk.vtkMath.Distance2BetweenPoints(dp,ap)) < 1E-6:
            LSA.append(j)
            #print 'section of LSA is', j

        else:
            #print 'non raggiunta'
            continue

section = []

cl_brachio1 = MC.ReadPolyData(file_centerline_brachio1)
cl_brachio1.GetCellData().SetActiveScalars('TractIds')
file_centerline_asceding = base+'cl_ascending.vtp'
M.Treshold(0, cl_brachio1, file_centerline_asceding)

section.append(STJ[-1])
section.append(int(MC.ReadPolyData(file_centerline_asceding).GetNumberOfPoints())-10) #brachio
section.append(LSA[0]) # LSA
section.append(LSA[0]+100) #10 cm away from LSA
section.append(LSA[0]+200) #20 cm away from LSA
section.append(index) #celiac
print 'section at level STJ, brachio, LSA, 10cm, 20cm, celiac',section
print 'area at level STJ, brachio, LSA, 10cm, 20 cm, celiac', area[section[0]], area[section[1]], area[section[2]], area[section[3]], area[section[4]], area[section[5]]


L = [] # L1, L2, L3 ,L
L1 = 0.0
for k in range(section[0],section[1]-1):
    
    p1 = centerline_arch.GetPoint(k)
    p2 = centerline_arch.GetPoint(k+1)
    
    distSquared = vtk.vtkMath.Distance2BetweenPoints(p1,p2)
    dist = math.sqrt(distSquared)
    L1 = L1 + dist

L.append(L1)

L2 = 0.0
for k in range(section[1],section[2]-1):
    
    p1 = centerline_arch.GetPoint(k)
    p2 = centerline_arch.GetPoint(k+1)
    
    distSquared = vtk.vtkMath.Distance2BetweenPoints(p1,p2)
    dist = math.sqrt(distSquared)
    L2 = L2 + dist

L.append(L2)

L3 = 0.0
for k in range(section[2],section[-1]-1):
    
    p1 = centerline_arch.GetPoint(k)
    p2 = centerline_arch.GetPoint(k+1)
    
    distSquared = vtk.vtkMath.Distance2BetweenPoints(p1,p2)
    dist = math.sqrt(distSquared)
    L3 = L3 + dist

L.append(L3)
L.append(L1+L2+L3)
print 'length, L1,L2,L3, L',L


# min e max diameter
def ComputeBranchSectionShape(branchSection, center):
    
    sizeRange = [0.0,0.0]
    branchSection.BuildCells()
    if branchSection.GetNumberOfCells() == 0:
        sizeRange[0] = sizeRange[1] = 0.0
        return
    
    #sectionPolygon = vtk.vtkPolygon()
    #sectionPolygon.SetInput(branchSection.GetCell(0))
    #sectionPolygon.Update()
    #sectionPolygon = sectionPolygon.GetOutput()
    sectionPolygon = branchSection.GetCell(0)
    
    numberOfSectionPolygonPoints = sectionPolygon.GetNumberOfPoints()
    
    minDistance = 0.0
    maxDistance = 0.0

    for i in range(numberOfSectionPolygonPoints):
        point = [0.0, 0.0, 0.0]
        sectionPolygon.GetPoints().GetPoint(i,point)
        distance = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point,center))
        
        if (distance > maxDistance):
            maxDistance = distance
        
        if (distance < minDistance):
            minDistance = distance
    
    
    sizeRange[0] = minDistance
    sizeRange[1] = maxDistance

  #sectionShape = minDistance / maxDistance;

    return sizeRange


d_min2 = []
d_max2 = []
min_line = vtk.vtkPolyData()
max_line = vtk.vtkPolyData()
min_points = vtk.vtkPoints()
max_points = vtk.vtkPoints()
appendSection = vtk.vtkAppendPolyData()

for j in section:


    section_read=MC.ReadPolyData(base+'_section_arch'+str(j)+'.vtp')
    
    d_max_array = vtk.vtkFloatArray()
    d_max_array.SetName('D_max')
    d_max_array.SetNumberOfTuples(section_read.GetNumberOfCells())
    d_max_array.SetNumberOfComponents(1)
    d_max_array.FillComponent(0,0)
    section_read.GetCellData().AddArray(d_max_array)
    
    d_min_array = vtk.vtkFloatArray()
    d_min_array.SetName('D_min')
    d_min_array.SetNumberOfTuples(section_read.GetNumberOfCells())
    d_min_array.SetNumberOfComponents(1)
    d_min_array.FillComponent(0,0)
    section_read.GetCellData().AddArray(d_min_array)
    
    area_array = vtk.vtkFloatArray()
    area_array.SetName('Area')
    area_array.SetNumberOfTuples(section_read.GetNumberOfCells())
    area_array.SetNumberOfComponents(1)
    area_array.FillComponent(0,0)
    section_read.GetCellData().AddArray(area_array)
    area_array.SetTuple1(0,area[j])
    
    
    shape_array = vtk.vtkFloatArray()
    shape_array.SetName('sectionShape')
    shape_array.SetNumberOfTuples(section_read.GetNumberOfCells())
    shape_array.SetNumberOfComponents(1)
    shape_array.FillComponent(0,0)
    section_read.GetCellData().AddArray(shape_array)
    
        
    dist = []
    #dist2 =[]
    temp_i = []
    temp_k = []
    for i in range(section_read.GetNumberOfPoints()):
        tmp = []
        #dist2 =[]
        for k in range(section_read.GetNumberOfPoints()):
            
            distSquared = vtk.vtkMath.Distance2BetweenPoints(section_read.GetPoint(i),section_read.GetPoint(k))
            tmp.append(math.sqrt(distSquared))
        
        dist_max = max(tmp)
        dist.append(dist_max)
        id_max = tmp.index(dist_max)
        temp_i.append(i)
        temp_k.append(id_max)
            #dist2.append(min(tmp))
            #print dist
    dist_max = max(dist)
    dist_min = min(dist)
    id_max = dist.index(dist_max)
    id_min = dist.index(dist_min)
    
    min_points.InsertNextPoint(section_read.GetPoint(temp_i[id_min]))
    min_points.InsertNextPoint(section_read.GetPoint(temp_k[id_min]))
    max_points.InsertNextPoint(section_read.GetPoint(temp_i[id_max]))
    max_points.InsertNextPoint(section_read.GetPoint(temp_k[id_max]))
    
    d_max_array.SetTuple1(0,dist_max)
    d_min_array.SetTuple1(0,dist_min)
    shape_array.SetTuple1(0,(dist_min/dist_max))
    
    
    appendSection.AddInput(section_read)
    #appendSection.Update()
    MC.WritePolyData(section_read,base+'_section_arch'+str(j)+'.vtp')
    
    d_min2.append(min(dist))
    d_max2.append(max(dist))

min_line.SetPoints(min_points)
max_line.SetPoints(max_points)
MC.WritePolyData(min_line,base+'_min_line.vtp')
MC.WritePolyData(max_line,base+'_max_line.vtp')
MC.WritePolyData(appendSection.GetOutput(),base+'_Sections.vtp')

print 'min & max diameter second method', d_min2,  '\n', d_max2