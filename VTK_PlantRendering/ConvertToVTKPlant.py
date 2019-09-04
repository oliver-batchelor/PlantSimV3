import vtk
from vtk.util.numpy_support import vtk_to_numpy
import geoplantrep as PG
import numpy as np


POT_SIZE = 0.08


class vtkPlantData():
    """Class containing functions needed to convert geoplantrep. into vtk data"""
    def __init__(self, plant_data_representation):
        self.plant_data = plant_data_representation

        self.StemPolydataList = []
        self.StemMapperList = []
        self.StemActorList = []
        self.MeshPolydataList = []
        self.MeshMapperList = []
        self.MeshActorList = []

        self.BackgroundOriPoses = []
        self.BackgroundPolydataList = []
        self.BackgroundMapperList = []
        self.BackgroundActorList = []
        self.vtkObjects = []

    def ConstructVTKStem(self):
        """Build vtk spline tubes around stem points"""
        for stem_idx, points in enumerate(self.plant_data.tubePntSets):
            rads = self.plant_data.tubeRadSets[stem_idx]
            cols = self.plant_data.tubeColSets[stem_idx]
            # Fit a spline to the points
            spline_points = vtk.vtkPoints()
            for point in points:
                spline_points.InsertNextPoint(point)

            spline = vtk.vtkParametricSpline()
            spline.SetPoints(spline_points)
            functionSource = vtk.vtkParametricFunctionSource()
            functionSource.SetParametricFunction(spline)
            functionSource.SetUResolution(10 * spline_points.GetNumberOfPoints())
            functionSource.Update()

            # Interpolate the scalars
            interpolatedRadius = vtk.vtkTupleInterpolator()
            interpolatedColor = vtk.vtkTupleInterpolator()
            interpolatedColor.SetInterpolationTypeToLinear()
            interpolatedColor.SetNumberOfComponents(3)
            interpolatedRadius.SetInterpolationTypeToLinear()
            interpolatedRadius.SetNumberOfComponents(1)
            for idx in range(len(points)):
                interpolatedRadius.AddTuple(idx, tuple([rads[idx]]))
                interpolatedColor.AddTuple(idx, cols[idx])

            # Generate the radius and color scalars
            tubeColors = vtk.vtkUnsignedCharArray()
            tubeColors.SetNumberOfComponents(3)
            tubeRadius = vtk.vtkDoubleArray()
            n = functionSource.GetOutput().GetNumberOfPoints()
            tubeColors.SetNumberOfTuples(n)
            tubeRadius.SetNumberOfTuples(n)
            tubeRadius.SetName("TubeRadius")
            tubeColors.SetName("Colors")
            cMin = interpolatedColor.GetMinimumT()
            cMax = interpolatedColor.GetMaximumT()
            tMin = interpolatedRadius.GetMinimumT()
            tMax = interpolatedRadius.GetMaximumT()
            r = [0.1]
            color_tuple = [1, 1, 1]
            for i in range(n):
                t = (tMax - tMin) / (n - 1) * i + tMin
                c = (cMax - cMin) / (n - 1) * i + cMin
                interpolatedRadius.InterpolateTuple(t, r)
                interpolatedColor.InterpolateTuple(c, color_tuple)
                tubeRadius.SetTuple1(i, r[0])
                tubeColors.SetTuple3(i, round(color_tuple[0]),
                                     round(color_tuple[1]),
                                     round(color_tuple[2]))

            # Add the scalars to the polydata
            tubePolyData = vtk.vtkPolyData()
            tubePolyData = functionSource.GetOutput()
            tubePolyData.GetPointData().AddArray(tubeRadius)
            tubePolyData.GetPointData().AddArray(tubeColors)
            tubePolyData.GetPointData().SetActiveScalars("TubeRadius")

            # Create the tubes and spheres
            tuber = vtk.vtkTubeFilter()
            tuber.SetNumberOfSides(round(np.random.rand()*4) + 6)
            tuber.SetInputData(tubePolyData)
            tuber.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
            tuber.CappingOn()
            tuber.Update()

            smoothing_filter = vtk.vtkSmoothPolyDataFilter()
            smoothing_filter.SetInputConnection(tuber.GetOutputPort())
            smoothing_filter.FeatureEdgeSmoothingOn()
            smoothing_filter.BoundarySmoothingOn()
            smoothing_filter.Update()

            self.vtkObjects.extend(
                [spline_points, spline, functionSource, interpolatedRadius, tubeColors, tubeRadius,
                 tubePolyData, tuber])
            self.StemPolydataList.append(smoothing_filter.GetOutput())

    def ConstructVTKMesh(self):
        """Construct polydata triangle mesh objects from data"""
        for set_indx, set_points in enumerate(self.plant_data.triMeshPntSets):
            mesh_points = vtk.vtkPoints()
            mesh_triangles = vtk.vtkCellArray()
            set_polydata = vtk.vtkPolyData()
            set_cols = vtk.vtkUnsignedCharArray()
            set_cols.SetNumberOfComponents(3)
            set_cols.SetName("Colors")

            for point_idx, point in enumerate(set_points):
                mesh_points.InsertNextPoint(point)
                set_cols.InsertNextTuple((round(255 * self.plant_data.triMeshColSets[set_indx][point_idx][0]),
                                          round(255 * self.plant_data.triMeshColSets[set_indx][point_idx][1]),
                                          round(255 * self.plant_data.triMeshColSets[set_indx][point_idx][2])))
            for indx_tuple in self.plant_data.triMeshPntIndxSets[set_indx]:
                triangle = vtk.vtkTriangle()
                triangle.GetPointIds().SetId(0, indx_tuple[0])
                triangle.GetPointIds().SetId(1, indx_tuple[1])
                triangle.GetPointIds().SetId(2, indx_tuple[2])
                mesh_triangles.InsertNextCell(triangle)

            set_polydata.SetPoints(mesh_points)
            set_polydata.SetPolys(mesh_triangles)
            set_polydata.GetPointData().SetScalars(set_cols)
            self.MeshPolydataList.append(set_polydata)
            self.vtkObjects.extend([mesh_points, mesh_triangles, set_cols])

    def CreatePolyDataMappers(self, disp_pot=False):
        """Sets up all the mappers to map stem polydata objects to the stemactors"""
        for stem_data in self.StemPolydataList:
            poly_mapper = vtk.vtkPolyDataMapper()
            poly_mapper.ScalarVisibilityOn()
            poly_mapper.SetScalarModeToUsePointFieldData()
            poly_mapper.SelectColorArray("Colors")
            poly_mapper.SetInputData(stem_data)
            poly_mapper.SetScalarRange(stem_data.GetScalarRange())
            self.StemMapperList.append(poly_mapper)
        for set_idx, mesh_data in enumerate(self.MeshPolydataList):
            poly_mapper = vtk.vtkPolyDataMapper()
            poly_mapper.SetInputData(mesh_data)
            poly_mapper.SetScalarModeToUsePointFieldData()
            poly_mapper.SetScalarRange(0, 255)
            poly_mapper.SelectColorArray("Colors")
            self.MeshMapperList.append(poly_mapper)
        if disp_pot:
            plant_bag = vtk.vtkCubeSource()
            plant_bag.SetCenter(0, -POT_SIZE / 2, 0)
            plant_bag.SetXLength(POT_SIZE)
            plant_bag.SetYLength(POT_SIZE)
            plant_bag.SetZLength(POT_SIZE)
            plant_bag.Update()
            plant_bag_mapper = vtk.vtkPolyDataMapper()
            plant_bag_mapper.SetInputData(plant_bag.GetOutput())
            self.BackgroundPolydataList.append(plant_bag)
            self.BackgroundMapperList.append(plant_bag_mapper)


    def CreatePolyDataActors(self):
        """Creates and adds actors for polydata mappers to renderer"""
        for idx, mapper in enumerate(self.StemMapperList):
            poly_actor = vtk.vtkActor()
            poly_actor.SetMapper(mapper)
            self.StemActorList.append(poly_actor)
        for idx, mapper in enumerate(self.MeshMapperList):
            poly_actor = vtk.vtkActor()
            poly_actor.SetMapper(mapper)
            self.MeshActorList.append(poly_actor)
        for idx, mapper in enumerate(self.BackgroundMapperList):
            poly_actor = vtk.vtkActor()
            poly_actor.SetMapper(mapper)
            self.BackgroundActorList.append(poly_actor)


    def SetActorPostions(self):
        """Set all actor origin positions on the fly based on lower stem segment vectors"""
        stem_end_positions = np.zeros((self.plant_data.numTubeSets, 3))
        for idx, actor in enumerate(self.StemActorList):
            actor_origin_pos = [0, 0, 0]
            # Sum all stem vectors of lower connections
            prev_seg_idx = self.plant_data.tubeConnIdxs[idx]
            while (prev_seg_idx != -1):
                actor_origin_pos += self.plant_data.tubePntSets[prev_seg_idx][-1]
                prev_seg_idx = self.plant_data.tubeConnIdxs[prev_seg_idx]
            actor.SetPosition(actor_origin_pos)
            stem_end_positions[idx] = actor_origin_pos + self.plant_data.tubePntSets[idx][-1]
        for idx, actor in enumerate(self.MeshActorList):
            leaf_origin_pos = stem_end_positions[self.plant_data.triMeshConnIdxs[idx]]
            actor.SetPosition(leaf_origin_pos)