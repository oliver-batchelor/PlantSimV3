import geoplantrep as PG
import numpy as np
import math
import vtk


LEAF_G_MEAN = 0.8
LEAF_G_VAR = 0.2
LEAF_R_MEAN = 0.4
LEAF_R_VAR = 0.1
LEAF_B_MEAN = 0.3
LEAF_B_VAR = 0.1

LEAF_MAX_SPLINE_SEGS = 5
LEAF_WIDTH_MEAN = 0.08
LEAF_WIDTH_VAR = 0.03
LEAF_LENGTH_MEAN = 0.14
LEAF_LENGTH_VAR = 0.04
LEAF_DEPTH_MEAN = 0.03
LEAF_DEPTH_VAR = 0.003
LEAF_X_ROT_VAR = 10
LEAF_Z_ROT_VAR = 10
LEAF_Y_ROT_VAR = 270

LEAF_UPMUL_MEAN = 1.2
LEAF_UPMUL_VAR = 0.6
LEAF_HORIZ_VAR = 0.06
LEAF_OUTRAD_MUL = 0.01

########## 2 ###########
MAX_LEAF_DEPTH = 0.15
LEAF_CURVE_MEAN = 0
LEAF_CURVE_RELMAX = 0.5
LEAF_Z_CURl_MUL = 0.12
LEAF_X_CURL_MUL = 0.04

LEAF_Z_VAR = 2
LEAF_NORMAL_VERTICAL_BIAS = 2
LEAF_KDTREE_POINTS = 50
MAX_BOUND_ITTS = 10
LEAF_SIZE_MULTIPLIER = 0.995


#TODO: make fully vtk free leaf generator
#TODO: make fully vtk free triangle intersection detector


def GenRandMeshLeaves_1(geo_plant_rep, end_seg_indxs ):
    """Generates randomised leaves given randomisation parameters"""

    for stem_seg_idx in end_seg_indxs:
        leafTriPoints = []
        leafTriIndices = []
        leafTriCols = []

        ################# Randomisations ##########################
        leaf_length = np.random.normal(loc=LEAF_LENGTH_MEAN, scale=LEAF_LENGTH_VAR)
        leaf_width = np.random.normal(loc=LEAF_WIDTH_MEAN, scale=LEAF_WIDTH_VAR)
        leaf_depth = np.random.normal(loc=LEAF_DEPTH_MEAN, scale=LEAF_DEPTH_VAR)

        up_bias = np.random.normal(loc=LEAF_UPMUL_MEAN, scale=LEAF_UPMUL_VAR)
        dx = np.random.normal(loc=LEAF_OUTRAD_MUL * np.sign(geo_plant_rep.tubePntSets[stem_seg_idx][-1][0]),
                              scale=LEAF_HORIZ_VAR)
        dz = np.random.normal(loc=LEAF_OUTRAD_MUL * np.sign(geo_plant_rep.tubePntSets[stem_seg_idx][-1][2]),
                              scale=LEAF_HORIZ_VAR)
        dy = np.sqrt(dx ** 2 + dz ** 2) * up_bias
        ###########################################################
        leaf_up_vec = [dx, dy, dz] / np.linalg.norm([dx, dy, dz], 2)

        point_spacing = leaf_width / 9
        z_res = max(math.trunc(leaf_length * 9 / leaf_width), 9)

        leafTriPoints.append([0.002, 0, 0.002])

        for z in range(1, z_res):
            z_prob_adj = 0.5 - abs(z_res / 2 - z) / z_res
            z_depth_mean = - 0.12 * (z / z_res) ** 2
            for x in range(-math.floor(9 / 2), math.floor(9 / 2) + 1):
                x_prob_adj = 0.5 - abs(x) / z_res

                ################# Randomisations ##########################
                rand_num = np.random.rand()
                ###########################################################

                probability_level = 0.5 * (z_prob_adj + x_prob_adj)
                if rand_num < probability_level:
                    x_depth_mean = z_depth_mean - 0.04 * abs(x / math.floor(z_res / 2))

                    ################# Randomisations ##########################
                    random_height = np.random.normal(loc=x_depth_mean, scale=LEAF_DEPTH_VAR)
                    ###########################################################

                    clipped_height = min(max(random_height, -leaf_depth), 0)
                    leafTriPoints.append([x * point_spacing, clipped_height, z * point_spacing]*leaf_up_vec*2)

        vtk_points = vtk.vtkPoints()
        for point in leafTriPoints:
            vtk_points.InsertNextPoint(point)
        polyData = vtk.vtkPolyData()
        polyData.SetPoints(vtk_points)
        mesh_generator = vtk.vtkDelaunay2D()
        mesh_generator.SetInputData(polyData)
        mesh_generator.Update()
        mesh_polydata = mesh_generator.GetOutput()

        for point_idx in range(mesh_polydata.GetNumberOfPoints()):
            ################# Randomisations ##########################
            leafTriCols.append((min(max(np.random.normal(loc=LEAF_R_MEAN, scale=LEAF_R_VAR), 0), 1),
                            min(max(np.random.normal(loc=LEAF_G_MEAN, scale=LEAF_G_VAR), 0), 1),
                            min(max(np.random.normal(loc=LEAF_B_MEAN, scale=LEAF_B_VAR), 0), 1)))
            ###########################################################

        triangle_array = mesh_polydata.GetPolys()
        triangle_array.InitTraversal()
        idList = vtk.vtkIdList()
        while (triangle_array.GetNextCell(idList)):
            leafTriIndices.append((idList.GetId(0),
                             idList.GetId(1),
                             idList.GetId(2)))

        geo_plant_rep.triMeshPntSets.append(leafTriPoints)
        geo_plant_rep.triMeshPntIndxSets.append(leafTriIndices)
        geo_plant_rep.triMeshColSets.append(leafTriCols)
        geo_plant_rep.triMeshConnIdxs.append(stem_seg_idx)
        geo_plant_rep.triMeshSetLabels.append("Leaf")
    geo_plant_rep.numMeshSets = len(geo_plant_rep.triMeshPntSets)


def GenRandMeshLeaves_2(geo_plant_rep, end_seg_indxs ):
    """Generates randomised leaves given randomisation parameters"""

    res = 2
    size_multiplier = 1

    for stem_seg_idx in end_seg_indxs:
        leafTriPoints = []
        leafTriIndices = []
        leafTriCols = []
        stemVec = geo_plant_rep.tubePntSets[stem_seg_idx][-1]

        z_axis = [np.random.normal(loc=stemVec[0], scale=LEAF_Z_VAR),
                  np.random.normal(loc=stemVec[1], scale=LEAF_Z_VAR),
                  np.random.normal(loc=stemVec[2], scale=LEAF_Z_VAR)]
        z_axis /= np.linalg.norm(z_axis, 2)
        arb_axis = [np.random.normal(), np.random.normal(), np.random.normal()]
        x_axis = np.cross(z_axis, arb_axis)
        x_axis /= np.linalg.norm(x_axis, 2)

        reduced_res = max(round(size_multiplier * res), 9)
        leaf_width = max(size_multiplier * np.random.normal(loc=LEAF_WIDTH_MEAN, scale=LEAF_WIDTH_VAR), 0.015)
        leaf_length = max(size_multiplier * np.random.normal(loc=LEAF_LENGTH_MEAN, scale=LEAF_LENGTH_VAR), 0.02)
        point_spacing = leaf_width / reduced_res
        z_res = max(math.trunc(leaf_length * reduced_res / leaf_width), 9)
        mesh_points = vtk.vtkPoints()
        mesh_polyData = vtk.vtkPolyData()
        mesh_grid = vtk.vtkDelaunay2D()
        mesh_points.InsertNextPoint(0, 0, 0)
        mesh_points.InsertNextPoint(-point_spacing, 0, point_spacing)
        mesh_points.InsertNextPoint(point_spacing, 0, point_spacing)
        clipped_height = 0

        for z in range(1, z_res):
            z_prob_adj = 0.5 - abs(z_res / 2 - z) / z_res
            z_depth_mean = - LEAF_Z_CURl_MUL * (z / z_res) ** 2
            for x in range(-math.floor(reduced_res / 2), math.floor(reduced_res / 2) + 1):
                x_prob_adj = 0.5 - abs(x) / reduced_res

                rand_num = np.random.rand()
                probability_level = LEAF_CURVE_MEAN + LEAF_CURVE_RELMAX * (z_prob_adj + x_prob_adj)
                if rand_num < probability_level:
                    x_depth_mean = z_depth_mean - LEAF_X_CURL_MUL * abs(x / math.floor(res / 2))

                    # place first point with normal distribution
                    random_height = np.random.normal(loc=x_depth_mean, scale=LEAF_DEPTH_VAR)
                    clipped_height = min(max(random_height, -MAX_LEAF_DEPTH), 0)
                    point_pos = [min(max(x * point_spacing, -2), 2), min(max(clipped_height, -1), 1), min(max(z * point_spacing, -2), 2)]
                    leafTriPoints.append(point_pos)
                    mesh_points.InsertNextPoint(point_pos[0], point_pos[1], point_pos[2])

        mesh_points.InsertNextPoint(0, clipped_height, (z_res + 1) * point_spacing)
        mesh_polyData.SetPoints(mesh_points)
        mesh_grid.SetInputData(mesh_polyData)
        mesh_grid.Update()
        leaf_transform = vtk.vtkTransform()
        leaf_transform.PostMultiply()
        leaf_transform.RotateWXYZ(np.random.rand()*360, 0.5-np.random.rand(), 0.5-np.random.rand(), 0.5-np.random.rand())
        leaf_transform.Update()
        transform_filter = vtk.vtkTransformPolyDataFilter()
        transform_filter.SetTransform(leaf_transform)
        transform_filter.SetInputConnection(mesh_grid.GetOutputPort())
        transform_filter.Update()
        mesh_polydata = transform_filter.GetOutput()

        for point_idx in range(mesh_polydata.GetNumberOfPoints()):
            ################# Randomisations ##########################
            leafTriCols.append((min(max(np.random.normal(loc=LEAF_R_MEAN, scale=LEAF_R_VAR), 0), 1),
                                min(max(np.random.normal(loc=LEAF_G_MEAN, scale=LEAF_G_VAR), 0), 1),
                                min(max(np.random.normal(loc=LEAF_B_MEAN, scale=LEAF_B_VAR), 0), 1)))
            ###########################################################

        triangle_array = mesh_polydata.GetPolys()
        triangle_array.InitTraversal()
        idList = vtk.vtkIdList()
        while (triangle_array.GetNextCell(idList)):
            leafTriIndices.append((idList.GetId(0),
                                   idList.GetId(1),
                                   idList.GetId(2)))

        geo_plant_rep.triMeshPntSets.append(leafTriPoints)
        geo_plant_rep.triMeshPntIndxSets.append(leafTriIndices)
        geo_plant_rep.triMeshColSets.append(leafTriCols)
        geo_plant_rep.triMeshConnIdxs.append(stem_seg_idx)
        geo_plant_rep.triMeshSetLabels.append("Leaf")
    geo_plant_rep.numMeshSets = len(geo_plant_rep.triMeshPntSets)