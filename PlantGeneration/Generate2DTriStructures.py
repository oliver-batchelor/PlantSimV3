import geoplantrep as PG
import numpy as np
import math
import vtk

LEAF_G_MEAN = 0.9
LEAF_G_VAR = 0.1
LEAF_R_MEAN = 0.5
LEAF_R_VAR = 0.05
LEAF_B_MEAN = 0.3
LEAF_B_VAR = 0.02

LEAF_LENGTH_MEAN = 0.08
LEAF_LENGTH_VAR = 0.02
LEAF_WIDTH_MEAN = 0.06
LEAF_WIDTH_VAR = 0.01
LEAFGRID_RES = 100e-6    #m

LEAF_LVEC_VAR = 0.2
LEAF_LVEC_VERT_MUL = 0.3
LEAF_UP_BIAS = 4
LEAFSPLNE_PNTCNT_MEAN = 12
LEAFSPLNE_PNTCNT_VAR = 4
LEAFSPLNE_PNT_SEP_VAR = 0.003
LEAFSPLNE_VEC_VAR = 0.003

LEAFCOL_MAX_ITTS = 80

#TODO: make fully vtk free leaf generator
#TODO: make fully vtk free triangle intersection detector


def PartCollisionCheck(leaf_tri_pnts, leaf_tri_indcs, part_tri_pnts, part_tri_indcs):
    """Scan through part given and find intersections with leaf triangles"""
    collision = False
    for leaf_tri in leaf_tri_indcs:
        for part_tri in part_tri_indcs:
            if vtk.vtkTriangle.TrianglesIntersect(leaf_tri_pnts[leaf_tri[0]], leaf_tri_pnts[leaf_tri[1]],
                                                  leaf_tri_pnts[leaf_tri[2]], part_tri_pnts[part_tri[0]],
                                                  part_tri_pnts[part_tri[1]], part_tri_pnts[part_tri[2]]):
                collision = True
                break
        if collision:
            break
    return collision


def ScanForIntersection(geo_plant_rep, leaf_tri_pnts, leaf_tri_indcs):
    """Scan through plant part list and Check for intersections"""

    collision = False
    for idx in range(geo_plant_rep.numMeshSets):
        if PartCollisionCheck(leaf_tri_pnts, leaf_tri_indcs,
                           geo_plant_rep.meshPntSets[idx], geo_plant_rep.triMeshPntIndxSets[idx]):
            collision = True
            break
    return collision


def GenSplinePoints(vec, offset, length, num_points):
    point_sep_mean = length / num_points
    spline_pnt_set = []
    for point_n in range(num_points):
        ################# Randomisations ##########################
        distance_to_prev = np.random.normal(loc=point_sep_mean, scale=LEAFSPLNE_PNT_SEP_VAR)
        vecpnt_offset = np.random.normal(loc=0, scale=LEAFSPLNE_VEC_VAR, size=3) + offset
        offset_mul = 1 - max(abs(0.5 - point_n/num_points), 0.01)
        spline_pnt_set.append(np.array(vec)*distance_to_prev*(point_n+1) + np.array(vecpnt_offset)*offset_mul)
        ###########################################################
    return spline_pnt_set


def GenRandLeaves(geo_plant_rep, end_seg_idxs):
    """Create fine grain leaf structures"""
    for stem_seg_idx in end_seg_idxs:
        leafTriPoints = [[0, 0, 0]]
        leafTriCols = []
        leafTriIndices = []

        stemend_pnts = geo_plant_rep.tubePntSets[stem_seg_idx]
        stemend_vec = stemend_pnts[-1] - stemend_pnts[-2]

        ################# Randomisations ##########################
        leaf_length = np.random.normal(loc=LEAF_LENGTH_MEAN, scale=LEAF_LENGTH_VAR)
        leaf_width = np.random.normal(loc=LEAF_WIDTH_MEAN, scale=LEAF_WIDTH_VAR)
        leaf_lvec = stemend_vec + np.random.normal(scale=LEAF_LVEC_VERT_MUL, size=3)
        leaf_hvec = np.cross(np.random.normal(size=3) + [0, LEAF_UP_BIAS, 0], leaf_lvec)
        n_splne_pts_l = max(round(np.random.normal(loc=LEAFSPLNE_PNTCNT_MEAN, scale=LEAFSPLNE_PNTCNT_VAR)), 3)
        n_splne_pts_c = max(round(np.random.normal(loc=LEAFSPLNE_PNTCNT_MEAN, scale=LEAFSPLNE_PNTCNT_VAR)), 3)
        n_splne_pts_r = max(round(np.random.normal(loc=LEAFSPLNE_PNTCNT_MEAN, scale=LEAFSPLNE_PNTCNT_VAR)), 3)
        ###########################################################
        leaf_lvec /= np.linalg.norm(leaf_lvec, 2)
        leaf_hvec /= np.linalg.norm(leaf_hvec, 2)

        leafTriPoints.extend(GenSplinePoints(leaf_lvec, [0, 0, 0], leaf_length, n_splne_pts_c))
        leafTriPoints.extend(GenSplinePoints(leaf_lvec, -0.5*leaf_width*leaf_hvec, 0.7*leaf_length, n_splne_pts_l))
        leafTriPoints.extend(GenSplinePoints(leaf_lvec, 0.5*leaf_width*leaf_hvec, 0.7*leaf_length, n_splne_pts_r))

        vtk_points = vtk.vtkPoints()
        for point in leafTriPoints:
            vtk_points.InsertNextPoint(point)
        polyData = vtk.vtkPolyData()
        polyData.SetPoints(vtk_points)
        mesh_generator = vtk.vtkDelaunay2D()
        mesh_generator.SetInputData(polyData)
        mesh_generator.Update()
        mesh_polydata = mesh_generator.GetOutput()

        for n in range(LEAFCOL_MAX_ITTS):
            triangle_array = mesh_polydata.GetPolys()
            triangle_array.InitTraversal()
            idList = vtk.vtkIdList()
            while triangle_array.GetNextCell(idList):
                leafTriIndices.append((idList.GetId(0),
                                       idList.GetId(1),
                                       idList.GetId(2)))

            if not ScanForIntersection(geo_plant_rep, leafTriPoints, leafTriIndices):
                break

        for point_idx in range(mesh_polydata.GetNumberOfPoints()):
            ################# Randomisations ##########################
            leafTriCols.append((min(max(np.random.normal(loc=LEAF_R_MEAN, scale=LEAF_R_VAR), 0), 1),
                                min(max(np.random.normal(loc=LEAF_G_MEAN, scale=LEAF_G_VAR), 0), 1),
                                min(max(np.random.normal(loc=LEAF_B_MEAN, scale=LEAF_B_VAR), 0), 1)))
            ###########################################################

        geo_plant_rep.triMeshPntSets.append(leafTriPoints)
        geo_plant_rep.triMeshPntIndxSets.append(leafTriIndices)
        geo_plant_rep.triMeshColSets.append(leafTriCols)
        geo_plant_rep.triMeshConnIdxs.append(stem_seg_idx)
        geo_plant_rep.triMeshSetLabels.append("Leaf")
    geo_plant_rep.numMeshSets = len(geo_plant_rep.triMeshPntSets)
