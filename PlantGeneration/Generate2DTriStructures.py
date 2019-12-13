import geoplantrep as PG
import numpy as np
import math
import vtk
import cv2

LEAF_G_MEAN = 0.7
LEAF_G_VAR = 0.02
LEAF_R_MEAN = 0.37
LEAF_R_VAR = 0.01
LEAF_B_MEAN = 0.28
LEAF_B_VAR = 0.01

LEAF_LENGTH_MEAN = 0.12
LEAF_LENGTH_VAR = 0.03
LEAF_WIDTH_MEAN = 0.07
LEAF_WIDTH_VAR = 0.02
LEAF_DEPTH_VAR = 0.002
LEAFGRID_RES = 2e-3    #m

LEAF_LVEC_VAR = 0.2
LEAF_LVEC_VERT_MUL = 0.3
LEAF_UP_BIAS = 100
LEAFSPLNE_PNTCNT_MEAN = 12
LEAFSPLNE_PNTCNT_VAR = 4
LEAFSPLNE_PNT_SEP_VAR = 0.003
LEAFSPLNE_VEC_VAR = 0.003

LEAFCOL_MAX_ITTS = 80

#TODO: make fully vtk free leaf generator
#TODO: make fully vtk free triangle intersection detector


def TriCollisionCheck(leaf_tri_pnts, leaf_tri_indcs, part_tri_pnts, part_tri_indcs):
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


def SplineTriCollisionCheck():
    """Scan tube point splines for intersection with triangles"""



def ScanForIntersection(geo_plant_rep, leaf_tri_pnts, leaf_tri_indcs):
    """Scan through plant part list and Check for intersections"""

    collision = False
    for idx in range(geo_plant_rep.numMeshSets):
        if TriCollisionCheck(leaf_tri_pnts, leaf_tri_indcs,
                           geo_plant_rep.meshPntSets[idx], geo_plant_rep.triMeshPntIndxSets[idx]):
            collision = True
            break
    if not collision:
        for idx in range(geo_plant_rep.numTubeSets):
            if SplineTriCollisionCheck():
                collision = True
                break

    return collision


def GenNoiseImg(shape, kernel_s, d_scale):
    """Creates noise image of specified size"""
    # Generate random smooth surface
    lowres_noise = np.random.normal(scale=d_scale*LEAF_DEPTH_VAR, size=kernel_s)
    lowres_noise -= lowres_noise[0, math.floor(kernel_s[0]/2)]

    return cv2.resize(lowres_noise, shape, interpolation=cv2.INTER_CUBIC)


def GenRandomLeaf(shape):
    """Cuts a leaf shape out of randomised image with gaussian curves"""
    leaf_img_w = round(shape[1] / LEAFGRID_RES)
    if leaf_img_w % 2 == 0:
        leaf_img_w += 1
    img_s = (round(shape[0] / LEAFGRID_RES), leaf_img_w)

    noise_img = np.zeros(img_s)
    for k in range(3, 13, 2):
        noise_img += GenNoiseImg(np.shape(noise_img)[::-1], (k, k), 3/k)

    top_curve = np.zeros((1, img_s[0]))
    bot_curve = np.zeros((1, img_s[0]))
    for s in range(3):
        T = (img_s[0] - 1) / (round(1 + s*np.random.rand())*math.pi)*np.random.normal(loc=1, scale=0.2)
        top_curve += [math.sin(i/T)/(s+1) for i in range(img_s[0])]
        T = (img_s[0] - 1) / (round(1 + s * np.random.rand()) * math.pi) * np.random.normal(loc=1, scale=0.2)
        bot_curve += [math.sin(i/T)/(s+1) for i in range(img_s[0])]

    mid_leaf_idx = math.floor(leaf_img_w/2)

    top_curve = top_curve - top_curve[0][0]
    top_curve *= mid_leaf_idx / max(0.001, np.max(top_curve))
    top_curve = (top_curve[0] + mid_leaf_idx + 1).astype(int)
    bot_curve = bot_curve - bot_curve[0][0]
    bot_curve *= mid_leaf_idx / max(0.001, np.max(bot_curve))
    bot_curve = (mid_leaf_idx - bot_curve[0]).astype(int)

    point_set = []
    for row_i, row in enumerate(noise_img):
        l_lim = bot_curve[row_i]
        u_lim = top_curve[row_i]
        for col_i, z in enumerate(row[l_lim:u_lim]):
            point_set.append([row_i*LEAFGRID_RES,
                              (l_lim + col_i - mid_leaf_idx)*LEAFGRID_RES,
                              noise_img[row_i, col_i + l_lim]])

    # noise_img += abs(np.min(noise_img))
    # noise_img /= np.max(noise_img)
    # cv2.imshow("test", noise_img)
    # cv2.waitKey(10000)

    return point_set


def GenRandLeaves(geo_plant_rep, end_seg_idxs):
    """Create fine grain leaf structures"""
    for stem_seg_idx in end_seg_idxs:
        stemend_pnts = geo_plant_rep.tubePntSets[stem_seg_idx]
        stemend_vec = stemend_pnts[-1] - stemend_pnts[-2]
        stemend_vec /= np.linalg.norm(stemend_vec, 2)

        ################# Randomisations ##########################
        leaf_length = np.random.normal(loc=LEAF_LENGTH_MEAN, scale=LEAF_LENGTH_VAR)
        leaf_width = np.random.normal(loc=LEAF_WIDTH_MEAN, scale=LEAF_WIDTH_VAR)
        leaf_xvec = stemend_vec + np.random.normal(scale=LEAF_LVEC_VERT_MUL, size=3)
        leaf_zvec = np.cross(np.random.normal(size=3) + [0, LEAF_UP_BIAS, 0], leaf_xvec)
        ###########################################################
        leaf_xvec /= np.linalg.norm(leaf_xvec, 2)
        leaf_zvec /= np.linalg.norm(leaf_zvec, 2)
        leaf_yvec = np.cross(leaf_zvec, leaf_xvec)
        leaf_yvec /= np.linalg.norm(leaf_yvec, 2)

        leaf_frame = np.transpose([leaf_xvec, leaf_yvec, leaf_zvec])

        leafTriPoints = GenRandomLeaf((leaf_length, leaf_width))
        leafTriCols = []
        leafTriIndices = []

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

        R = np.matmul(leaf_frame, [[1, 0, 0],
                                    [0, 0, 1],
                                    [0, -1, 0]])

        sy = math.sqrt(R[0, 0] * R[0, 0] + R[1, 0] * R[1, 0])
        singular = sy < 1e-6
        if not singular:
            x = math.atan2(R[2, 1], R[2, 2])
            y = math.atan2(-R[2, 0], sy)
            z = math.atan2(R[1, 0], R[0, 0])
        else:
            x = math.atan2(-R[1, 2], R[1, 1])
            y = math.atan2(-R[2, 0], sy)
            z = 0
        geo_plant_rep.triMeshOris.append([math.degrees(x), math.degrees(y), math.degrees(z)])

        geo_plant_rep.triMeshPntSets.append(leafTriPoints)
        geo_plant_rep.triMeshPntIndxSets.append(leafTriIndices)
        geo_plant_rep.triMeshColSets.append(leafTriCols)
        geo_plant_rep.triMeshConnIdxs.append(stem_seg_idx)
        geo_plant_rep.triMeshSetLabels.append("Leaf")
    geo_plant_rep.numMeshSets = len(geo_plant_rep.triMeshPntSets)
