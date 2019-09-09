import geoplantrep as PG
import numpy as np

STEMSEG_CNT_VAR = 3
STEMSEG_LENGTH_MEAN = 0.04
STEMSEG_LENGTH_VAR = 0.03
STEMSEG_UPMUL_MEAN = 2
STEMSEG_UPMUL_VAR = 0.6
STEMSEG_HORIZ_VAR = 3
STEMSEG_OUTRAD_MUL = 0.0005

STEMPNT_CNT_MEAN = 2
STEMPNT_CNT_VAR = 1
STEMPNT_SEP_VAR = 0.002
STEMPNT_VECOFFSET_VAR = 0.01

STEMTIPS_RAD = 0.001
STEMRAD_RAD_VAR = 0.0001
STEMRAD_MIN = 0.001
STEMRAD_MAX = 0.01

STEMCOL_R_MEAN = 100
STEMCOL_R_VAR = 12
STEMCOL_G_MEAN = 160
STEMCOL_G_VAR = 10
STEMCOL_B_MEAN = 50
STEMCOL_B_VAR = 5


def FindEndSegIndxs(geo_plant_rep):
    """Returns a list of indexes of stem segments which do not lead to another seg"""
    branch_end_mask = np.ones((geo_plant_rep.numTubeSets, 1))
    branch_end_idxs = []
    for seg_conn_idx in geo_plant_rep.tubeConnIdxs:
        branch_end_mask[seg_conn_idx] = 0
    for idx, end in enumerate(branch_end_mask):
        if end:
            branch_end_idxs.append(idx)
    return branch_end_idxs


def GenRandStemRadii(geo_plant_rep, mode=0):
    """Generates semi-randomised radii in a given mode"""
    for stemseg_idx, seg_points in enumerate(geo_plant_rep.tubePntSets):
        spline_rad_set = []

        for point_n in range(len(seg_points)):
            rad_mul = 1 + 0.25*np.sqrt(geo_plant_rep.numTubeSets - stemseg_idx) + 0.04*point_n
            ################# Randomisations ##########################
            pnt_rad = min(max(STEMTIPS_RAD*rad_mul, STEMRAD_MIN)  + rad_mul*np.random.normal(scale=STEMRAD_RAD_VAR), STEMRAD_MAX)
            ###########################################################

            spline_rad_set.append(pnt_rad)
        if stemseg_idx > 0:
            spline_rad_set[0] = geo_plant_rep.tubeRadSets[geo_plant_rep.tubeConnIdxs[stemseg_idx]][-1]*0.99
            spline_rad_set[1] = geo_plant_rep.tubeRadSets[geo_plant_rep.tubeConnIdxs[stemseg_idx]][-1]*1.1
            spline_rad_set[2] = geo_plant_rep.tubeRadSets[geo_plant_rep.tubeConnIdxs[stemseg_idx]][-1] * 1.2
        geo_plant_rep.tubeRadSets.append(spline_rad_set)


def GenRandStemCols(geo_plant_rep, mode=0):
    """Genrates semi-randomised colour profile in a given mode"""
    for stemseg_idx, seg_points in enumerate(geo_plant_rep.tubePntSets):
        spline_col_set = []

        for point_n in range(len(seg_points)):
            ################# Randomisations ##########################
            pnt_col = (np.random.normal(loc=STEMCOL_R_MEAN, scale=STEMCOL_R_VAR),
                       np.random.normal(loc=STEMCOL_G_MEAN, scale=STEMCOL_G_VAR),
                       np.random.normal(loc=STEMCOL_B_MEAN, scale=STEMCOL_B_VAR))
            ###########################################################

            spline_col_set.append(pnt_col)
        if stemseg_idx > 0:
            spline_col_set[0] = geo_plant_rep.tubeColSets[geo_plant_rep.tubeConnIdxs[stemseg_idx]][-2]
            spline_col_set[1] = geo_plant_rep.tubeColSets[geo_plant_rep.tubeConnIdxs[stemseg_idx]][-1]
        geo_plant_rep.tubeColSets.append(spline_col_set)


def GenRandSplineSeg(geo_plant_rep, seg_vec, mode=0):
    """Generates randomised spline segment for stem given randomised parameters"""

    ################# Randomisations ##########################
    up_bias = np.random.normal(loc=STEMSEG_UPMUL_MEAN, scale=STEMSEG_UPMUL_VAR)

    dx = np.random.normal(loc=STEMSEG_OUTRAD_MUL * np.sign(seg_vec[0]),
                          scale=STEMSEG_HORIZ_VAR)
    dz = np.random.normal(loc=STEMSEG_OUTRAD_MUL * np.sign(seg_vec[2]),
                          scale=STEMSEG_HORIZ_VAR)
    dy = np.sqrt(dx ** 2 + dz ** 2) * up_bias

    new_seg_vec = seg_vec + np.array([dx, dy, dz])*0.1/np.linalg.norm([dx, dy, dz], 2)
    new_seg_vec /= np.linalg.norm(new_seg_vec, 2)

    seg_length = np.random.normal(loc=STEMSEG_LENGTH_MEAN, scale=STEMSEG_LENGTH_VAR)
    num_points = max(round(np.random.normal(loc=STEMPNT_CNT_MEAN, scale=STEMPNT_CNT_VAR)), 1)
    ###########################################################

    point_sep_mean = seg_length / num_points
    spline_pnt_set = [-seg_vec*0.0005, [0, 0, 0], seg_vec*0.0005]

    for point_n in range(num_points):

        ################# Randomisations ##########################
        distance_to_prev = np.random.normal(loc=point_sep_mean, scale=STEMPNT_SEP_VAR)
        vecpnt_offset = [np.random.normal(loc=0, scale=STEMPNT_VECOFFSET_VAR),
                         np.random.normal(loc=0, scale=STEMPNT_VECOFFSET_VAR),
                         np.random.normal(loc=0, scale=STEMPNT_VECOFFSET_VAR)]
        ###########################################################

        spline_pnt_set.append(np.array(spline_pnt_set[-1]) +
                              np.array(new_seg_vec)*distance_to_prev +
                              np.array(vecpnt_offset))

    geo_plant_rep.tubePntSets.append(spline_pnt_set)
    geo_plant_rep.tubeSetLabels.append("Stem")


def GenRandSplineStem(geo_plant_rep, num_segs):
    """Generates randomised splines given randomisation parameters"""

    # generate random stem vectors

    ################# Randomisations ##########################
    num_segs_rand = max(round(np.random.normal(loc=num_segs, scale=STEMSEG_CNT_VAR)), 1)
    ###########################################################

    #print("NumberOfSegs: {}".format(num_segs_rand))
    prev_seg_idx = -1
    prev_end_pnt = [0, 0, 0]
    for seg_n in range(num_segs_rand):
        if prev_seg_idx == -1:
            prev_end_vec = [0, 1, 0]
        else:
            prev_end_vec = geo_plant_rep.tubePntSets[prev_seg_idx][-1] - geo_plant_rep.tubePntSets[prev_seg_idx][-2]
        prev_end_vec /= np.linalg.norm(prev_end_vec, 2)
        GenRandSplineSeg(geo_plant_rep, prev_end_vec)
        geo_plant_rep.tubeConnIdxs.append(prev_seg_idx)

        ################# Randomisations ##########################
        biasrand_seg_idx = min(max(np.random.normal(loc=2/3, scale=1/3), 0), 1)
        prev_seg_idx = round(biasrand_seg_idx * (len(geo_plant_rep.tubePntSets) - 1))
        ###########################################################

    geo_plant_rep.numTubeSets = num_segs_rand
    # generate stem radii
    GenRandStemRadii(geo_plant_rep)

    # generate stem colours
    GenRandStemCols(geo_plant_rep)