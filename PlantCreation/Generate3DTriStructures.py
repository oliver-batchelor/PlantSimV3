import geoplantrep as PG
import numpy as np


MEAN_FRUIT_LENGTH = 0.07
FRUIT_LENGTH_VAR = 0.02
MEAN_FRUIT_WIDTH = 0.015
FRUIT_WIDTH_VAR = 0.003
FRUIT_MAX_POINTS = 100
FRUIT_ROT_VAR = 20


def GenRandFruit(geo_plant_rep):
    """Generates randomised fruit given randomisation parameters"""

    # num_points = max(round(np.random.rand() * FRUIT_MAX_POINTS), 5)
    # fruit_length = np.random.normal(loc=MEAN_FRUIT_LENGTH, scale=FRUIT_LENGTH_VAR)
    # fruit_width = np.random.normal(loc=MEAN_FRUIT_WIDTH, scale=FRUIT_WIDTH_VAR)
    # self.fruitPoints = np.zeros((num_points, 3))
    # self.stemIDX = stem_idx
    # self.fruitPoints[0] = [0, 0, 0]
    #
    # for point in range(num_points - 1):
    #     p_x = fruit_width * (np.random.rand() - 0.5)
    #     p_z = fruit_width * (np.random.rand() - 0.5)
    #     p_y = -(fruit_length * point / num_points + 0.001*(np.random.rand() - 0.5))
    #     self.fruitPoints[point + 1] = [p_x, p_y, p_z]
    #
    # self.fruitRelRot = np.random.normal(scale=FRUIT_ROT_VAR, size=3)