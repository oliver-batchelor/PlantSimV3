import geoplantrep as PG
import DisplayVTKPlant as DP
import vtk
import numpy as np
from tqdm import tqdm
import cv2


INTERACT_WITH_PLANT = True


DIR_PATH = 'plant_samples_08_07_19/'
OUTPUT_DIR = '/local/plantData/plant_dataset_sample_3'

rep_plant = PG.PlantData()
rep_plant.LoadPlantFile(DIR_PATH + 'random_plant_V2_5.plant')

vtk_plant = DP.plantVTKDataDisplay(rep_plant)
vtk_plant.ConstructVTKStem()
vtk_plant.ConstructVTKMesh()

vtk_plant.InitRenderWindow( stereo_on=True, axes_on=False, bkgnd=[0, 0, 0], res_x=1200, res_y=600 )
vtk_plant.InitInteractor()
vtk_plant.InitLighting()
vtk_plant.CreatePolyDataMappers()
vtk_plant.CreatePolyDataActors()
vtk_plant.AddActors()

if INTERACT_WITH_PLANT:
    vtk_plant.RenderPlant()
else:
    vtk_plant.renderWindow.Render()
    vtk_plant.renderWindow.Modified()
    window_to_image.SetInput(vtk_plant.renderWindow)
    window_to_image.Modified()
    image_to_jpeg.SetFileName(OUTPUT_DIR + '/train/' + str(image_num) + '.jpg')
    image_to_jpeg.Modified()
    image_to_jpeg.Write()