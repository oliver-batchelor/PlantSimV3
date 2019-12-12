import geoplantrep as PG
import VTK_PlantRendering.ConvertToVTKPlant as CVVTK
import VTK_PlantRendering.DisplayVTKPlant as DP
import DataSetCreation.imageProccessing as IP
import VTK_PlantRendering.GenerateBackground as BG
import PlantGeneration.GenerateStemStructure as GStem
import PlantGeneration.Generate2DTriStructures as G2D
import PlantGeneration.Generate3DTriStructures as G3D
import vtk
import numpy as np
from tqdm import tqdm


LOCAL_HOME_PATH = '/local/'
OUTDIR_PATH = LOCAL_HOME_PATH + 'Dropbox/PlantSimData/plant_HM_dataset_16_10_19/'
BCKGND_IMAGE_DIR = LOCAL_HOME_PATH + 'Dropbox/PlantPhotos/'
LABEL = 'train/'

START_IND = 5000
IM_PER_PLANT = 50
INTERACT_WITH_PLANT = False

for plant_num in tqdm(range(100)):
    rep_plant = PG.PlantData()
    vtk_plant_list = []

    rep_plant = PG.PlantData()
    GStem.GenRandSplineStem(rep_plant, 30)
    #end_stem_indxs = GStem.FindEndSegIndxs(rep_plant)
    #G2D.GenRandLeaves(rep_plant, end_seg_idxs=end_stem_indxs)

    vtk_plant = CVVTK.vtkPlantData(rep_plant)
    vtk_plant.BuildComponents()
    vtk_plant.SetActorPoses()
    vtk_plant_list.append(vtk_plant)

    plant_display = DP.plantVTKDataDisplay(vtk_plant_list)
    plant_display.InitRenderWindow( axes_on=False, bkgnd=[0.0, 0.0, 0.0], res_x=600, res_y=600 )
    bkgnd_scene = BG.BackgroundScene()
    #bkgnd_scene.GeneratePlantPots(vtk_plant_list)
    plant_display.AddActors(bckgnd_actors=bkgnd_scene.BackgroundActorList)
    plant_display.InitInteractor()
    plant_display.InitLighting(mode=1, intensity_mul=2)

    # Add red spheres at node points
    annotation_actors = []
    annotation_mappers = []
    annotation_objects = []
    for plant in vtk_plant_list:
        for stem_seg_actor in plant.StemActorList:
            stem_base_node = vtk.vtkSphereSource()
            stem_base_node.SetCenter(stem_seg_actor.GetPosition())
            stem_base_node.SetRadius(0.002)
            stem_base_node.Update()
            annotation_objects.append(stem_base_node)
            stem_base_mapper = vtk.vtkPolyDataMapper()
            stem_base_mapper.SetInputData(stem_base_node.GetOutput())
            annotation_mappers.append(stem_base_mapper)
            stem_base_actor = vtk.vtkActor()
            stem_base_actor.SetMapper(stem_base_mapper)
            stem_base_actor.GetProperty().SetColor((1, 0, 0))
            stem_base_actor.VisibilityOff()
            annotation_actors.append(stem_base_actor)
            plant_display.renderer.AddActor(stem_base_actor)


    if INTERACT_WITH_PLANT:
        plant_display.RenderPlant()
        break
    else:
        plant_display.renderWindow.OffScreenRenderingOn()
        start_ind = IM_PER_PLANT*plant_num + START_IND
        index_a = np.arange(start_ind, start_ind + IM_PER_PLANT)
        for image_num in index_a:
            # Show plant and save image
            plant_display.SetPlantVisible(True)
            plant_display.SetActorsVisible(bkgnd_scene.BackgroundActorList, True)
            plant_display.SetActorsVisible(annotation_actors, False)

            plant_display.SaveCameraImage(OUTDIR_PATH, LABEL, str(image_num))

            # Hide plant, display annotations, save image
            plant_display.SetPlantVisible( False )
            plant_display.SetActorsVisible(bkgnd_scene.BackgroundActorList, False)
            plant_display.SetActorsVisible(annotation_actors, True)

            plant_display.SaveCameraImage(OUTDIR_PATH, LABEL, 'a' + str(image_num))

            IP.pyrDownSmple(OUTDIR_PATH + LABEL + 'a' + str(image_num) + '.jpg')
            IP.AddRandomBackground(OUTDIR_PATH + LABEL + str(image_num) + '.jpg', BCKGND_IMAGE_DIR)

            # Move camera
            plant_display.MoveCamera(1 - 2*np.random.rand(3), 0.1 - 0.2*np.random.rand(3), 0.1 - 0.2*np.random.rand(3))