import geoplantrep as PG
import VTK_PlantRendering.ConvertToVTKPlant as CVVTK
import VTK_PlantRendering.DisplayVTKPlant as DP
import DataSetCreation.imageProccessing as IP
import vtk
import numpy as np
from tqdm import tqdm


INTERACT_WITH_PLANT = True

LOCAL_HOME_PATH = '/local/'
DIR_PATH = LOCAL_HOME_PATH + ''
PLANT_FILENAME = ''
OUTPUT_DIR = LOCAL_HOME_PATH + ''


rep_plant = PG.PlantData()
rep_plant.LoadPlantFile(DIR_PATH + PLANT_FILENAME)

# Display Plants
vtk_plant = CVVTK.vtkPlantData(rep_plant)
vtk_plant.ConstructVTKStem()
vtk_plant.ConstructVTKMesh()
vtk_plant.CreatePolyDataMappers(disp_pot=True)
vtk_plant.CreatePolyDataActors()
vtk_plant.SetActorPostions()

plant_display = DP.plantVTKDataDisplay(vtk_plant)
plant_display.InitRenderWindow( axes_on=False, bkgnd=[0, 0, 0], res_x=600, res_y=600 )
plant_display.AddActors()
plant_display.InitInteractor()
plant_display.InitLighting()
plant_display.RenderPlant()


# Add red spheres at node points
annotation_actors = []
annotation_mappers = []
annotation_objects = []
for stem_seg_points in rep_plant.tubePntSets:
    stem_base_node = vtk.vtkSphereSource()
    stem_base_node.SetCenter(stem_seg_points[0])
    stem_base_node.SetRadius(0.001)
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
    vtk_plant.renderer.AddActor(stem_base_actor)

if INTERACT_WITH_PLANT:
    vtk_plant.RenderPlant()
else:
    vtk_plant.renderWindow.OffScreenRenderingOn()
    for image_num in tqdm(range(0, 5)):
        # Show plant and save image
        vtk_plant.SetPlantVisible(True)
        for ann_actor in annotation_actors:
            ann_actor.VisibilityOff()
        vtk_plant.SaveCameraImage(OUTPUT_DIR, 'test/', str(image_num))

        # Hide plant, display annotations, save image
        vtk_plant.SetPlantVisible( False )
        for ann_actor in annotation_actors:
            ann_actor.VisibilityOn()
        vtk_plant.SaveCameraImage(OUTPUT_DIR, 'test/', 'a' + str(image_num))

        IP.pyrDownSmple(OUTPUT_DIR + 'test/a' + str(image_num) + '.jpg')

        # Move camera
        vtk_plant.MoveCamera(1 - 2*np.random.rand(3), 0.1 - 0.2*np.random.rand(3), 0.1 - 0.2*np.random.rand(3))