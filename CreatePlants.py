import geoplantrep as PG
import VTK_PlantRendering.ConvertToVTKPlant as CVVTK
import VTK_PlantRendering.DisplayVTKPlant as DP
import VTK_PlantRendering.GenerateBackground as BG
import PlantGeneration.GenerateStemStructure as GStem
import PlantGeneration.Generate2DTriStructures as G2D
import PlantGeneration.Generate3DTriStructures as G3D


# TODO: Stem radius based on height/stem length
# TODO: make fully vtk free leaf generator
# TODO: Polydata noise addition and smoothing
# TODO: Leaf-stem merge for plant generation
# TODO: Surface texturing
# TODO: Surface reflectivity
# TODO: Shadows
# TODO: Plant movement, growth


LOCAL_HOME_PATH = '/local/'
OUTDIR_PATH = LOCAL_HOME_PATH + 'Dropbox/PlantSimData/V3Plants/'


DISPLAY_PLANT_SAMPLES = True
LOAD_PLANT = False
SAVE_PLANT = False

rep_plant = PG.PlantData()
vtk_plant_list = []

for plant_n in range(1):
    if LOAD_PLANT:
        rep_plant = PG.PlantData()
        rep_plant.LoadPlantFile(OUTDIR_PATH + '10_09_19_' + str(plant_n) + '.plant')
        vtk_plant = CVVTK.vtkPlantData(rep_plant)
        vtk_plant.BuildComponents()
        vtk_plant.SetActorPoses()
        vtk_plant_list.append(vtk_plant)
    else:
        for plant_c in range(1):
            for plant_r in range(1):
                # Generate Randomised plant using a combination of algorithms
                rep_plant = PG.PlantData()
                GStem.GenRandSplineStem(rep_plant, 50)
                end_stem_indxs = GStem.FindEndSegIndxs(rep_plant)
                G2D.GenRandLeaves(rep_plant, end_stem_indxs)
                #G3D.GenRandFruit(rep_plant)

                vtk_plant = CVVTK.vtkPlantData(rep_plant)
                vtk_plant.BuildComponents()
                vtk_plant.SetActorPoses(offset=[0.6 * plant_r, 0, 0.2 * plant_c])
                vtk_plant_list.append(vtk_plant)

    if DISPLAY_PLANT_SAMPLES:
        plant_display = DP.plantVTKDataDisplay(vtk_plant_list)
        plant_display.InitRenderWindow( stereo_on=False, axes_on=True, bkgnd=[0.0, 0.0, 0.0], res_x=1920, res_y=1080 )
        bkgnd_scene = BG.BackgroundScene()
        #bkgnd_scene.GeneratePlantPots(vtk_plant_list)
        plant_display.AddActors(bckgnd_actors=bkgnd_scene.BackgroundActorList)
        plant_display.InitInteractor()
        plant_display.InitLighting(mode=1, intensity_mul=0.8)
        plant_display.RenderPlant()
        plant_display.ClearRenderer()
        vtk_plant_list = []

    if SAVE_PLANT:
        rep_plant.SavePlantFile(OUTDIR_PATH + '10_09_19_' + str(plant_n) + ".plant")