import geoplantrep as PG
import VTK_PlantRendering.ConvertToVTKPlant as CVVTK
import VTK_PlantRendering.DisplayVTKPlant as DP
import PlantGeneration.GenerateStemStructure as GStem
import PlantGeneration.Generate2DTriStructures as G2D
import PlantGeneration.Generate3DTriStructures as G3D


# TODO: Leaf-stem merge for plant generation
# TODO: Camera pose and movement, saving stereo images with pose and doing a sweep of a plant/s
# TODO: Surface texturing
# TODO: Surface reflectivity
# TODO: Shadows
# TODO: Plant movement, growth, and node relative positioning
# TODO: leaf occlusion test
#TODO: make fully vtk free leaf generator


LOCAL_HOME_PATH = '/local/'
OUTDIR_PATH = LOCAL_HOME_PATH + 'Dropbox/PlantSimData/'


DISPLAY_PLANT_SAMPLES = True
LOAD_PLANT = True
SAVE_PLANT = False

rep_plant = PG.PlantData()

if LOAD_PLANT:
    rep_plant.LoadPlantFile('test_save.plant')
else:
    # Generate Randomised plant using a combination of algorithms
    GStem.GenRandSplineStem(rep_plant, 40)
    end_stem_indxs = GStem.FindEndSegIndxs(rep_plant)
    G2D.GenRandMeshLeaves_1(rep_plant, end_stem_indxs[::2])
    G2D.GenRandMeshLeaves_2(rep_plant, end_stem_indxs[1::2])
    #G3D.GenRandFruit(rep_plant)


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

if SAVE_PLANT:
    rep_plant.SavePlantFile('test_save.plant')