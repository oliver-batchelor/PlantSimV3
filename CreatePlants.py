import geoplantrep as PG
import VTK_PlantRendering.ConvertToVTKPlant as CVVTK
import VTK_PlantRendering.DisplayVTKPlant as DP
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
OUTDIR_PATH = LOCAL_HOME_PATH + 'Dropbox/PlantSimData/'


DISPLAY_PLANT_SAMPLES = True
LOAD_PLANT = False
SAVE_PLANT = False

rep_plant = PG.PlantData()

if LOAD_PLANT:
    rep_plant.LoadPlantFile('test_save.plant')
else:
    # Generate Randomised plant using a combination of algorithms
    GStem.GenRandSplineStem(rep_plant, 20)
    end_stem_indxs = GStem.FindEndSegIndxs(rep_plant)
    G2D.GenRandLeaves(rep_plant, end_stem_indxs)
    #G3D.GenRandFruit(rep_plant)


# Display Plants
vtk_plant = CVVTK.vtkPlantData(rep_plant)
vtk_plant.ConstructVTKStem()
vtk_plant.ConstructVTKMesh()
vtk_plant.CreatePolyDataMappers(disp_pot=True)
vtk_plant.CreatePolyDataActors()
vtk_plant.SetActorPostions()

plant_display = DP.plantVTKDataDisplay(vtk_plant)
plant_display.InitRenderWindow( axes_on=False, bkgnd=[0, 0, 0], res_x=1920, res_y=1080 )
plant_display.AddActors()
plant_display.InitInteractor()
plant_display.InitLighting()
plant_display.RenderPlant()

if SAVE_PLANT:
    rep_plant.SavePlantFile('test_save.plant')