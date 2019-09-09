import vtk
import VTK_PlantRendering.ConvertToVTKPlant as CVVTK
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np


POT_SIZE = 0.08

TOP_LIGHT_INTENS_MEAN = 0.3
BOT_LIGHT_INTENS_MEAN = 0.08
LIGHT_INTENS_VAR = 0.03

SHADOWS_RESOLUTION = 1000


class plantVTKDataDisplay():
    """Class controlling render environment and interactors"""
    def __init__(self, plant_list):
        self.vtk_plant_list = plant_list

        self.BackgroundOriPoses = []
        self.BackgroundPolydataList = []
        self.BackgroundMapperList = []
        self.BackgroundActorList = []

        self.renderWindow = None
        self.renderer = None
        self.bakerPass = None
        self.windowInteractor = None
        self.windowFilter = None
        self.camera = None
        self.writer = None


    def InitRenderPasses(self):
        """Initiate and add realism passes to renderer"""
        cameraP = vtk.vtkCameraPass()
        opaque = vtk.vtkOpaquePass()
        peeling = vtk.vtkDepthPeelingPass()
        peeling.SetMaximumNumberOfPeels(100)
        peeling.SetOcclusionRatio(0.1)

        translucent = vtk.vtkTranslucentPass()
        peeling.SetTranslucentPass(translucent)

        volume = vtk.vtkVolumetricPass()
        overlay = vtk.vtkOverlayPass()

        lights = vtk.vtkLightsPass()
        opaqueSequence = vtk.vtkSequencePass()

        passes2 = vtk.vtkRenderPassCollection()
        passes2.AddItem(lights)
        passes2.AddItem(opaque)
        opaqueSequence.SetPasses(passes2)

        opaqueCameraPass = vtk.vtkCameraPass()
        opaqueCameraPass.SetDelegatePass(opaqueSequence)

        shadowsBaker = vtk.vtkShadowMapBakerPass()
        shadowsBaker.SetOpaqueSequence(opaqueCameraPass)
        shadowsBaker.SetResolution(SHADOWS_RESOLUTION)

        shadows = vtk.vtkShadowMapPass()
        shadows.SetShadowMapBakerPass(shadowsBaker)
        shadows.SetOpaqueSequence(opaqueSequence)

        seq = vtk.vtkSequencePass()
        passes = vtk.vtkRenderPassCollection()
        passes.AddItem(shadowsBaker)
        passes.AddItem(shadows)
        passes.AddItem(lights)
        passes.AddItem(peeling)
        passes.AddItem(volume)
        passes.AddItem(overlay)
        seq.SetPasses(passes)
        cameraP.SetDelegatePass(seq)

        self.renderer.SetPass(cameraP)


    def InitRenderWindow(self, stereo_on=False, axes_on=False, bkgnd=[0.8, 0.8, 0.8], res_x=600, res_y=600):
        """Sets up the visualisation environment"""
        self.renderer = vtk.vtkOpenGLRenderer()
        self.renderWindow = vtk.vtkRenderWindow()
        self.renderWindow.SetStereoCapableWindow(stereo_on)
        self.renderWindow.SetStereoRender(stereo_on)
        self.renderWindow.SetStereoTypeToSplitViewportHorizontal()
        self.renderWindow.AddRenderer(self.renderer)
        self.renderWindow.SetMultiSamples(0)
        self.renderWindow.SetAlphaBitPlanes(1)

        self.renderer.SetBackground(bkgnd)
        self.renderWindow.SetSize(res_x, res_y)
        self.renderer.ResetCamera()
        self.renderer.GetActiveCamera().Zoom(1)

        self.InitRenderPasses()

        axes = vtk.vtkAxesActor()
        axes.SetConeRadius(0.15)
        axes.AxisLabelsOn()
        if axes_on:
            self.renderer.AddActor(axes)

        self.camera = self.renderer.GetActiveCamera()
        self.windowFilter = vtk.vtkWindowToImageFilter()
        self.windowFilter.SetInput(self.renderWindow)
        self.windowFilter.ReadFrontBufferOff()
        self.writer = vtk.vtkJPEGWriter()
        self.writer.SetInputConnection(self.windowFilter.GetOutputPort())


    def InitLighting(self, mode=0):
        """Sets up lights in scene- default is even overhead lighting (glasshouse)"""
        light_array = []

        if mode == 0:
            light_corner_pos = np.array([-2, 2.5, -2])
            intensity = TOP_LIGHT_INTENS_MEAN
            for n in range(2):
                for light_idx in range(9):
                    light_array.append(vtk.vtkLight())
                    light_array[-1].SetLightTypeToSceneLight()
                    light_array[-1].SetIntensity(np.random.normal(loc=intensity, scale=LIGHT_INTENS_VAR))
                    light_array[-1].SetPosition(light_corner_pos[0], light_corner_pos[1], light_corner_pos[2])
                    light_array[-1].SetPositional(True)
                    light_array[-1].SetConeAngle(36)
                    light_array[-1].SetFocalPoint(0, 2, 0)
                    light_array[-1].SetDiffuseColor(1, 1, 0.984)
                    light_array[-1].SetAmbientColor(1, 1, 0.984)
                    light_array[-1].SetSpecularColor(1, 1, 0.984)
                    light_corner_pos = light_corner_pos + [((light_idx + 1) % 3 == 0 and light_idx > 0) * 2, 0,
                                                           2 - ((light_idx + 1) % 3 == 0 and light_idx > 0) * 6]
                light_corner_pos = np.array([-2, -2, -2])
                intensity = BOT_LIGHT_INTENS_MEAN
        elif mode == 1:
            light_array.append(vtk.vtkLight())
            light_array[-1].SetLightTypeToSceneLight()
            light_array[-1].SetIntensity(1)
            light_array[-1].SetPosition(0, 5, 0)
            light_array[-1].SetPositional(True)
            light_array[-1].SetConeAngle(30)
            light_array[-1].SetFocalPoint(0, 4, 0)
            light_array[-1].SetDiffuseColor(1, 1, 0.984)
            light_array[-1].SetAmbientColor(1, 1, 0.984)
            light_array[-1].SetSpecularColor(1, 1, 0.984)

        self.renderer.RemoveAllLights()
        for light in light_array:
            self.renderer.AddLight(light)


    def InitBackground(self, disp_pots=False):
        """Create randomised background"""
        if disp_pots:
            for plant in self.vtk_plant_list:
                plant_bag = vtk.vtkCubeSource()
                plant_pos = plant.StemActorList[0].GetPosition()
                plant_bag.SetCenter(plant_pos[0], -POT_SIZE / 2 + plant_pos[1], plant_pos[2])
                plant_bag.SetXLength(POT_SIZE)
                plant_bag.SetYLength(POT_SIZE)
                plant_bag.SetZLength(POT_SIZE)
                plant_bag.Update()
                plant_bag_mapper = vtk.vtkPolyDataMapper()
                plant_bag_mapper.SetInputData(plant_bag.GetOutput())
                self.BackgroundPolydataList.append(plant_bag)
                self.BackgroundMapperList.append(plant_bag_mapper)
        for idx, mapper in enumerate(self.BackgroundMapperList):
            poly_actor = vtk.vtkActor()
            poly_actor.SetMapper(mapper)
            self.BackgroundActorList.append(poly_actor)


    def AddActors(self):
        """Adds all actors to renderer"""
        for plant in self.vtk_plant_list:
            for actor in plant.StemActorList:
                self.renderer.AddActor(actor)
            for actor in plant.MeshActorList:
                self.renderer.AddActor(actor)
        for background_obj in self.BackgroundActorList:
            self.renderer.AddActor(background_obj)


    def InitInteractor(self):
        """Sets up default VTK 3D interactor"""
        self.windowInteractor = vtk.vtkRenderWindowInteractor()
        self.windowInteractor.SetRenderWindow(self.renderWindow)
        self.windowInteractor.Initialize()


    def SetPlantVisible(self, show_actors):
        """Shows/Hides plant actors"""
        for actor in (self.vtk_plant_list[0].StemActorList + self.vtk_plant_list[0].MeshActorList):
            actor.SetVisibility(show_actors)
        self.renderWindow.Render()


    def SetBackgroundVisible(self, show_actors):
        """Shows/Hides plant actors"""
        for actor in self.vtk_plant_rep.BackgroundActorList:
            actor.SetVisibility(show_actors)
        self.renderWindow.Render()


    def RenderPlant(self):
        """Render and interact - blocking"""
        self.renderWindow.Render()
        self.windowInteractor.Start()


    def MoveCamera(self, position, focal_point, up_vector):
        """Moves camera position/rotation"""
        self.camera.SetPosition(position)
        self.camera.SetFocalPoint(focal_point)
        self.camera.SetViewUp(up_vector)


    def UpdateWIFilter(self):
        """Updates window to image buffer"""
        self.renderWindow.Render()
        self.renderWindow.Modified()
        self.windowFilter.Update()
        self.windowFilter.Modified()


    def GetWindowImage(self):
        """Convert window image to numpy array and return"""
        self.UpdateWIFilter()
        vtk_image = self.windowFilter.GetOutput()
        width, height, _ = vtk_image.GetDimensions()
        vtk_array = vtk_image.GetPointData().GetScalars()
        components = vtk_array.GetNumberOfComponents()

        return vtk_to_numpy(vtk_array).reshape(height, width, components)


    def SaveCameraImage(self, dataset_path, data_subdir, file_name, save_stereo_lr=False):
        """Saves the currently rendered image to file"""
        self.UpdateWIFilter()
        self.writer.SetFileName(dataset_path + data_subdir + file_name + '.jpg')
        self.writer.Modified()
        self.writer.Write()