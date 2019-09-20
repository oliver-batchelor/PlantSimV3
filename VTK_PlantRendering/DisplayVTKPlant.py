import vtk
import VTK_PlantRendering.ConvertToVTKPlant as CVVTK
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np

TOP_LIGHT_INTENS_MEAN = 0.3
BOT_LIGHT_INTENS_MEAN = 0.08
LIGHT_INTENS_VAR = 0.03

SHADOW_RENDER_RES = 4096

class plantVTKDataDisplay():
    """Class controlling render environment and interactors"""
    def __init__(self, plant_list):
        self.vtk_plant_list = plant_list

        self.renderWindow = None
        self.renderer = None
        self.bakerPass = None
        self.windowInteractor = None
        self.windowFilter = None
        self.camera = None
        self.writer = None


    def ClearRenderer(self):
        """Removes all actors from renderer"""
        self.renderer.RemoveAllViewProps()


    def InitRenderPasses(self):
        """Initiate and add realism passes to renderer"""
        cameraP = vtk.vtkCameraPass()
        opaque = vtk.vtkOpaquePass()
        peeling = vtk.vtkDepthPeelingPass()
        peeling.SetMaximumNumberOfPeels(1000)
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
        shadowsBaker.SetResolution(SHADOW_RENDER_RES)

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


    def InitLighting(self, mode=0, intensity_mul=1):
        """Sets up lights in scene- default is even overhead lighting (glasshouse)"""
        light_array = []

        if mode != 0:
            self.renderer.RemoveAllLights()
            if mode == 1:
                light_corner_pos = np.array([-2, 2.5, -2])
                intensity = TOP_LIGHT_INTENS_MEAN
                for n in range(2):
                    for light_idx in range(9):
                        light_array.append(vtk.vtkLight())
                        light_array[-1].SetLightTypeToSceneLight()
                        light_array[-1].SetIntensity(intensity_mul * np.random.normal(loc=intensity, scale=LIGHT_INTENS_VAR))
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
            elif mode == 2:
                light_array.append(vtk.vtkLight())
                light_array[-1].SetLightTypeToSceneLight()
                light_array[-1].SetIntensity(intensity_mul)
                light_array[-1].SetPosition(0, 5, 0)
                light_array[-1].SetPositional(True)
                light_array[-1].SetConeAngle(30)
                light_array[-1].SetFocalPoint(0, 4, 0)
                light_array[-1].SetDiffuseColor(1, 1, 0.984)
                light_array[-1].SetAmbientColor(1, 1, 0.984)
                light_array[-1].SetSpecularColor(1, 1, 0.984)


        for light in light_array:
            self.renderer.AddLight(light)


    def AddActors(self, bckgnd_actors=[]):
        """Adds all actors to renderer"""
        for plant in self.vtk_plant_list:
            for actor in plant.StemActorList:
                self.renderer.AddActor(actor)
            for actor in plant.MeshActorList:
                self.renderer.AddActor(actor)
        for bckgnd_a in bckgnd_actors:
            self.renderer.AddActor(bckgnd_a)


    def InitInteractor(self, observer=None):
        """Sets up default VTK 3D interactor"""
        self.windowInteractor = vtk.vtkRenderWindowInteractor()
        self.windowInteractor.SetRenderWindow(self.renderWindow)
        style = vtk.vtkInteractorStyleJoystickCamera()
        self.windowInteractor.SetInteractorStyle(style)
        self.windowInteractor.Initialize()
        if observer != None:
            self.windowInteractor.AddObserver(vtk.vtkCommand.MouseMoveEvent, observer, 1)


    def SetPlantVisible(self, show_actors):
        """Shows/Hides plant actors"""
        for actor in (self.vtk_plant_list[0].StemActorList + self.vtk_plant_list[0].MeshActorList):
            actor.SetVisibility(show_actors)
        self.renderWindow.Render()


    def SetActorsVisible(self, actors_list, show_actors):
        """Shows/Hides plant actors"""
        for actor in actors_list:
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