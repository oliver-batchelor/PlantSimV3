import vtk
import VTK_PlantRendering.ConvertToVTKPlant as CVVTK
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np


class plantVTKDataDisplay():
    """Class controlling render environment and interactors"""
    def __init__(self, vtk_represenation):
        self.vtk_plant_rep = vtk_represenation
        self.renderWindow = None
        self.renderer = None
        self.windowInteractor = None
        self.windowFilter = None
        self.camera = None
        self.writer = None


    def InitRenderWindow(self, stereo_on=False, axes_on=False, bkgnd=[0.8, 0.8, 0.8], res_x=600, res_y=600):
        """Sets up the visualisation environment"""
        self.renderer = vtk.vtkRenderer()
        self.renderWindow = vtk.vtkRenderWindow()
        self.renderWindow.SetStereoCapableWindow(stereo_on)
        self.renderWindow.SetStereoRender(stereo_on)
        self.renderWindow.SetStereoTypeToSplitViewportHorizontal()
        self.renderWindow.AddRenderer(self.renderer)

        self.renderer.SetBackground(bkgnd)
        self.renderWindow.SetSize(res_x, res_y)
        self.renderer.ResetCamera()
        self.renderer.GetActiveCamera().Zoom(1)

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


    def InitLighting(self):
        """Sets up lights in scene- default is even overhead lighting (glasshouse)"""
        light_array = []
        light_corner_pos = np.array([-2, 2, -2])
        intensity = 0.25
        for n in range(2):
            for light_idx in range(9):
                light_array.append(vtk.vtkLight())
                light_array[-1].SetLightTypeToSceneLight()
                light_array[-1].SetIntensity(intensity)
                light_array[-1].SetPosition(light_corner_pos[0], light_corner_pos[1], light_corner_pos[2])
                light_array[-1].SetPositional(True)
                light_array[-1].SetConeAngle(60)
                light_array[-1].SetFocalPoint(0, 0, 0)
                light_array[-1].SetDiffuseColor(1, 1, 0.984)
                light_array[-1].SetAmbientColor(1, 1, 0.984)
                light_array[-1].SetSpecularColor(1, 1, 0.984)
                light_corner_pos = light_corner_pos + [((light_idx + 1) % 3 == 0 and light_idx > 0) * 2, 0,
                                                       2 - ((light_idx + 1) % 3 == 0 and light_idx > 0) * 6]
            light_corner_pos = np.array([-2, -2, -2])
            intensity = 0.1

        self.renderer.RemoveAllLights()
        for light in light_array:
            self.renderer.AddLight(light)


    def AddActors(self):
        """Adds all actors to renderer"""
        for actor in self.vtk_plant_rep.StemActorList:
            self.renderer.AddActor(actor)
        for actor in self.vtk_plant_rep.MeshActorList:
            self.renderer.AddActor(actor)
        for actor in self.vtk_plant_rep.BackgroundActorList:
            self.renderer.AddActor(actor)


    def InitInteractor(self):
        """Sets up default VTK 3D interactor"""
        self.windowInteractor = vtk.vtkRenderWindowInteractor()
        self.windowInteractor.SetRenderWindow(self.renderWindow)
        self.windowInteractor.Initialize()


    def SetPlantVisible(self, show_actors):
        """Shows/Hides plant actors"""
        for actor in (self.vtk_plant_rep.StemActorList + self.vtk_plant_rep.MeshActorList):
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