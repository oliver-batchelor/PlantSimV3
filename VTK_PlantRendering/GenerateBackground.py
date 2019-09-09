import vtk
import numpy as np


POT_SIZE = 0.08
POT_TRANS_VAR = 0.006
POT_ROT_VAR = 2


class BackgroundScene():
    def __init__(self):
        self.BackgroundOriPoses = []
        self.BackgroundPolydataList = []
        self.BackgroundMapperList = []
        self.BackgroundActorList = []


    def GeneratePlantPots(self, plant_list):
        """Create randomised background"""
        for plant in plant_list:
            plant_bag = vtk.vtkCubeSource()
            plant_pos = plant.StemActorList[0].GetPosition()

            plant_bag.SetXLength(POT_SIZE)
            plant_bag.SetYLength(POT_SIZE)
            plant_bag.SetZLength(POT_SIZE)
            plant_bag.Update()
            plant_bag_mapper = vtk.vtkPolyDataMapper()
            plant_bag_mapper.SetInputData(plant_bag.GetOutput())

            center_offset = np.random.normal(loc=0, scale=POT_TRANS_VAR, size=3)
            center_rot = np.random.normal(loc=0, scale=POT_ROT_VAR, size=3)
            poly_actor = vtk.vtkActor()
            poly_actor.SetMapper(plant_bag_mapper)
            poly_actor.SetPosition(plant_pos[0] + center_offset[0],
                                   -POT_SIZE / 2 + plant_pos[1] + center_offset[1] + 0.006,
                                   plant_pos[2] + center_offset[2])
            poly_actor.SetOrientation(center_rot)
            self.BackgroundPolydataList.append(plant_bag)
            self.BackgroundMapperList.append(plant_bag_mapper)
            self.BackgroundActorList.append(poly_actor)