import pickle


class PlantData():
    """Class to represent geometric plant structure as generalised sets of tube-splines and triangle strip mesh"""
    def __init__(self):
        self.worldPos = []
        self.worldOri = []

        self.numTubeSets = 0
        self.tubePntSets = []
        self.tubeRadSets = []
        self.tubeColSets = []
        self.tubeConnIdxs = []
        self.tubeSetLabels = []

        self.numMeshSets = 0
        self.triMeshOris = []
        self.triMeshPntSets = []
        self.triMeshPntIndxSets = []
        self.triMeshColSets = []
        self.triMeshConnIdxs = []
        self.triMeshSetLabels = []

        self.textureMatrix = []
        self.textureMapSet = []

    def LoadPlantFile(self, filename):
        """Loads plant data file from disk"""
        with open(filename, 'rb') as input:
            plant = pickle.load(input)

            self.numTubeSets = len(plant.tubePntSets)
            self.tubePntSets = plant.tubePntSets
            self.tubeRadSets = plant.tubeRadSets
            self.tubeColSets = plant.tubeColSets
            self.tubeConnIdxs = plant.tubeConnIdxs
            self.tubeSetLabels = plant.tubeSetLabels

            self.numMeshSets = len(plant.triMeshPntSets)
            self.triMeshPntSets = plant.triMeshPntSets
            self.triMeshPntIndxSets = plant.triMeshPntIndxSets
            self.triMeshColSets = plant.triMeshColSets
            self.triMeshConnIdxs = plant.triMeshConnIdxs
            self.triMeshSetLabels = plant.triMeshSetLabels


    def SavePlantFile(self, filename):
        """Saves plant data file to disk"""
        with open(filename, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)