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
            self.tubePntSets = plant.tubePntSets
            self.tubeRadSets = plant.tubeRadSets
            self.tubeColSets = plant.tubeColSets
            self.triMeshPntSets = plant.triPntSets
            self.triMeshPntIndxSets = plant.triPntIndxSets
            self.triMeshColSets = plant.triColSets


    def SavePlantFile(self, filename):
        """Saves plant data file to disk"""
        with open(filename, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)