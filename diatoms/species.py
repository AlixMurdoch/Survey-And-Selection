import util

import json

class Species:
    def __init__(self, species):
        self.name = species
        self.constants = util.Constants()

        if '2' in self.name:
            self.readDiatomProperties()
        else:
            self.readAtomProperties()

    def readDiatomProperties(self):
        jsonFilePath = 'diatomData.json'
        with open(jsonFilePath, 'r') as jsonFile:
            jsonData = json.load(jsonFile)

            speciesProperties = jsonData[self.name]
            self.charge = speciesProperties['charge']
            self.noElectrons = speciesProperties['noElectrons']
            self.reducedMass = speciesProperties['reducedMass'] * self.constants.atomicMassUnits
            self.symmetryFactor = speciesProperties['symmetryFactor']
            self.spectralData = speciesProperties['spectralData']
            self.nMax = len(self.spectralData['Te'])

    def readAtomProperties(self):
        jsonFilePath = util.getSpectralDataDir() + '/atomData.json'
        with open(jsonFilePath, 'r') as jsonFile:
            jsonData = json.load(jsonFile)

            speciesProperties = jsonData[self.name]
            self.atomicNumber = speciesProperties['atomicNumber']
            self.charge = speciesProperties['charge']
            self.noElectrons = speciesProperties['noElectrons']
            self.Io = speciesProperties['Io'] * self.constants.eVToCm_1
            self.cores = []
            for core in speciesProperties['cores'].items():
                self.cores.append(core[1])
    

#He = Species('He')
