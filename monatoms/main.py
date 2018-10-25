import spectra
import thermo
import util
import species
import plotter

import matplotlib.pyplot as plt

constants = util.Constants()

#maximum principal quantum number to model
nMax = 40

#list of all species to model. Ions required for extrapolation
speciesList = ['C', 'C+', 'N', 'N+', 'N++', 'O', 'O+', 'O++', 'F+', 'F++', 'F+++', 'Ne++', 'Ne+++']

#list of species to generate thermodynamic properties for
useList = ['N+']#, 'O+', 'N', 'N+']

allSpecies = []
for speciesStr in speciesList:
    speciesObj = species.Species(speciesStr)
    allSpecies.append([speciesObj, spectra.readNISTSpectra(speciesObj)])
    
for useName in useList:
    #load species objects containing their atomic data
    use = species.Species(useName)

    #create class containing experimental spectra data
    NIST = spectra.readNISTSpectra(use)

    #create class containing all expected energy levels
    theory = spectra.calculateExpectedStates(use, nMax)

    #fill theory class with observed energy levels
    calcEnergy = spectra.CalcEnergy(nMax, use, NIST, theory, allSpecies)

    #fill gaps in observed spectra using extrapolation
    completeLevels, calculatedLevels = calcEnergy.populateTheory()

    #sort levels by energy
    #a filename may be specified to dump this data to
    completeLevels = spectra.sortSpectra(completeLevels, None)
    #calculatedLevels = spectra.sortSpectra(calculatedLevels, None)
    #NISTSorted = spectra.sortSpectra(NIST, None)

    #range of temperatures to produce thermodynamic properties for
    tempRange = range(200, 50100, 100)


    #choice of ionization potential lowering
    #may also specify 'DebyeHuckel' and a number density
    ionizationLowering = 500
    numberDensity = None
    #ionizationLowering = 'DebyeHuckel'
    #numberDensity = 3E23
    

    #setup thermodynamic object
    thermoObj = thermo.Thermo(use, completeLevels, tempRange, ionizationLowering, numberDensity)

    
    #plotter.plotCp(use, tempRange, thermoObj, ionizationLowering, True, False)

    #calculate and return thermodynamic property tables
    thermoProps = thermoObj.calcThermoPropsRange(tempRange)
    CpArr = thermoProps[0]
    HArr = thermoProps[1]
    SArr = thermoProps[2]

    #write properties to files for polynomial fitting
    fileName = 'polyFits/' + useName + '-output.csv'
    with open(fileName, 'w') as exportFile:
        exportFile.write('T, Cp, H, S\n')
        for i, temp in enumerate(tempRange):
            string = str(temp) + ', ' + str(CpArr[i]) + ', ' + str(HArr[i]) + ', ' + str(SArr[i]) + '\n'
            exportFile.write(string)
