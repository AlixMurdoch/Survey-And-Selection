import outputDat
import plotting
import CEAData


import csv
from scipy import optimize, integrate
import matplotlib.pyplot as plt
import numpy as np
import math

class CoefficientList:
    def __init__(self, species):
        self.name = species
        self.cpList = {}
        self.enthalpyList = {}
        self.entropyList = {}
        self.luaString = 'db[\'' + self.name + '\'].thermoCoeffs = {\n'

    def setCp(self, tempRange, coefficients):
        self.cpList[tempRange] = coefficients

    def setEnthalpy(self, tempRange, coefficients):
        self.enthalpyList[tempRange] = coefficients

    def setEntropy(self, tempRange, coefficients):
        self.entropyList[tempRange] = coefficients

speciesList = ['O', 'O+', 'O2', 'O2+', 'N', 'N+', 'N2', 'N2+', 'NO', 'NO+', 'e']
heatOfFormationList = {'O': 249.175, 'O+': 1568.787, 'O2': 0.00, 'O2+': 1171.828, \
                       'N': 462.68, 'N+': 1882.128, 'N2': 0.00, 'N2+': 1509.508, \
                       'NO': 91.271, 'NO+':990.810, 'e': 0.}

temps = {species: [] for i,species in enumerate(speciesList)}
CpList = {species: [] for i,species in enumerate(speciesList)}
enthalpyList = {species: [] for i,species in enumerate(speciesList)}
entropyList = {species: [] for i,species in enumerate(speciesList)}

R = 8.3144598

speciesList = ['O', 'O+', 'N', 'N+']

for species in speciesList:
    fileName = species + '-output.csv'
    heatofFormation = heatOfFormationList[species]

    with open(fileName, 'r') as readFile:
        readData = csv.reader(readFile, delimiter=',')
        next(readData)
        for row in readData:
            temps[species].append(float(row[0]))
            CpList[species].append(float(row[1]) / R)

            enthalpyList[species].append( (float(row[2])/1000. + heatofFormation) / (R * float(row[0])))
            entropyList[species].append(float(row[3]) / R)


regions = {'O':   [0, 6000, 14000, 25000, 50000],  
           'O+':  [0, 14000, 32000, 50000],        
           'N':   [0, 6000, 13000, 21000, 50000],  
           'N+':  [0, 6000, 24000, 38000, 50000],  
           'O2':  [0, 600, 4000, 10000, 50000],     
           'O2+': [0, 600, 4000, 14000, 50000],    
           'N2':  [0, 600, 10000, 20000, 50000],    
           'N2+': [0, 600, 5000, 12000, 50000],     
           'NO':  [0, 600, 4000, 15000, 50000],     
           'NO+': [0, 600, 12000, 20000, 50000],    
           'e':   [0, 50000]}                       
           

def CpFunc(x, a, b, c, d, e, f, g):
    y = a*x**-2. + b*x**-1 + c + d*x + e*x**2. + f*x**3. + g*x**4.
    return y

def makeEnthalpy(a, b, c, d, e, f, g):
    def enthalpy(x, intCoeff):
        y = -a*x**-2. + b*np.log(x)/x + c + d*x/2. + e*x**2./3. + f*x**3./4. + g*x**4./5. + intCoeff/x
        return y/x
    return enthalpy

def enthFunc(x, a, b, c, d, e, f, g, intCoeff):
    y = -a*x**-2. + b*np.log(x)/x + c + d*x/2. + e*x**2./3. + f*x**3./4. + g*x**4./5. + intCoeff/x
    return y/x

def makeEntropy(a, b, c, d, e, f, g):
    def entropy(x, intCoeff):
        y = (-a*x**-2.)/2. - b*x**-1. + c*np.log(x) + d*x + (e*x**2.)/2. + (f*x**3.)/3. + (g*x**4.)/4. + intCoeff
        return y
    return entropy

def entrFunc(x, a, b, c, d, e, f, g, intCoeff):
    y = (-a*x**-2.)/2. - b*x**-1. + c*np.log(x) + d*x + (e*x**2.)/2. + (f*x**3.)/3. + (g*x**4.)/4. + intCoeff
    return y

allSpeciesList = {}
for species in speciesList:
    speciesObj = CoefficientList(species)

    speciesRegions = sorted(regions[species])

    CpSpecies = CpList[species]

    enthalpySpecies = sorted(enthalpyList[species])
    entropySpecies = sorted(entropyList[species])

    for i in range(len(speciesRegions[:-1])):

        tempsMin = min([x for x in temps[species] if x >= speciesRegions[i]])
        tempsMax = max([x for x in temps[species] if x <= speciesRegions[i+1]])

        tempsMinIdx = temps[species].index(tempsMin)
        tempsMaxIdx = temps[species].index(tempsMax) + 1

        tempRegion = temps[species][tempsMinIdx:tempsMaxIdx]
        cpRegion = CpSpecies[tempsMinIdx:tempsMaxIdx]
        
        cpWeight = np.empty(len(cpRegion))
        cpWeight.fill(1000)
        cpWeight[0] = cpWeight[-1] = 1e-5

        p, e = optimize.curve_fit(CpFunc, tempRegion, cpRegion, sigma = cpWeight)
        
        speciesObj.setCp((tempsMin, tempsMax), p)

        ######### Enthalpy ###########
        lenTemps = len(temps[species][:tempsMaxIdx])

        if i == 0:
            boundaryTempEnth = temps[species][int(lenTemps)]
            boundaryEnth = enthalpySpecies[int(lenTemps)]
            Ep, Ee = optimize.curve_fit(makeEnthalpy(*p), boundaryTempEnth, boundaryEnth)
            prevEnth = enthFunc(tempsMax, *p, *Ep)
        else:
            
            boundaryEnth = prevEnth
            Ep, Ee = optimize.curve_fit(makeEnthalpy(*p), tempsMin, boundaryEnth)
            prevEnth = enthFunc(tempsMax, *p, *Ep)
        speciesObj.setEnthalpy((tempsMin, tempsMax), Ep)
        
        ######### Entropy ###########
        
        if i == 0:
            boundaryTempEntr = temps[species][int(lenTemps)]
            boundaryEntr = entropySpecies[int(lenTemps)]
            Sp, Se = optimize.curve_fit(makeEntropy(*p), boundaryTempEntr, boundaryEntr)
            prevEntr = entrFunc(tempsMax, *p, *Sp)
            
        else:
            boundaryEntr = prevEntr
            Sp, Se = optimize.curve_fit(makeEntropy(*p), tempsMin, boundaryEntr)
            prevEntr = entrFunc(tempsMax, *p, *Sp)
        
        speciesObj.setEntropy((tempsMin, tempsMax), Sp)
        
    allSpeciesList[species] = speciesObj


outputDat.outputDat(allSpeciesList)

CEA = CEAData.getCEA()
#CEA = None



plotSpecies = 'N+'

speciesList = ['O', 'O+', 'O2', 'O2+', 'N', 'N+', 'N2', 'N2+', 'NO', 'NO+', 'e']
CapTemps = []
CapList = {species: [] for i,species in enumerate(speciesList)}

with open('CapitelliCp.csv', 'r') as readFile:
    readData = csv.reader(readFile, delimiter=',')
    next(readData)
    for row in readData:
        CapTemps.append(float(row[0]))
        for i, species in enumerate(speciesList):
            CapList[species].append(float(row[i + 1]) / R)

if plotSpecies != '' and plotSpecies != None:
    speciesList = [plotSpecies]

for species in speciesList:

    
    plt.figure()
    plt.title(species + ' Cp')
    plt.title('Cp Polynomial Fit for ' + plotSpecies)
    Capitelli = [temps[species], CpList[species]]
    plotTemps = range(200, 50000, 100)
    plotting.plotCp(species, plotTemps, allSpeciesList, CEA, Capitelli)
    plt.ylim([min(Capitelli[1])*0.8, max(Capitelli[1])*1.2])
    plt.xlabel('Temperature (K)')
    plt.ylabel('Cp/R')
    
    plt.figure()
    plt.title(species + ' Enthalpy')
    plt.title('Enthalpy Polynomial Fit for ' + plotSpecies)
    #plt.title('Comparison of Enthalpy fits for N+')
    Capitelli = [temps[species], enthalpyList[species]]
    plotTemps = range(200, 50000, 100)
    plotting.plotEnthalpy(species, plotTemps, allSpeciesList, CEA, Capitelli)
    plt.xlabel('Temperature (K)')
    plt.ylabel('H/RT')
    plt.ylim([0.0, 0.08])
    plt.xlim([0, 50000])
    plt.legend(loc=1)

    plt.figure()
    plt.title(species + ' Entropy')
    plt.title('Entropy Polynomial Fit for ' + plotSpecies)
    Capitelli = [temps[species], entropyList[species]]
    plotTemps = range(200, 50000, 100)
    plotting.plotEntropy(species, plotTemps, allSpeciesList, CEA, Capitelli)
    plt.xlabel('Temperature (K)')
    plt.ylabel('S/R')
    
    plt.show()
    

    


