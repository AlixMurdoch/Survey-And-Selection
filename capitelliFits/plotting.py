import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import curve_fit


def CpFunc(x, a, b, c, d, e, f, g):
    y = a*x**-2. + b*x**-1 + c + d*x + e*x**2. + f*x**3. + g*x**4.
    return y

def enthFunc(x, a, b, c, d, e, f, g, intCoeff):
    y = -a*x**-2. + b*np.log(x)/x + c + d*x/2. + e*x**2./3. + f*x**3./4. + g*x**4./5. + intCoeff/x
    return y/x

def entrFunc(x, a, b, c, d, e, f, g, intCoeff):
    y = (-a*x**-2.)/2. - b*x**-1. + c*np.log(x) + d*x + (e*x**2.)/2. + (f*x**3.)/3. + (g*x**4.)/4. + intCoeff
    return y


def plotCpSpecies(allSpeciesList, plotSpecies, data):
    for species in allSpeciesList:
        speciesRegions = sorted(species.cpList.items())

        for key, _ in speciesRegions:
            tempRegion = [key[0], key[1]]

            tempsMin = min(tempRegion)
            tempsMax = max(tempRegion)
            xd = np.linspace(tempsMin, tempsMax, 50)

            p = species.cpList[(tempsMin, tempsMax)]
            #Ep = species.enthalpyList[(tempsMin, tempsMax)]

            if species == plotSpecies or plotSpecies == 'all':
                plt.plot(xd, CpFunc(xd, *p))
                plt.title(species)

        for series in data:
            plt.plot(series[0], series[1], 'kx')


def plotEnthSpecies(allSpeciesList, plotSpecies, data):
    for species in allSpeciesList:
        speciesRegions = sorted(species.cpList.items())

        for key, _ in speciesRegions:
            tempRegion = [key[0], key[1]]

            tempsMin = min(tempRegion)
            tempsMax = max(tempRegion)
            xd = np.linspace(tempsMin, tempsMax, 50)

            p = species.cpList[(tempsMin, tempsMax)]
            Ep = species.enthalpyList[(tempsMin, tempsMax)]

            if species == plotSpecies or plotSpecies == 'all':
                plt.plot(xd, enthFunc(xd, *p, *Ep))
                plt.title(species)

        for series in data:
            plt.plot(series[0], series[1], 'kx')
    

def getCoeffs(temp, tempList):
    for key, value in tempList:
        if temp < key[1] and temp >= key[0]:
            return value
    return 0

def plotCp(speciesStr, plotTemps, allSpeciesList, CEA, Capitelli):
    yplot = []
    species = allSpeciesList[speciesStr]
    for temp in plotTemps:
        CpCoeff = getCoeffs(temp, species.cpList.items())
        yplot.append(CpFunc(temp, *CpCoeff))

    plt.plot(plotTemps, yplot, 'k.', label='Data Series')
    plt.plot(Capitelli[0], Capitelli[1], 'r-', label = 'Capitelli (2005)')

    if CEA:
        plotTemps = range(300, 20000)
        yplot = []
        for temp in plotTemps:
            coeff = getCoeffs(temp, CEA.items())
            yplot.append(CpFunc(temp, *coeff[:7]))
        plt.plot(plotTemps, yplot, label = 'CEA (2002)')

    plt.grid()
    plt.legend(loc = 1)


def plotEnthalpy(speciesStr, plotTemps, allSpeciesList, CEA, Capitelli):
    if CEA:
        yplot = []
        plotTemps = range(300, 20000)
        for temp in plotTemps:
            coeff = getCoeffs(temp, CEA.items())
            yplot.append(enthFunc(temp, *coeff[:8]))
        plt.plot(plotTemps, yplot, '-.', label = 'CEA (2002)')
    plt.plot(Capitelli[0], Capitelli[1], 'k.', label = 'Capitelli (2005)')
    plt.grid()

    yplot = []
    species = allSpeciesList[speciesStr]
    for temp in plotTemps:
        CpCoeff = getCoeffs(temp, species.cpList.items())
        enthalpyCoeff = getCoeffs(temp, species.enthalpyList.items())
        yplot.append(enthFunc(temp, *CpCoeff, enthalpyCoeff))
    plt.plot(plotTemps, yplot, 'r-', label='This work')
    plt.legend(loc = 4) 

def plotEntropy(speciesStr, plotTemps, allSpeciesList, CEA, Capitelli):
    species = allSpeciesList[speciesStr]

    data = []
    tempAr = range(300, 50000, 100)
    for temp in tempAr:
        CpCoeff = getCoeffs(temp, species.cpList.items())
        entropyCoeff = getCoeffs(temp, species.entropyList.items())
        data.append(entrFunc(temp, *CpCoeff, entropyCoeff))
    plt.plot(tempAr, data, 'k.', label='Data series')
    
    yplot = []
    for temp in plotTemps:
        CpCoeff = getCoeffs(temp, species.cpList.items())
        entropyCoeff = getCoeffs(temp, species.entropyList.items())
        yplot.append(entrFunc(temp, *CpCoeff, entropyCoeff))
    plt.plot(plotTemps, yplot, 'r-', label='Polynomial fit')

    if CEA:
        yplot = []
        plotTemps = range(300, 20000)
        for temp in plotTemps:
            coeff = getCoeffs(temp, CEA.items())
            coeff = coeff[:7] + [coeff[-1]]
            yplot.append(entrFunc(temp, *coeff))
        plt.plot(plotTemps, yplot, label = 'CEA (2002)')
    plt.grid()
    plt.legend(loc = 4)
