import matplotlib.pyplot as plt
import csv


def plotCp(species, tempRange, userCp, ionizationLowering, Capitelli=True, CEA=False):
    
    Capitelli1000 = [[100, 200, 500, 700, 1000, 2000, 3000, 5000, 10000, 12000, 13000, 14000,
                    15000, 16000, 17000, 18000, 19000, 20000, 22000, 23000, 24000, 25000, 26000, 27000, 28000,
                    30000, 34000, 40000, 44000, 50000],
                [23.70, 22.74, 21.26, 21.04, 20.91, 20.83, 20.94, 21.80, 23.29, 24.35, 25.74, 28.15,
                    31.89, 37.19, 44.04, 52.19, 61.11, 70.00, 84.37, 88.54, 90.36, 89.94, 87.66, 84.00, 79.46,
                    69.37, 51.52, 35.92, 30.6, 26.30]]
    
    Capitelli500 = [[100, 1000, 10000, 15000, 20000, 22000, 25000, 30000, 40000],
                    [23.7, 20.91, 23.45, 43.49, 109.1, 114.7, 94.51, 54.89, 28.3]]
    
    plt.figure()
    
    CpPlot = userCp.calcCpRange(tempRange)
    plt.plot(tempRange, CpPlot, 'k.', label = 'My Results')

    if ionizationLowering == 250 and Capitelli:
        with open('capitelliCpData.csv', 'r') as capReadFile:
            capReadData = csv.reader(capReadFile, delimiter=',')
            speciesList = next(capReadData)[1:]
            capCpList = {species: [[],[]] for i,species in enumerate(speciesList)}
            for row in capReadData:
                for i, speciesInp in enumerate(speciesList):
                    capCpList[speciesInp][0].append(float(row[0]))
                    capCpList[speciesInp][1].append(float(row[i + 1]))
        plt.plot(capCpList[species.name][0], capCpList[species.name][1], 'rx', label = 'Capitelli')
    
    elif ionizationLowering == 500 and Capitelli:
        plt.plot(Capitelli500[0], Capitelli500[1], 'rx', label = 'Capitelli')

    elif ionizationLowering == 1000 and Capitelli:
        plt.plot(Capitelli1000[0], Capitelli1000[1], 'rx', label = 'Capitelli')

    plt.legend()
    plt.show()


if __name__ == '__main__':
    species = 'O'
    tempRange = range(200, 50000, 10)
    userCp = 3
    ionizationLowering = 250

    plotCp(species, tempRange, userCp, ionizationLowering)


