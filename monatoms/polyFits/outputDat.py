import numpy as np

def outputDat(allSpeciesList):
    with open('output.dat', 'w') as outputFile:
        for _, species in sorted(allSpeciesList.items()):
            outputFile.writelines('\n\n~~~~~~~' + species.name + '~~~~~~~\n')
            speciesRegions = sorted(species.cpList.items())

            for key, _ in speciesRegions:
                tempRegion = [key[0], key[1]]

                tempsMin = min(tempRegion)
                tempsMax = max(tempRegion)

                p = species.cpList[(tempsMin, tempsMax)]
                Ep = species.enthalpyList[(tempsMin, tempsMax)]
                Sp = species.entropyList[(tempsMin, tempsMax)]

                coeffStr = str(tempsMin) + '-' + str(tempsMax) 
                spaces = 16 - len(coeffStr)
                for i in range(spaces):
                    coeffStr += ' '
                coeffStr += '| \n'
                outputFile.writelines(coeffStr)

                def writeCoeff(coeff):
                    coeffStr = ''
                    if coeff > 0:
                        coeffStr += "       {:.9e}".format(coeff) + ',\n'
                    else:
                        coeffStr += "      {:.9e}".format(coeff) + ',\n'
                    outputFile.writelines(coeffStr)      

                for coeffSeries in [p, Ep, Sp]:
                    for coeff in coeffSeries:
                        writeCoeff(coeff)


