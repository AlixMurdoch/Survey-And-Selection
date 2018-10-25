import matplotlib.pyplot as plt
import csv

species = 'O'

fileName = species + '-output.csv'

R = 8.3144598
heatofFormation = 249.175

temps = []
CpList = []
enthalpyList = []
entropyList = []


with open(fileName, 'r') as readFile:
    readData = csv.reader(readFile, delimiter=',')
    next(readData)
    for row in readData:
        temps.append(float(row[0]))
        CpList.append(float(row[1]))# / R)
        enthalpyList.append( (float(row[2])))# + heatofFormation) / (R * float(row[0])))
        entropyList.append(float(row[3]))# / R)


plt.figure()
plt.plot(temps, CpList)
plt.title('Cp')

plt.figure()
plt.plot(temps, enthalpyList)
plt.title('enthalpy')

plt.figure()
plt.plot(temps, entropyList)
plt.title('entropy')

plt.show()

