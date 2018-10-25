import matplotlib.pyplot as plt

cells = [30, 60, 120]

plotCells = 120
species = ['N', 'O']
rR = ['2_5e-5', '1e-4']
rR = rR[0]

folderName = 'nonReac/'
folderName = 'reac/'

CEAFileName = species[0] + '_' +  species[1] + '-CEA.dat'
CapFileName = species[0] + '_' + species[1] + '-Cap.dat'

CEAFileName = 'CEA.dat'
CapFileName = 'Cap.dat'

CEA = {}
with open(folderName + CEAFileName, 'r') as CEAFile:
    for row in CEAFile:
        data = row.split(' ')
        ID = (data[0], int(data[2]))
        CEA[ID] = float(data[8])

Cap = {}
with open(folderName + CapFileName, 'r') as CapFile:
    for row in CapFile:
        data = row.split(' ')
        ID = (data[0], int(data[2]))
        Cap[ID] = float(data[8])


plt.xlabel('Velocity ($km/s$)')
plt.ylabel('Proportional Shock Standoff ($D/d$)')
#plt.title('Lobb-Sphere Shock Standoff Comparison Using\nNon-Reacting Monatomic Species')
plt.title('Lobb-Sphere Shock Standoff Comparison')

CEAPoints = [[],[]]
CapPoints = [[],[]]

for keys, results in CEA.items():
    uIdx = keys[0].index('u')
    vel = float(keys[0][uIdx+1:])
    
    if keys[0][:8] == 'pR2.5e-5':
        plt.figure(1)
    if keys[1] == plotCells or plotCells == None:
        CEAPoints[0].append(vel)
        CEAPoints[1].append(results)
    elif keys[1] == plotCells or plotCells == None:
        plt.plot(vel, results, 'g+')
    elif keys[1] == plotCells or plotCells == None:
        plt.plot(vel, results, 'r+')


for keys, results in Cap.items():
    uIdx = keys[0].index('u')
    vel = float(keys[0][uIdx+1:])
    if keys[0][:8] == 'pR2.5e-5':
        plt.figure(1)

    if keys[1] == plotCells or plotCells == None:
        CapPoints[0].append(vel)
        CapPoints[1].append(results)
    elif keys[1] == plotCells or plotCells == None:
        plt.plot(vel, results, 'go')
    elif keys[1] == plotCells or plotCells == None:
        plt.plot(vel, results, 'ro')

Lobb = [[5.1, 5.95, 5.55, 6.45, 4.8, 5.15], [0.0545, 0.0505, 0.0499, 0.048, 0.048, 0.0451]]
Zandar = [[8.66, 9.71], [0.036, 0.0345]]

plt.figure(1)
plt.plot(CEAPoints[0], CEAPoints[1], 'kx', markersize=7, label='CEA Coefficients')
plt.plot(CapPoints[0], CapPoints[1], 'r+', markersize=10, label='This work')
plt.plot(Lobb[0], Lobb[1], 'b^', markersize=7, label='Lobb Experimental Data')
plt.plot(Zandar[0], Zandar[1], 'bs', markersize=7, label='Zandar Experimental Data (Mean)')
plt.ylim([0.0, 0.13])
plt.ylim([0.0, 0.08])
plt.legend(loc = 1)
plt.xlim([2.0, 11.0])
plt.grid()
plt.show()
