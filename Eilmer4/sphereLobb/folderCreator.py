import os
import shutil

velocityList = xrange(12000, 15000, 1000)
rhoRList = ['2.5e-5']#'2.5e-5',

cwd = os.getcwd()
childDir = cwd + '/runsNewGeom/'
masterDir = cwd + '/master/'

for i in xrange(len(rhoRList)):
	rhoR = rhoRList[i]
	for j in xrange(len(velocityList)):
		velocity = velocityList[j]
		
		velocityKm = float(velocity)/1000.
		velocityStr = "%.1f" % velocityKm
		
		folderName = 'pR' + rhoR + '-' + 'u' + velocityStr

		print(folderName)

		intermedFolder = 'pR' + rhoR + '/'

		folderDir = childDir + intermedFolder + folderName

		shutil.copytree(masterDir, folderDir)
		
		#print(str(velocity))

		with open(folderDir + '/inf-prop.lua', 'w') as propertyFile:
			propertyFile.write('u_inf = ' + str(velocity) + ' -- m/s\n')
			propertyFile.write('rR = ' + rhoR + ' -- kg/m^2')
		
		with open(folderDir + '/run-tinaroo.qsub', 'w') as qsubFile:
			
			qsubString = '#PBS -S /bin/bash\n' +\
							 '#PBS -N ' + folderName + '\n' +\
							 '#PBS -A UQ-EAIT-MechMining\n' +\
							 '#PBS -l select=1:ncpus=4:mpiprocs=4:mem=5g\n' +\
							 '#PBS -l walltime=02:00:00\n\n' +\
							 'module purge\n' +\
							 'module load gnu\n' +\
							 'module load openmpi_ib/1.8.4\n\n' +\
							 'cd $PBS_O_WORKDIR\n'+\
							 'lua run-calculation-in-stages.lua\n'

			qsubFile.write(qsubString)
			




