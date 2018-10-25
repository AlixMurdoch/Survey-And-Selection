Eilmer4 code used to validate the difference in simulated results between the CEA
and Capitelli scripts

steps for producing results are as follows:


1/ 	specify the inflow conditions in the folderCreator.py file, lines 3 and 4
	These should be contained in a list, with the script creating all combinations 
	of these properties
	
2/	give your output directory a name on line 7

3/	choose either the CEA or Capitelli based coefficients by changing the filename
	to just <air-11sp.lua>

4/	check your meshing and timestep values are reasonable in the /master/lobb.lua file
	check your qsub file is ok, its much easier to fix this now than later!

5/	run folderCreator.py

6/	on tinaroo copy the runTinaroo.sh file to /ouputDirectory/density/runTinaroo.sh
	running this will run all run-tinaroo.qsub files, adding each job to the queue
	note there's a maximum of around 20 concurrent jobs you may run
	
7/ 	on a workstation you may wish to stage jobs to run in parallel one after the other
	in this case do the same as 6/ but using the runWorkstation.sh file
	note you should specify the maximum number of jobs you want to run concurrently.
	to max a CPU you should choose your processors number of vCores (threads) divided by 4

8/	good luck and have fun!