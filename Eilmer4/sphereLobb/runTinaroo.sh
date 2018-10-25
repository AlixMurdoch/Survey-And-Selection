#!/bin/bash

for dir in */; 
do
	echo "$dir"
	cd "$dir"
	qsub run-tinaroo.qsub
	cd ..

done

