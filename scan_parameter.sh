#!/bin/bash

MODES=( 1 2 3 4 5 6 7 8 )		 

FOLDBASE=KHI

for VV in ${MODES[@]}
do

	FOLD=${FOLDBASE}_${VV}
	
	# create directory 
	if [ ! -e ${FOLD} ]
	then
		mkdir ${FOLD}
	fi
	
	# copy BOUT.inp files		
	cp input.inp ${FOLD}
	cp a.out ${FOLD}
	cp jobRMHD.sh ${FOLD}
		
	# replace variable
	cd ${FOLD}
	sed -i s/mode/${VV}/ input.inp
#	cd ..

	# set output dir. & submit job
#	sed -i s/OUTPUT_DIR/${FOLD}/ jobRMHD.sh
	qsub jobRMHD.sh
#	sed -i s/${FOLD}/OUTPUT_DIR/ jobRMHD.sh
	cd ../
	# echo
	echo Running job with MODE Number = ${VV} to folder ${FOLD}

done

