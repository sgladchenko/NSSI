#!/bin/bash

# The home directory of the project on Sk.Pardus
cd "/home/Oleg.Kharlanov/TurboNLM.NSSI/output"

folders=*
if [ "$*" != "" ] ; then
	folders=$*
fi

echo "Scanning folder(s) ${folders} in $(pwd)..."

for taskDir in $folders
do
	if [ ${taskDir} == "pbslogs" ] || [ ! -d "${taskDir}" ]
	then 
		continue
	fi
	
	echo "Filtering output folder ${taskDir}..."
	
	if [ ! -d "../results/${taskDir}" ]
	then
		mkdir "../results/${taskDir}"
	fi
	
	# Copy plots and output logs
	cd ${taskDir}
	cp --parents ./*/*/plots/*.*    ../../results/${taskDir}
	cp --parents ./*/*/averages/*.* ../../results/${taskDir}
	cp --parents ./*/*/*.out        ../../results/${taskDir}
	cp --parents ./*/*/*.json       ../../results/${taskDir}
	cd ..
	
	pushd "../results/${taskDir}"
	zipName="../$(basename ${taskDir}).zip"
	echo "Zipping the folder to $zipName..."
	zip -r "$zipName" *
	popd
done


