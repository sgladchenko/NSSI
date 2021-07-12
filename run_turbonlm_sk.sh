#!/bin/bash
#PBS -l nodes=1:ppn=8
# commented PBS -l mem=2GB
#PBS -l walltime=60:00:00
#PBS -e output/pbslogs/tnlm.nssi.$PBS_JOBID.err
#PBS -o output/pbslogs/tnlm.nssi.$PBS_JOBID.out

export OMP_NUM_THREADS=8
executable='./nssi --numthreads=8'

# ------------- Choose libraries/MPI environment -------------------
module delete Compilers

# The most recent GCC available (not tested with MPI).
module load Compilers/GCC/8.3.0

# Python to generate plots and initial conditions
module load ScriptLang/python/3.6i

# ------------------------------------------------------------------

date
echo "Running TruboNLM.NSSI task " $PBS_JOBNAME " at " $(eval hostname)

cd "/home/Oleg.Kharlanov/TurboNLM.NSSI"

# Where to place output subfolders
outputRoot=./output
# Where to search for tasks
taskDirs='./tasks/mugscan_05Jul2021'


for taskDir in ${taskDirs}
do
	taskDir_base=$(basename ${taskDir})
	echo "Doing task ${taskDir_base}..."
	
	outputDir=${outputRoot}/${taskDir_base}
	
	# Create output directory
	if [ ! -d "${outputDir}" ]
	then
		mkdir "${outputDir}"
	fi
	
	# Allow for tasks in subdirectories, create output subdirs for them
	for subDir in ${taskDir}/*/
	do
		taskSubdir_base=$(basename ${subDir})
		echo "Found subtask ${taskSubdir_base}!"
		outputSubdir=${outputDir}/${taskSubdir_base}
		if [ ! -d "${outputSubdir}" ]
		then
			mkdir "${outputSubdir}"
		fi
		
		for paramsFile in ${subDir}*.json
		do
			noiseFile=${taskDir}/Noise.json
			
			# Search if the calculation with the same noise & params and, if one is found, don't do the same calculation again
			doCalc=true
			for sd in ${outputSubdir}/*/
			do
				#echo "Scanning directory ${sd} for completed calculations with params ${paramsFile}..."
				if cmp -s "${noiseFile}" "${sd}Noise.json" ; then 
					if cmp -s "${paramsFile}" "${sd}Parameters.json" ; then
						doCalc=false
						echo "Calculation with parameters ${paramsFile} has already been done in ${sd}, skipping recalculation."
					fi
				fi
			done
			
			if [ "${doCalc}" == "false" ] ; then
				continue
			fi
			
			# Prepare calculation params
			nssiParams="--root=${outputSubdir} --periodN_x=10 --periodN_z=10 --parameters=$paramsFile --noise=${noiseFile}"
			logFile=${outputSubdir}/tnlm.nssi.$(eval date +"%Y%m%d-%H%M%S").out
			echo "Submitting a calculation with parameters:   ${nssiParams}..."
			${executable} ${nssiParams} > ${logFile}
			# Now extract the output folder name from the log and hide the log in it
			outFolder=$(tail --lines=1 ${logFile})
			mv ${logFile} "${outFolder}"
			
			# Postprocess the data (e.g., draw plots, find averages, etc.)
			echo "Postprocessing data in output folder ${outFolder}..."
			echo "python Plots.py --periodN_x=1 --periodN_z=1 --dir=\"${outFolder}\""
			python Plots.py --periodN_x=1 --periodN_z=1 --dir="${outFolder}"
			
			printf "End of task ${taskDir_base}.${taskSubdir_base}!\n\n"
		done
		
	done
		
done

echo "Done!"
