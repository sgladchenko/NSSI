#!/bin/sh

# The home directory of the project on Sk.Pardus
cd "/home/Oleg.Kharlanov/TurboNLM.NSSI"

# determines if a given option is passed in the command line
function optionSupplied {
	for a in ${BASH_ARGV[*]} ; do
		if [ "$a" == "$1" ] 
		then
			echo "true"
			return
		fi
	done
	echo "false"
}

# ------------- Choose compiler/MPI environment -------------------
echo "Setting up the environment..."

module delete Compilers

# The most recent GCC available (not tested with MPI).
module load Compilers/GCC/8.3.0

# Python to generate plots and initial conditions
module load ScriptLang/python/3.7.9

# --------------- Compilation of the project [+enqueueing] --------

echo "1. Compiling TurboNLM.NSSI project..."
 
#make -f makefile_cpp.intel
make -f makefile_cpp.gcc Main
returncode_buildnssi=$?
 
printf "1. Done!\n\n"

# If the user has asked to generate the initial conditions, do that right now!
if [ "$(optionSupplied -generate)" == "true" ]
then
	echo "2. Generating initial conditions..."
	python ./MakeNoise.py --N_Noise=250 --sigma=0.01
	printf "2. Done!\n\n"
fi

printf "All done!\n\n"

# Check the compiler's return code(s) and, if the caller has asked to enqueue, do it!
if [ "$returncode_buildnssi" == "0" ] && [ "$(optionSupplied -E)" == "true" ]
then
    source "./Enqueue_Sk.sh" $*
fi
