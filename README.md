# NSSI
Computer simulations of the collective neutrino oscillations in the presence of scalar and pseudoscalar coupling coming from extensions of the Standard Model

## This branch: SkPardus
A working version for IBM Pardus cluster at Skoltech. Differs from master branch in the wrappers of std::filesystem functions, to handle lack of C++17 support, and in the shell scripts that perform compilation, enqueueing, postprocessing, and yielding of the calculation results. 
Wrappers and shell scripts © okharl, all other project units: © sgladchenko

### Usage note
Due to a definitive change of directory structure in the master branch _after_ the SkPardus branch had been physically created, the latter now retains the original folders. 
Please **do not change** the directory structure in place, as this branch is a snapshot of a working deployed code. Consider creating a new branch instead.