#$ -S /bin/sh -cwd
#$  -pe mpi.singlehost 8

test_mscreenCLI_prepare.bat
test_mscreenCLI_screening.bat
test_mscreenCLI_analysis.bat


