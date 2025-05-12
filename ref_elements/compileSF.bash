#!bin/bash
#################################################
#
# Author: Sohail Rathore
#
#################################################

##Load the MPI module file with compiler
#module load apps/gcc/openmpi-4.0.1

# Clear the object and module files
#
rm -r *".o" *".mod"              

## Build the object Files
##
mpif90 -fopenmp -c OrthoHeat_ALE.f90  -I $PARAFEM_DIR/include/mpi \
                                      -I ../src

mpif90 -fopenmp -c Solids_UPAS.f90 -I $PARAFEM_DIR/include/mpi \
                                   -I ../src

mpif90 -fopenmp -c Solids_Traction.f90 -I $PARAFEM_DIR/include/mpi \
                                       -I ../src
