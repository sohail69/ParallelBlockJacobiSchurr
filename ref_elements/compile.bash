#!bin/bash
#################################################
#
# Author: Sohail Rathore
#
#################################################

##Load the MPI module file with compiler
#module load apps/gcc/openmpi-4.0.1

clear
clear

# Clear the object and module files
#
rm -r *".o" *".mod"              

## Build the object Files
##
export PARAFEM_DIR="/home/sar-local/X_Software/parafem/"
#mpif90 -c Solids_UAS.f90  -fallow-argument-mismatch \
#                          -I $PARAFEM_DIR/include/mpi \
#                          -I ../src

mpif90 -c Solids_UPAS.f90 -fallow-argument-mismatch \
                          -I $PARAFEM_DIR/include/mpi \
                          -I ../src
