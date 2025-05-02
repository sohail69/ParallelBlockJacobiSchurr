#!/bin/bash
#
# Compile the modules necessary
# for the finite element schurr
# complement preconditioner
#
# Author : Sohail Rathore
# Date   : 15/02/2025
#

# Clear the object and module files
#
clear
rm -r *.o *.mod


## Build the object Files
##
export PARAFEM_DIR="/home/sar-local/X_Software/parafem/"
mpif90 -c Parallel_supplementary_Maths.f90    -fallow-argument-mismatch  \
                                              -I $PARAFEM_DIR/include/mpi
mpif90 -c Parallel_IO.f90                     -fallow-argument-mismatch  \
                                              -I $PARAFEM_DIR/include/mpi
mpif90 -c Parallel_ELEMENT_Colouring.f90      -fallow-argument-mismatch  \
                                              -I $PARAFEM_DIR/include/mpi
mpif90 -c Parallel_FEA_LinearSolvers.f90      -fallow-argument-mismatch  \
                                              -I $PARAFEM_DIR/include/mpi
mpif90 -c Parallel_BoundaryConditions.f90     -fallow-argument-mismatch  \
                                              -I $PARAFEM_DIR/include/mpi

mpif90 -c BlockJacobiPrecon.f90 -fallow-argument-mismatch
mpif90 -c TensorElement.f90
