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
rm -r *.o *.mod


## Build the object Files
##
mpif90 -c Parallel_supplementary_Maths.f90    -I $PARAFEM_DIR/include/mpi
mpif90 -c Parallel_IO.f90                     -I $PARAFEM_DIR/include/mpi
mpif90 -c Parallel_ELEMENT_Colouring.f90      -I $PARAFEM_DIR/include/mpi
mpif90 -c Parallel_FEA_LinearSolvers.f90      -I $PARAFEM_DIR/include/mpi
mpif90 -c Parallel_BoundaryConditions.f90     -I $PARAFEM_DIR/include/mpi

mpif90 -c BlockJacobiPrecon.f90
mpif90 -c TensorElement.f90
