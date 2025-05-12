#!/bin/bash
#
# Compile the tests for
# the finite element schurr
# complement preconditioner
#
# Author : Sohail Rathore
# Date   : 15/02/2025
#

clear
export PARAFEM_DIR="/home/sar-local/X_Software/parafem/"
cd src
sh ./compile.bash
cd ..
cd ref_elements
sh ./compile.bash
cd ..
rm -r main_SolidsUPAS

#
# Compile the C++ files
#
mpicxx -fopenmp -std=c++11 -O3 -m64 -c main_SolidsUPAS.cpp

#
# Compile the fortran Interface(s)
#
mpif90 -fopenmp -c InterfaceF.f90 -fallow-argument-mismatch   \
                                  -I src/                     \
                                  -I ref_elements/            \
                                  -I $PARAFEM_DIR/include/mpi

#
# Make the executables
#
mpif90 -fopenmp -o main_SolidsUPAS main_SolidsUPAS.o -fallow-argument-mismatch  \
                                            -I src/                    \
                                            -I $PARAFEM_DIR/include/mpi \
                                            InterfaceF.o               \
                                            src/Parallel_IO.o                       \
                                            src/Parallel_supplementary_Maths.o      \
                                            src/Parallel_ELEMENT_Colouring.o        \
                                            src/Parallel_FEA_LinearSolvers.o        \
                                            src/Parallel_BoundaryConditions.o       \
                                            ref_elements/OrthoHeat_ALE.o            \
                                            ref_elements/Solids_Traction.o          \
                                            ref_elements/Solids_UPAS.o              \
                                            $PARAFEM_DIR/lib/libParaFEM_mpi.5.0.3.a \
                                            -lgfortran -lstdc++
rm -r *.o
