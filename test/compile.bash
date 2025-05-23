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
cd src
sh ./compile.bash
cd ..
rm -r BJSTest0 BJSTest1

mpif90 -o BJSTest0    test/BlockJacobiPreconTest0.f90 -fallow-argument-mismatch \
                                                      -I src/                   \
                                                       src/BlockJacobiPrecon.o

mpif90 -o BJSTest1    test/BlockJacobiPreconTest1.f90 -fallow-argument-mismatch \
                                                      -I src/                   \
                                                       src/BlockJacobiPrecon.o

mpif90 -o HierElmTest test/HierElmTest.f90            -I src/                  \
                                                      ref_elements/Solids_UAS.o    \
                                                      src/TensorElement.o -llapack -lblas
                                                      
mpif90 -o main_SolidsUPAS main_SolidsUPAS.o -fallow-argument-mismatch \
                                            -I src/                   \
                                            src/BlockJacobiPrecon.o
