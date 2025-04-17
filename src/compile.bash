#!/bin/bash
#
# Compile the modules necessary
# for the finite element schurr
# complement preconditioner
#
# Author : Sohail Rathore
# Date   : 15/02/2025
#

clear
rm -r *.mod *.o

mpif90 -c BlockJacobiPrecon.f90 -fallow-argument-mismatch
