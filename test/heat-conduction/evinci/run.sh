#!/bin/bash

source ../../../setenv.sh

mpirun -n 1 pampa $1 -verbose -vtk -ksp_atol 1.0e-9 -ksp_rtol 1.0e-9
