#!/bin/bash

source ../../../setenv.sh

mpirun -n 1 pampa input.pmp -verbose -vtk -ksp_type cr -pc_type bjacobi -ksp_min_it 5 -ksp_atol 1.0e-6 -ksp_rtol 1.0e-6
