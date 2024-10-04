#!/bin/bash

source ../../../setenv.sh

mpirun -n 1 pampa input.pmp -verbose -vtk
