#!/bin/bash

PAMPA_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

source $PAMPA_DIR/setenv.sh

mpirun -n $1 pampa $2 -vtk -ksp_atol 1.0e-9 -ksp_rtol 1.0e-9
