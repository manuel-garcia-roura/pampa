#!/bin/bash

source ../../../setenv.sh

mpirun -n 1 pampa $1 -vtk
