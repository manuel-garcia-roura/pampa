#!/bin/bash

source ../../../setenv.sh

mpirun -n 1 pampa input.pmp -ksp_type cgs
