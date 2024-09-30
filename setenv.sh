#!/bin/bash

export PETSC_DIR=/lib/petscdir/petsc3.12
export SLEPC_DIR=/lib/slepcdir/slepc3.12
export PETSC_ARCH=x86_64-linux-gnu-real

export EIGEN_DIR=/usr/include/eigen3
export PAMPA_WITH_EIGEN=true

export METIS_DIR=/usr/lib/x86_64-linux-gnu
export PAMPA_WITH_METIS=false

export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$SLEPC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH

export PAMPA_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
export PATH=$PAMPA_DIR/bin:$PATH
export LD_LIBRARY_PATH=$PAMPA_DIR/lib:$LD_LIBRARY_PATH
