#!/bin/bash

set -e

export PETSC_DIR=/lib/petscdir/petsc3.12
export SLEPC_DIR=/lib/slepcdir/slepc3.12
export PETSC_ARCH=x86_64-linux-gnu-real

cd src

make clean
make

cd cxx
make clean
make
cd ..

cd ..
