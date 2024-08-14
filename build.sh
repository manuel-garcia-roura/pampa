#!/bin/bash

set -e

export PETSC_DIR=/lib/petscdir/petsc3.12
export SLEPC_DIR=/lib/slepcdir/slepc3.12
export PETSC_ARCH=x86_64-linux-gnu-real

rm -rf lib bin
mkdir lib bin

cd src
make clean
make
cd ..

cd src/cxx
make clean
make
cd ../..

cd src/c
make clean
make
cd ../..

cd src/f90
make clean
make
cd ../..

cd bin
ln -sf pampa-f90 pampa
cd ..
