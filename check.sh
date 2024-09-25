#!/bin/bash

set -e

source setenv.sh

np=2
petscopts="-ksp_type cg"
slepcopts="-eps_type krylovschur -st_type sinvert -st_pc_factor_shift_type NONZERO"

cd test

{

echo "Run slabs/reflected-diffusion..."
cd slabs/reflected-diffusion
mpirun -n $np pampa input.pmp $slepcopts
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run slabs/reflected-s2..."
cd slabs/reflected-s2
mpirun -n $np pampa input.pmp $slepcopts
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run slabs/reflected-s4..."
cd slabs/reflected-s4
mpirun -n $np pampa input.pmp $slepcopts
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run pwr-iaea-benchmark/cartesian-diffusion..."
cd pwr-iaea-benchmark/cartesian-diffusion
mpirun -n $np pampa input.pmp $slepcopts
mpirun -n $np pampa input_dd.pmp $slepcopts
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run pwr-iaea-benchmark/cartesian-diffusion-3d..."
cd pwr-iaea-benchmark/cartesian-diffusion-3d
mpirun -n $np pampa input.pmp $slepcopts
mpirun -n $np pampa input_dd.pmp $slepcopts
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run pwr-iaea-benchmark/cartesian-sn..."
cd pwr-iaea-benchmark/cartesian-sn
mpirun -n $np pampa input.pmp $slepcopts
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run pwr-iaea-benchmark/unstructured-diffusion..."
cd pwr-iaea-benchmark/unstructured-diffusion
mpirun -n $np pampa input.pmp $slepcopts
mpirun -n $np pampa input_dd.pmp $slepcopts
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run pwr-iaea-benchmark/unstructured-diffusion-3d..."
cd pwr-iaea-benchmark/unstructured-diffusion-3d
mpirun -n $np pampa input.pmp $slepcopts
mpirun -n $np pampa input_dd.pmp $slepcopts
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run pwr-iaea-benchmark/unstructured-sn..."
cd pwr-iaea-benchmark/unstructured-sn
mpirun -n $np pampa input.pmp $slepcopts
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run pwr-iaea-benchmark/cartesian-conduction-3d..."
cd pwr-iaea-benchmark/cartesian-conduction-3d
mpirun -n $np pampa input.pmp $petscopts
mpirun -n $np pampa input_dd.pmp $petscopts
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run pwr-iaea-benchmark/cartesian-precursors-3d..."
cd pwr-iaea-benchmark/cartesian-precursors-3d
mpirun -n $np pampa input.pmp
mpirun -n $np pampa input_dd.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

} > check.txt

diff check.txt check_ref.txt

cmp --silent check.txt check_ref.txt && echo -e "\033[0;32mOK!" || echo -e "\033[0;31mKO!"

rm -f check.txt

cd ..
