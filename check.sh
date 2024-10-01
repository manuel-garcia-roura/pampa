#!/bin/bash

set -e

source setenv.sh

np=2

cd test

{

printf "\nRun slabs/reflected-diffusion...\n"
cd slabs/reflected-diffusion
mpirun -n $np pampa input.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
printf "Done.\n"

printf "\nRun slabs/reflected-s2...\n"
cd slabs/reflected-s2
mpirun -n $np pampa input.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
printf "Done.\n"

printf "\nRun slabs/reflected-s4...\n"
cd slabs/reflected-s4
mpirun -n $np pampa input.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
printf "Done.\n"

printf "\nRun pwr-iaea-benchmark/cartesian-diffusion...\n"
cd pwr-iaea-benchmark/cartesian-diffusion
mpirun -n $np pampa input.pmp
printf "Rerun...\n"
mpirun -n $np pampa input_dd.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
printf "Done.\n"

printf "\nRun pwr-iaea-benchmark/cartesian-diffusion-3d...\n"
cd pwr-iaea-benchmark/cartesian-diffusion-3d
mpirun -n $np pampa input.pmp
printf "Rerun...\n"
mpirun -n $np pampa input_dd.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
printf "Done.\n"

printf "\nRun pwr-iaea-benchmark/cartesian-sn...\n"
cd pwr-iaea-benchmark/cartesian-sn
mpirun -n $np pampa input.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
printf "Done.\n"

printf "\nRun pwr-iaea-benchmark/unstructured-diffusion...\n"
cd pwr-iaea-benchmark/unstructured-diffusion
mpirun -n $np pampa input.pmp
printf "Rerun...\n"
mpirun -n $np pampa input_dd.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
printf "Done.\n"

printf "\nRun pwr-iaea-benchmark/unstructured-diffusion-3d...\n"
cd pwr-iaea-benchmark/unstructured-diffusion-3d
mpirun -n $np pampa input.pmp
printf "Rerun...\n"
mpirun -n $np pampa input_dd.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
printf "Done.\n"

printf "\nRun pwr-iaea-benchmark/unstructured-sn...\n"
cd pwr-iaea-benchmark/unstructured-sn
mpirun -n $np pampa input.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
printf "Done.\n"

printf "\nRun pwr-iaea-benchmark/cartesian-conduction-3d...\n"
cd pwr-iaea-benchmark/cartesian-conduction-3d
mpirun -n $np pampa input.pmp
printf "Rerun...\n"
mpirun -n $np pampa input_dd.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
printf "Done.\n"

printf "\nRun pwr-iaea-benchmark/cartesian-precursors-3d...\n"
cd pwr-iaea-benchmark/cartesian-precursors-3d
mpirun -n $np pampa input.pmp
printf "Rerun...\n"
mpirun -n $np pampa input_dd.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
printf "Done.\n"

} > check.txt

diff check.txt check_ref.txt

cmp --silent check.txt check_ref.txt && echo -e "\033[0;32mOK!" || echo -e "\033[0;31mKO!"

rm -f check.txt

cd ..
