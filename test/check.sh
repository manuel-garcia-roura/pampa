#!/bin/bash

set -e

np=2

{

echo "Run slabs/reflected-diffusion..."
cd slabs/reflected-diffusion
../../run.sh slepc $np input.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run slabs/reflected-s2..."
cd slabs/reflected-s2
../../run.sh slepc $np input.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run slabs/reflected-s4..."
cd slabs/reflected-s4
../../run.sh slepc $np input.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run pwr-iaea-benchmark/cartesian-diffusion..."
cd pwr-iaea-benchmark/cartesian-diffusion
../../run.sh slepc $np input.pmp
../../run.sh slepc $np input_dd.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run pwr-iaea-benchmark/cartesian-diffusion-3d..."
cd pwr-iaea-benchmark/cartesian-diffusion-3d
../../run.sh slepc $np input.pmp
../../run.sh slepc $np input_dd.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run pwr-iaea-benchmark/cartesian-sn..."
cd pwr-iaea-benchmark/cartesian-sn
../../run.sh slepc $np input.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run pwr-iaea-benchmark/unstructured-diffusion..."
cd pwr-iaea-benchmark/unstructured-diffusion
../../run.sh slepc $np input.pmp
../../run.sh slepc $np input_dd.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run pwr-iaea-benchmark/unstructured-diffusion-3d..."
cd pwr-iaea-benchmark/unstructured-diffusion-3d
../../run.sh slepc $np input.pmp
../../run.sh slepc $np input_dd.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run pwr-iaea-benchmark/unstructured-sn..."
cd pwr-iaea-benchmark/unstructured-sn
../../run.sh slepc $np input.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run pwr-iaea-benchmark/cartesian-conduction-3d..."
cd pwr-iaea-benchmark/cartesian-conduction-3d
../../run.sh petsc $np input.pmp
../../run.sh petsc $np input_dd.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

echo "Run pwr-iaea-benchmark/cartesian-precursors-3d..."
cd pwr-iaea-benchmark/cartesian-precursors-3d
../../run.sh petsc $np input.pmp
../../run.sh petsc $np input_dd.pmp
rm -f -R -- */
rm -f output_*.vtk
cd ../..
echo "Done."

} > check.txt

diff check.txt check_ref.txt

cmp --silent check.txt check_ref.txt && echo -e "\033[0;32mOK!" || echo -e "\033[0;31mKO!"

rm -f check.txt
