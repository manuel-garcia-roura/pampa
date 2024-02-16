#!/bin/bash

set -e

np=3

{

echo "slabs/reflected..."
cd slabs/reflected
../../run.sh $np input.pmp
rm -f -R -- */
rm -f output.vtk
cd ../..
echo "done."

echo "pwr-iaea-benchmark/cartesian-diffusion..."
cd pwr-iaea-benchmark/cartesian-diffusion
../../run.sh $np input.pmp
../../run.sh $np input_dd.pmp
rm -f -R -- */
rm -f output.vtk
cd ../..
echo "done."

echo "pwr-iaea-benchmark/cartesian-diffusion-3d..."
cd pwr-iaea-benchmark/cartesian-diffusion-3d
../../run.sh $np input.pmp
../../run.sh $np input_dd.pmp
rm -f -R -- */
rm -f output.vtk
cd ../..
echo "done."

echo "pwr-iaea-benchmark/cartesian-sn..."
cd pwr-iaea-benchmark/cartesian-sn
../../run.sh $np input.pmp
rm -f -R -- */
rm -f output.vtk
cd ../..
echo "done."

echo "pwr-iaea-benchmark/unstructured-diffusion..."
cd pwr-iaea-benchmark/unstructured-diffusion
../../run.sh $np input.pmp
../../run.sh $np input_dd.pmp
rm -f -R -- */
rm -f output.vtk
cd ../..
echo "done."

echo "pwr-iaea-benchmark/unstructured-diffusion-3d..."
cd pwr-iaea-benchmark/unstructured-diffusion-3d
../../run.sh $np input.pmp
../../run.sh $np input_dd.pmp
rm -f -R -- */
rm -f output.vtk
cd ../..
echo "done."

echo "pwr-iaea-benchmark/unstructured-sn..."
cd pwr-iaea-benchmark/unstructured-sn
../../run.sh $np input.pmp
rm -f -R -- */
rm -f output.vtk
cd ../..
echo "done."

} > check.txt

diff check.txt check_ref.txt

cmp --silent check.txt check_ref.txt && echo -e "\033[0;32mOK!" || echo -e "\033[0;31mKO!"

rm -f check.txt
