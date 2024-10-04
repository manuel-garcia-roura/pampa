#!/bin/bash

source ../../../setenv.sh

ksp_types=("cg" "groppcg" "pipecg" "pipecgrr" "pipelcg" "pipeprcg" "nash" "stcg" "gltr" "fcg" \
           "pipefcg" "gmres" "fgmres" "lgmres" "dgmres" "pgmres" "tcqmr" "bcgs" "ibcgs" "qmrcgs" \
           "fbcgs" "fbcgsr" "bcgsl" "pipebcgs" "tfqmr" "cr" "pipecr" "bicg" "minres" "symmlq" \
           "lcd" "gcr" "pipegcr")

pc_types=("none" "jacobi" "sor" "bjacobi" "mg" "asm" "gasm" "gamg" "telescope" "hmg")

for ksp in "${ksp_types[@]}"
do
   for pc in "${pc_types[@]}"
   do
      echo -n "$ksp $pc "; mpirun -n 4 pampa input.pmp -ksp_type "$ksp" -pc_type "$pc" | grep "Solution time:" | awk -F':' '{printf $2}' | tr -d ' ' || true
      echo ""
   done
done
