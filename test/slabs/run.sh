#!/bin/bash

solver=krylovschur

if [[ "$solver" == "krylovschur" ]]; then
   mpirun -n 1 ../../../bin/pampa \
      -eps_nev 1 -eps_conv_abs -eps_tol 1e-9 \
      -eps_type krylovschur \
      -st_type sinvert -st_pc_factor_shift_type NONZERO
else
   echo "Wrong solver!"
fi
