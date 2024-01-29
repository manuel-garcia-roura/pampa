#!/bin/bash

solver=krylovschur

if [[ "$solver" == "krylovschur" ]]; then
   mpirun -n 2 ../../../bin/pampa \
      -eps_nev 1 -eps_conv_abs -eps_tol 1e-9 \
      -eps_type krylovschur \
      -st_type sinvert -st_pc_factor_shift_type NONZERO
elif [[ "$solver" == "power" ]]; then
   mpirun -n 2 ../../../bin/pampa \
      -eps_nev 1 -eps_conv_abs -eps_tol 1e-9 \
      -eps_type power \
      -st_type sinvert -st_pc_factor_shift_type NONZERO
elif [[ "$solver" == "subspace" ]]; then
   mpirun -n 2 ../../../bin/pampa \
      -eps_nev 1 -eps_conv_abs -eps_tol 1e-9 \
      -eps_type subspace \
      -st_type sinvert -st_pc_factor_shift_type NONZERO
elif [[ "$solver" == "gd" ]]; then
   mpirun -n 2 ../../../bin/pampa \
      -eps_nev 1 -eps_conv_abs -eps_tol 1e-9 \
      -eps_type gd \
      -st_type precond -st_pc_factor_shift_type NONZERO
elif [[ "$solver" == "jd" ]]; then
   mpirun -n 2 ../../../bin/pampa \
      -eps_nev 1 -eps_conv_abs -eps_tol 1e-9 \
      -eps_type jd \
      -st_type precond -st_pc_type asm -st_ksp_type gmres
else
   echo "Wrong solver!"
fi
