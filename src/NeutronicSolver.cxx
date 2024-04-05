#include "NeutronicSolver.hxx"

/* Solve the linear system to get the solution: */
int NeutronicSolver::solve(int n, double dt) {
   
   /* Build the coefficient matrices: */
   PAMPA_CALL(buildMatrices(n, dt), "unable to build the coefficient matrices");
   
   /* Create the EPS context: */
   if (n == 0 && eps == 0) {
      PAMPA_CALL(petsc::create(eps, R, F), "unable to create the EPS context");
   }
   
   /* Create the KSP context: */
   if (n > 0 && ksp == 0) {
      PAMPA_CALL(petsc::create(ksp, R), "unable to create the KSP context");
   }
   
   /* Solve the linear system and get the solution: */
   PAMPA_CALL(getSolution(n), "unable to solve the linear system and get the solution");
   
   /* Calculate the thermal power and the production rate: */
   PAMPA_CALL(calculatePowerAndProductionRate(), "unable to calculate the thermal power");
   
   return 0;
   
}

/* Normalize the scalar flux: */
int NeutronicSolver::normalizeScalarFlux() {
   
   /* Get the array with the raw data: */
   PetscScalar* phi_data;
   PETSC_CALL(VecGetArray(phi, &phi_data));
   
   /* Get the current power: */
   double p0 = 0.0;
   for (int iphi = 0, i = 0; i < num_cells; i++) {
      const Material& mat = materials(cells.materials(i));
      for (int g = 0; g < num_energy_groups; g++)
         p0 += phi_data[iphi++] * mat.e_sigma_fission(g) * cells.volumes(i);
   }
   MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, &p0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
   
   /* Normalize the scalar flux and check for negative fluxes: */
   double f = power(0) / p0;
   for (int iphi = 0, i = 0; i < num_cells; i++) {
      for (int g = 0; g < num_energy_groups; g++) {
         phi_data[iphi] *= f;
         PAMPA_CHECK(phi_data[iphi] < 0.0, 1, "negative values in the scalar-flux solution");
         iphi++;
      }
   }
   
   /* Restore the array with the raw data: */
   PETSC_CALL(VecRestoreArray(phi, &phi_data));
   
   return 0;
   
}

/* Calculate the thermal power and the production rate: */
int NeutronicSolver::calculatePowerAndProductionRate() {
   
   /* Get the arrays with the raw data: */
   PetscScalar *phi_data, *q_data, *P_data, *S_data;;
   PETSC_CALL(VecGetArray(phi, &phi_data));
   PETSC_CALL(VecGetArray(q, &q_data));
   PETSC_CALL(VecGetArray(P, &P_data));
   PETSC_CALL(VecGetArray(S, &S_data));
   
   /* Get the current power: */
   power(0) = 0.0;
   for (int iphi = 0, i = 0; i < num_cells; i++) {
      const Material& mat = materials(cells.materials(i));
      q_data[i] = 0.0;
      P_data[i] = 0.0;
      for (int g = 0; g < num_energy_groups; g++) {
         q_data[i] += phi_data[iphi] * mat.e_sigma_fission(g) * cells.volumes(i);
         P_data[i] += phi_data[iphi] * mat.nu_sigma_fission(g) * cells.volumes(i) / keff;
         iphi++;
      }
      S_data[i] = mat.beta_total * P_data[i];
      power(0) += q_data[i];
   }
   MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, &(power(0)), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
   
   /* Restore the arrays with the raw data: */
   PETSC_CALL(VecRestoreArray(phi, &phi_data));
   PETSC_CALL(VecRestoreArray(q, &q_data));
   PETSC_CALL(VecRestoreArray(P, &P_data));
   PETSC_CALL(VecRestoreArray(S, &S_data));
   
   return 0;
   
}

/* Print the solution summary to standard output: */
int NeutronicSolver::printLog(int n) const {
   
   /* Print out the multiplication factor and the total thermal power: */
   if (mpi::rank == 0) {
      if (n == 0)
         std::cout << "keff = " << keff << std::endl;
      std::cout << "n = " << n << ": P = " << power(0) << std::endl;
   }
   
   return 0;
   
}
