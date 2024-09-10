#include "NeutronicSolver.hxx"

/* Solve the linear system to get the solution: */
int NeutronicSolver::solve(int n, double dt, double t) {
   
   /* Print progress: */
   mpi::print("Run '" + name + "' solver...", true);
   
   /* Build the coefficient matrices and the RHS vector: */
   PAMPA_CALL(buildMatrices(n, dt, t), "unable to build the coefficient matrices");
   
   /* Manage the EPS and KSP contexts: */
   if (n == 0) {
      if (eps != 0) {
         PAMPA_CALL(petsc::destroy(eps), "unable to destroy the EPS context");
      }
      PAMPA_CALL(petsc::create(eps, R, F), "unable to create the EPS context");
   }
   else {
      if (eps != 0) {
         PAMPA_CALL(petsc::destroy(eps), "unable to destroy the EPS context");
         PAMPA_CALL(petsc::destroy(F), "unable to destroy the F coefficient matrix");
      }
      if (ksp != 0) {
         PAMPA_CALL(petsc::destroy(ksp), "unable to destroy the KSP context");
      }
      PAMPA_CALL(petsc::create(ksp, R), "unable to create the KSP context");
   }
   
   /* Solve the linear system and get the solution: */
   PAMPA_CALL(getSolution(n), "unable to solve the linear system and get the solution");
   
   /* Calculate the thermal power and the production rate: */
   PAMPA_CALL(calculatePowerAndProductionRate(), "unable to calculate the thermal power");
   
   /* Print progress: */
   mpi::print("Done.", true);
   
   return 0;
   
}

/* Normalize the scalar flux: */
int NeutronicSolver::normalizeScalarFlux() {
   
   /* Get the arrays with the raw data: */
   PetscScalar *T_data, *phi_data;
   PETSC_CALL(VecGetArray(T, &T_data));
   PETSC_CALL(VecGetArray(phi, &phi_data));
   
   /* Get the current power: */
   double p0 = 0.0;
   for (int iphi = 0, i = 0; i < num_cells; i++) {
      const Material* mat = materials(cells->materials(i));
      for (int g = 0; g < num_energy_groups; g++)
         p0 += phi_data[iphi++] * mat->sigmaKappaFission(g, T_data[i]) * cells->volumes(i);
   }
   MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, &p0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
   
   /* Normalize the scalar flux and check for negative fluxes: */
   double f = power / p0;
   for (int iphi = 0, i = 0; i < num_cells; i++) {
      for (int g = 0; g < num_energy_groups; g++) {
         phi_data[iphi] *= f;
         PAMPA_CHECK(phi_data[iphi] < 0.0, 1, "negative values in the scalar-flux solution");
         iphi++;
      }
   }
   
   /* Restore the arrays with the raw data: */
   PETSC_CALL(VecRestoreArray(T, &T_data));
   PETSC_CALL(VecRestoreArray(phi, &phi_data));
   
   return 0;
   
}

/* Calculate the thermal power and the production rate: */
int NeutronicSolver::calculatePowerAndProductionRate() {
   
   /* Get the arrays with the raw data: */
   PetscScalar *T_data, *S_data, *phi_data, *q_data, *P_data;
   PETSC_CALL(VecGetArray(T, &T_data));
   PETSC_CALL(VecGetArray(S, &S_data));
   PETSC_CALL(VecGetArray(phi, &phi_data));
   PETSC_CALL(VecGetArray(q, &q_data));
   PETSC_CALL(VecGetArray(P, &P_data));
   
   /* Get the current power: */
   power = 0.0;
   for (int iphi = 0, i = 0; i < num_cells; i++) {
      const Material* mat = materials(cells->materials(i));
      q_data[i] = 0.0;
      P_data[i] = 0.0;
      for (int g = 0; g < num_energy_groups; g++) {
         q_data[i] += phi_data[iphi] * mat->sigmaKappaFission(g, T_data[i]) * cells->volumes(i);
         P_data[i] += phi_data[iphi] * mat->sigmaNuFission(g, T_data[i]) * cells->volumes(i) / keff;
         iphi++;
      }
      S_data[i] = mat->beta() * P_data[i];
      power += q_data[i];
   }
   MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, &power, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
   
   /* Restore the arrays with the raw data: */
   PETSC_CALL(VecRestoreArray(T, &T_data));
   PETSC_CALL(VecRestoreArray(S, &S_data));
   PETSC_CALL(VecRestoreArray(phi, &phi_data));
   PETSC_CALL(VecRestoreArray(q, &q_data));
   PETSC_CALL(VecRestoreArray(P, &P_data));
   
   return 0;
   
}

/* Print the solution summary to standard output: */
int NeutronicSolver::printLog(int n) const {
   
   /* Print out the multiplication factor and the total thermal power: */
   if (n == 0) mpi::print("keff", keff);
   mpi::print("P", power);
   
   return 0;
   
}
