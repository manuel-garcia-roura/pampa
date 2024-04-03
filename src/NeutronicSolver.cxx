#include "NeutronicSolver.hxx"

/* Solve the linear system to get the solution: */
int NeutronicSolver::solve(int n, double dt) {
   
   /* Get a random neutron source: */
   if (n > 0) {
      PAMPA_CALL(petsc::random(q), "unable to initialize the neutron source");
      PAMPA_CALL(petsc::normalize(q, 0.0), "unable to normalize the neutron source");
   }
   
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
   
   /* Calculate the thermal power: */
   if (n > 0) {
      PAMPA_CALL(calculatePower(), "unable to calculate the thermal power");
   }
   
   return 0;
   
}

/* Normalize the scalar flux: */
int NeutronicSolver::normalizeScalarFlux() {
   
   /* Get the array for the scalar flux: */
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
   
   /* Restore the array for the scalar flux: */
   PETSC_CALL(VecRestoreArray(phi, &phi_data));
   
   return 0;
   
}

/* Calculate the thermal power: */
int NeutronicSolver::calculatePower() {
   
   /* Get the array for the scalar flux: */
   PetscScalar* phi_data;
   PETSC_CALL(VecGetArray(phi, &phi_data));
   
   /* Get the current power: */
   power(0) = 0.0;
   for (int iphi = 0, i = 0; i < num_cells; i++) {
      const Material& mat = materials(cells.materials(i));
      for (int g = 0; g < num_energy_groups; g++)
         power(0) += phi_data[iphi++] * mat.e_sigma_fission(g) * cells.volumes(i);
   }
   MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, &(power(0)), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
   
   /* Restore the array for the scalar flux: */
   PETSC_CALL(VecRestoreArray(phi, &phi_data));
   
   return 0;
   
}

/* Print the solution summary to standard output: */
int NeutronicSolver::printLog(int n) const {
   
   /* Print out the multiplication factor and the total thermal power: */
   if (mpi::rank == 0) {
      if (n == 0)
         std::cout << "keff = " << keff << std::endl;
      else
         std::cout << "P = " << power(0) << std::endl;
   }
   
   return 0;
   
}
