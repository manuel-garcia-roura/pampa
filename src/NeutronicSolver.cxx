#include "NeutronicSolver.hxx"

/* Solve the eigensystem to get the solution: */
int NeutronicSolver::solve(int n, double dt) {
   
   /* Check the time step: */
   PAMPA_CHECK(n > 0 || dt > 0.0, 1, "neutronic solvers only implemented for steady-state cases");
   
   /* Solve the eigensystem: */
   PAMPA_CALL(petsc::solve(eps), "unable to solve the eigensystem");
   
   /* Get the solution after solving the eigensystem: */
   PAMPA_CALL(getSolution(), "unable to get the solution");
   
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

/* Print the solution summary to standard output: */
int NeutronicSolver::printLog() const {
   
   /* Print out the multiplication factor: */
   if (mpi::rank == 0)
      std::cout << "keff = " << keff << std::endl;
   
   return 0;
   
}
