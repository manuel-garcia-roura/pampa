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

/* Output the solution: */
int NeutronicSolver::output(const std::string& filename) {
   
   /* Print out the multiplication factor: */
   if (mpi::rank == 0)
      std::cout << "keff = " << keff << std::endl;
   
   /* Write the solution to a plain-text file in .vtk format: */
   PAMPA_CALL(writeVTK(mpi::get_path(filename)), "unable to output the solution in .vtk format");
   
   /* Write the solution to a binary file in PETSc format: */
   PAMPA_CALL(writePETSc("flux.ptc"), "unable to output the solution in PETSc format");
   
   return 0;
   
}

/* Normalize the scalar flux: */
int NeutronicSolver::normalizeScalarFlux() {
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   const Cells& cells = mesh->getCells();
   
   /* Get the array for the scalar flux: */
   PetscScalar* data_phi;
   PETSC_CALL(VecGetArray(phi, &data_phi));
   
   /* Normalize the scalar flux: */
   /* TODO: normalize correctly with the power. */
   double vol = 0.0;
   for (int i = 0; i < num_cells; i++)
      if (materials(cells.materials(i)).nu_sigma_fission(1) > 0.0)
         vol += cells.volumes(i);
   MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
   double sum = 0.0;
   for (int iphi = 0, i = 0; i < num_cells; i++)
      for (int g = 0; g < num_energy_groups; g++)
         sum += data_phi[iphi++] * materials(cells.materials(i)).nu_sigma_fission(g) * 
            cells.volumes(i);
   MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
   double f = vol / sum;
   for (int iphi = 0, i = 0; i < num_cells; i++)
      for (int g = 0; g < num_energy_groups; g++)
         data_phi[iphi++] *= f;
   
   /* Check for negative fluxes: */
   for (int iphi = 0, i = 0; i < num_cells; i++) {
      for (int g = 0; g < num_energy_groups; g++) {
         PAMPA_CHECK(data_phi[iphi++] < 0.0, 1, "negative values in the scalar-flux solution");
      }
   }
   
   /* Restore the array for the scalar flux: */
   PETSC_CALL(VecRestoreArray(phi, &data_phi));
   
   return 0;
   
}
