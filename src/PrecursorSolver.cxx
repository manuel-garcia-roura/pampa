#include "PrecursorSolver.hxx"

/* Solve the linear system to get the solution: */
int PrecursorSolver::solve(int n, double dt) {
   
   /* Get a random production rate: */
   PAMPA_CALL(petsc::random(P), "unable to initialize the production rate");
   PAMPA_CALL(petsc::normalize(P, 1.0), "unable to normalize the production rate");
   
   /* Get the arrays for the precursor population and the production rate: */
   PetscScalar *data_C, *data_P;
   PETSC_CALL(VecGetArray(C, &data_C));
   PETSC_CALL(VecGetArray(P, &data_P));
   
   /* Get the steady-state or transient precursor population: */
   if (n == 0) {
      for (int i = 0; i < num_cells; i++) {
         const Material& mat = materials(cells.materials(i));
         for (int g = 0; g < num_precursor_groups; g++)
            data_C[index(i, g)] = (mat.beta(g)/mat.lambda(g)) * data_P[i];
      }
   }
   else {
      for (int i = 0; i < num_cells; i++) {
         const Material& mat = materials(cells.materials(i));
         for (int g = 0; g < num_precursor_groups; g++) {
            data_C[index(i, g)] += dt * mat.beta(g) * data_P[i];
            data_C[index(i, g)] /= 1.0 + dt*mat.lambda(g);
         }
      }
   }
   
   /* Restore the arrays for the precursor population and the production rate: */
   PETSC_CALL(VecRestoreArray(P, &data_P));
   PETSC_CALL(VecRestoreArray(C, &data_C));
   
   return 0;
   
}

/* Check the material data: */
int PrecursorSolver::checkMaterials() const {
   
   /* Check the materials: */
   for (int i = 0; i < materials.size(); i++) {
      PAMPA_CHECK(materials(i).num_precursor_groups != num_precursor_groups, 1, 
         "wrong number of delayed-neutron precursor groups");
      PAMPA_CHECK(materials(i).lambda.empty(), 2, "missing precursor decay constants");
      PAMPA_CHECK(materials(i).beta.empty(), 3, "missing precursor fractions");
   }
   
   return 0;
   
}

/* Build the solution and source vectors: */
int PrecursorSolver::build() {
   
   /* Create the precursor-population vector: */
   PAMPA_CALL(petsc::create(C, num_cells*num_precursor_groups, 
      num_cells_global*num_precursor_groups, vectors), 
      "unable to create the precursor-population vector");
   
   /* Create the production-rate vector: */
   PAMPA_CALL(petsc::create(P, num_cells, num_cells_global, vectors), 
      "unable to create the production-rate vector");
   
   return 0;
   
}

/* Print the solution summary to standard output: */
int PrecursorSolver::printLog() const {
   
   /* Print out the minimum and maximum precursor populations: */
   double Cmin, Cmax;
   PETSC_CALL(VecMin(C, NULL, &Cmin));
   PETSC_CALL(VecMax(C, NULL, &Cmax));
   if (mpi::rank == 0)
      std::cout << "Cmin = " << Cmin << ", Cmax = " << Cmax << std::endl;
   
   return 0;
   
}

/* Write the solution to a plain-text file in .vtk format: */
int PrecursorSolver::writeVTK(const std::string& filename) const {
   
   /* Write the precursor population in .vtk format: */
   PAMPA_CALL(vtk::write(filename, "precursor", C, num_cells, num_precursor_groups), 
      "unable to write the precursor population");
   
   return 0;
   
}

/* Write the solution to a binary file in PETSc format: */
int PrecursorSolver::writePETSc() const {
   
   /* Write the precursor population in PETSc format: */
   PAMPA_CALL(petsc::write("precursors.ptc", C), "unable to write the precursor population");
   
   return 0;
   
}
