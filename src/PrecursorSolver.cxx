#include "PrecursorSolver.hxx"

/* Solve the linear system to get the solution: */
int PrecursorSolver::solve(int n, double dt) {
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   const Cells& cells = mesh->getCells();
   
   /* Get a random production rate: */
   PAMPA_CALL(petsc::normalize(P, 1.0, true), "unable to normalize the production rate");
   
   /* Get the arrays for the precursor population and the production rate: */
   PetscScalar *data_C, *data_P;
   PETSC_CALL(VecGetArray(C, &data_C));
   PETSC_CALL(VecGetArray(P, &data_P));
   
   /* Calculate the precursor population for each cell i and group g: */
   for (int i = 0; i < num_cells; i++) {
      for (int g = 0; g < num_precursor_groups; g++) {
         
         /* Get the material for cell i: */
         const Material& mat = materials(cells.materials(i));
         
         /* Get the steady-state precursor population or integrate it: */
         if (n == 0) {
            data_C[index(i, g)] = (mat.beta(g)/mat.lambda(g)) * data_P[i];
         }
         else {
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

/* Output the solution: */
int PrecursorSolver::output(const std::string& filename) {
   
   /* Get the number of cells: */
   int num_cells = mesh->getNumCells();
   
   /* Write the mesh: */
   PAMPA_CALL(mesh->writeVTK(filename), "unable to write the mesh");
   
   /* Write the precursor population to a .vtk file: */
   PAMPA_CALL(vtk::write(filename, "precursor", C, num_cells, num_precursor_groups), 
      "unable to write the precursor population");
   
   /* Write the precursor population to a PETSc binary file: */
   PAMPA_CALL(petsc::write("precursors.ptc", C), "unable to output the solution in PETSc format");
   
   return 0;
   
}

/* Check the material data: */
int PrecursorSolver::checkMaterials() {
   
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
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   int num_cells_global = mesh->getNumCellsGlobal();
   
   /* Create the precursor-population vector: */
   PAMPA_CALL(petsc::create(C, num_cells*num_precursor_groups, 
      num_cells_global*num_precursor_groups, vectors), 
      "unable to create the precursor-population vector");
   
   /* Create the production-rate vector: */
   PAMPA_CALL(petsc::create(P, num_cells, num_cells_global, vectors), 
      "unable to create the production-rate vector");
   
   return 0;
   
}
