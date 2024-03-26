#include "PrecursorSolver.hxx"

/* Initialize: */
int PrecursorSolver::initialize(int argc, char* argv[]) {
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   int num_cells_global = mesh->getNumCellsGlobal();
   
   /* Check the material data: */
   PAMPA_CALL(checkMaterials(), "wrong material data");
   
   /* Create the precursor-population vector: */
   PAMPA_CALL(petsc::create_vector(C, num_cells*num_precursor_groups, 
      num_cells_global*num_precursor_groups), "unable to create the precursor-population vector");
   
   /* Create the production-rate vector: */
   PAMPA_CALL(petsc::create_vector(P, num_cells, num_cells_global), 
      "unable to create the production-rate vector");
   
   return 0;
   
}

/* Solve the linear system to get the solution: */
int PrecursorSolver::solve(int n, double dt) {
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   const Cells& cells = mesh->getCells();
   
   /* Get a random production rate: */
   PAMPA_CALL(petsc::normalize_vector(P, 1.0, true), "unable to normalize the production rate");
   
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
            data_C[index(i, g, num_precursor_groups)] = (mat.beta(g)/mat.lambda(g)) * data_P[i];
         }
         else {
            data_C[index(i, g, num_precursor_groups)] += dt * mat.beta(g) * data_P[i];
            data_C[index(i, g, num_precursor_groups)] /= 1.0 + dt*mat.lambda(g);
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
   
   /* Write the solution to a plain-text file in .vtk format: */
   PAMPA_CALL(writeVTK(mpi::get_path(filename)), "unable to output the solution in .vtk format");
   
   /* Write the solution to a binary file in PETSc format: */
   PetscBool write, flag;
   PETSC_CALL(PetscOptionsGetBool(NULL, NULL, "-petsc_write_solution", &write, &flag));
   if (flag && write) {
      PAMPA_CALL(petsc::write("precursors.ptc", C), "unable to write the solution");
   }
   
   return 0;
   
}

/* Finalize: */
int PrecursorSolver::finalize() {
   
   /* Destroy the production-rate vector: */
   PETSC_CALL(VecDestroy(&P));
   
   /* Destroy the precursor-population vector: */
   PETSC_CALL(VecDestroy(&C));
   
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

/* Write the solution to a plain-text file in .vtk format: */
int PrecursorSolver::writeVTK(const std::string& filename) const {
   
   /* Get the number of cells: */
   int num_cells = mesh->getNumCells();
   
   /* Get the array for the precursor population: */
   PetscScalar* data_C;
   PETSC_CALL(VecGetArray(C, &data_C));
   
   /* Write the mesh: */
   PAMPA_CALL(mesh->writeVTK(filename), "unable to write the mesh");
   
   /* Open the output file: */
   std::ofstream file(filename, std::ios_base::app);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Write the precursor population: */
   for (int g = 0; g < num_precursor_groups; g++) {
      file << "SCALARS precursor_" << (g+1) << " double 1" << std::endl;
      file << "LOOKUP_TABLE default" << std::endl;
      for (int i = 0; i < num_cells; i++)
         file << data_C[index(i, g, num_precursor_groups)] << std::endl;
      file << std::endl;
   }
   
   /* Restore the array for the precursor population: */
   PETSC_CALL(VecRestoreArray(C, &data_C));
   
   return 0;
   
}
