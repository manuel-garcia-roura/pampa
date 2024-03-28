#include "Solver.hxx"

/* Initialize: */
int Solver::initialize() {
   
   /* Check the material data: */
   PAMPA_CALL(checkMaterials(), "wrong material data");
   
   /* Build the matrices, vectors and solver contexts: */
   PAMPA_CALL(build(), "unable to build the solver");
   
   return 0;
   
}

/* Output the solution: */
int Solver::output(const std::string& filename) const {
   
   /* Print the solution summary to standard output: */
   PAMPA_CALL(printLog(), "unable to print the solution summary to standard output");
   
   /* Write the mesh in .vtk format: */
   PAMPA_CALL(mesh->writeVTK(mpi::get_path(filename)), "unable to write the mesh in .vtk format");
   
   /* Write the solution in .vtk format: */
   PAMPA_CALL(writeVTK(mpi::get_path(filename)), "unable to write the solution in .vtk format");
   
   /* Write the solution in PETSc format: */
   PAMPA_CALL(writePETSc(), "unable to write the solution in PETSc format");
   
   return 0;
   
}

/* Finalize: */
int Solver::finalize() {
   
   /* Destroy the EPS context: */
   PETSC_CALL(EPSDestroy(&eps));
   
   /* Destroy the KSP context: */
   PETSC_CALL(KSPDestroy(&ksp));
   
   /* Destroy the PETSc vectors: */
   for (int i = 0; i < vectors.size(); i++) {
      PETSC_CALL(VecDestroy(vectors(i)));
   }
   
   /* Destroy the PETSc matrices: */
   for (int i = 0; i < matrices.size(); i++) {
      PETSC_CALL(MatDestroy(matrices(i)));
   }
   
   return 0;
   
}
