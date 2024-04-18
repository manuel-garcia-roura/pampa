#include "PhysicsSolver.hxx"

/* Initialize: */
int PhysicsSolver::initialize(bool transient) {
   
   /* Check the material data: */
   PAMPA_CALL(checkMaterials(transient), "wrong material data");
   
   /* Build the matrices, vectors and solver contexts: */
   PAMPA_CALL(build(), "unable to build the solver");
   
   return 0;
   
}

/* Output the solution: */
int PhysicsSolver::output(const std::string& filename, int n, bool write_mesh) const {
   
   /* Print the solution summary to standard output: */
   PAMPA_CALL(printLog(n), "unable to print the solution summary to standard output");
   
   /* Write the mesh in .vtk format: */
   if (write_mesh) {
      PAMPA_CALL(mesh->writeVTK(filename), "unable to write the mesh in .vtk format");
   }
   
   /* Write the solution in .vtk format: */
   PAMPA_CALL(writeVTK(filename), "unable to write the solution in .vtk format");
   
   /* Write the solution in PETSc format: */
   PAMPA_CALL(writePETSc(), "unable to write the solution in PETSc format");
   
   return 0;
   
}

/* Finalize: */
int PhysicsSolver::finalize() {
   
   /* Destroy the EPS context: */
   PAMPA_CALL(petsc::destroy(eps), "unable to destroy the EPS context");
   
   /* Destroy the KSP context: */
   PAMPA_CALL(petsc::destroy(ksp), "unable to destroy the KSP context");
   
   /* Destroy the PETSc vectors: */
   for (int i = 0; i < vectors.size(); i++) {
      PAMPA_CALL(petsc::destroy(*vectors(i)), "unable to destroy the PETSc vectors");
   }
   
   /* Destroy the PETSc matrices: */
   for (int i = 0; i < matrices.size(); i++) {
      PAMPA_CALL(petsc::destroy(*matrices(i)), "unable to destroy the PETSc matrices");
   }
   
   return 0;
   
}
