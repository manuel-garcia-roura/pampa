#include "Solver.hxx"

/* Initialize: */
int Solver::initialize() {
   
   /* Check the material data: */
   PAMPA_CALL(checkMaterials(), "wrong material data");
   
   /* Build the matrices, vectors and solver contexts: */
   PAMPA_CALL(build(), "unable to build the solver");
   
   return 0;
   
}

/* Finalize: */
int Solver::finalize() {
   
   /* Destroy the EPS context: */
   PETSC_CALL(EPSDestroy(&eps));
   
   /* Destroy the KSP context: */
   PETSC_CALL(KSPDestroy(&ksp));
   
   /* Destroy the PETSc vectors: */
   for (int i = 0; i < vectors.size(); i++)
      PETSC_CALL(VecDestroy(vectors(i)));
   
   /* Destroy the PETSc matrices: */
   for (int i = 0; i < matrices.size(); i++)
      PETSC_CALL(MatDestroy(matrices(i)));
   
   return 0;
   
}
