#include "petsc.hxx"

/* Create, preallocate and set up a matrix: */
int petsc::create_matrix(Mat& M, int n, int m) {
   
   /* Create the matrix: */
   PETSC_CALL(MatCreate(MPI_COMM_WORLD, &M));
   
   /* Preallocate the matrix: */
   PETSC_CALL(MatSetSizes(M, PETSC_DECIDE, PETSC_DECIDE, n, n));
   PETSC_CALL(MatSeqAIJSetPreallocation(M, m, NULL));
   PETSC_CALL(MatMPIAIJSetPreallocation(M, m, NULL, m, NULL));
   
   /* Set up the matrix: */
   PETSC_CALL(MatSetFromOptions(M));
   PETSC_CALL(MatSetUp(M));
   
   return 0;
   
}
