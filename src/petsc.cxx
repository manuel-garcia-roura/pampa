#include "petsc.hxx"

/* Create, preallocate and set up a matrix: */
int petsc::create_matrix(Mat& M, int nl, int ng, int m) {
   
   /* Create the matrix: */
   PETSC_CALL(MatCreate(MPI_COMM_WORLD, &M));
   
   /* Set the matrix options: */
   PETSC_CALL(MatSetFromOptions(M));
   
   /* Preallocate the matrix: */
   PETSC_CALL(MatSetSizes(M, nl, nl, ng, ng));
   PETSC_CALL(MatSeqAIJSetPreallocation(M, m, NULL));
   PETSC_CALL(MatMPIAIJSetPreallocation(M, m, NULL, m, NULL));
   
   /* Set up the matrix: */
   PETSC_CALL(MatSetUp(M));
   
   return 0;
   
}
