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

/* Create a vector from a matrix: */
int petsc::create_vector(Vec& v, const Mat& M, bool random) {
   
   /* Create the vector: */
   PETSC_CALL(MatCreateVecs(M, NULL, &v));
   
   /* Initialize with random values: */
   if (random)
      PETSC_CALL(VecSetRandom(v, NULL));
   
   return 0;
   
}

/* Create a vector from its dimensions: */
int petsc::create_vector(Vec& v, int nl, int ng, bool random) {
   
   /* Create the vector: */
   PETSC_CALL(VecCreate(MPI_COMM_WORLD, &v));
   PETSC_CALL(VecSetSizes(v, nl, ng));
   
   /* Set the vector options: */
   PETSC_CALL(VecSetFromOptions(v));
   
   /* Initialize with random values: */
   if (random)
      PETSC_CALL(VecSetRandom(v, NULL));
   
   return 0;
   
}

/* Normalize a vector by its 1-norm: */
int petsc::normalize_vector(Vec& v, double x, bool random) {
   
   /* Initialize with random values: */
   if (random)
      PETSC_CALL(VecSetRandom(v, NULL));
   
   /* Get the 1-norm: */
   double x0;
   PETSC_CALL(VecNorm(v, NORM_1, &x0));
   
   /* Scale the vector to get the 1-norm right: */
   PETSC_CALL(VecScale(v, x/x0));
   
   return 0;
   
}

/* Normalize a vector by its 1-norm and add it to another vector: */
int petsc::normalize_vector(Vec& v, double x, Vec& v0, bool random) {
   
   /* Initialize with random values: */
   if (random)
      PETSC_CALL(VecSetRandom(v, NULL));
   
   /* Get the 1-norm: */
   double x0;
   PETSC_CALL(VecNorm(v, NORM_1, &x0));
   
   /* Scale the vector to get the 1-norm right and add it to the other vector: */
   PETSC_CALL(VecAYPX(v, x/x0, v0));
   
   return 0;
   
}

/* Write a vector to a binary file: */
int petsc::write(const std::string& filename, const Vec& v) {
   
   /* Create the viewer: */
   PetscViewer viewer;
   PETSC_CALL(PetscViewerBinaryOpen(MPI_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &viewer));
   
   /* Write the vector: */
   PETSC_CALL(VecView(v, viewer));
   
   /* Destroy the viewer: */
   PETSC_CALL(PetscViewerDestroy(&viewer));
   
   return 0;
   
}
