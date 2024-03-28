#include "petsc.hxx"

/* Create, preallocate and set up a matrix: */
int petsc::create(Mat& M, int nl, int ng, int m, Array1D<Mat*>& matrices) {
   
   /* Create the matrix: */
   PETSC_CALL(MatCreate(MPI_COMM_WORLD, &M));
   matrices.pushBack(&M);
   
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
int petsc::create(Vec& v, const Mat& M, Array1D<Vec*>& vectors, bool random) {
   
   /* Create the vector: */
   PETSC_CALL(MatCreateVecs(M, NULL, &v));
   vectors.pushBack(&v);
   
   /* Initialize with random values: */
   if (random)
      PETSC_CALL(VecSetRandom(v, NULL));
   
   return 0;
   
}

/* Create a vector from its dimensions: */
int petsc::create(Vec& v, int nl, int ng, Array1D<Vec*>& vectors, bool random) {
   
   /* Create the vector: */
   PETSC_CALL(VecCreate(MPI_COMM_WORLD, &v));
   PETSC_CALL(VecSetSizes(v, nl, ng));
   vectors.pushBack(&v);
   
   /* Set the vector options: */
   PETSC_CALL(VecSetFromOptions(v));
   
   /* Initialize with random values: */
   if (random)
      PETSC_CALL(VecSetRandom(v, NULL));
   
   return 0;
   
}

/* Create a KSP context: */
int petsc::create(KSP& ksp, const Mat& A) {
   
   /* Create the KSP context: */
   PETSC_CALL(KSPCreate(MPI_COMM_WORLD, &ksp));
   
   /* Set the matrix: */
   PETSC_CALL(KSPSetOperators(ksp, A, A));
   
   /* Set the KSP options: */
   PETSC_CALL(KSPSetFromOptions(ksp));
   
   return 0;
   
}

/* Create an EPS context: */
int petsc::create(EPS& eps, const Mat& A, const Mat& B) {
   
   /* Create the EPS context: */
   PETSC_CALL(EPSCreate(MPI_COMM_WORLD, &eps));
   
   /* Set the matrices: */
   PETSC_CALL(EPSSetOperators(eps, A, B));
   
   /* Set the EPS options: */
   PETSC_CALL(EPSSetFromOptions(eps));
   
   return 0;
   
}

/* Normalize a vector by its 1-norm: */
int petsc::normalize(Vec& v, double x, bool random) {
   
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
int petsc::normalize(Vec& v, double x, const Vec& v0, bool random) {
   
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

/* Solve a linear system: */
int petsc::solve(KSP& ksp, const Vec& b, Vec& x) {
   
   /* Solve the linear system: */
   double t1 = MPI_Wtime();
   PETSC_CALL(KSPSolve(ksp, b, x));
   double t2 = MPI_Wtime();
   
   /* Get the KSP information: */
   KSPType ksp_type;
   PetscInt ksp_num_iterations;
   PetscReal ksp_residual_norm;
   PETSC_CALL(KSPGetType(ksp, &ksp_type));
   PETSC_CALL(KSPGetTotalIterations(ksp, &ksp_num_iterations));
   PETSC_CALL(KSPGetResidualNorm(ksp, &ksp_residual_norm));
   
   /* Print out the solver information: */
   PetscBool print, flag;
   PETSC_CALL(PetscOptionsGetBool(NULL, NULL, "-petsc_print_solver_info", &print, &flag));
   if (flag && print && mpi::rank == 0) {
      std::cout << "Elapsed time: " << t2-t1 << std::endl;
      std::cout << "KSP type: " << ksp_type << std::endl;
      std::cout << "Number of KSP iterations: " << ksp_num_iterations << std::endl;
      std::cout << "KSP residual norm: " << ksp_residual_norm << std::endl;
   }
   
   return 0;
   
}

/* Solve an eigensystem: */
int petsc::solve(EPS& eps) {
   
   /* Solve the eigensystem: */
   double t1 = MPI_Wtime();
   PETSC_CALL(EPSSolve(eps));
   double t2 = MPI_Wtime();
   
   /* Get the EPS information: */
   EPSType eps_type;
   PetscInt eps_num_iterations;
   PETSC_CALL(EPSGetType(eps, &eps_type));
   PETSC_CALL(EPSGetIterationNumber(eps, &eps_num_iterations));
   
   /* Get the KSP information: */
   ST st;
   KSP ksp;
   KSPType ksp_type;
   PetscInt ksp_num_iterations;
   PetscReal ksp_residual_norm;
   PETSC_CALL(EPSGetST(eps, &st));
   PETSC_CALL(STGetKSP(st, &ksp));
   PETSC_CALL(KSPGetType(ksp, &ksp_type));
   PETSC_CALL(KSPGetTotalIterations(ksp, &ksp_num_iterations));
   PETSC_CALL(KSPGetResidualNorm(ksp, &ksp_residual_norm));
   
   /* Print out the solver information: */
   PetscBool print, flag;
   PETSC_CALL(PetscOptionsGetBool(NULL, NULL, "-petsc_print_solver_info", &print, &flag));
   if (flag && print && mpi::rank == 0) {
      std::cout << "Elapsed time: " << t2-t1 << std::endl;
      std::cout << "EPS type: " << eps_type << std::endl;
      std::cout << "Number of EPS iterations: " << eps_num_iterations << std::endl;
      std::cout << "KSP type: " << ksp_type << std::endl;
      std::cout << "Number of KSP iterations: " << ksp_num_iterations << std::endl;
      std::cout << "KSP residual norm: " << ksp_residual_norm << std::endl;
   }
   
   return 0;
   
}

/* Write a solution vector to a PETSc binary file: */
int petsc::write(const std::string& filename, const Vec& v) {
   
   /* Check if the PETSc output is switched on: */
   PetscBool write, flag;
   PETSC_CALL(PetscOptionsGetBool(NULL, NULL, "-petsc_write_solution", &write, &flag));
   if (flag && write) {
      
      /* Create the viewer: */
      PetscViewer viewer;
      PETSC_CALL(PetscViewerBinaryOpen(MPI_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &viewer));
      
      /* Write the vector: */
      PETSC_CALL(VecView(v, viewer));
      
      /* Destroy the viewer: */
      PETSC_CALL(PetscViewerDestroy(&viewer));
      
   }
   
   return 0;
   
}
