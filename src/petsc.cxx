#include "petsc.hxx"

/* Random number context: */
namespace petsc {
   PetscRandom rctx = 0;
}

/* Initialize: */
int petsc::initialize(int argc, char* argv[]) {
   
   /* Initialize PETSc: */
   static char petsc_help[] = "Solver for the linear system A*x = b.\n";
   PETSC_CALL(PetscInitialize(&argc, &argv, (char*)0, petsc_help));
   
   /* Initialize SLEPc: */
   static char slepc_help[] = "Solver for the generalized eigensystem A*x = lambda*B*x.\n";
   PETSC_CALL(SlepcInitialize(&argc, &argv, (char*)0, slepc_help));
   
   /* Get the switch for verbose output: */
   PetscBool verbose, flag;
   PETSC_CALL(PetscOptionsGetBool(nullptr, nullptr, "-verbose", &verbose, &flag));
   mpi::verbose = flag && verbose;
   
   return 0;
   
}

/* Create, preallocate and set up a matrix: */
int petsc::create(Mat& M, int nl, int ng, int m, Array1D<Mat*>& matrices) {
   
   /* Create the matrix: */
   PETSC_CALL(MatCreate(MPI_COMM_WORLD, &M));
   matrices.pushBack(&M);
   
   /* Set the matrix options: */
   PETSC_CALL(MatSetFromOptions(M));
   
   /* Preallocate the matrix: */
   PETSC_CALL(MatSetSizes(M, nl, nl, ng, ng));
   PETSC_CALL(MatSeqAIJSetPreallocation(M, m, nullptr));
   PETSC_CALL(MatMPIAIJSetPreallocation(M, m, nullptr, m, nullptr));
   
   /* Set up the matrix: */
   PETSC_CALL(MatSetUp(M));
   
   return 0;
   
}

/* Create a vector from a matrix: */
int petsc::create(Vec& v, const Mat& M, Array1D<Vec*>& vectors) {
   
   /* Create the vector: */
   PETSC_CALL(MatCreateVecs(M, nullptr, &v));
   vectors.pushBack(&v);
   
   /* Set all elements to zero: */
   PETSC_CALL(VecZeroEntries(v));
   
   return 0;
   
}

/* Create a vector from its dimensions: */
int petsc::create(Vec& v, int nl, int ng, Array1D<Vec*>& vectors) {
   
   /* Create the vector: */
   PETSC_CALL(VecCreate(MPI_COMM_WORLD, &v));
   PETSC_CALL(VecSetSizes(v, nl, ng));
   vectors.pushBack(&v);
   
   /* Set the vector options: */
   PETSC_CALL(VecSetFromOptions(v));
   
   /* Set all elements to zero: */
   PETSC_CALL(VecZeroEntries(v));
   
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
   
   /* Set the initial guess to non-zero: */
   PETSC_CALL(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
   
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

/* Destroy a matrix: */
int petsc::destroy(Mat& M) {
   
   /* Destroy the matrix and reset it to 0: */
   PETSC_CALL(MatDestroy(&M));
   M = 0;
   
   return 0;
   
}

/* Destroy a vector: */
int petsc::destroy(Vec& v) {
   
   /* Destroy the vector and reset it to 0: */
   PETSC_CALL(VecDestroy(&v));
   v = 0;
   
   return 0;
   
}

/* Destroy a KSP context: */
int petsc::destroy(KSP& ksp) {
   
   /* Destroy the KSP context and reset it to 0: */
   PETSC_CALL(KSPDestroy(&ksp));
   ksp = 0;
   
   return 0;
   
}

/* Destroy an EPS context: */
int petsc::destroy(EPS& eps) {
   
   /* Destroy the EPS context and reset it to 0: */
   PETSC_CALL(EPSDestroy(&eps));
   eps = 0;
   
   return 0;
   
}

/* Initialize a vector with random values: */
int petsc::random(Vec& v) {
   
   /* Create the random number context: */
   if (rctx == 0) {
      PETSC_CALL(PetscRandomCreate(MPI_COMM_WORLD, &rctx));
   }
   
   /* Set the random values: */
   PETSC_CALL(VecSetRandom(v, rctx));
   
   return 0;
   
}

/* Normalize a vector by its 1-norm: */
int petsc::normalize(Vec& v, double x) {
   
   /* Get the 1-norm: */
   PetscScalar norm;
   PETSC_CALL(VecNorm(v, NORM_1, &norm));
   
   /* Scale the vector to get the 1-norm right: */
   if (fabs(norm) > 0.0) {
      PETSC_CALL(VecScale(v, x/norm));
   }
   else {
      PETSC_CALL(VecSet(v, 0.0));
   }
   
   return 0;
   
}

/* Get the difference between two vectors using a p-norm: */
int petsc::difference(const Vec& v1, const Vec& v2, double p, double& eps, bool relative) {
   
   /* Get the vector size: */
   PetscInt n;
   PETSC_CALL(VecGetLocalSize(v1, &n));
   
   /* Get the arrays with the raw data: */
   PetscScalar *v1_data, *v2_data;
   PETSC_CALL(VecGetArray(v1, &v1_data));
   PETSC_CALL(VecGetArray(v2, &v2_data));
   
   /* Get the p-norm of the difference: */
   eps = 0.0;
   for (int i = 0; i < n; i++)
      eps += std::pow(v1_data[i]-v2_data[i], p);
   MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, &eps, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
   eps = std::pow(eps, 1.0/p);
   
   /* Restore the arrays with the raw data: */
   PETSC_CALL(VecRestoreArray(v1, &v1_data));
   PETSC_CALL(VecRestoreArray(v2, &v2_data));
   
   /* Normalize the difference with the 2-norm of the first vector: */
   if (relative) {
      PetscScalar norm;
      PETSC_CALL(VecNorm(v1, NORM_2, &norm));
      eps /= norm;
   }
   
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
   PetscScalar ksp_residual_norm;
   PETSC_CALL(KSPGetType(ksp, &ksp_type));
   PETSC_CALL(KSPGetTotalIterations(ksp, &ksp_num_iterations));
   PETSC_CALL(KSPGetResidualNorm(ksp, &ksp_residual_norm));
   
   /* Print out the solver information: */
   PetscBool print, flag;
   PETSC_CALL(PetscOptionsGetBool(nullptr, nullptr, "-petsc_print_solver_info", &print, &flag));
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
   PetscScalar ksp_residual_norm;
   PETSC_CALL(EPSGetST(eps, &st));
   PETSC_CALL(STGetKSP(st, &ksp));
   PETSC_CALL(KSPGetType(ksp, &ksp_type));
   PETSC_CALL(KSPGetTotalIterations(ksp, &ksp_num_iterations));
   PETSC_CALL(KSPGetResidualNorm(ksp, &ksp_residual_norm));
   
   /* Print out the solver information: */
   PetscBool print, flag;
   PETSC_CALL(PetscOptionsGetBool(nullptr, nullptr, "-petsc_print_solver_info", &print, &flag));
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
   PETSC_CALL(PetscOptionsGetBool(nullptr, nullptr, "-petsc_write_solution", &write, &flag));
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

/* Finalize: */
int petsc::finalize() {
   
   /* Destroy the random number context: */
   PETSC_CALL(PetscRandomDestroy(&rctx));
   
   /* Finalize SLEPc: */
   PETSC_CALL(SlepcFinalize());
   
   /* Finalize PETSCc: */
   PETSC_CALL(PetscFinalize());
   
   return 0;
   
}
