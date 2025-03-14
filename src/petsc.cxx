#include "petsc.hxx"

/* The petsc namespace: */
namespace petsc {
   
   /* Switch for verbose output: */
   bool verbose = false;
   
   /* Switch to write the solution in PETSc format: */
   bool dump = false;
   
   /* Random number context: */
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
   PAMPA_CHECK(get("-petsc_verbose", verbose), "unable to get the 'petsc_verbose' switch");
   
   /* Get the switch to write the solution in PETSc format: */
   PAMPA_CHECK(get("-petsc_dump", dump), "unable to get the 'petsc_dump' switch");
   
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

/* Get a switch from the command-line arguments: */
int petsc::get(const std::string& name, bool& on) {
   
   /* Get the switch: */
   PetscBool value, present;
   PETSC_CALL(PetscOptionsGetBool(nullptr, nullptr, name.c_str(), &value, &present));
   on = present && value;
   
   return 0;
   
}

/* Get an int value from the command-line arguments: */
int petsc::get(const std::string& name, int& x) {
   
   /* Get the value: */
   PetscInt value;
   PetscBool present;
   PETSC_CALL(PetscOptionsGetInt(nullptr, nullptr, name.c_str(), &value, &present));
   if (present) x = value;
   
   return 0;
   
}

/* Set an option: */
int petsc::set(const std::string& name, const std::string& value) {
   
   /* Set the option in the default global database: */
   std::string option = "-" + name;
   PETSC_CALL(PetscOptionsSetValue(NULL, option.c_str(), value.c_str()));
   
   return 0;
   
}

/* Create, preallocate and set up a matrix: */
int petsc::create(Mat& M, int nl, int ng, int m, Array1D<Mat*>& matrices, bool seq) {
   
   /* Create the matrix: */
   MPI_Comm comm = seq ? MPI_COMM_SELF : MPI_COMM_WORLD;
   PETSC_CALL(MatCreate(comm, &M));
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

/* Create a vector from its dimensions: */
int petsc::create(Vec& v, int nl, int ng, Array1D<Vec*>& vectors, bool seq) {
   
   /* Create the vector: */
   MPI_Comm comm = seq ? MPI_COMM_SELF : MPI_COMM_WORLD;
   PETSC_CALL(VecCreate(comm, &v));
   PETSC_CALL(VecSetSizes(v, nl, ng));
   vectors.pushBack(&v);
   
   /* Set the vector options: */
   PETSC_CALL(VecSetFromOptions(v));
   
   /* Set all elements to zero: */
   PETSC_CALL(VecZeroEntries(v));
   
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

/* Create a KSP context: */
int petsc::create(KSP& ksp, const Mat& A, bool seq) {
   
   /* Create the KSP context: */
   MPI_Comm comm = seq ? MPI_COMM_SELF : MPI_COMM_WORLD;
   PETSC_CALL(KSPCreate(comm, &ksp));
   
   /* Set the matrix: */
   PETSC_CALL(KSPSetOperators(ksp, A, A));
   
   /* Get the PC object: */
   PC pc;
   PETSC_CALL(KSPGetPC(ksp, &pc));
   
   /* Set the default options: */
   PETSC_CALL(KSPSetType(ksp, KSPCR));
   PETSC_CALL(PCSetType(pc, PCBJACOBI));
   
   /* Set the KSP options: */
   PETSC_CALL(KSPSetFromOptions(ksp));
   
   /* Set the initial guess to non-zero: */
   PETSC_CALL(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
   
   return 0;
   
}

/* Create an EPS context: */
int petsc::create(EPS& eps, const Mat& A, const Mat& B, bool seq) {
   
   /* Create the EPS context: */
   MPI_Comm comm = seq ? MPI_COMM_SELF : MPI_COMM_WORLD;
   PETSC_CALL(EPSCreate(comm, &eps));
   
   /* Set the matrices: */
   PETSC_CALL(EPSSetOperators(eps, A, B));
   
   /* Get the ST, KSP and PC objects: */
   ST st;
   KSP ksp;
   PC pc;
   PETSC_CALL(EPSGetST(eps, &st));
   PETSC_CALL(STGetKSP(st, &ksp));
   PETSC_CALL(KSPGetPC(ksp, &pc));
   
   /* Set the default options: */
   PETSC_CALL(EPSSetType(eps, EPSKRYLOVSCHUR));
   PETSC_CALL(STSetType(st, STSINVERT));
   PETSC_CALL(KSPSetType(ksp, KSPPREONLY));
   PETSC_CALL(PCSetType(pc, PCLU));
   PETSC_CALL(PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO));
   
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

/* Initialize a vector to a single value: */
int petsc::set(Vec& v, double x) {
   
   /* Set the vector values: */
   PETSC_CALL(VecSet(v, x));
   
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
   PETSC_CALL(VecScale(v, x/norm));
   
   return 0;
   
}

/* Normalize a vector by its 1-norm and add it to another vector: */
int petsc::normalize(Vec& v, double x, const Vec& v0) {
   
   /* Get the 1-norm: */
   PetscScalar norm;
   PETSC_CALL(VecNorm(v, NORM_1, &norm));
   
   /* Scale the vector to get the 1-norm right and add it to the other vector: */
   PETSC_CALL(VecAYPX(v, x/norm, v0));
   
   return 0;
   
}

/* Copy the values from a vector to a raw array: */
int petsc::copy(const Vec* v, double* x) {
   
   /* Get the vector size: */
   PetscInt n;
   PETSC_CALL(VecGetLocalSize(*v, &n));
   
   /* Get the array with the raw data: */
   PetscScalar* v_data;
   PETSC_CALL(VecGetArray(*v, &v_data));
   
   /* Copy the values: */
   for (int i = 0; i < n; i++)
      x[i] = v_data[i];
   
   /* Restore the array with the raw data: */
   PETSC_CALL(VecRestoreArray(*v, &v_data));
   
   return 0;
   
}

/* Copy the values from a raw array to a vector: */
int petsc::copy(const double* x, Vec* v) {
   
   /* Get the vector size: */
   PetscInt n;
   PETSC_CALL(VecGetLocalSize(*v, &n));
   
   /* Get the array with the raw data: */
   PetscScalar* v_data;
   PETSC_CALL(VecGetArray(*v, &v_data));
   
   /* Copy the values: */
   for (int i = 0; i < n; i++)
      v_data[i] = x[i];
   
   /* Restore the array with the raw data: */
   PETSC_CALL(VecRestoreArray(*v, &v_data));
   
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
   
   /* Print out the solver information: */
   if (verbose) {
      
      /* Get the PC object: */
      PC pc;
      PETSC_CALL(KSPGetPC(ksp, &pc));
      
      /* Get the KSP information: */
      KSPType ksp_type;
      PetscInt ksp_num_iterations;
      PetscScalar ksp_residual_norm;
      KSPConvergedReason ksp_converged_reason;
      PETSC_CALL(KSPGetType(ksp, &ksp_type));
      PETSC_CALL(KSPGetTotalIterations(ksp, &ksp_num_iterations));
      PETSC_CALL(KSPGetResidualNorm(ksp, &ksp_residual_norm));
      PETSC_CALL(KSPGetConvergedReason(ksp, &ksp_converged_reason));
      
      /* Get the PC information: */
      PCType pc_type;
      PETSC_CALL(PCGetType(pc, &pc_type));
      
      /* Print out the KSP information: */
      output::print("PETSc solver information:");
      output::indent();
      output::print("Elapsed time", t2-t1, true, 3);
      output::print("KSP type: " + std::string(ksp_type) + ".");
      output::print("Number of KSP iterations: " + std::to_string(ksp_num_iterations) + ".");
      output::print("KSP residual norm", ksp_residual_norm, true, 3);
      output::print("KSP convergence reason: " + std::to_string(ksp_converged_reason) + ".");
      output::print("PC type: " + std::string(pc_type) + ".");
      output::outdent();
      
   }
   
   return 0;
   
}

/* Solve an eigensystem: */
int petsc::solve(EPS& eps) {
   
   /* Solve the eigensystem: */
   double t1 = MPI_Wtime();
   PETSC_CALL(EPSSolve(eps));
   double t2 = MPI_Wtime();
   
   /* Print out the solver information: */
   if (verbose) {
      
      /* Get the ST, KSP and PC objects: */
      ST st;
      KSP ksp;
      PC pc;
      PETSC_CALL(EPSGetST(eps, &st));
      PETSC_CALL(STGetKSP(st, &ksp));
      PETSC_CALL(KSPGetPC(ksp, &pc));
      
      /* Get the EPS information: */
      EPSType eps_type;
      PetscInt eps_num_iterations;
      PETSC_CALL(EPSGetType(eps, &eps_type));
      PETSC_CALL(EPSGetIterationNumber(eps, &eps_num_iterations));
      
      /* Get the ST information: */
      EPSType st_type;
      PETSC_CALL(STGetType(st, &st_type));
      
      /* Get the KSP information: */
      KSPType ksp_type;
      PetscInt ksp_num_iterations;
      PetscScalar ksp_residual_norm;
      KSPConvergedReason ksp_converged_reason;
      PETSC_CALL(KSPGetType(ksp, &ksp_type));
      PETSC_CALL(KSPGetTotalIterations(ksp, &ksp_num_iterations));
      PETSC_CALL(KSPGetResidualNorm(ksp, &ksp_residual_norm));
      PETSC_CALL(KSPGetConvergedReason(ksp, &ksp_converged_reason));
      
      /* Get the PC information: */
      PCType pc_type;
      PETSC_CALL(PCGetType(pc, &pc_type));
      
      /* Print out the EPS and KSP information: */
      output::print("SLEPc solver information:");
      output::indent();
      output::print("Elapsed time", t2-t1, true, 3);
      output::print("EPS type: " + std::string(eps_type) + ".");
      output::print("Number of EPS iterations: " + std::to_string(eps_num_iterations) + ".");
      output::print("ST type: " + std::string(st_type) + ".");
      output::print("KSP type: " + std::string(ksp_type) + ".");
      output::print("Number of KSP iterations: " + std::to_string(ksp_num_iterations) + ".");
      output::print("KSP residual norm", ksp_residual_norm, true, 3);
      output::print("KSP convergence reason: " + std::to_string(ksp_converged_reason) + ".");
      output::print("PC type: " + std::string(pc_type) + ".");
      output::outdent();
      
   }
   
   return 0;
   
}

/* Write a solution vector to a PETSc binary file: */
int petsc::write(const std::string& prefix, int n, const Vec& v) {
   
   /* Write the solution: */
   if (dump) {
      
      /* Create the viewer: */
      PetscViewer viewer;
      std::string filename = prefix + "_" + std::to_string(n) + ".ptc";
      PETSC_CALL(PetscViewerBinaryOpen(MPI_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &viewer));
      
      /* Write the vector: */
      PETSC_CALL(VecView(v, viewer));
      
      /* Destroy the viewer: */
      PETSC_CALL(PetscViewerDestroy(&viewer));
      
   }
   
   return 0;
   
}
