#include "Solver.hxx"

/* Initialize: */
int Solver::initialize(int argc, char* argv[]) {
   
   /* Check the materials: */
   for (int i = 0; i < materials.size(); i++) {
      PAMPA_CHECK(materials(i).method.type != method.type, 1, "wrong transport method");
      PAMPA_CHECK(materials(i).method.num_groups != method.num_groups, 2, 
         "wrong number of energy groups");
   }
   
   /* Initialize SLEPc: */
   static char help[] = "Solver for the generalized eigensystem R*x = (1/keff)*F*x.\n";
   PETSC_CALL(SlepcInitialize(&argc, &argv, (char*)0, help));
   
   /* Build the coefficient matrices and solution vectors: */
   PAMPA_CALL(build(), "unable to build the solver");
   
   /* Create the EPS context: */
   PETSC_CALL(EPSCreate(MPI_COMM_WORLD, &eps));
   PETSC_CALL(EPSSetOperators(eps, R, F));
   PETSC_CALL(EPSSetFromOptions(eps));
   
   /* Get the initial condition, if given: */
   PetscBool flag;
   char filename[PETSC_MAX_PATH_LEN];
   PETSC_CALL(PetscOptionsGetString(NULL, NULL, "-petsc_initial_condition", filename, 
      PETSC_MAX_PATH_LEN, &flag));
   if (flag) {
      PetscViewer viewer;
      Vec x0;
      PETSC_CALL(PetscViewerBinaryOpen(MPI_COMM_WORLD, filename, FILE_MODE_READ, &viewer));
      PETSC_CALL(VecCreate(MPI_COMM_WORLD, &x0));
      PETSC_CALL(VecLoad(x0, viewer));
      PETSC_CALL(EPSSetInitialSpace(eps, 1, &x0));
      PETSC_CALL(VecDestroy(&x0));
      PETSC_CALL(PetscViewerDestroy(&viewer));
   }
   
   return 0;
   
}

/* Solve the eigensystem to get the neutron flux and the multiplication factor: */
int Solver::solve() {
   
   /* Solve the eigensystem: */
   double t1 = MPI_Wtime();
   PETSC_CALL(EPSSolve(eps));
   double t2 = MPI_Wtime();
   
   /* Get the solver information: */
   ST st;
   KSP ksp;
   EPSType eps_type;
   PetscInt num_eigenvalues, max_eps_iterations, num_eps_iterations, num_ksp_iterations;
   PetscReal eps_tol;
   PETSC_CALL(EPSGetST(eps, &st));
   PETSC_CALL(STGetKSP(st, &ksp));
   PETSC_CALL(EPSGetType(eps, &eps_type));
   PETSC_CALL(EPSGetTolerances(eps, &eps_tol, &max_eps_iterations));
   PETSC_CALL(EPSGetIterationNumber(eps, &num_eps_iterations));
   PETSC_CALL(KSPGetTotalIterations(ksp, &num_ksp_iterations));
   
   /* Print out the solver information: */
   PetscBool print, flag;
   PETSC_CALL(PetscOptionsGetBool(NULL, NULL, "-petsc_print_info", &print, &flag));
   if (flag && print && mpi::rank == 0) {
      std::cout << "Elapsed time: " << t2-t1 << std::endl;
      std::cout << "Solution method: " << eps_type << std::endl;
      std::cout << "EPS convergence tolerance: " << eps_tol << std::endl;
      std::cout << "Maximum number of EPS iterations: " << max_eps_iterations << std::endl;
      std::cout << "Number of EPS iterations: " << num_eps_iterations << std::endl;
      std::cout << "Number of KSP iterations: " << num_ksp_iterations << std::endl;
   }
   
   /* Get the solution after solving the eigensystem: */
   PAMPA_CALL(getSolution(), "unable to get the solution");
   
   return 0;
   
}

/* Output the solution: */
int Solver::output(const std::string& filename) {
   
   /* Print out the multiplication factor: */
   if (mpi::rank == 0)
      std::cout << "keff = " << keff << std::endl;
   
   /* Write to a rank directory in parallel runs: */
   std::string path;
   if (mpi::size > 1) {
      std::string dir = std::to_string(mpi::rank);
      PAMPA_CALL(utils::create(dir), "unable to create the output directory");
      path = dir + "/" + filename;
   }
   else
      path = filename;
   
   /* Write the solution to a plain-text file in .vtk format: */
   PAMPA_CALL(writeVTK(path), "unable to output the solution in .vtk format");
   
   /* Write the solution to a binary file in PETSc format: */
   PetscBool write, flag;
   PETSC_CALL(PetscOptionsGetBool(NULL, NULL, "-petsc_write_solution", &write, &flag));
   if (flag && write) {
      PAMPA_CALL(writePETSc("flux.ptc"), "unable to output the solution in PETSc format");
   }
   
   return 0;
   
}

/* Finalize: */
int Solver::finalize() {
   
   /* Destroy the EPS context: */
   PETSC_CALL(EPSDestroy(&eps));
   
   /* Destroy the solution vectors: */
   PAMPA_CALL(destroyVectors(), "unable to destroy the solution vectors");
   
   /* Destroy the coefficient matrices: */
   PETSC_CALL(MatDestroy(&R));
   PETSC_CALL(MatDestroy(&F));
   
   /* Finalize SLEPc: */
   PETSC_CALL(SlepcFinalize());
   
   return 0;
   
}

/* Normalize the scalar flux: */
int Solver::normalizeScalarFlux() {
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   const Cells& cells = mesh->getCells();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Get the array for the scalar flux: */
   PetscScalar* data_phi;
   PETSC_CALL(VecGetArray(phi, &data_phi));
   
   /* Normalize the scalar flux (TODO: normalize correctly with the power!): */
   double vol = 0.0;
   for (int i = 0; i < num_cells; i++)
      if (materials(cells.materials(i)).nu_sigma_fission(1) > 0.0)
         vol += cells.volumes(i);
   MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
   double sum = 0.0;
   for (int iphi = 0, i = 0; i < num_cells; i++)
      for (int g = 0; g < num_groups; g++)
         sum += data_phi[iphi++] * materials(cells.materials(i)).nu_sigma_fission(g) * 
            cells.volumes(i);
   MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
   double f = vol / sum;
   for (int iphi = 0, i = 0; i < num_cells; i++)
      for (int g = 0; g < num_groups; g++)
         data_phi[iphi++] *= f;
   
   /* Restore the array for the scalar flux: */
   PETSC_CALL(VecRestoreArray(phi, &data_phi));
   
   return 0;
   
}
