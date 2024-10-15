#include "Driver.hxx"

/* Initialize the calculation: */
int Driver::initialize(int argc, char* argv[], Array1D<double>& dt) {
   
   /* Get the input file name: */
   PAMPA_CHECK(argc < 2, "missing input file");
   std::string filename(argv[1]);
   
   /* Initialize MPI: */
   PAMPA_CHECK(mpi::initialize(argc, argv), "unable to initialize MPI");
   
   /* Initialize PETSc and SLEPc: */
   PAMPA_CHECK(petsc::initialize(argc, argv), "unable to initialize PETSc and SLEPc");
   
   /* Initialize the terminal output: */
   PAMPA_CHECK(output::initialize(), "unable to initialize the terminal output");
   
   /* Initialize the .vtk output: */
   PAMPA_CHECK(vtk::initialize(), "unable to initialize the .vtk output");
   
   /* Print info: */
   output::print("\nInitialize...");
   output::indent();
   
   /* Read the main input file: */
   Parser parser;
   output::print("Parse the input file...", true);
   PAMPA_CHECK(parser.read(filename, &mesh, &mesh_nodal, materials, solvers, dt), 
      "unable to parse " + filename);
   output::print("Done.", true);
   
   /* Initialize the solver: */
   bool transient = !(dt.empty());
   output::print("Initialize the solver...", true);
   PAMPA_CHECK(getMainSolver(), "unable to get the main solver");
   PAMPA_CHECK(solver->initialize(transient), "unable to initialize the solver");
   output::print("Done.", true);
   
   /* Print info: */
   output::outdent();
   output::print("Done.");
   
   return 0;
   
}

/* Get the solution: */
int Driver::solve(int n, double dt, double t) {
   
   /* Print info: */
   output::print("\n--------------------------------");
   if (n == 0)
      output::print("\nSolve steady state...\n");
   else
      output::print("\nSolve time step " + std::to_string(n) + "...\n");
   
   /* Get the solution: */
   double t1 = MPI_Wtime();
   PAMPA_CHECK(solver->solve(n, dt, t), "unable to get the solution");
   double t2 = MPI_Wtime();
   
   /* Output the solution: */
   output::print("", true && !(petsc::verbose));
   output::print("Solution time", t2-t1, true, 3, true);
   if (n > 0) {
      output::print("Time step size", dt, true, 3);
      output::print("Physical time", t, true, 3);
   }
   PAMPA_CHECK(solver->output(mpi::get_path(), n), "unable to output the solution");
   
   /* Print info: */
   output::print("\nDone.");
   
   return 0;
   
}

/* Finalize the calculation: */
int Driver::finalize() {
   
   /* Print info: */
   output::print("\n--------------------------------");
   output::print("\nFinalize...");
   output::indent();
   
   /* Finalize the solver: */
   output::print("Finalize the solver...", true);
   PAMPA_CHECK(solver->finalize(), "unable to finalize the solver");
   output::print("Done.", true);
   
   /* Finalize PETSc and SLEPc: */
   PAMPA_CHECK(petsc::finalize(), "unable to finalize PETSc and SLEPc");
   
   /* Finalize MPI: */
   PAMPA_CHECK(mpi::finalize(), "unable to finalize MPI");
   
   /* Free the meshes, materials and solvers: */
   utils::free(&mesh);
   utils::free(&mesh_nodal);
   for (int i = 0; i < materials.size(); i++)
      utils::free(&materials(i));
   for (int i = 0; i < solvers.size(); i++)
      utils::free(&solvers(i));
   
   /* Print info: */
   output::outdent();
   output::print("Done.\n");
   
   return 0;
   
}

/* Get the main solver: */
int Driver::getMainSolver() {
   
   /* Check if there's at least one solver: */
   PAMPA_CHECK(solvers.empty(), "no solvers defined");
   
   /* Look for the main solver: */
   if (solvers.size() == 1)
      solver = solvers(0);
   else {
      PAMPA_CHECK(utils::find("main", solvers, &solver), "unable to find the main solver");
   }
   
   return 0;
   
}
