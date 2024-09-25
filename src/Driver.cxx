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
   
   /* Read the main input file: */
   Parser parser;
   PAMPA_CHECK(parser.read(filename, &mesh, &mesh_nodal, materials, solvers, dt), 
      "unable to parse " + filename);
   
   /* Get the main solver: */
   PAMPA_CHECK(getMainSolver(), "unable to get the main solver");
   
   /* Initialize the solver: */
   bool transient = !(dt.empty());
   PAMPA_CHECK(solver->initialize(transient), "unable to initialize the solver");
   
   return 0;
   
}

/* Get the solution: */
int Driver::solve(int n, double dt, double t) {
   
   /* Print the time-step number: */
   output::print("--------------------------------");
   output::print("n = " + std::to_string(n) + ":");
   
   /* Get the solution: */
   double t1 = MPI_Wtime();
   PAMPA_CHECK(solver->solve(n, dt, t), "unable to get the solution");
   double t2 = MPI_Wtime();
   output::print("Solution time: " + std::to_string(t2-t1));
   
   /* Output the solution: */
   std::string filename = "output_" + std::to_string(n) + ".vtk";
   PAMPA_CHECK(solver->output(mpi::get_path(filename), n), "unable to output the solution");
   
   return 0;
   
}

/* Finalize the calculation: */
int Driver::finalize() {
   
   /* Finalize the solver: */
   PAMPA_CHECK(solver->finalize(), "unable to finalize the solver");
   
   /* Finalize the terminal output: */
   PAMPA_CHECK(output::finalize(), "unable to finalize the terminal output");
   
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
   
   /* Print the last line: */
   output::print("--------------------------------");
   
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
