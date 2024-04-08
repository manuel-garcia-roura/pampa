#include <string>
#include <iostream>

#include "Parser.hxx"
#include "Mesh.hxx"
#include "Material.hxx"
#include "Solver.hxx"
#include "mpi.hxx"
#include "petsc.hxx"
#include "utils.hxx"

/* Initialize: */
int WARN_UNUSED initialize(int argc, char* argv[]) {
   
   /* Initialize MPI: */
   PAMPA_CALL(mpi::initialize(argc, argv), "unable to initialize MPI");
   
   /* Initialize PETSc and SLEPc: */
   PAMPA_CALL(petsc::initialize(argc, argv), "unable to initialize PETSc and SLEPc");
   
   return 0;
   
}

/* Get the main solver: */
int WARN_UNUSED get_main_solver(Array1D<Solver*>& solvers, Solver** solver) {
   
   /* Check if there's at least one solver: */
   PAMPA_CHECK(solvers.empty(), 1, "no solvers defined");
   
   /* Look for the main solver: */
   if (solvers.size() == 1)
      *solver = solvers(0);
   else {
      PAMPA_CALL(utils::find("main", solvers, solver), "unable to find main solver");
   }
   
   return 0;
   
}

/* Finalize: */
int WARN_UNUSED finalize(Mesh* mesh, Array1D<Solver*>& solvers) {
   
   /* Finalize PETSc and SLEPc: */
   PAMPA_CALL(petsc::finalize(), "unable to finalize PETSc and SLEPc");
   
   /* Finalize MPI: */
   PAMPA_CALL(mpi::finalize(), "unable to finalize MPI");
   
   /* Free the mesh and the solvers: */
   delete mesh;
   for (int i = 0; i < solvers.size(); i++)
      delete solvers(i);
   
   return 0;
   
}

/* The main function: */
int main(int argc, char* argv[]) {
   
   /* Main data structures: */
   Parser parser;
   Mesh* mesh = NULL;
   Array1D<Material> materials;
   Array1D<Solver*> solvers;
   Solver* solver;
   Array1D<double> dt;
   
   /* Get the input file name: */
   PAMPA_CHECK(argc < 2, 1, "missing input file");
   std::string filename(argv[1]);
   
   /* Initialize the program: */
   PAMPA_CALL(initialize(argc, argv), "unable to initialize the program");
   
   /* Read the main input file: */
   PAMPA_CALL(parser.read(filename, &mesh, materials, solvers, dt), "unable to parse " + filename);
   
   /* Get the main solver: */
   PAMPA_CALL(get_main_solver(solvers, &solver), "unable to get the main solver");
   
   /* Initialize the solver: */
   PAMPA_CALL(solver->initialize(), "unable to initialize the solver");
   
   /* Get the steady-state solution: */
   PAMPA_CALL(solver->solve(), "unable to get the steady-state solution");
   
   /* Output the steady-state solution: */
   PAMPA_CALL(solver->output(mpi::get_path("output_0.vtk")), "unable to output the solution");
   
   /* Run the time-stepping loop: */
   for (int n = 0; n < dt.size(); n++) {
      
      /* Get the transient solution: */
      PAMPA_CALL(solver->solve(n+1, dt(n)), "unable to get the transient solution");
      
      /* Output the transient solution: */
      PAMPA_CALL(solver->output(mpi::get_path("output_" + std::to_string(n+1) + ".vtk"), n+1), 
         "unable to output the solution");
      
   }
   
   /* Finalize the solver: */
   PAMPA_CALL(solver->finalize(), "unable to finalize the solver");
   
   /* Finalize the program: */
   PAMPA_CALL(finalize(mesh, solvers), "unable to finalize the program");
   
   return 0;
   
}
