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

/* Finalize: */
int WARN_UNUSED finalize(Mesh* mesh, Solver* solver) {
   
   /* Finalize PETSc and SLEPc: */
   PAMPA_CALL(petsc::finalize(), "unable to finalize PETSc and SLEPc");
   
   /* Finalize MPI: */
   PAMPA_CALL(mpi::finalize(), "unable to finalize MPI");
   
   /* Free the mesh and the solver: */
   delete mesh;
   delete solver;
   
   return 0;
   
}

/* The main function: */
int main(int argc, char* argv[]) {
   
   /* Main objects: */
   Parser parser;
   Mesh* mesh = NULL;
   Array1D<Material> materials;
   Solver* solver = NULL;
   Array1D<double> dt;
   
   /* Get the input file name: */
   PAMPA_CHECK(argc < 2, 1, "missing input file");
   std::string filename(argv[1]);
   
   /* Initialize the program: */
   PAMPA_CALL(initialize(argc, argv), "unable to initialize the program");
   
   /* Read the main input file: */
   PAMPA_CALL(parser.read(filename, &mesh, materials, &solver, dt), "unable to parse " + filename);
   
   /* Initialize the solver: */
   PAMPA_CALL(solver->initialize(), "unable to initialize the solver");
   
   /* Get the steady-state solution: */
   PAMPA_CALL(solver->solve(), "unable to solve the linear system");
   
   /* Output the steady-state solution: */
   PAMPA_CALL(solver->output("output_0.vtk"), "unable to output the solution");
   
   /* Run the time-stepping loop: */
   for (int n = 0; n < dt.size(); n++) {
      
      /* Get the transient solution: */
      PAMPA_CALL(solver->solve(n+1, dt(n)), "unable to solve the linear system");
      
      /* Output the transient solution: */
      filename = "output_" + std::to_string(n+1) + ".vtk";
      PAMPA_CALL(solver->output(filename, n+1), "unable to output the solution");
      
   }
   
   /* Finalize the solver: */
   PAMPA_CALL(solver->finalize(), "unable to finalize the solver");
   
   /* Finalize the program: */
   PAMPA_CALL(finalize(mesh, solver), "unable to finalize the program");
   
   return 0;
   
}
