#include <string>
#include <iostream>

#include <petsc.h>
#include <slepc.h>

#include "Parser.hxx"
#include "Mesh.hxx"
#include "Material.hxx"
#include "Solver.hxx"
#include "mpi.hxx"
#include "utils.hxx"

/* Initialize: */
int WARN_UNUSED initialize(int argc, char* argv[]) {
   
   /* Initialize MPI: */
   PAMPA_CALL(mpi::initialize(argc, argv), "unable to initialize MPI");
   
   /* Initialize PETSc: */
   static char petsc_help[] = "Solver for the linear system A*x = b.\n";
   PETSC_CALL(PetscInitialize(&argc, &argv, (char*)0, petsc_help));
   
   /* Initialize SLEPc: */
   static char slepc_help[] = "Solver for the generalized eigensystem A*x = lambda*B*x.\n";
   PETSC_CALL(SlepcInitialize(&argc, &argv, (char*)0, slepc_help));
   
   return 0;
   
}

/* Finalize: */
int WARN_UNUSED finalize() {
   
   /* Finalize SLEPc: */
   PETSC_CALL(SlepcFinalize());
   
   /* Finalize PETSCc: */
   PETSC_CALL(PetscFinalize());
   
   /* Finalize MPI: */
   PAMPA_CALL(mpi::finalize(), "unable to finalize MPI");
   
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
   PAMPA_CALL(solver->initialize(argc, argv), "unable to initialize the solver");
   
   /* Solve the linear system to get the steady-state solution: */
   PAMPA_CALL(solver->solve(), "unable to solve the linear system");
   
   /* Output the steady-state solution: */
   PAMPA_CALL(solver->output("output_0.vtk"), "unable to output the solution");
   
   /* Run the time-stepping loop: */
   for (int i = 0; i < dt.size(); i++) {
      
      /* Solve the linear system to get the transient solution: */
      PAMPA_CALL(solver->solve(i+1, dt(i)), "unable to solve the linear system");
      
      /* Output the transient solution: */
      PAMPA_CALL(solver->output("output_" + std::to_string(i+1) + ".vtk"), 
         "unable to output the solution");
      
   }
   
   /* Finalize the solver: */
   PAMPA_CALL(solver->finalize(), "unable to finalize the solver");
   
   /* Finalize the program: */
   PAMPA_CALL(finalize(), "unable to finalize the program");
   
   /* Free the mesh and the solver: */
   /* TODO: this should be done somewhere else. */
   delete mesh;
   delete solver;
   
   return 0;
   
}
