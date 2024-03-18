#include <string>
#include <iostream>

#include "Parser.hxx"
#include "Mesh.hxx"
#include "Material.hxx"
#include "Solver.hxx"
#include "DiffusionSolver.hxx"
#include "SNSolver.hxx"
#include "mpi.hxx"
#include "utils.hxx"

/* The main function: */
int main(int argc, char* argv[]) {
   
   /* Main objects (parser, mesh, materials, solver): */
   Parser parser;
   Mesh* mesh = NULL;
   Array1D<Material> materials;
   TransportMethod method;
   GradientScheme gradient;
   Solver* solver = NULL;
   
   /* Get the input file name: */
   PAMPA_CHECK(argc < 2, 1, "missing input file");
   std::string filename(argv[1]);
   
   /* Initialize MPI: */
   PAMPA_CALL(mpi::initialize(argc, argv), "unable to initialize MPI");
   
   /* Read the main input file: */
   PAMPA_CALL(parser.read(filename, &mesh, materials, method, gradient), 
      "unable to parse " + filename);
   
   /* Build the mesh: */
   PAMPA_CALL(mesh->build(), "unable to build the mesh");
   
   /* Partition the mesh and swap the meshes: */
   if (mpi::size > 1 && !mesh->isPartitioned()) {
      Mesh* submesh = NULL;
      PAMPA_CALL(mesh->partition(&submesh), "unable to partition the mesh");
      delete mesh;
      mesh = submesh;
   }
   
   /* Create the solver: */
   /* TODO: this should be done somewhere else. */
   switch (method.type) {
      case TM::DIFFUSION : {
         solver = new DiffusionSolver(mesh, materials, method);
         break;
      }
      case TM::SN : {
         solver = new SNSolver(mesh, materials, method, gradient);
         break;
      }
   }
   
   /* Initialize the solver: */
   PAMPA_CALL(solver->initialize(argc, argv), "unable to initialize the solver");
   
   /* Solve the eigensystem to get the flux and the multiplication factor: */
   PAMPA_CALL(solver->solve(), "unable to solve the eigensystem");
   
   /* Output the solution: */
   PAMPA_CALL(solver->output("output.vtk"), "unable to output the solution");
   
   /* Finalize the solver: */
   PAMPA_CALL(solver->finalize(), "unable to initialize the solver");
   
   /* Finalize MPI: */
   PAMPA_CALL(mpi::finalize(), "unable to finalize MPI");
   
   /* Free the mesh and the solver: */
   /* TODO: this should be done somewhere else. */
   delete mesh;
   delete solver;
   
   return 0;
   
}
