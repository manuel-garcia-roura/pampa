#include <string>
#include <iostream>

#include "Parser.hxx"
#include "Model.hxx"
#include "Solver.hxx"
#include "utils.hxx"

/* The main function: */
int main(int argc, char* argv[]) {
   
   /* Main data: */
   Parser parser;
   Model model;
   Solver solver;
   
   /* Main input file: */
   std::string filename = "input.pmp";
   
   /* Read the main input file: */
   PAMPA_CALL(parser.read(filename, model), "unable to parse " + filename);
   
   /* Build the mesh: */
   PAMPA_CALL((model.mesh)->build(), "unable to build the mesh");
   
   /* Write the mesh: */
   PAMPA_CALL((model.mesh)->write("mesh.vtk"), "unable to write the mesh");
   
   /* Set the model materials: */
   (model.mesh)->setModelMaterials(&(model.materials));
   
   /* Initialize the solver: */
   PAMPA_CALL(solver.initialize(argc, argv, model), "unable to initialize the solver");
   
   /* Solve the eigensystem to get the flux and the multiplication factor: */
   PAMPA_CALL(solver.solve(), "unable to solve the eigensystem");
   
   /* Finalize the solver: */
   PAMPA_CALL(solver.finalize(), "unable to initialize the solver");
   
   return 0;
   
};
