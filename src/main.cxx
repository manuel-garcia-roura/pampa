#include <string>
#include <iostream>

#include "Parser.hxx"
#include "Model.hxx"
#include "utils.hxx"

/* The main function: */
int main(int argc, char* argv[]) {
   
   /* Main data: */
   Parser parser;
   Model model;
   
   /* Main input file: */
   std::string filename = "input.pmp";
   
   /* Read the main input file: */
   PAMPA_CALL(parser.read(filename, model), "unable to parse " + filename);
   
   /* Build the mesh: */
   PAMPA_CALL((model.mesh)->build(), "unable to build the mesh");
   
   /* Write the mesh: */
   PAMPA_CALL((model.mesh)->write("mesh.vtk"), "unable to write the mesh");
   
   /* Initialize the solver: */
   model.solver.initialize(argc, argv);
   
   /* Finalize the solver: */
   model.solver.finalize();
   
   return 0;
   
};
