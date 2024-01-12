#include <string>
#include <iostream>

#include "Parser.hxx"
#include "Model.hxx"
#include "Solver.hxx"

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
   
   /* Build the model: */
   PAMPA_CALL(model.build(), "unable to build the model");
   
   /* Initialize the solver: */
   PAMPA_CALL(solver.initialize(argc, argv, model), "unable to initialize the solver");
   
   /* Solve the eigensystem to get the flux and the multiplication factor: */
   PAMPA_CALL(solver.solve(), "unable to solve the eigensystem");
   
   /* Output the solution: */
   PAMPA_CALL(solver.output("output.vtk", model), "unable to output the solution");
   
   /* Finalize the solver: */
   PAMPA_CALL(solver.finalize(), "unable to initialize the solver");
   
   return 0;
   
}
