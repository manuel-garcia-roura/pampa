#include <string>
#include <iostream>

#include "Parser.hxx"
#include "Model.hxx"

/* The main function: */
int main(int argc, char* argv[]) {
   
   std::string filename = "input.pmp";
   
   Parser parser;
   Model model;
   
   if (!parser.read(filename, model)) {
      std::cout << "Error: unable to parse " << filename << "!\n";
      return 1;
   }
   
   if (!((model.mesh)->build())) {
      std::cout << "Error: unable to build the mesh!" << std::endl;
      return 1;
   }
   
   if (!((model.mesh)->write("mesh.vtk"))) {
      std::cout << "Error: unable to write the mesh!" << std::endl;
      return 1;
   }
   
   return 0;
   
};
