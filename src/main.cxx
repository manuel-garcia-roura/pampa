#include <iostream>

#include "Parser.hxx"
#include "Config.hxx"
#include "Mesh.hxx"
#include "Material.hxx"

int main() {
   
   Parser parser;
   Config config;
   Mesh *mesh;
   std::vector<Material> materials;
   
   parser.read("input.pmp", config, &mesh, materials);
   
   if (!(mesh->build())) std::cout << "Error building the mesh!" << std::endl;
   if (!(mesh->write("mesh.vtk"))) std::cout << "Error writing the mesh!" << std::endl;
   
};
