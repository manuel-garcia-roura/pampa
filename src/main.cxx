#include <iostream>

#include "Parser.hxx"
#include "Config.hxx"
#include "Mesh.hxx"
#include "Material.hxx"
#include "CartesianMesh.hxx"
#include "UnstructuredExtrudedMesh.hxx"

int main() {
   
   Parser parser;
   Config config;
   Mesh mesh;
   std::vector<Material> materials;
   
   parser.read("input.pmp", config, mesh, materials);
   
   CartesianMesh cartesian_mesh;
   UnstructuredExtrudedMesh unstructured_extruded_mesh;
   
   if (!cartesian_mesh.read("cartesian_mesh.pmp"))
      std::cout << "Error reading cartesian_mesh.pmp!" << std::endl;
   if (!cartesian_mesh.build())
      std::cout << "Error building cartesian_mesh!" << std::endl;
   if (!cartesian_mesh.write("cartesian_mesh.vtk"))
      std::cout << "Error writing cartesian_mesh.vtk!" << std::endl;
   
   if (!unstructured_extruded_mesh.read("unstructured_extruded_mesh.pmp"))
      std::cout << "Error reading unstructured_extruded_mesh.pmp!" << std::endl;
   if (!unstructured_extruded_mesh.build())
      std::cout << "Error building unstructured_extruded_mesh!" << std::endl;
   if (!unstructured_extruded_mesh.write("unstructured_extruded_mesh.vtk"))
      std::cout << "Error writing unstructured_extruded_mesh.vtk!" << std::endl;
   
}
