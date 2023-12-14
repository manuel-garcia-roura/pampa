#include <iostream>

#include "CartesianMesh.hxx"
#include "UnstructuredExtrudedMesh.hxx"

int main() {
   
   std::cout << "Run..." << std::endl;
   
   CartesianMesh cartesian_mesh;
   UnstructuredExtrudedMesh unstructured_extruded_mesh;
   
   if (!cartesian_mesh.read("cartesian_mesh.dat"))
      std::cout << "Error reading cartesian_mesh.dat!" << std::endl;
   if (!cartesian_mesh.build())
      std::cout << "Error building cartesian_mesh!" << std::endl;
   if (!cartesian_mesh.write("cartesian_mesh.vtk"))
      std::cout << "Error writing cartesian_mesh.vtk!" << std::endl;
   
   std::cout << "----------------------------------------------------" << std::endl;
   
   if (!unstructured_extruded_mesh.read("unstructured_extruded_mesh.dat"))
      std::cout << "Error reading unstructured_extruded_mesh.dat!" << std::endl;
   if (!unstructured_extruded_mesh.build())
      std::cout << "Error building unstructured_extruded_mesh!" << std::endl;
   if (!unstructured_extruded_mesh.write("unstructured_extruded_mesh.vtk"))
      std::cout << "Error writing unstructured_extruded_mesh.vtk!" << std::endl;
   
   std::cout << "Done." << std::endl;
   
}
