#pragma once

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "Mesh.hxx"

/* The UnstructuredExtrudedMesh class */
class UnstructuredExtrudedMesh : public Mesh {
   
   private:
      
      /* Mesh dimensions in the xy-plane and in z: */
      int num_xy_points, num_xy_cells, nz;
      
      /* Mesh points in the xy-plane: */
      std::vector<std::array<double, 2>> xy_points;
      
      /* Mesh cells in the xy-plane: */
      std::vector<std::vector<int>> xy_cells;
      
      /* Mesh spacing in z: */
      std::vector<double> dz;
   
   public:
      
      /* The UnstructuredExtrudedMesh constructor: */
      UnstructuredExtrudedMesh();
      
      /* The UnstructuredExtrudedMesh destructor: */
      ~UnstructuredExtrudedMesh();
      
      /* Read the mesh from a plain-text file: */
      bool read(const std::string &filename);
      
      /* Build the mesh: */
      bool build();
   
};
