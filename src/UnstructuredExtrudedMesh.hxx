#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "Mesh.hxx"
#include "utils.hxx"
#include "math.hxx"

/* The UnstructuredExtrudedMesh class: */
class UnstructuredExtrudedMesh : public Mesh {
   
   private:
      
      /* Mesh dimensions in the xy-plane and in z: */
      int num_xy_points, num_xy_cells, nz;
      
      /* Mesh points in the xy-plane: */
      std::vector<std::vector<double>> xy_points;
      
      /* Mesh cells in the xy-plane: */
      std::vector<std::vector<int>> xy_cells;
      
      /* Mesh spacing in z: */
      std::vector<double> dz;
   
   public:
      
      /* The UnstructuredExtrudedMesh constructor: */
      UnstructuredExtrudedMesh();
      
      /* The UnstructuredExtrudedMesh destructor: */
      ~UnstructuredExtrudedMesh();
      
      /* Read the mesh from a plain-text input file: */
      bool read(const std::string &filename);
      
      /* Build the mesh: */
      bool build();
   
};
