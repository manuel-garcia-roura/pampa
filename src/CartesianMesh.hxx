#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "Mesh.hxx"
#include "math.hxx"
#include "utils.hxx"

/* The CartesianMesh class: */
class CartesianMesh : public Mesh {
   
   private:
      
      /* Mesh dimensions in x, y and z: */
      int nx, ny, nz;
      
      /* Mesh spacing in x, y and z: */
      std::vector<double> dx, dy, dz;
      
      /* Boundary conditions in x, y and z (-x, +x, -y, +y, -z, +z): */
      int bc_x[2], bc_y[2], bc_z[2];
   
   public:
      
      /* The CartesianMesh constructor: */
      CartesianMesh();
      
      /* The CartesianMesh destructor: */
      ~CartesianMesh();
      
      /* Read the mesh from a plain-text input file: */
      int read(const std::string &filename);
      
      /* Build the mesh: */
      int build();
   
};
