#pragma once

#include "Mesh.hxx"

/* The CartesianMesh class: */
class CartesianMesh : public Mesh {
   
   private:
      
      /* Mesh dimensions in x, y and z: */
      int nx = 1, ny = 0, nz = 0;
      
      /* Mesh spacing in x, y and z: */
      Array1D<double> dx{1, 1.0}, dy{1, 0.0}, dz{1, 0.0};
   
   public:
      
      /* The CartesianMesh constructor: */
      CartesianMesh() {partitioned = false;}
      
      /* The CartesianMesh destructor: */
      ~CartesianMesh() {}
      
      /* Read the mesh from a plain-text input file: */
      int read(const std::string& filename);
      
      /* Build the mesh: */
      int build();
   
};
