#pragma once

#include "Mesh.hxx"

/* The UnstructuredExtrudedMesh class: */
class UnstructuredExtrudedMesh : public Mesh {
   
   private:
      
      /* Mesh dimensions in the xy-plane and in z: */
      int num_xy_points = 0, num_xy_cells = 0, num_xy_boundaries = 0, nz = 0;
      
      /* Mesh points in the xy-plane: */
      Array2D<double> xy_points;
      
      /* Mesh cells in the xy-plane: */
      Vector2D<int> xy_cells;
      
      /* Mesh boundary points in the xy-plane: */
      Vector2D<int> xy_boundaries;
      
      /* Mesh spacing in z: */
      Array1D<double> dz{1, 0.0};
   
   public:
      
      /* The UnstructuredExtrudedMesh constructor: */
      UnstructuredExtrudedMesh() {partitioned = false;}
      
      /* The UnstructuredExtrudedMesh destructor: */
      ~UnstructuredExtrudedMesh() {}
      
      /* Read the mesh from a plain-text input file: */
      int read(const std::string& filename);
      
      /* Build the mesh: */
      int build();
   
};
