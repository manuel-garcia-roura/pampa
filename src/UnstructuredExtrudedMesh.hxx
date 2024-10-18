#pragma once

#include "Mesh.hxx"

/* The UnstructuredExtrudedMesh class: */
class UnstructuredExtrudedMesh : public Mesh {
   
   private:
      
      /* Dimensions in the xy-plane and in z: */
      int num_xy_points = 0, num_xy_cells = 0, num_xy_boundaries = 0, nz = 0;
      
      /* Points in the xy-plane: */
      Array2D<double> xy_points;
      
      /* Cells in the xy-plane: */
      Vector2D<int> xy_cells;
      
      /* Boundary points in the xy-plane: */
      Vector2D<int> xy_boundary_points;
      
      /* Boundary names in the xy-plane: */
      Array1D<std::string> xy_boundary_names;
      
      /* Default boundary in the xy-plane: */
      int xy_default_boundary = -1;
      
      /* Spacing in z: */
      Array1D<double> dz{1, 0.0};
   
   public:
      
      /* The UnstructuredExtrudedMesh constructor: */
      UnstructuredExtrudedMesh() {partitioned = false;}
      
      /* The UnstructuredExtrudedMesh destructor: */
      ~UnstructuredExtrudedMesh() {}
      
      /* Read the mesh from a plain-text input file: */
      int WARN_UNUSED read(const std::string& filename);
      
      /* Build the mesh: */
      int WARN_UNUSED build();
   
};
