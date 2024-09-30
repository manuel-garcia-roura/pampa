#pragma once

#include "petsc.hxx"
#include "utils.hxx"

/* The vtk namespace: */
namespace vtk {
   
   /* Switch for .vtk output: */
   extern bool on;
   
   /* Output interval in time steps: */
   extern int dn;
   
   /* Initialize: */
   int WARN_UNUSED initialize();
   
   /* Write a mesh to a .vtk file: */
   int WARN_UNUSED write(const std::string& prefix, int n, const Array2D<double>& points, 
      int num_points, const Vector2D<int>& cells, int num_cells, const Array1D<int>& materials, 
      const Array1D<int>& nodal_indices);
   
   /* Write a solution vector to a .vtk file: */
   int WARN_UNUSED write(const std::string& prefix, int n, const std::string& name, const Vec& v, 
      int num_cells, int num_groups = 1, int num_directions = 1);
   
}
