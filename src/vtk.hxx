#pragma once

#include <string>

#include "petsc.hxx"
#include "utils.hxx"

/* The vtk namespace: */
namespace vtk {
   
   /* Write a solution vector to a .vtk file: */
   int WARN_UNUSED write(const std::string& filename, const std::string& name, const Vec& v, 
      int num_cells, int num_groups = 1, int num_directions = 1);
   
}
