#pragma once

#include <petsc.h>

#include "utils.hxx"

/* The petsc namespace: */
namespace petsc {
   
   /* Create, preallocate and set up a matrix: */
   int WARN_UNUSED create_matrix(Mat& M, int nl, int ng, int m);
   
   /* Create a vector from a matrix: */
   int WARN_UNUSED create_vector(Vec& v, const Mat& M, bool random = false);
   
   /* Create a vector from the dimensions: */
   int WARN_UNUSED create_vector(Vec& v, int nl, int ng, bool random = false);
   
   /* Write a vector to a binary file: */
   int WARN_UNUSED write(const std::string& filename, const Vec& v);
   
}
