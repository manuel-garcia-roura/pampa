#pragma once

#include <petsc.h>

#include "utils.hxx"

/* The petsc namespace: */
namespace petsc {
   
   /* Create, preallocate and set up a matrix: */
   int WARN_UNUSED create_matrix(Mat& M, int nl, int ng, int m);
   
   /* Create a vector from a matrix: */
   int WARN_UNUSED create_vector(Vec& v, const Mat& M, bool random = false);
   
   /* Create a vector from its dimensions: */
   int WARN_UNUSED create_vector(Vec& v, int nl, int ng, bool random = false);
   
   /* Normalize a vector by its 1-norm: */
   int WARN_UNUSED normalize_vector(Vec& v, double x, bool random = false);
   
   /* Normalize a vector by its 1-norm and add it to another vector: */
   int WARN_UNUSED normalize_vector(Vec& v, double x, Vec& v0, bool random = false);
   
   /* Write a vector to a binary file: */
   int WARN_UNUSED write(const std::string& filename, const Vec& v);
   
}
