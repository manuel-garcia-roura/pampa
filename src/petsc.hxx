#pragma once

#include <petsc.h>

#include "utils.hxx"

/* The petsc namespace: */
namespace petsc {
   
   /* Create, preallocate and set up a matrix: */
   int create_matrix(Mat& M, int nl, int ng, int m);
   
}
