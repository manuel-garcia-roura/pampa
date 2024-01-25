#pragma once

#include <slepceps.h>

#include "utils.hxx"

/* The petsc namespace: */
namespace petsc {
   
   /* Create, preallocate and set up a matrix: */
   int create_matrix(Mat& M, int n, int m);
   
}
