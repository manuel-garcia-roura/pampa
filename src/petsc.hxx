#pragma once

#include <string>

#include <petsc.h>
#include <petscksp.h>

#include <slepc.h>
#include <slepceps.h>

#include "mpi.hxx"
#include "utils.hxx"

/* The petsc namespace: */
namespace petsc {
   
   /* Create, preallocate and set up a matrix: */
   int WARN_UNUSED create(Mat& M, int nl, int ng, int m, Array1D<Mat*>& matrices);
   
   /* Create a vector from a matrix: */
   int WARN_UNUSED create(Vec& v, const Mat& M, Array1D<Vec*>& vectors, bool random = false);
   
   /* Create a vector from its dimensions: */
   int WARN_UNUSED create(Vec& v, int nl, int ng, Array1D<Vec*>& vectors, bool random = false);
   
   /* Create a KSP context: */
   int WARN_UNUSED create(KSP& ksp, const Mat& A);
   
   /* Create an EPS context: */
   int WARN_UNUSED create(EPS& eps, const Mat& A, const Mat& B);
   
   /* Normalize a vector by its 1-norm: */
   int WARN_UNUSED normalize(Vec& v, double x, bool random = false);
   
   /* Normalize a vector by its 1-norm and add it to another vector: */
   int WARN_UNUSED normalize(Vec& v, double x, const Vec& v0, bool random = false);
   
   /* Solve a linear system: */
   int WARN_UNUSED solve(KSP& ksp, const Vec& b, Vec& x);
   
   /* Solve an eigensystem: */
   int WARN_UNUSED solve(EPS& eps);
   
   /* Write a solution vector to a PETSc binary file: */
   int WARN_UNUSED write(const std::string& filename, const Vec& v);
   
}
