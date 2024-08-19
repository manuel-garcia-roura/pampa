#pragma once

#include <string>
#include <cmath>

#include <petsc.h>
#include <petscksp.h>

#include <slepc.h>
#include <slepceps.h>

#include "mpi.hxx"
#include "utils.hxx"

/* The petsc namespace: */
namespace petsc {
   
   /* Random number context: */
   extern PetscRandom rctx;
   
   /* Initialize: */
   int WARN_UNUSED initialize(int argc, char* argv[]);
   
   /* Create, preallocate and set up a matrix: */
   int WARN_UNUSED create(Mat& M, int nl, int ng, int m, Array1D<Mat*>& matrices);
   
   /* Create a vector from a matrix: */
   int WARN_UNUSED create(Vec& v, const Mat& M, Array1D<Vec*>& vectors);
   
   /* Create a vector from its dimensions: */
   int WARN_UNUSED create(Vec& v, int nl, int ng, Array1D<Vec*>& vectors);
   
   /* Create a KSP context: */
   int WARN_UNUSED create(KSP& ksp, const Mat& A);
   
   /* Create an EPS context: */
   int WARN_UNUSED create(EPS& eps, const Mat& A, const Mat& B);
   
   /* Destroy a matrix: */
   int WARN_UNUSED destroy(Mat& M);
   
   /* Destroy a vector: */
   int WARN_UNUSED destroy(Vec& v);
   
   /* Destroy a KSP context: */
   int WARN_UNUSED destroy(KSP& ksp);
   
   /* Destroy an EPS context: */
   int WARN_UNUSED destroy(EPS& eps);
   
   /* Initialize a vector to a single value: */
   int WARN_UNUSED set(Vec& v, double x);
   
   /* Initialize a vector with random values: */
   int WARN_UNUSED random(Vec& v);
   
   /* Normalize a vector by its 1-norm: */
   int WARN_UNUSED normalize(Vec& v, double x);
   
   /* Normalize a vector by its 1-norm and add it to another vector: */
   int WARN_UNUSED normalize(Vec& v, double x, const Vec& v0);
   
   /* Copy the values from a vector to a raw array: */
   int WARN_UNUSED copy(const Vec* v, double* x);
   
   /* Copy the values from a raw array to a vector: */
   int WARN_UNUSED copy(const double* x, Vec* v);
   
   /* Get the difference between two vectors using a p-norm: */
   int WARN_UNUSED difference(const Vec& v1, const Vec& v2, double p, double& eps, bool relative);
   
   /* Solve a linear system: */
   int WARN_UNUSED solve(KSP& ksp, const Vec& b, Vec& x);
   
   /* Solve an eigensystem: */
   int WARN_UNUSED solve(EPS& eps);
   
   /* Write a solution vector to a PETSc binary file: */
   int WARN_UNUSED write(const std::string& filename, const Vec& v);
   
   /* Finalize: */
   int WARN_UNUSED finalize();
   
}
