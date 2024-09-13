#pragma once

#include "utils.hxx"

/* The math namespace: */
namespace math {
   
   /* Get the area of a polygon with n sides: */
   double area(const Array2D<double>& pts, const int* ids, int n);
   
   /* Get the centroid of a polygon with n sides: */
   void centroid(double* c, const Array2D<double>& pts, const int* ids, int n, double a = -1.0);
   
   /* Get the distance between two points in n dimensions: */
   double distance(const Array2D<double>& pts, int i1, int i2, int n);
   
   /* Get the midpoint between two points in n dimensions: */
   void midpoint(double* p, const Array2D<double>& pts, int i1, int i2, int n);
   
   /* Get the normal between two points in 2 dimensions: */
   void normal(double* n, const Array2D<double>& pts, int i1, int i2);
   
   /* Get the distance between two points in n dimensions: */
   double distance(const double* p1, const double* p2, int n);
   
   /* Get the dot product of two vectors in n dimensions: */
   double dot_product(const double* v1, const double* v2, int n);
   
   /* Get the L2 norm of a vector in n dimensions: */
   double l2_norm(const double* v, int n);
   
   /* Get the L2 norm squared of a vector in n dimensions: */
   double l2_norm_2(const double* v, int n);
   
   /* Get the difference between two points in n dimensions: */
   void subtract(double* dp, const double* p1, const double* p2, int n);
   
   /* Perform a SAXPY operation (v = a*x + y) in n dimensions: */
   void saxpy(double* v, double a, const double* x, const double* y, int n);
   
   /* Get the surface leakage factor for two centroids and a normal in 3 dimensions: */
   double surface_leakage_factor(const double* p1, const double* p2, const double* n);
   
}
