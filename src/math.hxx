#pragma once

#include <cmath>
#include <vector>

#include "utils.hxx"

/* The math namespace: */
namespace math {
   
   /* Get the area of a polygon: */
   double get_area(const Array2D<double>& pts, const std::vector<int>& ids);
   
   /* Get the centroid of a polygon: */
   std::vector<double> get_centroid(const Array2D<double>& pts, const std::vector<int>& ids, 
      double a);
   
   /* Get the distance between two points in n dimensions: */
   double get_distance(const Array2D<double>& pts, int i1, int i2, int n);
   
   /* Get the midpoint between two points in n dimensions: */
   std::vector<double> get_midpoint(const Array2D<double>& pts, int i1, int i2, int n);
   
   /* Get the normal between two points in 2 dimensions: */
   std::vector<double> get_normal(const Array2D<double>& pts, int i1, int i2);
   
   /* Extrude an edge to form a face: */
   std::vector<int> extrude_edge(const std::vector<int>& ids, int i1, int n);
   
   /* Get the distance between two points in n dimensions: */
   double get_distance(const double* p1, const double* p2, int n);
   
   /* Get the dot product of two vectors (x = v1 * v2) in n dimensions: */
   double dot_product(const double* v1, const double* v2, int n);
   
   /* Get the L2 norm of a vector in n dimensions: */
   double l2_norm(const double* v, int n);
   
   /* Get the L2 norm squared of a vector in n dimensions: */
   double l2_norm_2(const double* v, int n);
   
   /* Perform a SAXPY operation (v = a*x + y) in n dimensions: */
   void saxpy(double* v, double a, const double* x, const double* y, int n);
   
   /* Get the surface leakage factor for two centroids and a normal in 3 dimensions: */
   double surface_leakage_factor(const double* p1, const double* p2, const double* n);
   
}
