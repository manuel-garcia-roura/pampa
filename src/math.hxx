#pragma once

#include <cmath>
#include <vector>

/* The math namespace: */
namespace math {
   
   /* Calculate the area of a polygon: */
   double get_area(const std::vector<std::vector<double>>& pts, const std::vector<int>& ids);
   
   /* Calculate the centroid of a polygon: */
   std::vector<double> get_centroid(const std::vector<std::vector<double>>& pts, 
      const std::vector<int>& ids, double a);
   
   /* Calculate the distance between two points in n dimensions: */
   double get_distance(const std::vector<std::vector<double>>& pts, int i1, int i2, int n);
   
   /* Calculate the midpoint between two points in n dimensions: */
   std::vector<double> get_midpoint(const std::vector<std::vector<double>>& pts, int i1, int i2, 
      int n);
   
   /* Calculate the normal between two points in 2 dimensions: */
   std::vector<double> get_normal(const std::vector<std::vector<double>>& pts, int i1, int i2);
   
   /* Extrude an edge to form a face: */
   std::vector<int> extrude_edge(const std::vector<int>& ids, int i1, int n);
   
   /* Subtract two vectors (v = v1 - v2) in n dimensions: */
   std::vector<double> subtract(const std::vector<double>& v1, const std::vector<double>& v2, 
      int n);
   
   /* Get the distance between two points in n dimensions: */
   double get_distance(const double* p1, const double* p2, int n);
   
   /* Get the dot product of two vectors (x = v1 * v2) in n dimensions: */
   double dot_product(const std::vector<double>& v1, const std::vector<double>& v2, int n);
   
   /* Get the L2 norm of a vector in n dimensions: */
   double l2_norm(const std::vector<double>& v, int n);
   
   /* Get the L2 norm squared of a vector in n dimensions: */
   double l2_norm_2(const std::vector<double>& v, int n);
   
   /* Perform a SAXPY operation (v = a*x + y) in n dimensions: */
   std::vector<double> saxpy(double a, const std::vector<double>& x, const std::vector<double>& y, 
      int n);
   
   /* Get the surface leakage factor for two centroids and a normal in 3 dimensions: */
   double surface_leakage_factor(const std::vector<double>& p1, const std::vector<double>& p2, 
      const std::vector<double>& n);
   
}
