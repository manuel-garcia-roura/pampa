#pragma once

#include <cmath>
#include <vector>

/* The math namespace: */
namespace math {
   
   /* Calculate the area of a polygon: */
   double get_area(const std::vector<std::vector<double>> &pts, const std::vector<int> &ids);
   
   /* Calculate the centroid of a polygon: */
   std::vector<double> get_centroid(const std::vector<std::vector<double>> &pts, 
      const std::vector<int> &ids, double a);
   
   /* Calculate the distance between two points in n dimensions: */
   double get_distance(const std::vector<std::vector<double>> &pts, int i1, int i2, int n);
   
   /* Calculate the midpoint between two points in n dimensions: */
   std::vector<double> get_midpoint(const std::vector<std::vector<double>> &pts, int i1, int i2, 
      int n);
   
   /* Calculate the normal between two points in 2 dimensions: */
   std::vector<double> get_normal(const std::vector<std::vector<double>> &pts, int i1, int i2);
   
   /* Extrude an edge to form a face: */
   std::vector<int> extrude_edge(const std::vector<int> &ids, int i1, int n);
   
}
