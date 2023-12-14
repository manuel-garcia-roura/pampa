#include "math.hxx"

/* Calculate the area of a polygon: */
double math::get_area(const std::vector<std::vector<double>> &pts, const std::vector<int> &ids) {
   
   /* Calculate the area with the shoelace formula: */
   int n = ids.size();
   double a = 0.0;
   for (int i = 0; i < n; i++) {
      const std::vector<double> &p1 = pts[ids[i]];
      const std::vector<double> &p2 = pts[ids[(i+1)%n]];
      a += p1[0]*p2[1] - p2[0]*p1[1];
   }
   a *= 0.5;
   
   return a;
   
};

/* Calculate the centroid of a polygon: */
std::vector<double> math::get_centroid(const std::vector<std::vector<double>> &pts, 
   const std::vector<int> &ids, double a) {
   
   /* Calculate the area if not given: */
   if (a < 0.0)
      a = get_area(pts, ids);
   
   /* Calculate the centroid: */
   int n = ids.size();
   std::vector<double> p0(2);
   for (int i = 0; i < n; i++) {
      const std::vector<double> &p1 = pts[ids[i]];
      const std::vector<double> &p2 = pts[ids[(i+1)%n]];
      double da = p1[0]*p2[1] - p2[0]*p1[1];
      p0[0] += (p1[0]+p2[0]) * da;
      p0[1] += (p1[1]+p2[1]) * da;
   }
   p0[0] *= (1.0/(6.0*a));
   p0[1] *= (1.0/(6.0*a));
   
   return p0;
   
};

/* Calculate the distance between two points in n dimensions: */
double math::get_distance(const std::vector<std::vector<double>> &pts, int i1, int i2, int n) {
   
   /* Calculate the distance: */
   const std::vector<double> &p1 = pts[i1];
   const std::vector<double> &p2 = pts[i2];
   double d = 0.0;
   for (int i = 0; i < n; i++)
      d += std::pow(p2[i]-p1[i], 2);
   d = sqrt(d);
   
   return d;
   
};

/* Calculate the midpoint between two points in n dimensions: */
std::vector<double> math::get_midpoint(const std::vector<std::vector<double>> &pts, int i1, int i2, 
   int n) {
   
   /* Calculate the midpoint: */
   const std::vector<double> &p1 = pts[i1];
   const std::vector<double> &p2 = pts[i2];
   std::vector<double> p0(n);
   for (int i = 0; i < n; i++)
      p0[i] = 0.5 * (p1[i]+p2[i]);
   
   return p0;
   
};

/* Calculate the normal between two points in 2 dimensions: */
std::vector<double> math::get_normal(const std::vector<std::vector<double>> &pts, int i1, int i2) {
   
   /* Calculate the normal as (dy, -dx): */
   /* Note: (dy, -dx) should have the correct orientation, as opposed to (-dy, dx). */
   const std::vector<double> &p1 = pts[i1];
   const std::vector<double> &p2 = pts[i2];
   std::vector<double> n(2);
   n[0] = p2[1] - p1[1];
   n[1] = p1[0] - p2[0];
   
   /* Get a unit normal: */
   double norm = sqrt(n[0]*n[0]+n[1]*n[1]);
   n[0] /= norm;
   n[1] /= norm;
   
   return n;
   
};

/* Extrude an edge to form a face: */
std::vector<int> math::extrude_edge(const std::vector<int> &ids, int i1, int n) {
   
   /* Get the (counterclockwise oriented) extruded indexes: */
   return std::vector<int>{ids[i1], ids[(i1+1)%n], ids[n+(i1+1)%n], ids[n+i1]};
   
};
