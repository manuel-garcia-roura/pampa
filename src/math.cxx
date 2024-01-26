#include "math.hxx"

/* Calculate the area of a polygon: */
double math::get_area(const std::vector<std::vector<double>>& pts, const std::vector<int>& ids) {
   
   /* Calculate the area with the shoelace formula: */
   int n = ids.size();
   double a = 0.0;
   for (int i = 0; i < n; i++) {
      const std::vector<double>& p1 = pts[ids[i]];
      const std::vector<double>& p2 = pts[ids[(i+1)%n]];
      a += p1[0]*p2[1] - p2[0]*p1[1];
   }
   a *= 0.5;
   
   return a;
   
}

/* Calculate the centroid of a polygon: */
std::vector<double> math::get_centroid(const std::vector<std::vector<double>>& pts, 
   const std::vector<int>& ids, double a) {
   
   /* Calculate the area if not given: */
   if (a < 0.0) a = get_area(pts, ids);
   
   /* Calculate the centroid: */
   int n = ids.size();
   std::vector<double> p0(2);
   for (int i = 0; i < n; i++) {
      const std::vector<double>& p1 = pts[ids[i]];
      const std::vector<double>& p2 = pts[ids[(i+1)%n]];
      double da = p1[0]*p2[1] - p2[0]*p1[1];
      p0[0] += (p1[0]+p2[0]) * da;
      p0[1] += (p1[1]+p2[1]) * da;
   }
   p0[0] *= (1.0/(6.0*a));
   p0[1] *= (1.0/(6.0*a));
   
   return p0;
   
}

/* Calculate the distance between two points in n dimensions: */
double math::get_distance(const std::vector<std::vector<double>>& pts, int i1, int i2, int n) {
   
   /* Calculate the distance: */
   const std::vector<double>& p1 = pts[i1];
   const std::vector<double>& p2 = pts[i2];
   double d = 0.0;
   for (int i = 0; i < n; i++)
      d += std::pow(p2[i]-p1[i], 2);
   d = sqrt(d);
   
   return d;
   
}

/* Calculate the midpoint between two points in n dimensions: */
std::vector<double> math::get_midpoint(const std::vector<std::vector<double>>& pts, int i1, int i2, 
   int n) {
   
   /* Calculate the midpoint: */
   const std::vector<double>& p1 = pts[i1];
   const std::vector<double>& p2 = pts[i2];
   std::vector<double> p0(n);
   for (int i = 0; i < n; i++)
      p0[i] = 0.5 * (p1[i]+p2[i]);
   
   return p0;
   
}

/* Calculate the normal between two points in 2 dimensions: */
std::vector<double> math::get_normal(const std::vector<std::vector<double>>& pts, int i1, int i2) {
   
   /* Calculate the normal as (dy, -dx): */
   /* Note: (dy, -dx) should have the correct orientation, as opposed to (-dy, dx). */
   const std::vector<double>& p1 = pts[i1];
   const std::vector<double>& p2 = pts[i2];
   std::vector<double> n(2);
   n[0] = p2[1] - p1[1];
   n[1] = p1[0] - p2[0];
   
   /* Get a unit normal: */
   double norm = sqrt(n[0]*n[0]+n[1]*n[1]);
   n[0] /= norm;
   n[1] /= norm;
   
   return n;
   
}

/* Extrude an edge to form a face: */
std::vector<int> math::extrude_edge(const std::vector<int>& ids, int i1, int n) {
   
   /* Get the (counterclockwise oriented) extruded indexes: */
   std::vector<int> ids3d{ids[i1], ids[(i1+1)%n], ids[n+(i1+1)%n], ids[n+i1]};
   
   return ids3d;
   
}

/* Subtract two vectors (v = v1 - v2) in n dimensions: */
std::vector<double> math::subtract(const std::vector<double>& v1, const std::vector<double>& v2, 
   int n) {
   
   /* Get the vector difference: */
   std::vector<double> v(n);
   for (int i = 0; i < n; i++)
      v[i] = v1[i] - v2[i];
   
   return v;
   
}

/* Get the dot product of two vectors (x = v1 * v2) in n dimensions: */
double math::dot_product(const std::vector<double>& v1, const std::vector<double>& v2, int n) {
   
   /* Get the dot product: */
   double x = 0.0;
   for (int i = 0; i < n; i++)
      x += v1[i] * v2[i];
   
   return x;
   
}

/* Get the L2 norm of a vector in n dimensions: */
double math::l2_norm(const std::vector<double>& v, int n) {
   
   /* Get the L2 norm: */
   double x = 0.0;
   for (int i = 0; i < n; i++)
      x += v[i] * v[i];
   x = sqrt(x);
   
   return x;
   
}

/* Get the L2 norm squared of a vector in n dimensions: */
double math::l2_norm_2(const std::vector<double>& v, int n) {
   
   /* Get the L2 norm squared: */
   double x = 0.0;
   for (int i = 0; i < n; i++)
      x += v[i] * v[i];
   
   return x;
   
}

/* Perform a SAXPY operation (v = a*x + y) in n dimensions: */
std::vector<double> math::saxpy(double a, const std::vector<double>& x, 
   const std::vector<double>& y, int n) {
   
   /* Get the SAXPY vector: */
   std::vector<double> v(n);
   for (int i = 0; i < n; i++)
      v[i] = a*x[i] + y[i];
   
   return v;
   
}

/* Get the surface leakage factor for two centroids and a normal in 3 dimensions: */
double math::surface_leakage_factor(const std::vector<double>& p1, const std::vector<double>& p2, 
   const std::vector<double>& n) {
   
   /* Get the difference between the two centroids: */
   std::vector<double> dp = math::subtract(p2, p1, 3);
   
   /* Get the surface leakage factor w = (p2-p1)*n / |p2-p1|^2: */
   double w = math::dot_product(dp, n, 3) / math::l2_norm_2(dp, 3);
   
   return w;
   
}
