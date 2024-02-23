#include "math.hxx"

/* Get the area of a polygon with n sides: */
double math::area(const Array2D<double>& pts, const int* ids, int n) {
   
   /* Get the area with the shoelace formula: */
   double a = 0.0;
   for (int i = 0; i < n; i++) {
      const double* p1 = pts(ids[i]);
      const double* p2 = pts(ids[(i+1)%n]);
      a += p1[0]*p2[1] - p2[0]*p1[1];
   }
   a *= 0.5;
   
   return a;
   
}

/* Get the centroid of a polygon with n sides: */
void math::centroid(double* c, const Array2D<double>& pts, const int* ids, int n, double a) {
   
   /* Get the area if not given: */
   if (a < 0.0) a = area(pts, ids, n);
   
   /* Get the centroid: */
   for (int i = 0; i < n; i++) {
      const double* p1 = pts(ids[i]);
      const double* p2 = pts(ids[(i+1)%n]);
      double da = p1[0]*p2[1] - p2[0]*p1[1];
      c[0] += (p1[0]+p2[0]) * da;
      c[1] += (p1[1]+p2[1]) * da;
   }
   c[0] *= (1.0/(6.0*a));
   c[1] *= (1.0/(6.0*a));
   
}

/* Get the distance between two points in n dimensions: */
double math::distance(const Array2D<double>& pts, int i1, int i2, int n) {
   
   /* Get the distance: */
   const double* p1 = pts(i1);
   const double* p2 = pts(i2);
   double d = 0.0;
   for (int i = 0; i < n; i++)
      d += std::pow(p2[i]-p1[i], 2);
   d = sqrt(d);
   
   return d;
   
}

/* Get the midpoint between two points in n dimensions: */
void math::midpoint(double* p, const Array2D<double>& pts, int i1, int i2, int n) {
   
   /* Get the midpoint: */
   const double* p1 = pts(i1);
   const double* p2 = pts(i2);
   for (int i = 0; i < n; i++)
      p[i] = 0.5 * (p1[i]+p2[i]);
   
}

/* Get the normal between two points in 2 dimensions: */
void math::normal(double* n, const Array2D<double>& pts, int i1, int i2) {
   
   /* Get the normal as (dy, -dx): */
   /* Note: (dy, -dx) should have the correct orientation, as opposed to (-dy, dx). */
   const double* p1 = pts(i1);
   const double* p2 = pts(i2);
   n[0] = p2[1] - p1[1];
   n[1] = p1[0] - p2[0];
   
   /* Get a unit normal: */
   double norm = sqrt(n[0]*n[0]+n[1]*n[1]);
   n[0] /= norm;
   n[1] /= norm;
   
}

/* Get the distance between two points in n dimensions: */
double math::distance(const double* p1, const double* p2, int n) {
   
   /* Get the distance: */
   double d = 0.0;
   for (int i = 0; i < n; i++)
      d += std::pow(p2[i]-p1[i], 2);
   d = sqrt(d);
   
   return d;
   
}

/* Get the dot product of two vectors in n dimensions: */
double math::dot_product(const double* v1, const double* v2, int n) {
   
   /* Get the dot product: */
   double x = 0.0;
   for (int i = 0; i < n; i++)
      x += v1[i] * v2[i];
   
   return x;
   
}

/* Get the L2 norm of a vector in n dimensions: */
double math::l2_norm(const double* v, int n) {
   
   /* Get the L2 norm: */
   double x = 0.0;
   for (int i = 0; i < n; i++)
      x += v[i] * v[i];
   x = sqrt(x);
   
   return x;
   
}

/* Get the L2 norm squared of a vector in n dimensions: */
double math::l2_norm_2(const double* v, int n) {
   
   /* Get the L2 norm squared: */
   double x = 0.0;
   for (int i = 0; i < n; i++)
      x += v[i] * v[i];
   
   return x;
   
}

/* Get the difference between two points in n dimensions: */
void math::subtract(double* dp, const double* p1, const double* p2, int n) {
   
   /* Get the difference: */
   for (int i = 0; i < n; i++)
      dp[i] = p1[i] - p2[i];
   
}

/* Perform a SAXPY operation (v = a*x + y) in n dimensions: */
void math::saxpy(double* v, double a, const double* x, const double* y, int n) {
   
   /* Get the SAXPY vector: */
   for (int i = 0; i < n; i++)
      v[i] = a*x[i] + y[i];
   
}

/* Get the surface leakage factor for two centroids and a normal in 3 dimensions: */
double math::surface_leakage_factor(const double* p1, const double* p2, const double* n) {
   
   /* Get the difference between the two centroids: */
   double dp[3];
   subtract(dp, p2, p1, 3);
   
   /* Get the surface leakage factor w = (p2-p1)*n / |p2-p1|^2: */
   double w = math::dot_product(dp, n, 3) / math::l2_norm_2(dp, 3);
   
   return w;
   
}

double math::determinant(const std::vector<std::vector<double>> vect) {
    if(vect.size() != vect[0].size()) {
        throw std::runtime_error("Matrix is not quadratic");
    } 
    int dimension = vect.size();

    if(dimension == 0) {
        return 1;
    }

    if(dimension == 1) {
        return vect[0][0];
    }

    //Formula for 2x2-matrix
    if(dimension == 2) {
        return vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
    }

    double result = 0;
    int sign = 1;
    for(int i = 0; i < dimension; i++) {

        //Submatrix
        std::vector<std::vector<double>> subVect(dimension - 1, std::vector<double> (dimension - 1));
        for(int m = 1; m < dimension; m++) {
            int z = 0;
            for(int n = 0; n < dimension; n++) {
                if(n != i) {
                    subVect[m-1][z] = vect[m][n];
                    z++;
                }
            }
        }

        //recursive call
        result = result + sign * vect[0][i] * determinant(subVect);
        sign = -sign;
    }

    return result;
}

std::vector<std::vector<double>> math::transpose(const std::vector<std::vector<double>> matrix1) {

    //Transpose-matrix: height = width(matrix), width = height(matrix)
    std::vector<std::vector<double>> solution(matrix1[0].size(), std::vector<double> (matrix1.size()));

    //Filling solution-matrix
    for(size_t i = 0; i < matrix1.size(); i++) {
        for(size_t j = 0; j < matrix1[0].size(); j++) {
            solution[j][i] = matrix1[i][j];
        }
    }
    return solution;
}

std::vector<std::vector<double>> math::cofactor(const std::vector<std::vector<double>> vect) {
    if(vect.size() != vect[0].size()) {
        throw std::runtime_error("Matrix is not quadratic");
    } 

    std::vector<std::vector<double>> solution(vect.size(), std::vector<double> (vect.size()));
    std::vector<std::vector<double>> subVect(vect.size() - 1, std::vector<double> (vect.size() - 1));

    for(std::size_t i = 0; i < vect.size(); i++) {
        for(std::size_t j = 0; j < vect[0].size(); j++) {

            int p = 0;
            for(size_t x = 0; x < vect.size(); x++) {
                if(x == i) {
                    continue;
                }
                int q = 0;

                for(size_t y = 0; y < vect.size(); y++) {
                    if(y == j) {
                        continue;
                    }

                    subVect[p][q] = vect[x][y];
                    q++;
                }
                p++;
            }
            solution[i][j] = pow(-1, i + j) * determinant(subVect);
        }
    }
    return solution;
}

std::vector<std::vector<double>> math::inverse(const std::vector<std::vector<double>> vect) {
    if(determinant(vect) == 0) {
        throw std::runtime_error("Determinant is 0");
    }

    double d = 1.0/determinant(vect);
    std::vector<std::vector<double>> solution(vect.size(), std::vector<double> (vect.size()));

    for(size_t i = 0; i < vect.size(); i++) {
        for(size_t j = 0; j < vect.size(); j++) {
            solution[i][j] = vect[i][j]; 
        }
    }

    solution = transpose(cofactor(solution));

    for(size_t i = 0; i < vect.size(); i++) {
        for(size_t j = 0; j < vect.size(); j++) {
            solution[i][j] *= d;
        }
    }

    return solution;
}
