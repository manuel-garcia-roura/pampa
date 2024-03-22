#pragma once

#include <cmath>

#include "math.hxx"
#include "utils.hxx"

#define TOL 1.0e-6

/* The AngularQuadratureSet class: */
class AngularQuadratureSet {
   
   private:
      
      /* Order (N) of the SN method: */
      int order = -1;
      
      /* Number of discrete directions: */
      int num_directions = -1;
      
      /* Discrete directions: */
      Array2D<double> directions;
      
      /* Weights: */
      Array1D<double> weights;
      
      /* Cartesian axes to reflect directions: */
      Array2D<double> axes;
      
      /* Reflected directions with respect to the (+/-)x, (+/-)y and (+/-)z normals: */
      Array2D<int> reflected_directions;
   
   public:
      
      /* The AngularQuadratureSet default constructor: */
      AngularQuadratureSet() {}
      
      /* The AngularQuadratureSet constructor: */
      AngularQuadratureSet(int order) : order(order) {}
      
      /* The AngularQuadratureSet destructor: */
      ~AngularQuadratureSet() {}
      
      /* Get the number of discrete directions: */
      int getNumDirections() const {return num_directions;}
      
      /* Get the discrete directions: */
      const Array2D<double>& getDirections() const {return directions;}
      
      /* Get the weights: */
      const Array1D<double>& getWeights() const {return weights;}
      
      /* Get the Cartesian axes to reflect directions: */
      const Array2D<double>& getAxes() const {return axes;}
      
      /* Get the reflected directions with respect to the (+/-)x, (+/-)y and (+/-)z normals: */
      const Array2D<int>& getReflectedDirections() const {return reflected_directions;}
      
      /* Build the angular quadrature set: */
      int WARN_UNUSED build();
   
};
