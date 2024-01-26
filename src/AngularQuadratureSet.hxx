#pragma once

#include <cmath>
#include <vector>
#include <array>

#include "math.hxx"
#include "utils.hxx"

#define TOL 1.0e-6

/* The AngularQuadratureSet class: */
class AngularQuadratureSet {
   
   /* Order (N) of the SN method: */
   int order = -1;
   
   /* Number of discrete directions: */
   int num_directions = -1;
   
   /* Discrete directions: */
   std::vector<std::vector<double>> directions;
   
   /* Weights: */
   std::vector<double> weights;
   
   /* Reflected directions with respect to the (+/-)x, (+/-)y and (+/-)z normals: */
   std::vector<std::array<int, 3>> reflected_directions;
   
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
      const std::vector<std::vector<double>>& getDirections() const {return directions;}
      
      /* Get the weights: */
      const std::vector<double>& getWeights() const {return weights;}
      
      /* Get the reflected directions with respect to the (+/-)x, (+/-)y and (+/-)z normals: */
      const std::vector<std::array<int, 3>>& getReflectedDirections() const 
         {return reflected_directions;}
      
      /* Build the angular quadrature set: */
      int build();
   
};
