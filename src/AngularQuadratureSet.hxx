#pragma once

#include <cmath>
#include <vector>

#include "utils.hxx"

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
      
      /* Build the angular quadrature set: */
      int build();
   
};
