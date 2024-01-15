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
      
      /* The AngularQuadratureSet constructor: */
      AngularQuadratureSet() {}
      
      /* The AngularQuadratureSet destructor: */
      ~AngularQuadratureSet() {}
      
      /* Set the angular quadrature order: */
      void setOrder(int order) {this->order = order;}
      
      /* Build the angular quadrature set: */
      int build();
   
};
