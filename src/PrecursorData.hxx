#pragma once

#include "input.hxx"
#include "utils.hxx"

/* The PrecursorData class: */
class PrecursorData {
   
   public:
      
      /* Number of delayed-neutron precursor groups: */
      int num_precursor_groups = -1;
      
      /* Precursor decay constants: */
      Array1D<double> lambda;
      
      /* Precursor fractions: */
      Array1D<double> beta;
      double beta_total = 0.0;
      
      /* The PrecursorData constructor: */
      PrecursorData() {}
      
      /* The PrecursorData destructor: */
      ~PrecursorData() {}
      
      /* Read the precursor data from a plain-text input file: */
      int WARN_UNUSED read(std::ifstream& file);
      
      /* Check the precursor data to use it in a solver: */
      int WARN_UNUSED check(int num_precursor_groups) const;
   
};
