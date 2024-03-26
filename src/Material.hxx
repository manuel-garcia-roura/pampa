#pragma once

#include "utils.hxx"

/* The Material class: */
class Material {
   
   public:
      
      /* Number of energy groups: */
      int num_groups = -1;
      
      /* Cross sections: */
      Array1D<double> sigma_total;
      Array1D<double> nu_sigma_fission;
      Array2D<double> sigma_scattering;
      
      /* Diffusion coefficients: */
      Array1D<double> diffusion_coefficient;
      
      /* Fission spectrum: */
      Array1D<double> chi;
      
      /* Number of delayed-neutron precursor groups: */
      int num_precursor_groups;
      
      /* Precursor decay constants: */
      Array1D<double> lambda;
      
      /* Precursor fractions: */
      Array1D<double> beta;
      
      /* Thermal properties: */
      double rho = -1.0, cp = -1.0, k = -1.0;
      
      /* The Material constructor: */
      Material() {}
      
      /* The Material destructor: */
      ~Material() {}
      
      /* Read the material from a plain-text input file: */
      int WARN_UNUSED read(const std::string& filename);
   
};
