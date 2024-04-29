#pragma once

#include "utils.hxx"

/* The NuclearData class: */
class NuclearData {
   
   public:
      
      /* Number of energy groups: */
      int num_energy_groups = -1;
      
      /* Cross sections: */
      Array1D<double> sigma_total, nu_sigma_fission, kappa_sigma_fission, sigma_transport;
      Array2D<double> sigma_scattering;
      
      /* Diffusion coefficients: */
      Array1D<double> diffusion_coefficient;
      
      /* Fission spectrum: */
      Array1D<double> chi_prompt, chi_delayed, chi_eff;
      
      /* Neutron velocity: */
      Array1D<double> velocity;
      
      /* Default neutron yield and fission energy: */
      double nu = 2.4355, kappa = 3.2e-11;
      
      /* The NuclearData constructor: */
      NuclearData() {}
      
      /* The NuclearData destructor: */
      ~NuclearData() {}
      
      /* Read the nuclear data from a plain-text input file: */
      int WARN_UNUSED read(std::ifstream& file);
   
};
