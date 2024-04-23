#pragma once

#include <string>
#include <cmath>

#include "utils.hxx"

/* The TH::Properties enum: */
namespace TH {
   enum Properties {CONSTANT, GRAPHITE_H_451, GRAPHITE_MATRIX_A3_27};
}

/* The Material class: */
class Material {
   
   public:
      
      /* Material name: */
      const std::string name = "";
      
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
      
      /* Number of delayed-neutron precursor groups: */
      int num_precursor_groups = -1;
      
      /* Precursor decay constants: */
      Array1D<double> lambda;
      
      /* Precursor fractions: */
      Array1D<double> beta;
      double beta_total = 0.0;
      
      /* Thermal properties: */
      TH::Properties thermal_properties = TH::CONSTANT;
      double k0 = -1.0, rho0 = -1.0, cp0 = -1.0;
      
      /* The Material constructor: */
      Material(const std::string& name) : name(name) {}
      
      /* The Material destructor: */
      ~Material() {}
      
      /* Read the material from a plain-text input file: */
      int WARN_UNUSED read(const std::string& filename);
      
      /* Read the material from a plain-text input file: */
      int WARN_UNUSED read(std::ifstream& file);
      
      /* Get the thermal conductivity: */
      double k(double T) const;
      
      /* Get the density: */
      double rho(double T) const;
      
      /* Get the specific heat capacity: */
      double cp(double T) const;
   
};
