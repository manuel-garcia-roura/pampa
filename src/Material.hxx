#pragma once

#include <string>

#include "PrecursorData.hxx"
#include "ThermalProperties.hxx"
#include "ConstantProperties.hxx"
#include "GraphiteProperties.hxx"
#include "GraphiteMatrixProperties.hxx"
#include "utils.hxx"

/* The Material class: */
class Material {
   
   private:
      
      /* Precursor data: */
      PrecursorData* precursor_data = nullptr;
      
      /* Thermal properties: */
      ThermalProperties* thermal_properties = nullptr;
   
   public:
      
      /* Material name: */
      const std::string name;
      
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
      
      /* The Material constructor: */
      Material(const std::string& name) : name(name) {}
      
      /* The Material destructor: */
      ~Material() {if (precursor_data != nullptr) delete precursor_data; if (thermal_properties != nullptr) delete thermal_properties;}
      
      /* Read the material from a plain-text input file: */
      int WARN_UNUSED read(const std::string& filename);
      
      /* Read the material from a plain-text input file: */
      int WARN_UNUSED read(std::ifstream& file);
      
      /* Check the precursor data: */
      int WARN_UNUSED checkPrecursorData(int num_precursor_groups) const;
      
      /* Check if the material has precursor data: */
      bool hasPrecursorData() const {return precursor_data != nullptr;}
      
      /* Check if the material has thermal properties: */
      bool hasThermalProperties() const {return thermal_properties != nullptr;}
      
      /* Check if the material has constant thermal properties: */
      bool hasConstantThermalProperties() const {return thermal_properties->constant;}
      
      /* Get a precursor decay constant: */
      double lambda(int g) const {return precursor_data->lambda(g);}
      
      /* Get a precursor fraction: */
      double beta(int g) const {return precursor_data->beta(g);}
      
      /* Get the total precursor fraction: */
      double beta() const {return hasPrecursorData() ? precursor_data->beta_total : 0.0;}
      
      /* Get the thermal conductivity: */
      double k(double T) const {return thermal_properties->k(T);}
      
      /* Get the density: */
      double rho(double T) const {return thermal_properties->rho(T);}
      
      /* Get the specific heat capacity: */
      double cp(double T) const {return thermal_properties->cp(T);}
   
};
