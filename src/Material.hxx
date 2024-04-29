#pragma once

#include <string>

#include "NuclearData.hxx"
#include "PrecursorData.hxx"
#include "ThermalProperties.hxx"
#include "ConstantProperties.hxx"
#include "GraphiteProperties.hxx"
#include "GraphiteMatrixProperties.hxx"
#include "utils.hxx"

/* The Material class: */
class Material {
   
   private:
      
      /* Nuclear data: */
      NuclearData* nuclear_data = nullptr;
      
      /* Precursor data: */
      PrecursorData* precursor_data = nullptr;
      
      /* Thermal properties: */
      ThermalProperties* thermal_properties = nullptr;
   
   public:
      
      /* Material name: */
      const std::string name;
      
      /* The Material constructor: */
      Material(const std::string& name) : name(name) {}
      
      /* The Material destructor: */
      ~Material() {
         utils::free(&nuclear_data);
         utils::free(&precursor_data);
         utils::free(&thermal_properties);
      }
      
      /* Read the material from a plain-text input file: */
      int WARN_UNUSED read(const std::string& filename);
      
      /* Read the material from a plain-text input file: */
      int WARN_UNUSED read(std::ifstream& file);
      
      /* Check if the material has nuclear data: */
      bool hasNuclearData() const {return nuclear_data != nullptr;}
      
      /* Check if the material has precursor data: */
      bool hasPrecursorData() const {return precursor_data != nullptr;}
      
      /* Check if the material has thermal properties: */
      bool hasThermalProperties() const {return thermal_properties != nullptr;}
      
      /* Check if the material has constant thermal properties: */
      bool hasConstantThermalProperties() const {return thermal_properties->constant;}
      
      /* Check the nuclear data: */
      int WARN_UNUSED checkNuclearData(int num_energy_groups, bool diffusion, bool transient) const;
      
      /* Check the precursor data: */
      int WARN_UNUSED checkPrecursorData(int num_precursor_groups) const;
      
      /* Get a total cross section: */
      double sigmaTotal(int g) const {return nuclear_data->sigma_total(g);}
      
      /* Get a nu-fission cross section: */
      double sigmaNuFission(int g) const {return nuclear_data->nu_sigma_fission(g);}
      
      /* Get a kappa-fission cross section: */
      double sigmaKappaFission(int g) const {return nuclear_data->kappa_sigma_fission(g);}
      
      /* Get a transport cross section: */
      double sigmaTransport(int g) const {return nuclear_data->sigma_transport(g);}
      
      /* Get a scattering cross section: */
      double sigmaScattering(int g, int g2) const {return nuclear_data->sigma_scattering(g, g2);}
      
      /* Get a diffusion coefficient: */
      double diffusionCoefficient(int g) const {return nuclear_data->diffusion_coefficient(g);}
      
      /* Get a prompt-fission spectrum point: */
      double chiPrompt(int g) const {return nuclear_data->chi_prompt(g);}
      
      /* Get a delayed-fission spectrum point: */
      double chiDelayed(int g) const {return nuclear_data->chi_delayed(g);}
      
      /* Get an effective-fission spectrum point: */
      double chiEffective(int g) const {return nuclear_data->chi_eff(g);}
      
      /* Get a neutron velocity: */
      double neutronVelocity(int g) const {return nuclear_data->velocity(g);}
      
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
