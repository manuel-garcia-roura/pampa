#pragma once

#include "NuclearData.hxx"

/* The ConstantNuclearData class: */
class ConstantNuclearData : public NuclearData {
   
   public:
      
      /* Number of energy groups: */
      int num_energy_groups = -1;
      
      /* Cross sections: */
      Array1D<double> sigma_total, nu_sigma_fission, kappa_sigma_fission, sigma_transport;
      Array2D<double> sigma_scattering;
      
      /* Diffusion coefficients: */
      Array1D<double> diffusion_coefficient;
      
      /* Fission spectrum: */
      Array1D<double> chi_prompt, chi_delayed, chi_effective;
      
      /* Neutron velocity: */
      Array1D<double> velocity;
      
      /* Default neutron yield and fission energy: */
      double nu = 2.4355, kappa = 3.2e-11;
      
      /* The ConstantNuclearData constructor: */
      ConstantNuclearData() {}
      
      /* The ConstantNuclearData destructor: */
      ~ConstantNuclearData() {}
      
      /* Read the nuclear data from a plain-text input file: */
      int WARN_UNUSED read(std::ifstream& file);
      
      /* Check the nuclear data after reading it: */
      int WARN_UNUSED check(double beta_total);
      
      /* Check the nuclear data to use it in a solver: */
      int WARN_UNUSED check(int num_energy_groups, bool diffusion, bool transient) const;
      
      /* Get a total cross section: */
      double sigmaTotal(int g, double T) {return sigma_total(g);}
      
      /* Get a nu-fission cross section: */
      double sigmaNuFission(int g, double T) {return nu_sigma_fission(g);}
      
      /* Get a kappa-fission cross section: */
      double sigmaKappaFission(int g, double T) {return kappa_sigma_fission(g);}
      
      /* Get a transport cross section: */
      double sigmaTransport(int g, double T) {return sigma_transport(g);}
      
      /* Get a scattering cross section: */
      double sigmaScattering(int g, int g2, double T) {return sigma_scattering(g, g2);}
      
      /* Get a diffusion coefficient: */
      double diffusionCoefficient(int g, double T) {return diffusion_coefficient(g);}
      
      /* Get a prompt-fission spectrum point: */
      double chiPrompt(int g, double T) {return chi_prompt(g);}
      
      /* Get a delayed-fission spectrum point: */
      double chiDelayed(int g, double T) {return chi_delayed(g);}
      
      /* Get an effective-fission spectrum point: */
      double chiEffective(int g, double T) {return chi_effective(g);}
      
      /* Get a neutron velocity: */
      double neutronVelocity(int g, double T) {return velocity(g);}
   
};
