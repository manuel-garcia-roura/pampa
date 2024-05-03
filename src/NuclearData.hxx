#pragma once

#include "utils.hxx"

/* The NuclearData class: */
class NuclearData {
   
   public:
      
      /* The NuclearData constructor: */
      NuclearData() {}
      
      /* The NuclearData destructor: */
      virtual ~NuclearData() {}
      
      /* Read the nuclear data from a plain-text input file: */
      virtual int WARN_UNUSED read(std::ifstream& file) {PAMPA_CHECK_VIRTUAL}
      
      /* Check the nuclear data after reading it: */
      virtual int WARN_UNUSED check(double beta_total) {PAMPA_CHECK_VIRTUAL}
      
      /* Check the nuclear data to use it in a solver: */
      virtual int WARN_UNUSED check(int num_energy_groups, bool diffusion, bool transient) const 
         {PAMPA_CHECK_VIRTUAL}
      
      /* Get a total cross section: */
      virtual double sigmaTotal(int g, double T) {PAMPA_CHECK_VIRTUAL}
      
      /* Get a nu-fission cross section: */
      virtual double sigmaNuFission(int g, double T) {PAMPA_CHECK_VIRTUAL}
      
      /* Get a kappa-fission cross section: */
      virtual double sigmaKappaFission(int g, double T) {PAMPA_CHECK_VIRTUAL}
      
      /* Get a transport cross section: */
      virtual double sigmaTransport(int g, double T) {PAMPA_CHECK_VIRTUAL}
      
      /* Get a scattering cross section: */
      virtual double sigmaScattering(int g, int g2, double T) {PAMPA_CHECK_VIRTUAL}
      
      /* Get a diffusion coefficient: */
      virtual double diffusionCoefficient(int g, double T) {PAMPA_CHECK_VIRTUAL}
      
      /* Get a prompt-fission spectrum point: */
      virtual double chiPrompt(int g, double T) {PAMPA_CHECK_VIRTUAL}
      
      /* Get a delayed-fission spectrum point: */
      virtual double chiDelayed(int g, double T) {PAMPA_CHECK_VIRTUAL}
      
      /* Get an effective-fission spectrum point: */
      virtual double chiEffective(int g, double T) {PAMPA_CHECK_VIRTUAL}
      
      /* Get a neutron velocity: */
      virtual double neutronVelocity(int g, double T) {PAMPA_CHECK_VIRTUAL}
   
};
