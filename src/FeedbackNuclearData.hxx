#pragma once

#include <cmath>

#include "NuclearData.hxx"
#include "ConstantNuclearData.hxx"

#define TOL 1.0e-6

/* The FeedbackNuclearData class: */
class FeedbackNuclearData : public NuclearData {
   
   public:
      
      /* Interpolation parameters: */
      int i1 = -1, i2 = -1;
      double T0 = -1.0, f = -1.0;
      
      /* Reference temperature: */
      Array1D<double> Tref;
      
      /* Constant nuclear data instances: */
      Array1D<ConstantNuclearData*> nuclear_data;
      
      /* The FeedbackNuclearData constructor: */
      FeedbackNuclearData() {}
      
      /* The FeedbackNuclearData destructor: */
      ~FeedbackNuclearData() {
         
         /* Free all the nuclear data instances: */
         for (int i = 0; i < nuclear_data.size(); i++)
            utils::free(&nuclear_data(i));
         
      }
      
      /* Read the nuclear data from a plain-text input file: */
      int WARN_UNUSED read(std::ifstream& file);
      
      /* Check the nuclear data after reading it: */
      int WARN_UNUSED check(double beta_total) {
         
         /* Check the nuclear data for all the nuclear data instances: */
         for (int i = 0; i < nuclear_data.size(); i++)
            PAMPA_CHECK(nuclear_data(i)->check(beta_total), "wrong nuclear data");
         
         return 0;
         
      }
      
      /* Check the nuclear data to use it in a solver: */
      int WARN_UNUSED check(int num_energy_groups, bool diffusion, bool transient) const {
         
         /* Check the nuclear data for all the nuclear data instances: */
         for (int i = 0; i < nuclear_data.size(); i++)
            PAMPA_CHECK(nuclear_data(i)->check(num_energy_groups, diffusion, transient), 
               "wrong nuclear data");
         
         return 0;
         
      }
      
      /* Get the interpolation parameters for a given temperature: */
      void interpolate(double T) {
         
         /* Check if the temperature has changed: */
         if (fabs(T-T0) > TOL) {
            
            /* Use the first value: */
            if (T < Tref(0)) {
               i1 = 0;
               i2 = i1;
               f = 0.0;
               return;
            }
            
            /* Use the last value: */
            if (T > Tref.back()) {
               i1 = Tref.size() - 1;
               i2 = i1;
               f = 0.0;
               return;
            }
            
            /* Use a pair of values: */
            i2 = Tref.lowerBound(T);
            i1 = i2 - 1;
            f = (T-Tref(i1)) / (Tref(i2)-Tref(i1));
            
         }
         
      }
      
      /* Get a total cross section: */
      double sigmaTotal(int g, double T) 
         {interpolate(T); return (1.0-f)*nuclear_data(i1)->sigma_total(g) + 
            f*nuclear_data(i2)->sigma_total(g);}
      
      /* Get a nu-fission cross section: */
      double sigmaNuFission(int g, double T) 
         {interpolate(T); return (1.0-f)*nuclear_data(i1)->nu_sigma_fission(g) + 
            f*nuclear_data(i2)->nu_sigma_fission(g);}
      
      /* Get a kappa-fission cross section: */
      double sigmaKappaFission(int g, double T) 
         {interpolate(T); return (1.0-f)*nuclear_data(i1)->kappa_sigma_fission(g) + 
            f*nuclear_data(i2)->kappa_sigma_fission(g);}
      
      /* Get a transport cross section: */
      double sigmaTransport(int g, double T) 
         {interpolate(T); return (1.0-f)*nuclear_data(i1)->sigma_transport(g) + 
            f*nuclear_data(i2)->sigma_transport(g);}
      
      /* Get a scattering cross section: */
      double sigmaScattering(int g, int g2, double T) 
         {interpolate(T); return (1.0-f)*nuclear_data(i1)->sigma_scattering(g, g2) + 
            f*nuclear_data(i2)->sigma_scattering(g, g2);}
      
      /* Get a diffusion coefficient: */
      double diffusionCoefficient(int g, double T) 
         {interpolate(T); return (1.0-f)*nuclear_data(i1)->diffusion_coefficient(g) + 
            f*nuclear_data(i2)->diffusion_coefficient(g);}
      
      /* Get a prompt-fission spectrum point: */
      double chiPrompt(int g, double T) 
         {interpolate(T); return (1.0-f)*nuclear_data(i1)->chi_prompt(g) + 
            f*nuclear_data(i2)->chi_prompt(g);}
      
      /* Get a delayed-fission spectrum point: */
      double chiDelayed(int g, double T) 
         {interpolate(T); return (1.0-f)*nuclear_data(i1)->chi_delayed(g) + 
            f*nuclear_data(i2)->chi_delayed(g);}
      
      /* Get an effective-fission spectrum point: */
      double chiEffective(int g, double T) 
         {interpolate(T); return (1.0-f)*nuclear_data(i1)->chi_effective(g) + 
            f*nuclear_data(i2)->chi_effective(g);}
      
      /* Get a neutron velocity: */
      double neutronVelocity(int g, double T) 
         {interpolate(T); return (1.0-f)*nuclear_data(i1)->velocity(g) + 
            f*nuclear_data(i2)->velocity(g);}
   
};
