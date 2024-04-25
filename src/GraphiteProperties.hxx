#pragma once

#include "ThermalProperties.hxx"

/* The GraphiteProperties class: */
class GraphiteProperties : public ThermalProperties {
   
   private:
      
      /* Coefficients for the thermal conductivity [NEA (2018)]: */
      double k0 = 1.69214e2, k1 = 1.24890e-1, k2 = 3.28248e-5;
      
      /* Constant density [NEA (2018)]: */
      double rho0 = 1.850e-3;
      
      /* Coefficients for the specific heat capacity [Butland and Maddison (1973)]: */
      double cp0 = 0.54212, cp1 = 2.42667e-6, cpm1 = 90.2725, cpm2 = 43449.3, cpm3 = 1.59309e7, 
         cpm4 = 1.43688e9;
      
   
   public:
      
      /* The GraphiteProperties constructor: */
      GraphiteProperties() : ThermalProperties(false) {}
      
      /* The GraphiteProperties destructor: */
      ~GraphiteProperties() {}
      
      /* Get the thermal conductivity for H-451 graphite: */
      double k(double T) const {
         
         /* Check the validity range: */
         if (T < 500.0) T = 500.0;
         if (T > 1800.0) T = 1800.0;
         
         /* Return the thermal conductivity in W/(cm^2*K) [NEA (2018)]: */
         return 1.0e-2 * (k0-k1*T+k2*T*T);
         
      }
      
      /* Get the density for H-451 graphite: */
      double rho(double T) const {
         
         /* Return the density in kg/cm^3 [NEA (2018)]: */
         return rho0;
         
      }
      
      /* Get the specific heat capacity for H-451 graphite: */
      double cp(double T) const {
         
         /* Check the validity range: */
         if (T < 200.0) T = 200.0;
         if (T > 3500.0) T = 3500.0;
         
         /* Return the specific heat capacity in J/(kg*K) [Butland and Maddison (1973)]: */
         return 4184.0 * (cp0-cp1*T-cpm1/T-cpm2/(T*T)+cpm3/(T*T*T)-cpm4/(T*T*T*T));
         
      }
   
};
