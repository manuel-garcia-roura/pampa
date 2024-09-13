#pragma once

#include "ThermalProperties.hxx"

/* The GraphiteMatrixProperties class: */
class GraphiteMatrixProperties : public ThermalProperties {
   
   private:
      
      /* Coefficients for the thermal conductivity [Gontard and Nabielek (1990)]: */
      double k100 = 47.4, alpha = 9.7556e-4, delta = -6.0360e-4;
      
      /* Constant density [NEA (2018)]: */
      double rho0 = 1.700e-3;
      
      /* Coefficients for the specific heat capacity [Butland and Maddison (1973)]: */
      double cp0 = 0.54212, cp1 = 2.42667e-6, cpm1 = 90.2725, cpm2 = 43449.3, cpm3 = 1.59309e7, 
         cpm4 = 1.43688e9;
      
   
   public:
      
      /* The GraphiteMatrixProperties constructor: */
      GraphiteMatrixProperties() {constant = false;}
      
      /* The GraphiteMatrixProperties destructor: */
      ~GraphiteMatrixProperties() {}
      
      /* Get the thermal conductivity for grade A3-27 graphite heat-treated at 1800°C: */
      double k(double T) const {
         
         /* Return the thermal conductivity in W/(cm^2*K) [Gontard and Nabielek (1990)]: */
         return 1.0e-2 * k100 * (1.0-alpha*(T-100.0)*std::exp(delta*T));
         
      }
      
      /* Get the density for grade A3-27 graphite heat-treated at 1800°C: */
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
