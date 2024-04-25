#pragma once

#include "ThermalProperties.hxx"

/* The ConstantProperties class: */
class ConstantProperties : public ThermalProperties {
   
   private:
      
      /* Constant thermal conductivity: */
      double k0;
      
      /* Constant density: */
      double rho0;
      
      /* Constant specific heat capacity: */
      double cp0;
   
   public:
      
      /* The ConstantProperties constructor: */
      ConstantProperties(double k0, double rho0, double cp0) : ThermalProperties(true), k0(k0), 
         rho0(rho0), cp0(cp0) {}
      
      /* The ConstantProperties destructor: */
      ~ConstantProperties() {}
      
      /* Get the thermal conductivity: */
      double k(double T) const {return k0;}
      
      /* Get the density: */
      double rho(double T) const {return rho0;}
      
      /* Get the specific heat capacity: */
      double cp(double T) const {return cp0;}
   
};
