#pragma once

#include "utils.hxx"

/* The ThermalProperties class: */
class ThermalProperties {
   
   public:
      
      /* Switches to manage the thermal properties: */
      bool constant = true;
      
      /* The ThermalProperties constructor: */
      ThermalProperties() {}
      
      /* The ThermalProperties destructor: */
      virtual ~ThermalProperties() {}
      
      /* Get the thermal conductivity: */
      virtual double k(double T) const {PAMPA_CHECK_VIRTUAL}
      
      /* Get the density: */
      virtual double rho(double T) const {PAMPA_CHECK_VIRTUAL}
      
      /* Get the specific heat capacity: */
      virtual double cp(double T) const {PAMPA_CHECK_VIRTUAL}
   
};
