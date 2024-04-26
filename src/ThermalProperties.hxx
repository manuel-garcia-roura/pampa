#pragma once

#include "utils.hxx"

/* The ThermalProperties class: */
class ThermalProperties {
   
   public:
      
      /* Switch to determine if the thermal properties are constant or temperature-dependent: */
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
