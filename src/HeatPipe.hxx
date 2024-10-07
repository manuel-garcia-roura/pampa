#pragma once

#include "utils.hxx"

/* The HeatPipe class: */
class HeatPipe {
   
   private:
      
      /* Number of heat pipes: */
      int n = 1;
      
      /* Condenser-side diameter, length and area: */
      double Dc = -1.0, Lc = -1.0, Ac = -1.0;
      
      /* Relaxation factor: */
      double w = -1.0;
      
      /* Condenser-side heat-transfer coefficient and temperature: */
      Function hc, Tc;
      
      /* Previous iteration used to relax the heat-pipe temperature: */
      double T0 = -1.0;
   
   public:
      
      /* The HeatPipe constructor: */
      HeatPipe(int n, double Dc, double Lc, double w, Function& hc, Function& Tc) : n(n), Dc(Dc), 
         Lc(Lc), Ac(n*M_PI*Dc*Lc), w(w), hc(hc), Tc(Tc) {}
      
      /* The HeatPipe destructor: */
      ~HeatPipe() {}
      
      /* Calculate the heat-pipe temperature for a given heat flow: */
      double calculateTemperature(double q, double t) {
         
         /* Get the heat-pipe temperature: */
         double T = Tc(t) + q/(Ac*hc(t));
         
         /* Relax the temperature: */
         if (T0 > 0.0) T = w*T + (1.0-w)*T0;
         T0 = T;
         
         return T;
         
      }
   
};
