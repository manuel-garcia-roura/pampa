#pragma once

#include "utils.hxx"

/* The HeatPipe class: */
class HeatPipe {
   
   private:
      
      /* Boundary condition: */
      BoundaryCondition* bc = nullptr;
      
      /* Number of heat pipes: */
      int n = 1;
      
      /* Condenser-side diameter, length and area: */
      double Dc = -1.0, Lc = -1.0, Ac = -1.0;
      
      /* Relaxation factor: */
      double w = -1.0;
      
      /* Condenser-side heat-transfer coefficient and temperature: */
      Function hc, Tc;
      
      /* Total heat source: */
      double q = -1.0;
   
   public:
      
      /* The HeatPipe constructor: */
      HeatPipe(BoundaryCondition* bc, int n, double Dc, double Lc, double w, Function& hc, 
         Function& Tc, double q) : bc(bc), n(n), Dc(Dc), Lc(Lc), Ac(n*M_PI*Dc*Lc), w(w), hc(hc), 
         Tc(Tc), q(q) {}
      
      /* The HeatPipe destructor: */
      ~HeatPipe() {}
      
      /* Set the heat source: */
      void setHeatSource(double q) {this->q = q;}
      
      /* Calculate the heat-pipe temperature: */
      void calculateTemperature(double t) {
         
         /* Get the heat-pipe temperature: */
         double T = Tc(t) + q/(Ac*hc(t));
         
         /* Relax the temperature and set it in the boundary condition: */
         T = w*T + (1.0-w)*bc->f(1)(0.0);
         bc->f(1) = Function(T);
         
      }
   
};
