#include "pampa.hxx"

/* Initialize: */
int initialize(int argc, char* argv[], Array1D<double>& dt) 
   {PAMPA_CALL(pampa.initialize(argc, argv, dt), "unable to initialize the calculation"); return 0;}

/* Solve: */
int solve(int n, double dt, double t) 
   {PAMPA_CALL(pampa.solve(n, dt, t), "unable to get the solution"); return 0;}

/* Finalize: */
int finalize() 
   {PAMPA_CALL(pampa.finalize(), "unable to finalize the calculation"); return 0;}
