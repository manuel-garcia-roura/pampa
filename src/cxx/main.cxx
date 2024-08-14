#include "Pampa.hxx"

/* The main function: */
int main(int argc, char* argv[]) {
   
   /* Calculation driver: */
   Pampa pampa;
   
   /* Time discretization: */
   Array1D<double> dt;
   
   /* Initialize the calculation: */
   PAMPA_CALL(pampa.initialize(argc, argv, dt), "unable to initialize the calculation");
   
   /* Run the time-stepping loop: */
   double t = 0.0;
   for (int n = 0; n < dt.size()+1; n++) {
      double dtn = (n == 0) ? 0.0 : dt(n-1);
      t += dtn;
      PAMPA_CALL(pampa.solve(n, dtn, t), "unable to get the solution");
   }
   
   /* Finalize the calculation: */
   PAMPA_CALL(pampa.finalize(), "unable to finalize the calculation");
   
   return 0;
   
}
