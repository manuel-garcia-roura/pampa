#include "pampa.hxx"

/* The main function: */
int main(int argc, char* argv[]) {
   
   /* Error signal: */
   int error;
   
   /* Time discretization: */
   double* dt;
   int ndt;
   
   /* Initialize: */
   pampa_initialize(argc, argv, &dt, &ndt, &error);
   
   /* Run the time-stepping loop: */
   double t = 0.0;
   for (int n = 0; n < ndt+1; n++) {
      double dtn = (n == 0) ? 0.0 : dt[n-1];
      t += dtn;
      pampa_solve(n, dtn, t, &error);
   }
   
   /* Finalize: */
   pampa_finalize(&dt, &error);
   
   return 0;
   
}
