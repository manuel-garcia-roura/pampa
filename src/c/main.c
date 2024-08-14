#include <stdio.h>

#include "pampa.hxx"

/* The main function: */
int main(int argc, char* argv[]) {
   
   /* Error signal: */
   int error;
   
   /* Time discretization: */
   double* dt;
   int ndt;
   
   /* Initialize the calculation: */
   pampa_initialize(argc, argv, &dt, &ndt, &error);
   if (error > 0) {printf("Error in pampa_initialize().\n"); return 1;}
   
   /* Run the time-stepping loop: */
   double t = 0.0;
   for (int n = 0; n < ndt+1; n++) {
      double dtn = (n == 0) ? 0.0 : dt[n-1];
      t += dtn;
      pampa_solve(n, dtn, t, &error);
      if (error > 0) {printf("Error in pampa_solve().\n"); return 1;}
   }
   
   /* Finalize the calculation: */
   pampa_finalize(&dt, &error);
   if (error > 0) {printf("Error in pampa_finalize().\n"); return 1;}
   
   return 0;
   
}
