#include "pampa.hxx"

/* The main function: */
int main(int argc, char* argv[]) {
   
   /* Time discretization: */
   double* dt;
   int ndt;
   
   /* Initialize: */
   initialize(argc, argv, &dt, &ndt);
   
   /* Run the time-stepping loop: */
   double t = 0.0;
   for (int n = 0; n < ndt+1; n++) {
      double dtn = (n == 0) ? 0.0 : dt[n-1];
      t += dtn;
      solve(n, dtn, t);
   }
   
   /* Finalize: */
   finalize(&dt);
   
   return 0;
   
}
