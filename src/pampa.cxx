#include "Pampa.hxx"

extern "C" {

/* Calculation driver: */
Pampa pampa;

/* Initialize: */
int initialize(int argc, char* argv[], double** dt, int* ndt) {
   
   Array1D<double> dt0;
   PAMPA_CALL(pampa.initialize(argc, argv, dt0), "unable to initialize the calculation");
   
   *ndt = dt0.size();
   *dt = new double[*ndt];
   for (int i = 0; i < *ndt; i++)
      (*dt)[i] = dt0(i);
   
   return 0;
   
}

/* Solve: */
int solve(int n, double dt, double t) {
   
   PAMPA_CALL(pampa.solve(n, dt, t), "unable to get the solution");
   
   return 0;
   
}

/* Finalize: */
int finalize(double** dt) {
   
   PAMPA_CALL(pampa.finalize(), "unable to finalize the calculation");
   
   delete[] *dt;
   
   return 0;
   
}

}
