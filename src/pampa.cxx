#include <stdio.h>

#include "Pampa.hxx"

extern "C" {

/* Calculation driver: */
Pampa pampa;

/* Initialize the calculation: */
int pampa_initialize(int argc, char* argv[], double** dt, int* ndt, int* error) {
   
   /* Copy the command line arguments to a regular C array: */
   /* Note: this is needed to avoid MPI errors when this function is called from Fortran90. */
   char** argv_tmp = new char*[argc];
   for (int i = 0; i < argc; i++)
      argv_tmp[i] = strdup(argv[i]);
   
   /* Initialize the calculation driver and get the time discretization: */
   Array1D<double> dt0;
   *error = pampa.initialize(argc, argv_tmp, dt0);
   if (*error > 0) {printf("Error in pampa.initialize().\n"); return 1;}
   
   /* Free the temporary arguments array: */
   for (int i = 0; i < argc; i++)
      free(argv_tmp[i]);
   free(argv_tmp);
   
   /* Allocate and copy the time-discretization array: */
   *ndt = dt0.size();
   *dt = new double[*ndt];
   for (int i = 0; i < *ndt; i++)
      (*dt)[i] = dt0(i);
   
   return 0;
   
}

/* Get the solution: */
int pampa_solve(int n, double dt, double t, int* error) {
   
   /* Get the solution from the calculation driver: */
   *error = pampa.solve(n, dt, t);
   if (*error > 0) {printf("Error in pampa.solve().\n"); return 1;}
   
   return 0;
   
}

/* Finalize the calculation: */
int pampa_finalize(double** dt, int* error) {
   
   /* Finalize the calculation driver: */
   *error = pampa.finalize();
   if (*error > 0) {printf("Error in pampa.finalize().\n"); return 1;}
   
   /* Free the time-discretization array: */
   delete[] *dt;
   
   return 0;
   
}

}
