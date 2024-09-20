#include <stdio.h>

#include "Driver.hxx"

extern "C" {

/* Calculation driver: */
Driver pampa;

/* Initialize the calculation: */
void pampa_initialize(int argc, char* argv[], double** dt, int* ndt, int* error) {
   
   /* Copy the command line arguments to a regular C array: */
   /* Note: this is needed to avoid MPI errors when this function is called from Fortran90. */
   char** argv_tmp = new char*[argc];
   for (int i = 0; i < argc; i++)
      argv_tmp[i] = strdup(argv[i]);
   
   /* Initialize the calculation driver and get the time discretization: */
   Array1D<double> dt0;
   *error = pampa.initialize(argc, argv_tmp, dt0);
   if (*error > 0) {printf("Error in pampa.initialize().\n"); return;}
   
   /* Free the temporary arguments array: */
   for (int i = 0; i < argc; i++)
      free(argv_tmp[i]);
   delete[] argv_tmp;
   
   /* Allocate and copy the time-discretization array: */
   *ndt = dt0.size();
   *dt = new double[*ndt];
   for (int i = 0; i < *ndt; i++)
      (*dt)[i] = dt0(i);
   
}

/* Initialize a steady-state calculation: */
void pampa_initialize_steady_state(int argc, char* argv[], int* error) {
   
   /* Copy the command line arguments to a regular C array: */
   /* Note: this is needed to avoid MPI errors when this function is called from Fortran90. */
   char** argv_tmp = new char*[argc];
   for (int i = 0; i < argc; i++)
      argv_tmp[i] = strdup(argv[i]);
   
   /* Initialize the calculation driver and get the time discretization: */
   Array1D<double> dt0;
   *error = pampa.initialize(argc, argv_tmp, dt0);
   if (*error > 0) {printf("Error in pampa.initialize().\n"); return;}
   
   /* Free the temporary arguments array: */
   for (int i = 0; i < argc; i++)
      free(argv_tmp[i]);
   delete[] argv_tmp;
   
}

/* Get the solution: */
void pampa_solve(int n, double dt, double t, int* error) {
   
   /* Get the solution from the calculation driver: */
   *error = pampa.solve(n, dt, t);
   if (*error > 0) {printf("Error in pampa.solve().\n"); return;}
   
}

/* Get the steady-state solution: */
void pampa_solve_steady_state(int* error) {
   
   /* Get the solution: */
   pampa_solve(0, 0.0, 0.0, error);
   if (*error > 0) {printf("Error in pampa_solve().\n"); return;}
   
}

/* Finalize the calculation: */
void pampa_finalize(double** dt, int* error) {
   
   /* Finalize the calculation driver: */
   *error = pampa.finalize();
   if (*error > 0) {printf("Error in pampa.finalize().\n"); return;}
   
   /* Free the time-discretization array: */
   delete[] *dt;
   
}

/* Finalize a steady-state calculation: */
void pampa_finalize_steady_state(int* error) {
   
   /* Finalize the calculation driver: */
   *error = pampa.finalize();
   if (*error > 0) {printf("Error in pampa.finalize().\n"); return;}
   
}

/* Get the values for a given field: */
void pampa_get_field(double* v, const char name[], int* error) {
   
   /* Get the field from the calculation driver: */
   *error = pampa.getField(v, std::string(name));
   if (*error > 0) {printf("Error in pampa.getField().\n"); return;}
   
}

/* Set the values for a given field: */
void pampa_set_field(const double* v, const char name[], int* error) {
   
   /* Set the field to the calculation driver: */
   *error = pampa.setField(v, std::string(name));
   if (*error > 0) {printf("Error in pampa.setField().\n"); return;}
   
}

}
