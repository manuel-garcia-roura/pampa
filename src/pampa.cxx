#include "Pampa.hxx"

extern "C" {

/* Calculation driver: */
Pampa pampa;

/* Initialize: */
int pampa_initialize(int argc, char* argv[], double** dt, int* ndt, int* error) {
   
    int i;
    char** new_argv;

    /* Allocate memory for new_argv array */
    new_argv = new char*[argc];
    if (new_argv == NULL) {
        fprintf(stderr, "Failed to allocate memory for argv\n");
        exit(EXIT_FAILURE);
    }

    /* Copy each argument to new_argv */
    for (i = 0; i < argc; i++) {
        new_argv[i] = strdup(argv[i]);
        if (new_argv[i] == NULL) {
            fprintf(stderr, "Failed to allocate memory for argv[%d]\n", i);
            exit(EXIT_FAILURE);
        }
    }
   
   Array1D<double> dt0;
   *error = pampa.initialize(argc, new_argv, dt0);
   
    /* Clean up */
    for (i = 0; i < argc; i++) {
        free(new_argv[i]);
    }
    free(new_argv);
   
   *ndt = dt0.size();
   *dt = new double[*ndt];
   for (int i = 0; i < *ndt; i++)
      (*dt)[i] = dt0(i);
   
   return 0;
   
}

/* Solve: */
int pampa_solve(int n, double dt, double t, int* error) {
   
   *error = pampa.solve(n, dt, t);
   
   return 0;
   
}

/* Finalize: */
int pampa_finalize(double** dt, int* error) {
   
   *error = pampa.finalize();
   
   delete[] *dt;
   
   return 0;
   
}

}
