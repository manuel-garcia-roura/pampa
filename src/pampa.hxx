#ifdef __cplusplus
extern "C" {
#endif

/* Initialize the calculation: */
int pampa_initialize(int argc, char* argv[], double** dt, int* ndt, int* error);

/* Get the solution: */
int pampa_solve(int n, double dt, double t, int* error);

/* Finalize the calculation: */
int pampa_finalize(double** dt, int* error);

#ifdef __cplusplus
}
#endif
