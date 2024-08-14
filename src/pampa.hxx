#ifdef __cplusplus
extern "C" {
#endif

/* Initialize: */
int pampa_initialize(int argc, char* argv[], double** dt, int* ndt, int* error);

/* Solve: */
int pampa_solve(int n, double dt, double t, int* error);

/* Finalize: */
int pampa_finalize(double** dt, int* error);

#ifdef __cplusplus
}
#endif
