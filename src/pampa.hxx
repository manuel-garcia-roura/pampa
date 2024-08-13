#ifdef __cplusplus
extern "C" {
#endif

/* Initialize: */
int initialize(int argc, char* argv[], double** dt, int* ndt);

/* Solve: */
int solve(int n, double dt, double t);

/* Finalize: */
int finalize(double** dt);

#ifdef __cplusplus
}
#endif
