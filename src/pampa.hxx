#ifdef __cplusplus
extern "C" {
#endif

/* Initialize the calculation: */
void pampa_initialize(int argc, char* argv[], double** dt, int* ndt, int* error);

/* Get the solution: */
void pampa_solve(int n, double dt, double t, int* error);

/* Finalize the calculation: */
void pampa_finalize(double** dt, int* error);

/* Get the values for a given field: */
void pampa_get_field(double* v, const char name[], int* error);

/* Set the values for a given field: */
void pampa_set_field(const double* v, const char name[], int* error);

#ifdef __cplusplus
}
#endif
