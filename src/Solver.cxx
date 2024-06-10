#include "Solver.hxx"

/* Read the solver from a plain-text input file: */
int Solver::read(const std::string& filename, Array1D<Solver*>& solvers) {
   
   /* Open the input file: */
   std::ifstream file(filename, std::ios_base::in);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Read the solver: */
   PAMPA_CALL(read(file, solvers), "unable to read the solver from " + filename);
   
   return 0;
   
}
