#include "Solver.hxx"

/* Read the solver from a plain-text input file: */
int Solver::read(const std::string& filename, Array1D<Solver*>& solvers) {
   
   /* Open the input file: */
   std::ifstream file(filename, std::ios_base::in);
   PAMPA_CHECK(!file.is_open(), "unable to open " + filename);
   
   /* Read the solver: */
   PAMPA_CHECK(read(file, solvers), "unable to read the solver from " + filename);
   
   return 0;
   
}

/* Get the values for a given field: */
int Solver::getField(double* v, const std::string& name) const {
   
   /* Find the field: */
   Field field;
   PAMPA_CHECK(utils::find(name, fields, field), "unable to find field '" + name + "'");
   
   /* Copy the values: */
   PAMPA_CHECK(petsc::copy(field.vec, v), "unable to copy the values for field '" + name + "'");
   
   return 0;
   
}

/* Set the values for a given field: */
int Solver::setField(const double* v, const std::string& name) {
   
   /* Find the field: */
   Field field;
   PAMPA_CHECK(utils::find(name, fields, field), "unable to find field '" + name + "'");
   
   /* Copy the values: */
   PAMPA_CHECK(petsc::copy(v, field.vec), "unable to copy the values for field '" + name + "'");
   
   return 0;
   
}
