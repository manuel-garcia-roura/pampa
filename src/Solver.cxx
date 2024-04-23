#include "Solver.hxx"

/* Add a boundary condition: */
int Solver::addBoundaryCondition(const BoundaryCondition& bc, int l) {
   
   /* Initialize the boundary-condition array if not done yet: */
   if (bcs.empty()) bcs.resize(1);
   
   /* Check the boundary-condition index: */
   PAMPA_CHECK(l != bcs.size(), 1, "wrong boundary-condition index");
   
   /* Keep the boundary condition: */
   bcs.pushBack(bc);
   
   return 0;
   
}

/* Read the solver from a plain-text input file: */
int Solver::read(const std::string& filename, Array1D<Solver*>& solvers) {
   
   /* Open the input file: */
   std::ifstream file(filename);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Read the solver: */
   PAMPA_CALL(read(file, solvers), "unable to read the solver from " + filename);
   
   return 0;
   
}

/* Find a solver from its name: */
int utils::find(const std::string& name, const Array1D<Solver*>& solvers, Solver** solver) {
   
   /* Find the solver: */
   bool found = false;
   for (int i = 0; i < solvers.size(); i++) {
      if (solvers(i)->getName() == name) {
         PAMPA_CHECK(found, 1, "duplicated solver '" + name + "'");
         *solver = solvers(i);
         found = true;
      }
   }
   PAMPA_CHECK(!found, 2, "undefined solver '" + name + "'");
   
   return 0;
   
}
