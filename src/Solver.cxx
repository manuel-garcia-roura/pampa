#include "Solver.hxx"

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
