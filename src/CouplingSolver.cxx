#include "CouplingSolver.hxx"

/* Initialize: */
int CouplingSolver::initialize() {
   
   /* Initialize all the solvers: */
   for (int i = 0; i < solvers.size(); i++) {
      PAMPA_CALL(solvers(i)->initialize(), "unable to initialize the solver");
   }
   
   return 0;
   
}

/* Get the solution: */
int CouplingSolver::solve(int n, double dt) {
   
   /* Get the solution from all the solvers: */
   for (int i = 0; i < solvers.size(); i++) {
      PAMPA_CALL(solvers(i)->solve(n, dt), "unable to get the solution from the solver");
   }
   
   return 0;
   
}

/* Output the solution: */
int CouplingSolver::output(const std::string& filename, int n) const {
   
   /* Output the solution from all the solvers: */
   for (int i = 0; i < solvers.size(); i++) {
      PAMPA_CALL(solvers(i)->output(filename, n), "unable to output the solution from the solver");
   }
   
   return 0;
   
}

/* Finalize: */
int CouplingSolver::finalize() {
   
   /* Finalize all the solvers: */
   for (int i = 0; i < solvers.size(); i++) {
      PAMPA_CALL(solvers(i)->finalize(), "unable to finalize the solver");
   }
   
   return 0;
   
}
