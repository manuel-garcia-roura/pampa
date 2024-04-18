#pragma once

#include "Solver.hxx"

/* The CouplingSolver class: */
class CouplingSolver : public Solver {
   
   private:
      
      /* Coupled solvers: */
      const Array1D<Solver*> solvers;
      
      /* Switch to use implicit coupling: */
      bool implicit = false;
      
      /* Convergence p-norm and tolerance for implicit coupling: */
      double tol = -1.0, p = -1.0;
   
   public:
      
      /* The CouplingSolver constructor: */
      CouplingSolver(const std::string& name, const Mesh* mesh, const Array1D<Solver*> solvers, 
         bool implicit = false, double tol = -1.0, double p = -1.0) : Solver(name, mesh), 
         solvers(solvers), implicit(implicit), tol(tol), p(p) {}
      
      /* The CouplingSolver destructor: */
      virtual ~CouplingSolver() {}
      
      /* Initialize: */
      int WARN_UNUSED initialize(bool transient = false);
      
      /* Get the solution: */
      int WARN_UNUSED solve(int n = 0, double dt = 0.0);
      
      /* Output the solution: */
      int WARN_UNUSED output(const std::string& filename, int n = 0, bool write_mesh = true) const;
      
      /* Finalize: */
      int WARN_UNUSED finalize();
   
};
