#pragma once

#include "Solver.hxx"

/* The CouplingSolver class: */
class CouplingSolver : public Solver {
   
   private:
      
      /* Coupled solvers: */
      const Array1D<Solver*> solvers;
   
   public:
      
      /* The CouplingSolver constructor: */
      CouplingSolver(const std::string& name, const Mesh* mesh, const Array1D<Solver*> solvers) : 
         Solver(name, mesh), solvers(solvers) {}
      
      /* The CouplingSolver destructor: */
      virtual ~CouplingSolver() {}
      
      /* Initialize: */
      int WARN_UNUSED initialize();
      
      /* Get the solution: */
      int WARN_UNUSED solve(int n = 0, double dt = 0.0);
      
      /* Output the solution: */
      int WARN_UNUSED output(const std::string& filename, int n = 0, bool write_mesh = true) const;
      
      /* Finalize: */
      int WARN_UNUSED finalize();
   
};
