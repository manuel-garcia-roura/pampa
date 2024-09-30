#pragma once

#include "Solver.hxx"

/* The CouplingSolver class: */
class CouplingSolver : public Solver {
   
   private:
      
      /* Coupled solvers: */
      Array1D<Solver*> coupled_solvers;
      
      /* Switch to use implicit coupling: */
      bool implicit = false;
   
   public:
      
      /* The CouplingSolver constructor: */
      CouplingSolver(const std::string& name, const Mesh* mesh) : Solver(name, mesh) {}
      
      /* The CouplingSolver destructor: */
      ~CouplingSolver() {}
      
      /* Read the solver from a plain-text input file: */
      int WARN_UNUSED read(std::ifstream& file, Array1D<Solver*>& solvers);
      
      /* Initialize: */
      int WARN_UNUSED initialize(bool transient = false);
      
      /* Get the solution: */
      int WARN_UNUSED solve(int n = 0, double dt = 0.0, double t = 0.0);
      
      /* Output the solution: */
      int WARN_UNUSED output(const std::string& path, int n = 0, bool write_mesh = true) const;
      
      /* Finalize: */
      int WARN_UNUSED finalize();
   
};
