#pragma once

#include "Solver.hxx"

/* The HeatConductionSolver class: */
class HeatConductionSolver : public Solver {
   
   private:
      
      /* Coefficient matrix for the linear system A*x = b: */
      Mat A = 0;
      
      /* Temperature (solution vector): */
      Vec T = 0;
      
      /* Source term from Dirichlet boundary conditions (right-hand side): */
      Vec qbc = 0;
      
      /* Volumetric heat source (right-hand side): */
      Vec q = 0;
      
      /* Previous time step: */
      double dt0 = 0.0;
      
      /* Check the material data: */
      int WARN_UNUSED checkMaterials();
      
      /* Build the coefficient matrix, the solution and RHS vectors, and the KSP context: */
      int WARN_UNUSED build();
      
      /* Build the coefficient matrix and the RHS vector: */
      int WARN_UNUSED buildMatrix();
      
      /* Set the time-derivative terms: */
      int WARN_UNUSED setTimeDerivative(double dt);
   
   public:
      
      /* The HeatConductionSolver constructor: */
      HeatConductionSolver(const Mesh* mesh, const Array1D<Material>& materials) : 
         Solver(mesh, materials) {}
      
      /* The HeatConductionSolver destructor: */
      ~HeatConductionSolver() {}
      
      /* Solve the linear system to get the solution: */
      int WARN_UNUSED solve(int n = 0, double dt = 0.0);
      
      /* Output the solution: */
      int WARN_UNUSED output(const std::string& filename);
   
};
