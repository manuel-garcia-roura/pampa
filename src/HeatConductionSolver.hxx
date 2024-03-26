#pragma once

#include <petscksp.h>

#include "Solver.hxx"

/* The HeatConductionSolver class: */
class HeatConductionSolver : public Solver {
   
   private:
      
      /* Coefficient matrix for the linear system A*x = b: */
      Mat A;
      
      /* Temperature (solution vector): */
      Vec T;
      
      /* Source term from Dirichlet boundary conditions (right-hand side): */
      Vec qbc;
      
      /* Volumetric heat source (right-hand side): */
      Vec q;
      
      /* Krylov Subspace Solver (KSP) context: */
      KSP ksp;
      
      /* Previous time step: */
      double dt0 = 0.0;
      
      /* Check the material data: */
      int WARN_UNUSED checkMaterials();
      
      /* Build the coefficient matrix and the solution and RHS vectors: */
      int WARN_UNUSED build();
      
      /* Build the coefficient matrix and the RHS vector: */
      int WARN_UNUSED buildMatrix();
      
      /* Set the time-derivative terms: */
      int WARN_UNUSED setTimeDerivative(double dt);
      
      /* Write the solution to a plain-text file in .vtk format: */
      int WARN_UNUSED writeVTK(const std::string& filename) const;
   
   public:
      
      /* The HeatConductionSolver constructor: */
      HeatConductionSolver(const Mesh* mesh, const Array1D<Material>& materials) : 
         Solver(mesh, materials) {}
      
      /* The HeatConductionSolver destructor: */
      ~HeatConductionSolver() {}
      
      /* Initialize: */
      int WARN_UNUSED initialize(int argc, char* argv[]);
      
      /* Solve the linear system to get the solution: */
      int WARN_UNUSED solve(int i = 0, double dt = 0.0);
      
      /* Output the solution: */
      int WARN_UNUSED output(const std::string& filename);
      
      /* Finalize: */
      int WARN_UNUSED finalize();
   
};
