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
      
      /* Volumetric heat source (right-hand side): */
      Vec q;
      
      /* Krylov Subspace Solver (KSP) context: */
      KSP ksp;
      
      /* Check the material data: */
      int WARN_UNUSED checkMaterials();
      
      /* Build the coefficient matrix and the solution and RHS vectors: */
      int WARN_UNUSED build();
      
      /* Build the coefficient matrix and the RHS vector: */
      int WARN_UNUSED buildMatrix();
      
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
      int WARN_UNUSED solve();
      
      /* Output the solution: */
      int WARN_UNUSED output(const std::string& filename);
      
      /* Finalize: */
      int WARN_UNUSED finalize();
   
};
