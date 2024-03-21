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
      int checkMaterials();
      
      /* Build the coefficient matrix and the solution and RHS vectors: */
      int build();
      
      /* Build the coefficient matrix and the RHS vector: */
      int buildMatrix();
      
      /* Write the solution to a plain-text file in .vtk format: */
      int writeVTK(const std::string& filename) const;
   
   public:
      
      /* The HeatConductionSolver constructor: */
      HeatConductionSolver(const Mesh* mesh, const Array1D<Material>& materials) : 
         Solver(mesh, materials) {}
      
      /* The HeatConductionSolver destructor: */
      ~HeatConductionSolver() {}
      
      /* Initialize: */
      int initialize(int argc, char* argv[]);
      
      /* Solve the linear system to get the solution: */
      int solve();
      
      /* Output the solution: */
      int output(const std::string& filename);
      
      /* Finalize: */
      int finalize();
   
};
