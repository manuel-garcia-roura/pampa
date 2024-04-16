#pragma once

#include "PhysicsSolver.hxx"

/* The HeatConductionSolver class: */
class HeatConductionSolver : public PhysicsSolver {
   
   private:
      
      /* Coefficient matrix for the linear system A*x = b: */
      Mat A = 0;
      
      /* Right-hand side for the linear system A*x = b: */
      Vec b = 0;
      
      /* Volumetric heat source: */
      Vec q = 0;
      
      /* Temperature: */
      Vec T = 0, T0 = 0;
      
      /* Fixed temperatures for given materials: */
      Array1D<double> fixed_temperatures;
      
      /* Check the material data: */
      int WARN_UNUSED checkMaterials() const;
      
      /* Build the coefficient matrix, and the solution and RHS vectors: */
      int WARN_UNUSED build();
      
      /* Build the coefficient matrix and the RHS vector: */
      int WARN_UNUSED buildMatrix(int n, double dt);
      
      /* Print the solution summary to standard output: */
      int WARN_UNUSED printLog(int n = 0) const;
      
      /* Write the solution to a plain-text file in .vtk format: */
      int WARN_UNUSED writeVTK(const std::string& filename) const;
      
      /* Write the solution to a binary file in PETSc format: */
      int WARN_UNUSED writePETSc() const;
   
   public:
      
      /* The HeatConductionSolver constructor: */
      HeatConductionSolver(const Mesh* mesh, const Array1D<Material>& materials) : 
         PhysicsSolver("conduction", mesh, materials), fixed_temperatures{materials.size(), -1.0} {}
      
      /* The HeatConductionSolver destructor: */
      ~HeatConductionSolver() {}
      
      /* Add a fixed temperature for a given material: */
      int WARN_UNUSED addFixedTemperature(int mat, double x);
      
      /* Solve the linear system to get the solution: */
      int WARN_UNUSED solve(int n = 0, double dt = 0.0);
   
};
