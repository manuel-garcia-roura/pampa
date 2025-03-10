#pragma once

#include "NeutronicSolver.hxx"

/* The DiffusionSolver class: */
class DiffusionSolver : public NeutronicSolver {
   
   private:
      
      /* Build the coefficient matrices and the RHS vector: */
      int WARN_UNUSED buildMatrices(int n, double dt, double t);
      
      /* Solve the linear system and get the solution: */
      int WARN_UNUSED getSolution(int n = 0);
      
      /* Check the material data: */
      int WARN_UNUSED checkMaterials(bool transient = false);
      
      /* Build the coefficient matrices and the solution vector: */
      int WARN_UNUSED build();
      
      /* Write the solution to a plain-text file in .vtk format: */
      int WARN_UNUSED writeVTK(const std::string& path, int n = 0) const;
      
      /* Write the solution to a binary file in PETSc format: */
      int WARN_UNUSED writePETSc(int n = 0) const;
   
   public:
      
      /* The DiffusionSolver constructor: */
      DiffusionSolver(const Mesh* mesh, const Array1D<Material*>& materials) : 
         NeutronicSolver("diffusion", mesh, materials) {}
      
      /* The DiffusionSolver destructor: */
      ~DiffusionSolver() {}
      
      /* Read the solver from a plain-text input file: */
      int WARN_UNUSED read(std::ifstream& file, Array1D<Solver*>& solvers);
   
};
