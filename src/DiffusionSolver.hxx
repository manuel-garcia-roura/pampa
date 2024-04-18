#pragma once

#include "NeutronicSolver.hxx"

/* The DiffusionSolver class: */
class DiffusionSolver : public NeutronicSolver {
   
   private:
      
      /* Build the coefficient matrices and the RHS vector: */
      int WARN_UNUSED buildMatrices(int n, double dt);
      
      /* Solve the linear system and get the solution: */
      int WARN_UNUSED getSolution(int n = 0);
      
      /* Check the material data: */
      int WARN_UNUSED checkMaterials(bool transient = false) const;
      
      /* Build the coefficient matrices and the solution vector: */
      int WARN_UNUSED build();
      
      /* Write the solution to a plain-text file in .vtk format: */
      int WARN_UNUSED writeVTK(const std::string& filename) const;
      
      /* Write the solution to a binary file in PETSc format: */
      int WARN_UNUSED writePETSc() const;
   
   public:
      
      /* The DiffusionSolver constructor: */
      DiffusionSolver(const Mesh* mesh, const Array1D<Material>& materials, 
         int num_energy_groups) : NeutronicSolver("diffusion", mesh, materials, num_energy_groups) 
         {}
      
      /* The DiffusionSolver destructor: */
      ~DiffusionSolver() {}
   
};
