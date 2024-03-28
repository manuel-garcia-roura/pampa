#pragma once

#include "NeutronicSolver.hxx"

/* The DiffusionSolver class: */
class DiffusionSolver : public NeutronicSolver {
   
   private:
      
      /* Check the material data: */
      int WARN_UNUSED checkMaterials();
      
      /* Build the coefficient matrices, the solution vector and the EPS context: */
      int WARN_UNUSED build();
      
      /* Build the coefficient matrices: */
      int WARN_UNUSED buildMatrices();
      
      /* Get the solution after solving the eigensystem: */
      int WARN_UNUSED getSolution();
      
      /* Write the solution to a plain-text file in .vtk format: */
      int WARN_UNUSED writeVTK(const std::string& filename) const;
      
      /* Write the solution to a binary file in PETSc format: */
      int WARN_UNUSED writePETSc(const std::string& filename) const;
   
   public:
      
      /* The DiffusionSolver constructor: */
      DiffusionSolver(const Mesh* mesh, const Array1D<Material>& materials, 
         int num_energy_groups) : NeutronicSolver(mesh, materials, num_energy_groups) {}
      
      /* The DiffusionSolver destructor: */
      ~DiffusionSolver() {}
   
};
