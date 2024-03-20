#pragma once

#include "NeutronicSolver.hxx"

/* The DiffusionSolver class: */
class DiffusionSolver : public NeutronicSolver {
   
   private:
      
      /* Check the material data: */
      int checkMaterials();
      
      /* Build the coefficient matrices and solution vectors: */
      int build();
      
      /* Build the coefficient matrices: */
      int buildMatrices();
      
      /* Build the solution vectors: */
      int buildVectors();
      
      /* Get the solution after solving the eigensystem: */
      int getSolution();
      
      /* Write the solution to a plain-text file in .vtk format: */
      int writeVTK(const std::string& filename) const;
      
      /* Write the solution to a binary file in PETSc format: */
      int writePETSc(const std::string& filename) const;
      
      /* Destroy the solution vectors: */
      int destroyVectors();
   
   public:
      
      /* The DiffusionSolver constructor: */
      DiffusionSolver(const Mesh* mesh, const Array1D<Material>& materials, int num_groups) : 
         NeutronicSolver(mesh, materials, num_groups) {}
      
      /* The DiffusionSolver destructor: */
      ~DiffusionSolver() {}
   
};
