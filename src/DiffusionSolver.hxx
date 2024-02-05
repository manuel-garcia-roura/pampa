#pragma once

#include <cmath>
#include <fstream>
#include <iostream>

#include <slepceps.h>

#include "Solver.hxx"
#include "Mesh.hxx"
#include "Material.hxx"
#include "mpi.hxx"
#include "petsc.hxx"
#include "math.hxx"
#include "utils.hxx"

/* The DiffusionSolver class: */
class DiffusionSolver : public Solver {
   
   private:
      
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
      DiffusionSolver(const Mesh* mesh, const Array1D<Material>& materials, 
         const TransportMethod& method) : Solver(mesh, materials, method) {}
      
      /* The DiffusionSolver destructor: */
      ~DiffusionSolver() {}
   
};
