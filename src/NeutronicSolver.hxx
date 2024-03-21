#pragma once

#include <slepceps.h>

#include "Solver.hxx"

/* The NeutronicSolver class: */
class NeutronicSolver : public Solver {
   
   protected:
      
      /* Number of energy groups: */
      int num_groups = -1;
      
      /* Coefficient matrices for the generalized eigensystem R*x = (1/keff)*F*x: */
      Mat R, F;
      
      /* Scalar neutron flux (eigenvector): */
      Vec phi;
      
      /* Multiplication factor (eigenvalue): */
      double keff = -1.0;
      
      /* Eigenvalue Problem Solver (EPS) context: */
      EPS eps;
      
      /* Get the flat index for cell i and group g: */
      int index(int i, int g, int ng) const {return i*ng + g;}
      
      /* Get the flat index for cell i, group g and direction m: */
      int index(int i, int g, int m, int ng, int nm) const {return i*nm*ng + g*nm + m;}
      
      /* Check the material data: */
      virtual int checkMaterials() 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Build the coefficient matrices and the solution vectors: */
      virtual int build() 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Get the solution after solving the eigensystem: */
      virtual int getSolution() 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Write the solution to a plain-text file in .vtk format: */
      virtual int writeVTK(const std::string& filename) const 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Write the solution to a binary file in PETSc format: */
      virtual int writePETSc(const std::string& filename) const 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Destroy the solution vectors: */
      virtual int destroyVectors() 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Normalize the scalar flux: */
      int normalizeScalarFlux();
   
   public:
      
      /* The NeutronicSolver constructor: */
      NeutronicSolver(const Mesh* mesh, const Array1D<Material>& materials, int num_groups) : 
         Solver(mesh, materials), num_groups(num_groups) {}
      
      /* The NeutronicSolver destructor: */
      ~NeutronicSolver() {}
      
      /* Initialize: */
      int initialize(int argc, char* argv[]);
      
      /* Solve the linear system to get the solution: */
      int solve();
      
      /* Output the solution: */
      int output(const std::string& filename);
      
      /* Finalize: */
      int finalize();
   
};
