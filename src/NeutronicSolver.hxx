#pragma once

#include "Solver.hxx"

/* The NeutronicSolver class: */
class NeutronicSolver : public Solver {
   
   protected:
      
      /* Number of energy groups: */
      int num_energy_groups = -1;
      
      /* Coefficient matrices for the generalized eigensystem R*x = (1/keff)*F*x: */
      Mat R = 0, F = 0;
      
      /* Scalar neutron flux (eigenvector): */
      Vec phi = 0;
      
      /* Multiplication factor (eigenvalue): */
      double keff = -1.0;
      
      /* Get the flat index for cell i and group g: */
      int index(int i, int g) const {return i*num_energy_groups + g;}
      
      /* Get the solution after solving the eigensystem: */
      virtual int WARN_UNUSED getSolution() 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Write the solution to a plain-text file in .vtk format: */
      virtual int WARN_UNUSED writeVTK(const std::string& filename) const 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Write the solution to a binary file in PETSc format: */
      virtual int WARN_UNUSED writePETSc(const std::string& filename) const 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Normalize the scalar flux: */
      int WARN_UNUSED normalizeScalarFlux();
   
   public:
      
      /* The NeutronicSolver constructor: */
      NeutronicSolver(const Mesh* mesh, const Array1D<Material>& materials, 
         int num_energy_groups) : Solver(mesh, materials), num_energy_groups(num_energy_groups) {}
      
      /* The NeutronicSolver destructor: */
      virtual ~NeutronicSolver() {}
      
      /* Solve the eigensystem to get the solution: */
      int WARN_UNUSED solve(int n = 0, double dt = 0.0);
      
      /* Output the solution: */
      int WARN_UNUSED output(const std::string& filename);
   
};
