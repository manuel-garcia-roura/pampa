#pragma once

#include "PhysicsSolver.hxx"

/* The PrecursorSolver class: */
class PrecursorSolver : public PhysicsSolver {
   
   private:
      
      /* Number of delayed-neutron precursor groups: */
      int num_precursor_groups = -1;
      
      /* Production rate: */
      Vec P = 0;
      
      /* Precursor population: */
      Vec C = 0, C0 = 0;
      
      /* Delayed neutron source: */
      Vec S = 0;
      
      /* Get the flat index for cell i and group g: */
      int index(int i, int g) const {return i*num_precursor_groups + g;}
      
      /* Check the material data: */
      int WARN_UNUSED checkMaterials(bool transient = false) const;
      
      /* Build the solution and source vectors: */
      int WARN_UNUSED build();
      
      /* Print the solution summary to standard output: */
      int WARN_UNUSED printLog(int n = 0) const;
      
      /* Write the solution to a plain-text file in .vtk format: */
      int WARN_UNUSED writeVTK(const std::string& filename) const;
      
      /* Write the solution to a binary file in PETSc format: */
      int WARN_UNUSED writePETSc() const;
   
   public:
      
      /* The PrecursorSolver constructor: */
      PrecursorSolver(const Mesh* mesh, const Array1D<Material>& materials, 
         int num_precursor_groups) : PhysicsSolver("precursors", mesh, materials), 
         num_precursor_groups(num_precursor_groups) {}
      
      /* The PrecursorSolver destructor: */
      ~PrecursorSolver() {}
      
      /* Solve the linear system to get the solution: */
      int WARN_UNUSED solve(int n = 0, double dt = 0.0);
   
};
