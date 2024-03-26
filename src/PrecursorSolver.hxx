#pragma once

#include <petsc.h>

#include "Solver.hxx"

/* The PrecursorSolver class: */
class PrecursorSolver : public Solver {
   
   private:
      
      /* Number of delayed-neutron precursor groups: */
      int num_precursor_groups = -1;
      
      /* Precursor population (solution vector): */
      Vec C;
      
      /* Production rate (source term): */
      Vec P;
      
      /* Get the flat index for cell i and group g: */
      int index(int i, int g, int ng) const {return i*ng + g;}
      
      /* Check the material data: */
      int WARN_UNUSED checkMaterials();
      
      /* Write the solution to a plain-text file in .vtk format: */
      int WARN_UNUSED writeVTK(const std::string& filename) const;
   
   public:
      
      /* The PrecursorSolver constructor: */
      PrecursorSolver(const Mesh* mesh, const Array1D<Material>& materials, 
         int num_precursor_groups) : Solver(mesh, materials), 
         num_precursor_groups(num_precursor_groups) {}
      
      /* The PrecursorSolver destructor: */
      ~PrecursorSolver() {}
      
      /* Initialize: */
      int WARN_UNUSED initialize(int argc, char* argv[]);
      
      /* Solve the linear system to get the solution: */
      int WARN_UNUSED solve(int n = 0, double dt = 0.0);
      
      /* Output the solution: */
      int WARN_UNUSED output(const std::string& filename);
      
      /* Finalize: */
      int WARN_UNUSED finalize();
   
};
