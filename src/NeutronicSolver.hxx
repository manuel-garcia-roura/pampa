#pragma once

#include "Solver.hxx"

/* The NeutronicSolver class: */
class NeutronicSolver : public Solver {
   
   protected:
      
      /* Number of energy groups: */
      int num_energy_groups = -1;
      
      /* Neutron velocity for each energy group: */
      Array1D<double> v;
      
      /* Coefficient matrices for the eigen- (R*x = (1/keff)*F*x) or linear (R*x = q) system: */
      Mat R = 0, F = 0;
      
      /* Scalar neutron flux (eigenvector): */
      Vec phi = 0;
      
      /* Neutron source (right-hand side): */
      Vec q = 0;
      
      /* Multiplication factor (eigenvalue): */
      double keff = -1.0;
      
      /* Get the flat index for cell i and group g: */
      int index(int i, int g) const {return i*num_energy_groups + g;}
      
      /* Build the coefficient matrices: */
      virtual int WARN_UNUSED buildMatrices(int n, double dt) {PAMPA_CHECK_VIRTUAL}
      
      /* Solve the linear system and get the solution: */
      virtual int WARN_UNUSED getSolution(int n = 0) {PAMPA_CHECK_VIRTUAL}
      
      /* Normalize the scalar flux: */
      int WARN_UNUSED normalizeScalarFlux();
      
      /* Calculate the thermal power: */
      int WARN_UNUSED calculatePower();
      
      /* Print the solution summary to standard output: */
      int WARN_UNUSED printLog(int n = 0) const;
   
   public:
      
      /* The NeutronicSolver constructor: */
      NeutronicSolver(const Mesh* mesh, const Array1D<Material>& materials, 
         int num_energy_groups) : Solver(mesh, materials), num_energy_groups(num_energy_groups) {}
      
      /* The NeutronicSolver destructor: */
      virtual ~NeutronicSolver() {}
      
      /* Solve the linear system to get the solution: */
      int WARN_UNUSED solve(int n = 0, double dt = 0.0);
   
};
