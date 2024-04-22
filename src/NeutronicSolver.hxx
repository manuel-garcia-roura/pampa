#pragma once

#include "PhysicsSolver.hxx"

/* The NeutronicSolver class: */
class NeutronicSolver : public PhysicsSolver {
   
   protected:
      
      /* Number of energy groups: */
      int num_energy_groups = -1;
      
      /* Coefficient matrices for the eigen- (R*x = (1/keff)*F*x) or linear (R*x = b) system: */
      Mat R = 0, F = 0;
      
      /* Right-hand side for the linear (R*x = b) system: */
      Vec b = 0;
      
      /* Delayed neutron source: */
      Vec S = 0;
      
      /* Multiplication factor: */
      double keff = -1.0;
      
      /* Scalar neutron flux: */
      Vec phi = 0, phi0 = 0;
      
      /* Thermal power: */
      Vec q = 0;
      
      /* Production rate: */
      Vec P = 0;
      
      /* Get the flat index for cell i and group g: */
      int index(int i, int g) const {return i*num_energy_groups + g;}
      
      /* Build the coefficient matrices and the RHS vector: */
      virtual int WARN_UNUSED buildMatrices(int n, double dt) {PAMPA_CHECK_VIRTUAL}
      
      /* Solve the linear system and get the solution: */
      virtual int WARN_UNUSED getSolution(int n = 0) {PAMPA_CHECK_VIRTUAL}
      
      /* Normalize the scalar flux: */
      int WARN_UNUSED normalizeScalarFlux();
      
      /* Calculate the thermal power and the production rate: */
      int WARN_UNUSED calculatePowerAndProductionRate();
      
      /* Print the solution summary to standard output: */
      int WARN_UNUSED printLog(int n = 0) const;
   
   public:
      
      /* The NeutronicSolver constructor: */
      NeutronicSolver(const std::string& name, const Mesh* mesh, 
         const Array1D<Material>& materials, int num_energy_groups) : 
         PhysicsSolver(name, mesh, materials), num_energy_groups(num_energy_groups) {}
      
      /* The NeutronicSolver destructor: */
      virtual ~NeutronicSolver() {}
      
      /* Solve the linear system to get the solution: */
      int WARN_UNUSED solve(int n = 0, double dt = 0.0, double t = 0.0);
   
};
