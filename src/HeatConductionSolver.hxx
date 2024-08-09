#pragma once

#include "PhysicsSolver.hxx"

/* The HeatConductionSolver class: */
class HeatConductionSolver : public PhysicsSolver {
   
   private:
      
      /* Meshes (main mesh and submeshes for microscale problems): */
      const Array1D<Mesh*> meshes;
      
      /* Boundary conditions (1-based indexed): */
      Array1D<BoundaryCondition> bcs;
      
      /* Fixed temperatures for given materials: */
      Array1D<Function> fixed_temperatures;
      
      /* Total thermal power: */
      Function power;
      
      /* Convergence tolerance and p-norm for nonlinear problems: */
      double tol = 1.0, p = 2.0;
      
      /* Switch for nonlinear problems where the thermal properties depend on the temperature: */
      bool nonlinear = false;
      
      /* Switches for multiscale problems: */
      bool multiscale = false, microscale = false;
      
      /* Solvers for microscale problems: */
      Array1D<HeatConductionSolver> solvers;
      
      /* Coefficient matrix for the linear system A*x = b: */
      Mat A = 0;
      
      /* Right-hand side for the linear system A*x = b: */
      Vec b = 0;
      
      /* Volumetric heat source: */
      Vec q = 0, qfull = 0;
      
      /* Effective volumetric heat source from microscale problems: */
      Array1D<double> s, s0, smicro, smacro, fb;
      
      /* Temperature: */
      Vec T = 0, Tprev = 0, T0 = 0, Tfull = 0;
      
      /* Check the material data: */
      int WARN_UNUSED checkMaterials(bool transient = false);
      
      /* Build the coefficient matrix, and the solution and RHS vectors: */
      int WARN_UNUSED build();
      
      /* Build the coefficient matrix and the RHS vector: */
      int WARN_UNUSED buildMatrix(int n, double dt, double t);
      
      /* Print the solution summary to standard output: */
      int WARN_UNUSED printLog(int n = 0) const;
      
      /* Write the solution to a plain-text file in .vtk format: */
      int WARN_UNUSED writeVTK(const std::string& filename) const;
      
      /* Write the solution to a binary file in PETSc format: */
      int WARN_UNUSED writePETSc(int n = 0) const;
      
      /* Get the effective thermal properties: */
      int WARN_UNUSED getEffectiveThermalProperties(double& k, double& rho, double& cp, 
         double& s, double& smicro, double& fb, double t, int& ic, PetscScalar* qfull_data, 
         PetscScalar* Tfull_data) const;
      
      /* Get the volume for fuel materials: */
      double getFuelVolume() const;
   
   public:
      
      /* The HeatConductionSolver constructor: */
      HeatConductionSolver(const Array1D<Mesh*>& meshes, const Array1D<Material*>& materials) : 
         PhysicsSolver("conduction", meshes(0), materials), meshes(meshes), 
         fixed_temperatures{materials.size()}, multiscale(meshes.size() > 1) {}
      
      /* The HeatConductionSolver destructor: */
      ~HeatConductionSolver() {}
      
      /* Read the solver from a plain-text input file: */
      int WARN_UNUSED read(std::ifstream& file, Array1D<Solver*>& solvers);
      
      /* Solve the linear system to get the solution: */
      int WARN_UNUSED solve(int n = 0, double dt = 0.0, double t = 0.0);
   
};
