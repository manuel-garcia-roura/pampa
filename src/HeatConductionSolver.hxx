#pragma once

#include "PhysicsSolver.hxx"

/* The HeatConductionSolver class: */
class HeatConductionSolver : public PhysicsSolver {
   
   private:
      
      /* Coarse mesh for nodal-level input and output: */
      const Mesh* mesh_nodal = nullptr;
      
      /* Boundary conditions (1-based indexed): */
      Array1D<BoundaryCondition> bcs;
      
      /* Boundary-condition indices for materials: */
      Array1D<int> bcmat_indices;
      
      /* Total thermal power: */
      Function power;
      
      /* Convergence tolerance and p-norm for nonlinear problems: */
      double tol = 1.0, p = 2.0;
      
      /* Switch for nonlinear problems where the thermal properties depend on the temperature: */
      bool nonlinear = false;
      
      /* Coefficient matrix for the linear system A*x = b: */
      Mat A = 0;
      
      /* Right-hand side for the linear system A*x = b: */
      Vec b = 0;
      
      /* Volumetric heat source: */
      Vec q = 0;
      
      /* Temperature: */
      Vec T = 0, Tprev = 0, T0 = 0;
      
      /* Nodal volumetric heat source: */
      Vec qnodal = 0;
      
      /* Nodal temperatures: */
      Array1D<Vec> Tnodal;
      
      /* Check the material data: */
      int WARN_UNUSED checkMaterials(bool transient = false);
      
      /* Build the coefficient matrix, and the solution and RHS vectors: */
      int WARN_UNUSED build();
      
      /* Initialize the volumetric heat source: */
      int WARN_UNUSED initializeHeatSource();
      
      /* Calculate the volumetric heat source from the nodal power: */
      int WARN_UNUSED calculateHeatSource();
      
      /* Calculate the nodal temperatures: */
      int WARN_UNUSED calculateNodalTemperatures();
      
      /* Build the coefficient matrix and the RHS vector: */
      int WARN_UNUSED buildMatrix(int n, double dt, double t);
      
      /* Print the solution summary to standard output: */
      int WARN_UNUSED printLog(int n = 0) const;
      
      /* Write the solution to a plain-text file in .vtk format: */
      int WARN_UNUSED writeVTK(const std::string& filename) const;
      
      /* Write the solution to a binary file in PETSc format: */
      int WARN_UNUSED writePETSc(int n = 0) const;
   
   public:
      
      /* The HeatConductionSolver constructor: */
      HeatConductionSolver(const Mesh* mesh_nodal, const Array1D<Material*>& materials) : 
         PhysicsSolver("conduction", materials), mesh_nodal(mesh_nodal), 
         bcmat_indices{materials.size(), -1} {}
      
      /* The HeatConductionSolver destructor: */
      ~HeatConductionSolver() {}
      
      /* Read the solver from a plain-text input file: */
      int WARN_UNUSED read(std::ifstream& file, Array1D<Solver*>& solvers);
      
      /* Solve the linear system to get the solution: */
      int WARN_UNUSED solve(int n = 0, double dt = 0.0, double t = 0.0);
   
};
