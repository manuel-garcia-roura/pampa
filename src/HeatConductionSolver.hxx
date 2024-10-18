#pragma once

#include "PhysicsSolver.hxx"
#include "HeatPipe.hxx"

/* The HeatPipeBoundaryConditions struct: */
struct HeatPipeBoundaryConditions {
   
   /* Boundary conditions: */
   Array1D<BoundaryCondition*> bcs;
   
   /* Heat pipes: */
   Array1D<HeatPipe> heat_pipes;
   
   /* Heat-pipe heat flows: */
   Array1D<double> q;
   
   /* Heat-pipe temperatures: */
   Array1D<double> T;
   
};

/* The HeatConductionSolver class: */
class HeatConductionSolver : public PhysicsSolver {
   
   private:
      
      /* Coarse mesh for nodal-level input and output: */
      const Mesh* mesh_nodal = nullptr;
      
      /* Boundary conditions (1-based indexed): */
      Array1D<BoundaryCondition> bcs;
      
      /* Boundary-condition indices for materials: */
      Array1D<int> mat_bc_indices;
      
      /* Total thermal power: */
      Function power;
      
      /* Number of heat pipes: */
      int num_heat_pipes = 0;
      
      /* Heat-pipe boundary conditions: */
      HeatPipeBoundaryConditions hp_bcs;
      
      /* Heat-pipe indices for boundary conditions: */
      Array1D<int> bc_hp_indices;
      
      /* Switch for nonlinear problems where the thermal properties depend on the temperature: */
      bool nonlinear = false;
      
      /* Coefficient matrix for the linear system A*x = b: */
      Mat A = 0;
      
      /* Right-hand side for the linear system A*x = b: */
      Vec b = 0;
      
      /* Volumetric heat source: */
      Vec q = 0;
      
      /* Temperature: */
      Vec T = 0, T0 = 0, T_prev = 0;
      ConvergenceError dT = ConvergenceError("temperature", NORM_2, false, 1.0);
      
      /* Nodal volumetric heat source: */
      Vec qnodal = 0;
      
      /* Nodal temperatures: */
      Array1D<Vec> Tnodal, Tnodal_max;
      
      /* Total heat flows: */
      double qin = -1.0, qout = -1.0;
      
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
      
      /* Calculate the total heat flows in and out of the system: */
      int WARN_UNUSED calculateHeatFlows(double t);
      
      /* Build the coefficient matrix and the RHS vector: */
      int WARN_UNUSED buildMatrix(int n, double dt, double t);
      
      /* Print the solution summary to standard output: */
      int WARN_UNUSED printLog(int n = 0) const;
      
      /* Write the solution to a plain-text file in .vtk format: */
      int WARN_UNUSED writeVTK(const std::string& path, int n = 0) const;
      
      /* Write the solution to a binary file in PETSc format: */
      int WARN_UNUSED writePETSc(int n = 0) const;
   
   public:
      
      /* The HeatConductionSolver constructor: */
      HeatConductionSolver(const Mesh* mesh, const Mesh* mesh_nodal, 
         const Array1D<Material*>& materials) : PhysicsSolver("conduction", mesh, materials), 
         mesh_nodal(mesh_nodal), mat_bc_indices{materials.size(), -1} {}
      
      /* The HeatConductionSolver destructor: */
      ~HeatConductionSolver() {}
      
      /* Read the solver from a plain-text input file: */
      int WARN_UNUSED read(std::ifstream& file, Array1D<Solver*>& solvers);
      
      /* Solve the linear system to get the solution: */
      int WARN_UNUSED solve(int n = 0, double dt = 0.0, double t = 0.0);
   
};
