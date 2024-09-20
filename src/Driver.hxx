#pragma once

#include "Parser.hxx"
#include "Mesh.hxx"
#include "Material.hxx"
#include "Solver.hxx"
#include "petsc.hxx"
#include "output.hxx"
#include "mpi.hxx"
#include "utils.hxx"

/* The Driver class: */
class Driver {
   
   private:
      
      /* Model mesh: */
      Mesh* mesh = nullptr;
      
      /* Coarse mesh for nodal-level input and output: */
      Mesh* mesh_nodal = nullptr;
      
      /* Model materials: */
      Array1D<Material*> materials;
      
      /* Calculation solvers and main solver to drive the calculation: */
      Array1D<Solver*> solvers;
      Solver* solver;
      
      /* Get the main solver: */
      int WARN_UNUSED getMainSolver();
   
   public:
      
      /* The Driver constructor: */
      Driver() {}
      
      /* The Driver destructor: */
      ~Driver() {}
      
      /* Initialize the calculation: */
      int WARN_UNUSED initialize(int argc, char* argv[], Array1D<double>& dt);
      
      /* Get the solution: */
      int WARN_UNUSED solve(int n = 0, double dt = 0.0, double t = 0.0);
      
      /* Finalize the calculation: */
      int WARN_UNUSED finalize();
      
      /* Get the values for a given field: */
      int WARN_UNUSED getField(double* v, const std::string& name) const 
         {PAMPA_CHECK(solver->getField(v, name), "unable to get the field"); return 0;}
      
      /* Set the values for a given field: */
      int WARN_UNUSED setField(const double* v, const std::string& name) 
         {PAMPA_CHECK(solver->setField(v, name), "unable to set the field"); return 0;}
   
};
