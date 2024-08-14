#pragma once

#include <string>
#include <iostream>

#include "Parser.hxx"
#include "Mesh.hxx"
#include "Material.hxx"
#include "Solver.hxx"
#include "mpi.hxx"
#include "petsc.hxx"
#include "utils.hxx"

/* The Pampa class: */
class Pampa {
   
   private:
      
      /* Model mesh: */
      Mesh* mesh = nullptr;
      
      /* Model materials: */
      Array1D<Material*> materials;
      
      /* Calculation solvers and main solver to drive the calculation: */
      Array1D<Solver*> solvers;
      Solver* solver;
      
      /* Get the main solver: */
      int WARN_UNUSED getMainSolver();
   
   public:
      
      /* The Pampa constructor: */
      Pampa() {}
      
      /* The Pampa destructor: */
      ~Pampa() {}
      
      /* Initialize the calculation: */
      int WARN_UNUSED initialize(int argc, char* argv[], Array1D<double>& dt);
      
      /* Get the solution: */
      int WARN_UNUSED solve(int n = 0, double dt = 0.0, double t = 0.0);
      
      /* Finalize the calculation: */
      int WARN_UNUSED finalize();
   
};
