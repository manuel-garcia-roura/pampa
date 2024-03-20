#pragma once

#include <petscksp.h>

#include "Solver.hxx"

/* The HeatConductionSolver class: */
class HeatConductionSolver : public Solver {
   
   public:
      
      /* The HeatConductionSolver constructor: */
      HeatConductionSolver(const Mesh* mesh, const Array1D<Material>& materials) : 
         Solver(mesh, materials) {}
      
      /* The HeatConductionSolver destructor: */
      ~HeatConductionSolver() {}
      
      /* Initialize: */
      int initialize(int argc, char* argv[]) {return 0;}
      
      /* Solve the eigensystem to get the neutron flux and the multiplication factor: */
      int solve() {return 0;}
      
      /* Output the solution: */
      int output(const std::string& filename) {return 0;}
      
      /* Finalize: */
      int finalize() {return 0;}
   
};
