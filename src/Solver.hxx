#pragma once

#include <string>
#include <fstream>
#include <iostream>

#include "Mesh.hxx"
#include "vtk.hxx"
#include "utils.hxx"

/* The Solver class: */
class Solver {
   
   protected:
      
      /* Solver name: */
      const std::string name = "";
      
      /* Mesh: */
      const Mesh* mesh;
      
      /* Total power: */
      Array1D<double> power{1, 1.0};
      
      /* Boundary conditions (1-based indexed): */
      Array1D<BoundaryCondition> bcs;
   
   public:
      
      /* The Solver constructor: */
      Solver(const std::string& name, const Mesh* mesh) : name(name), mesh(mesh) {}
      
      /* The Solver destructor: */
      virtual ~Solver() {}
      
      /* Set the total power: */
      void setPower(const Array1D<double>& power) {this->power = power;}
      
      /* Add a boundary condition: */
      int WARN_UNUSED addBoundaryCondition(const BoundaryCondition& bc, int l);
      
      /* Get the solver name: */
      const std::string& getName() const {return name;}
      
      /* Initialize: */
      virtual int WARN_UNUSED initialize() {PAMPA_CHECK_VIRTUAL}
      
      /* Get the solution: */
      virtual int WARN_UNUSED solve(int n = 0, double dt = 0.0) {PAMPA_CHECK_VIRTUAL}
      
      /* Output the solution: */
      virtual int WARN_UNUSED output(const std::string& filename, int n = 0, 
         bool write_mesh = true) const {PAMPA_CHECK_VIRTUAL}
      
      /* Finalize: */
      virtual int WARN_UNUSED finalize() {PAMPA_CHECK_VIRTUAL}
   
};

/* Find a solver from its name: */
namespace utils {
   int WARN_UNUSED find(const std::string& name, const Array1D<Solver*>& solvers, Solver** solver);
}
