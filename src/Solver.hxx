#pragma once

#include <string>
#include <fstream>
#include <iostream>

#include "vtk.hxx"
#include "utils.hxx"

/* The Solver class: */
class Solver {
   
   protected:
      
      /* Solver name: */
      const std::string name = "";
      
      /* Total power: */
      Array1D<double> power{1, 1.0};
   
   public:
      
      /* The Solver constructor: */
      Solver(const std::string& name) : name(name) {}
      
      /* The Solver destructor: */
      virtual ~Solver() {}
      
      /* Set the total power: */
      void setPower(const Array1D<double>& power) {this->power = power;}
      
      /* Get the solver name: */
      const std::string& getName() const {return name;}
      
      /* Initialize: */
      virtual int WARN_UNUSED initialize() {PAMPA_CHECK_VIRTUAL}
      
      /* Get the solution: */
      virtual int WARN_UNUSED solve(int n = 0, double dt = 0.0) {PAMPA_CHECK_VIRTUAL}
      
      /* Output the solution: */
      virtual int WARN_UNUSED output(const std::string& filename, int n = 0) const 
         {PAMPA_CHECK_VIRTUAL}
      
      /* Finalize: */
      virtual int WARN_UNUSED finalize() {PAMPA_CHECK_VIRTUAL}
   
};

/* Find a solver from its name: */
namespace utils {
   int WARN_UNUSED find(const std::string& name, const Array1D<Solver*>& solvers, Solver** solver);
}
