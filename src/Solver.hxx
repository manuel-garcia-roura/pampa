#pragma once

#include <iostream>

#include <slepceps.h>

#include "utils.hxx"

/* The Solver class: */
class Solver {
   
   private:
      
      /* Coefficient matrices for the generalized eigensystem R*phi = (1/k)*F*phi: */
      Mat R, F;
      
      /* Eigenvalue Problem Solver (EPS) context: */
      EPS eps;
   
   public:
      
      /* The Solver constructor: */
      Solver();
      
      /* The Solver destructor: */
      ~Solver();
      
      /* Initialize: */
      int initialize(int argc, char* argv[]);
      
      /* Finalize: */
      int finalize();
   
};
