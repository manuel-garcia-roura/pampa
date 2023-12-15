#pragma once

#include <slepceps.h>

/* The Solver class: */
class Solver {
   
   private:
      
      /* Coefficient matrices for the generalized eigensystem Rx = kFx: */
      Mat R, F;
      
      /* Eigenvalue Problem Solver (EPS): */
      EPS eps;
   
   public:
      
      /* The Solver constructor: */
      Solver();
      
      /* The Solver destructor: */
      ~Solver();
   
};
