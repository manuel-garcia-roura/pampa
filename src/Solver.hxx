#pragma once

#include <vector>
#include <iostream>

#include <slepceps.h>

#include "Model.hxx"
#include "Mesh.hxx"
#include "math.hxx"
#include "utils.hxx"

/* The Solver class: */
class Solver {
   
   private:
      
      /* Coefficient matrices for the generalized eigensystem R*phi = (1/keff)*F*phi: */
      Mat R, F;
      
      /* Neutron flux (eigenvector): */
      Vec phi;
      
      /* Multiplication factor (eigenvalue): */
      double keff;
      
      /* Eigenvalue Problem Solver (EPS) context: */
      EPS eps;
      
      /* Build the coefficient matrices: */
      int build_matrices(const Model &model);
   
   public:
      
      /* The Solver constructor: */
      Solver();
      
      /* The Solver destructor: */
      ~Solver();
      
      /* Initialize: */
      int initialize(int argc, char* argv[], const Model &model);
      
      /* Finalize: */
      int finalize();
      
      /* Solve the eigensystem to get the neutron flux and the multiplication factor: */
      int solve();
   
};
