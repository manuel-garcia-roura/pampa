#pragma once

#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

#include <slepceps.h>

#include "Model.hxx"
#include "Mesh.hxx"
#include "mpi.hxx"
#include "math.hxx"
#include "utils.hxx"

/* The Solver class: */
class Solver {
   
   private:
      
      /* Coefficient matrices for the generalized eigensystem R*x = (1/keff)*F*x: */
      Mat R, F;
      
      /* Neutron flux (eigenvector): */
      Vec phi_mpi, psi_mpi, phi_seq, psi_seq;
      
      /* Multiplication factor (eigenvalue): */
      double keff;
      
      /* Eigenvalue Problem Solver (EPS) context: */
      EPS eps;
      
      /* Build the coefficient matrices for the diffusion method: */
      int buildDiffusionMatrices(const Model& model);
      
      /* Build the coefficient matrices for the SN method: */
      int buildSNMatrices(const Model& model);
      
      /* Gather the solution from all ranks to the master rank: */
      int gatherSolution(const Model& model);
      
      /* Calculate the scalar flux: */
      int calculateScalarFlux(const Model& model);
      
      /* Normalize the flux: */
      int normalizeFlux(const Model& model);
   
   public:
      
      /* The Solver constructor: */
      Solver() {}
      
      /* The Solver destructor: */
      ~Solver() {}
      
      /* Initialize: */
      int initialize(int argc, char* argv[], const Model& model);
      
      /* Solve the eigensystem to get the neutron flux and the multiplication factor: */
      int solve();
      
      /* Output the solution: */
      int output(const std::string& filename, const Model& model);
      
      /* Finalize: */
      int finalize(const Model& model);
   
};
