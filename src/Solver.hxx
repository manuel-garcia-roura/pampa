#pragma once

#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

#include <slepceps.h>

#include "Mesh.hxx"
#include "Material.hxx"
#include "AngularQuadratureSet.hxx"
#include "mpi.hxx"
#include "math.hxx"
#include "utils.hxx"

/* The Solver class: */
class Solver {
   
   private:
      
      /* Mesh: */
      const Mesh* mesh;
      
      /* Materials: */
      const std::vector<Material>& materials;
      
      /* Transport method: */
      TransportMethod method;
      
      /* Angular quadrature set: */
      AngularQuadratureSet quadrature;
      
      /* Coefficient matrices for the generalized eigensystem R*x = (1/keff)*F*x: */
      Mat R, F;
      
      /* Neutron flux (eigenvector): */
      Vec phi_mpi, psi_mpi, phi_seq, psi_seq;
      
      /* Multiplication factor (eigenvalue): */
      double keff;
      
      /* Eigenvalue Problem Solver (EPS) context: */
      EPS eps;
      
      /* Build the coefficient matrices for the diffusion method: */
      int buildDiffusionMatrices();
      
      /* Build the coefficient matrices for the SN method: */
      int buildSNMatrices();
      
      /* Gather the solution from all ranks to the master rank: */
      int gatherSolution();
      
      /* Calculate the scalar flux: */
      int calculateScalarFlux();
      
      /* Normalize the flux: */
      int normalizeFlux();
   
   public:
      
      /* The Solver constructor: */
      Solver(const Mesh* mesh, const std::vector<Material>& materials, 
         const TransportMethod& method) : mesh(mesh), materials(materials), method(method) {}
      
      /* The Solver destructor: */
      ~Solver() {}
      
      /* Initialize: */
      int initialize(int argc, char* argv[]);
      
      /* Solve the eigensystem to get the neutron flux and the multiplication factor: */
      int solve();
      
      /* Output the solution: */
      int output(const std::string& filename);
      
      /* Finalize: */
      int finalize();
   
};
