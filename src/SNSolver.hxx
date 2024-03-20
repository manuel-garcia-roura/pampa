#pragma once

#include <cmath>
#include <fstream>
#include <iostream>

#include <slepceps.h>
#include <Eigen/Dense>

#include "Solver.hxx"
#include "Mesh.hxx"
#include "Material.hxx"
#include "AngularQuadratureSet.hxx"
#include "mpi.hxx"
#include "petsc.hxx"
#include "math.hxx"
#include "utils.hxx"

/* The SNSolver class: */
class SNSolver : public Solver {
   
   private:
      
      /* Gradient discretization and interpolation scheme: */
      GradientScheme gradient;
      
      /* Angular quadrature set: */
      AngularQuadratureSet quadrature;
      
      /* Angular neutron flux (eigenvector): */
      Vec psi;
      
      /* Cell-to-cell coupling coefficients for the gradient-discretization scheme: */
      Vector3D<double> grad_coefs, grad_coefs_bc;
      
      /* Mapping from cell indices to boundary-cell indices: */
      Array1D<int> ic_to_ibc;
      
      /* Check the material data: */
      int checkMaterials();
      
      /* Build the coefficient matrices and solution vectors: */
      int build();
      
      /* Get the mapping and the number of faces for boundary cells: */
      int getBoundaryCells(Array1D<int>& num_faces_bc);
      
      /* Build the coefficients for the Gauss gradient-discretization scheme: */
      int buildGaussGradientScheme(Vector3D<double>& coefs, bool bc);
      
      /* Build the coefficients for the least-squares gradient-discretization scheme: */
      int buildLSGradientScheme(Vector3D<double>& coefs, bool bc);
      
      /* Build the coefficient matrices: */
      int buildMatrices();
      
      /* Build the solution vectors: */
      int buildVectors();
      
      /* Get the solution after solving the eigensystem: */
      int getSolution();
      
      /* Calculate the scalar flux: */
      int calculateScalarFlux();
      
      /* Normalize the angular flux: */
      int normalizeAngularFlux();
      
      /* Write the solution to a plain-text file in .vtk format: */
      int writeVTK(const std::string& filename) const;
      
      /* Write the solution to a binary file in PETSc format: */
      int writePETSc(const std::string& filename) const;
      
      /* Destroy the solution vectors: */
      int destroyVectors();
   
   public:
      
      /* The SNSolver constructor: */
      SNSolver(const Mesh* mesh, const Array1D<Material>& materials, const TransportMethod& method, 
         const GradientScheme& gradient) : gradient(gradient), Solver(mesh, materials, method) {}
      
      /* The SNSolver destructor: */
      ~SNSolver() {}
   
};
