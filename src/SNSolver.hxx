#pragma once

#include <Eigen/Dense>

#include "NeutronicSolver.hxx"
#include "AngularQuadratureSet.hxx"

/* The SNSolver class: */
class SNSolver : public NeutronicSolver {
   
   private:
      
      /* Order (N) of the SN method: */
      int order = -1;
      
      /* Weight between upwind and linear schemes for face interpolation: */
      /* Note: 1.0 corresponds to pure upwind and 0.0 to pure linear interpolation. */
      double face_interpolation_delta;
      
      /* Switch to use the least-squares gradient for boundary interpolation: */
      bool boundary_interpolation_ls;
      
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
      
      /* Build the coefficient matrices and the solution vectors: */
      int build();
      
      /* Get the mapping and the number of faces for boundary cells: */
      int getBoundaryCells(Array1D<int>& num_faces_bc);
      
      /* Build the coefficients for the Gauss gradient-discretization scheme: */
      int buildGaussGradientScheme(Vector3D<double>& coefs, bool bc);
      
      /* Build the coefficients for the least-squares gradient-discretization scheme: */
      int buildLSGradientScheme(Vector3D<double>& coefs, bool bc);
      
      /* Build the coefficient matrices: */
      int buildMatrices();
      
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
      SNSolver(const Mesh* mesh, const Array1D<Material>& materials, int num_groups, int order, 
         double face_interpolation_delta, bool boundary_interpolation_ls) : 
         NeutronicSolver(mesh, materials, num_groups), order(order), 
         face_interpolation_delta(face_interpolation_delta), 
         boundary_interpolation_ls(boundary_interpolation_ls) {}
      
      /* The SNSolver destructor: */
      ~SNSolver() {}
   
};
