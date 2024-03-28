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
      
      /* Number of directions: */
      int num_directions = -1;
      
      /* Angular quadrature set: */
      AngularQuadratureSet quadrature;
      
      /* Angular neutron flux (eigenvector): */
      Vec psi = 0;
      
      /* Cell-to-cell coupling coefficients for the gradient-discretization scheme: */
      Vector3D<double> grad_coefs, grad_coefs_bc;
      
      /* Mapping from cell indices to boundary-cell indices: */
      Array1D<int> ic_to_ibc;
      
      /* Get the flat index for cell i, group g and direction m: */
      int index(int i, int g, int m) const 
         {return i*num_directions*num_energy_groups + g*num_directions + m;}
      
      /* Check the material data: */
      int WARN_UNUSED checkMaterials() const;
      
      /* Build the coefficient matrices, the solution vectors and the EPS context: */
      int WARN_UNUSED build();
      
      /* Get the mapping and the number of faces for boundary cells: */
      int WARN_UNUSED getBoundaryCells(Array1D<int>& num_faces_bc);
      
      /* Build the coefficients for the Gauss gradient-discretization scheme: */
      int WARN_UNUSED buildGaussGradientScheme(Vector3D<double>& coefs, bool bc);
      
      /* Build the coefficients for the least-squares gradient-discretization scheme: */
      int WARN_UNUSED buildLSGradientScheme(Vector3D<double>& coefs, bool bc);
      
      /* Build the coefficient matrices: */
      int WARN_UNUSED buildMatrices();
      
      /* Get the solution after solving the eigensystem: */
      int WARN_UNUSED getSolution();
      
      /* Calculate the scalar flux: */
      int WARN_UNUSED calculateScalarFlux();
      
      /* Normalize the angular flux: */
      int WARN_UNUSED normalizeAngularFlux();
      
      /* Write the solution to a plain-text file in .vtk format: */
      int WARN_UNUSED writeVTK(const std::string& filename) const;
      
      /* Write the solution to a binary file in PETSc format: */
      int WARN_UNUSED writePETSc() const;
   
   public:
      
      /* The SNSolver constructor: */
      SNSolver(const Mesh* mesh, const Array1D<Material>& materials, int num_energy_groups, 
         int order, double face_interpolation_delta, bool boundary_interpolation_ls) : 
         NeutronicSolver(mesh, materials, num_energy_groups), order(order), 
         face_interpolation_delta(face_interpolation_delta), 
         boundary_interpolation_ls(boundary_interpolation_ls) {}
      
      /* The SNSolver destructor: */
      ~SNSolver() {}
   
};
