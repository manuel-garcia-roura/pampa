#pragma once

#include <cmath>

#include "NeutronicSolver.hxx"
#include "AngularQuadratureSet.hxx"

/* The SNSolver class: */
class SNSolver : public NeutronicSolver {
   
   private:
      
      /* Order (N) of the SN method: */
      int order = -1;
      
      /* Weight between upwind and linear schemes for face interpolation: */
      /* Note: 1.0 corresponds to pure upwind and 0.0 to pure linear interpolation. */
      double face_interpolation_delta = 0.1;
      
      /* Switch to use the least-squares gradient for boundary interpolation: */
      bool boundary_interpolation_ls = false;
      
      /* Number of directions: */
      int num_directions = -1;
      
      /* Angular quadrature set: */
      AngularQuadratureSet quadrature;
      
      /* Angular neutron flux: */
      Vec psi = 0, psi0 = 0;
      
      /* Cell-to-cell coupling coefficients for the gradient-discretization scheme: */
      Vector3D<double> grad_coefs, grad_coefs_bc;
      
      /* Mapping from cell indices to boundary-cell indices: */
      Array1D<int> ic_to_ibc;
      
      /* Get the flat index for cell i, group g and direction m: */
      int index(int i, int g, int m) const 
         {return i*num_directions*num_energy_groups + g*num_directions + m;}
      
      /* Get the mapping and the number of faces for boundary cells: */
      int WARN_UNUSED getBoundaryCells(Array1D<int>& num_faces_bc);
      
      /* Build the coefficients for the Gauss gradient-discretization scheme: */
      int WARN_UNUSED buildGaussGradientScheme(Vector3D<double>& coefs, bool bc);
      
      /* Build the coefficients for the least-squares gradient-discretization scheme: */
      #ifdef WITH_EIGEN
      int WARN_UNUSED buildLSGradientScheme(Vector3D<double>& coefs, bool bc);
      #endif
      
      /* Calculate the scalar flux: */
      int WARN_UNUSED calculateScalarFlux();
      
      /* Normalize the angular flux: */
      int WARN_UNUSED normalizeAngularFlux();
      
      /* Build the coefficient matrices and the RHS vector: */
      int WARN_UNUSED buildMatrices(int n, double dt, double t);
      
      /* Solve the linear system and get the solution: */
      int WARN_UNUSED getSolution(int n = 0);
      
      /* Check the material data: */
      int WARN_UNUSED checkMaterials(bool transient = false);
      
      /* Build the coefficient matrices and the solution vectors: */
      int WARN_UNUSED build();
      
      /* Write the solution to a plain-text file in .vtk format: */
      int WARN_UNUSED writeVTK(const std::string& filename) const;
      
      /* Write the solution to a binary file in PETSc format: */
      int WARN_UNUSED writePETSc(int n = 0) const;
   
   public:
      
      /* The SNSolver constructor: */
      SNSolver(const Mesh* mesh, const Array1D<Material*>& materials) : 
         NeutronicSolver("sn", mesh, materials) {}
      
      /* The SNSolver destructor: */
      ~SNSolver() {}
      
      /* Read the solver from a plain-text input file: */
      int WARN_UNUSED read(std::ifstream& file, Array1D<Solver*>& solvers);
   
};
