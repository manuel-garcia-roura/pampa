#pragma once

#include <cmath>
#include <fstream>
#include <iostream>

#include <slepceps.h>

#include "Mesh.hxx"
#include "Material.hxx"
#include "mpi.hxx"
#include "petsc.hxx"
#include "math.hxx"
#include "utils.hxx"

/* The Solver class: */
class Solver {
   
   protected:
      
      /* Mesh: */
      const Mesh* mesh;
      
      /* Materials: */
      const Array1D<Material>& materials;
      
      /* Transport method: */
      const TransportMethod method;
      
      /* Coefficient matrices for the generalized eigensystem R*x = (1/keff)*F*x: */
      Mat R, F;
      
      /* Scalar neutron flux (eigenvector): */
      Vec phi;
      
      /* Multiplication factor (eigenvalue): */
      double keff;
      
      /* Eigenvalue Problem Solver (EPS) context: */
      EPS eps;
      
      /* Get the flat index for cell i and group g: */
      int index(int i, int g, int ng) const {return i*ng + g;}
      
      /* Get the flat index for cell i, group g and direction m: */
      int index(int i, int g, int m, int ng, int nm) const {return i*nm*ng + g*nm + m;}
      
      /* Build the coefficient matrices and solution vectors: */
      virtual int build() 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Get the solution after solving the eigensystem: */
      virtual int getSolution() 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Write the solution to a plain-text file in .vtk format: */
      virtual int writeVTK(const std::string& filename) const 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Write the solution to a binary file in PETSc format: */
      virtual int writePETSc(const std::string& filename) const 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Destroy the solution vectors: */
      virtual int destroyVectors() 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Normalize the scalar flux: */
      int normalizeScalarFlux();
   
   public:
      
      /* The Solver constructor: */
      Solver(const Mesh* mesh, const Array1D<Material>& materials, const TransportMethod& method) : 
         mesh(mesh), materials(materials), method(method) {}
      
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
