#pragma once

#include <cmath>
#include <string>
#include <fstream>
#include <iostream>

#include "Mesh.hxx"
#include "Material.hxx"
#include "mpi.hxx"
#include "petsc.hxx"
#include "math.hxx"
#include "vtk.hxx"
#include "utils.hxx"

/* The Solver class: */
class Solver {
   
   protected:
      
      /* Mesh: */
      const Mesh* mesh;
      
      /* Materials: */
      const Array1D<Material>& materials;
      
      /* Pointers to the PETSc matrices: */
      Array1D<Mat*> matrices;
      
      /* Pointers to the PETSc vectors: */
      Array1D<Vec*> vectors;
      
      /* Krylov Subspace Solver (KSP) context: */
      KSP ksp = 0;
      
      /* Eigenvalue Problem Solver (EPS) context: */
      EPS eps = 0;
      
      /* Total power: */
      Array1D<double> power;
      
      /* Check the material data: */
      virtual int WARN_UNUSED checkMaterials() 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Build the matrices, vectors and solver contexts: */
      virtual int WARN_UNUSED build() 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
   
   public:
      
      /* The Solver constructor: */
      Solver(const Mesh* mesh, const Array1D<Material>& materials) : mesh(mesh), 
         materials(materials) {}
      
      /* The Solver destructor: */
      virtual ~Solver() {}
      
      /* Set the total power: */
      void setPower(const Array1D<double>& power) {this->power = power;}
      
      /* Initialize: */
      int WARN_UNUSED initialize();
      
      /* Get the solution: */
      virtual int WARN_UNUSED solve(int n = 0, double dt = 0.0) 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Output the solution: */
      virtual int WARN_UNUSED output(const std::string& filename) 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Finalize: */
      int WARN_UNUSED finalize();
   
};
