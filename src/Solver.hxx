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
      virtual int WARN_UNUSED checkMaterials() {PAMPA_CHECK_VIRTUAL}
      
      /* Build the matrices, vectors and solver contexts: */
      virtual int WARN_UNUSED build() {PAMPA_CHECK_VIRTUAL}
      
      /* Print the solution summary to standard output: */
      virtual int WARN_UNUSED printLog() const {PAMPA_CHECK_VIRTUAL}
      
      /* Write the solution to a plain-text file in .vtk format: */
      virtual int WARN_UNUSED writeVTK(const std::string& filename) const {PAMPA_CHECK_VIRTUAL}
      
      /* Write the solution to a binary file in PETSc format: */
      virtual int WARN_UNUSED writePETSc() const {PAMPA_CHECK_VIRTUAL}
   
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
      virtual int WARN_UNUSED solve(int n = 0, double dt = 0.0) {PAMPA_CHECK_VIRTUAL}
      
      /* Output the solution: */
      int WARN_UNUSED output(const std::string& filename);
      
      /* Finalize: */
      int WARN_UNUSED finalize();
   
};
