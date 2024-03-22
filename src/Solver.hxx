#pragma once

#include <cmath>
#include <fstream>
#include <iostream>

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
   
   public:
      
      /* The Solver constructor: */
      Solver(const Mesh* mesh, const Array1D<Material>& materials) : mesh(mesh), 
         materials(materials) {}
      
      /* The Solver destructor: */
      virtual ~Solver() {}
      
      /* Initialize: */
      virtual int WARN_UNUSED initialize(int argc, char* argv[]) 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Solve the linear system to get the solution: */
      virtual int WARN_UNUSED solve() 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Output the solution: */
      virtual int WARN_UNUSED output(const std::string& filename) 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Finalize: */
      virtual int WARN_UNUSED finalize() 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
   
};
