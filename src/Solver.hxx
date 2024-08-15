#pragma once

#include <string>
#include <fstream>
#include <iostream>

#include "Mesh.hxx"
#include "petsc.hxx"
#include "vtk.hxx"
#include "mpi.hxx"
#include "math.hxx"
#include "utils.hxx"

/* The Field struct: */
struct Field {
   
   /* Field name: */
   std::string name;
   
   /* Pointer to the PETSc vector in the solver: */
   Vec* vec = nullptr;
   
   /* Input/output flags: */
   bool input = false, output = false;
   
   /* Previous iteration used to evaluate convergence: */
   Vec* vec0 = nullptr;
   
};

/* The Solver class: */
class Solver {
   
   public:
      
      /* Solver name: */
      const std::string name;
   
   protected:
      
      /* Mesh: */
      const Mesh* mesh = nullptr;
      
      /* Input/output fields: */
      Array1D<Field> fields;
      
      /* Index of the previous time step: */
      int n0 = -1;
   
   public:
      
      /* The Solver constructor: */
      Solver(const std::string& name, const Mesh* mesh) : name(name), mesh(mesh) {}
      
      /* The Solver destructor: */
      virtual ~Solver() {}
      
      /* Get the input/output fields: */
      Array1D<Field>& getFields() {return fields;}
      
      /* Read the solver from a plain-text input file: */
      int WARN_UNUSED read(const std::string& filename, Array1D<Solver*>& solvers);
      
      /* Read the solver from a plain-text input file: */
      virtual int WARN_UNUSED read(std::ifstream& file, Array1D<Solver*>& solvers) 
         {PAMPA_CHECK_VIRTUAL}
      
      /* Initialize: */
      virtual int WARN_UNUSED initialize(bool transient = false) {PAMPA_CHECK_VIRTUAL}
      
      /* Get the solution: */
      virtual int WARN_UNUSED solve(int n = 0, double dt = 0.0, double t = 0.0) 
         {PAMPA_CHECK_VIRTUAL}
      
      /* Output the solution: */
      virtual int WARN_UNUSED output(const std::string& filename, int n = 0, 
         bool write_mesh = true) const {PAMPA_CHECK_VIRTUAL}
      
      /* Finalize: */
      virtual int WARN_UNUSED finalize() {PAMPA_CHECK_VIRTUAL}
      
      /* Get the values for a given field: */
      int WARN_UNUSED getField(double* v, const std::string& name) const;
      
      /* Set the values for a given field: */
      int WARN_UNUSED setField(const double* v, const std::string& name);
   
};
