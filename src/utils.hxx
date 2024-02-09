#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "Array1D.hxx"
#include "Array2D.hxx"
#include "Array3D.hxx"
#include "Vector2D.hxx"

/* Check for errors with a condition: */
#define PAMPA_CHECK(condition, error, message) { \
   if (condition) { \
      std::cout << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << "(): "; \
      std::cout << "error (" << error << "): " << message << "." << std::endl; \
      return 1; \
   } \
}

/* Check for errors in general calls: */
#define PAMPA_CALL(function, message) { \
   int error = function; \
   if (error) { \
      std::cout << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << "(): "; \
      std::cout << "error (" << error << "): " << message << "." << std::endl; \
      return 1; \
   } \
}

/* Check for errors in MPI calls: */
#define MPI_CALL(function) { \
   int error = function; \
   if (error) { \
      std::cout << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << "(): "; \
      std::cout << "MPI error (" << error << ")." << std::endl; \
      return 1; \
   } \
}

/* Check for errors in PETSc/SLEPc calls: */
#define PETSC_CALL(function) { \
   int error = function; \
   if (error) { \
      std::cout << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << "(): "; \
      std::cout << "PETSc/SLEPc error (" << error << ")." << std::endl; \
      return 1; \
   } \
}

/* The BC::Type enum: */
namespace BC {
   enum Type {VACUUM, REFLECTIVE, ROBIN};
}

/* The BoundaryCondition struct: */
struct BoundaryCondition {
   
   /* Boundary condition type: */
   BC::Type type = BC::REFLECTIVE;
   
   /* Albedo factor for Robin boundary conditions: */
   double a = 0.0;
   
};

/* The TM::Type enum: */
namespace TM {
   enum Type {DIFFUSION, SN};
}

/* The TransportMethod struct: */
struct TransportMethod {
   
   /* Transport method type: */
   TM::Type type = TM::DIFFUSION;
   
   /* Order (N) of the SN method: */
   int order = -1;
   
   /* Number of energy groups: */
   int num_groups = -1;
   
};

/* The utils namespace: */
namespace utils {
   
   /* Trim and remove tabs and double spaces from a string: */
   void clean(std::string& s);
   
   /* Get the next line from a file stream: */
   std::vector<std::string> get_next_line(std::ifstream& file);
   
   /* Read an array with n elements of type double from a file stream: */
   int read(Array1D<double>& v, int n, std::ifstream& file);
   
   /* Read an array with n elements of type int from a file stream: */
   int read(Array1D<int>& v, int n, std::ifstream& file);
   
   /* Read an array with (n, m) elements of type int from a file stream: */
   int read(Array2D<double>& v, int n, int m, std::ifstream& file);
   
   /* Read a vector with n rows of type int and unknown size from a file stream: */
   int read(Vector2D<int>& v, int n, std::ifstream& file);
   
   /* Read a boundary condition from a line: */
   int read(BoundaryCondition& bc, const std::vector<std::string>& line, int& i);
   
}
