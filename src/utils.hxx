#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits.h>
#include <float.h>
#include <sys/stat.h>

#ifdef WITH_METIS
#include <metis.h>
#endif

#include "Array1D.hxx"
#include "Array2D.hxx"
#include "Array3D.hxx"
#include "Vector2D.hxx"
#include "Vector3D.hxx"
#include "Function.hxx"

/* Check that error signals returned by functions are not ignored: */
#define WARN_UNUSED __attribute__((warn_unused_result))

/* Check for errors with a condition: */
#define PAMPA_CHECK(condition, error, message) { \
   if (condition) { \
      std::cout << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << "(): "; \
      std::cout << "error (" << error << "): " << message << "." << std::endl; \
      return 1; \
   } \
}

/* Check for calls to virtual functions: */
#define PAMPA_CHECK_VIRTUAL { \
   PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1; \
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

/* Check for errors in METIS calls: */
#ifdef WITH_METIS
#define METIS_CALL(function) { \
   int error = function; \
   if (error != METIS_OK) { \
      std::cout << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << "(): "; \
      if (error == METIS_ERROR_INPUT) \
         std::cout << "METIS error (METIS_ERROR_INPUT, " << error << ")." << std::endl; \
      else if (error == METIS_ERROR_MEMORY) \
         std::cout << "METIS error (METIS_ERROR_MEMORY, " << error << ")." << std::endl; \
      else if (error == METIS_ERROR) \
         std::cout << "METIS error (METIS_ERROR, " << error << ")." << std::endl; \
      else \
         std::cout << "METIS error (" << error << ")." << std::endl; \
      return 1; \
   } \
}
#endif

/* The BC::Type enum: */
namespace BC {
   enum Type {VACUUM, REFLECTIVE, ROBIN, DIRICHLET};
}

/* The BoundaryCondition struct: */
struct BoundaryCondition {
   
   /* Boundary condition type: */
   BC::Type type = BC::REFLECTIVE;
   
   /* Albedo factor for Robin boundary conditions: */
   double a = 0.0;
   
   /* Fixed value for Dirichlet boundary conditions: */
   Function x = 0.0;
   
};

/* The utils namespace: */
namespace utils {
   
   /* Trim and remove tabs and double spaces from a string: */
   void clean(std::string& s);
   
   /* Get the next line from a file stream: */
   std::vector<std::string> get_next_line(std::ifstream& file);
   
   /* Read an array with n elements of type double from a file stream: */
   int WARN_UNUSED read(Array1D<double>& v, unsigned int n, std::ifstream& file);
   
   /* Read an array with n elements of type int from a file stream: */
   int WARN_UNUSED read(Array1D<int>& v, unsigned int n, std::ifstream& file);
   
   /* Read an array with (n, m) elements of type int from a file stream: */
   int WARN_UNUSED read(Array2D<double>& v, unsigned int n, unsigned int m, std::ifstream& file);
   
   /* Read a vector with n rows of type double and total size nt from a file stream: */
   int WARN_UNUSED read(Vector2D<double>& v, unsigned int n, unsigned int nt, std::ifstream& file);
   
   /* Read a vector with n rows of type int and total size nt from a file stream: */
   int WARN_UNUSED read(Vector2D<int>& v, unsigned int n, unsigned int nt, std::ifstream& file);
   
   /* Read a vector with (n1, n2, n3) elements of type double from a file stream: */
   int WARN_UNUSED read(Vector3D<double>& v, unsigned int n1, const Array1D<int>& n2, 
      unsigned int n3, std::ifstream& file);
   
   /* Read a bool value from a string: */
   int WARN_UNUSED read(bool& q, const std::string& s);
   
   /* Read an int value from a string: */
   int WARN_UNUSED read(int& x, double x1, double x2, const std::string& s);
   
   /* Read a double value from a string: */
   int WARN_UNUSED read(double& x, double x1, double x2, const std::string& s);
   
   /* Read a function from a file stream: */
   int WARN_UNUSED read(Function& f, const std::vector<std::string>& line, unsigned int& i, 
      std::ifstream& file);
   
   /* Read a boundary condition from a line: */
   int WARN_UNUSED read(BoundaryCondition& bc, const std::vector<std::string>& line, 
      unsigned int& i, std::ifstream& file);
   
   /* Remove a directory: */
   void remove_directory(const std::string& dir);
   
   /* Create a directory: */
   void create_directory(const std::string& dir);
   
   /* Free an object: */
   template <typename T>
   void free(T** ptr) {
      
      /* Delete the pointer, if allocated: */
      if (*ptr != nullptr) delete *ptr;
      
      /* Set the pointer to nullptr: */
      *ptr = nullptr;
      
   }
   
   /* Find an object by its name: */
   template <typename T>
   int find(const std::string& name, const Array1D<T*>& v, T** x) {
      
      /* Find the object: */
      bool found = false;
      for (int i = 0; i < v.size(); i++) {
         if (v(i)->name == name) {
            PAMPA_CHECK(found, 1, "duplicated name '" + name + "'");
            *x = v(i);
            found = true;
         }
      }
      PAMPA_CHECK(!found, 2, "undefined name '" + name + "'");
      
      return 0;
      
   }
   
}
