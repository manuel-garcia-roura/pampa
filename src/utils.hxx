#pragma once

#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <limits.h>
#include <float.h>
#include <sys/stat.h>

#include <mpi.h>

#include <petsc.h>
#include <petscksp.h>
#include <slepc.h>
#include <slepceps.h>

#include "Array1D.hxx"
#include "Array2D.hxx"
#include "Array3D.hxx"
#include "Vector2D.hxx"
#include "Vector3D.hxx"
#include "Function.hxx"

/* Output precision for .vtk files: */
#define VTK_PRECISION 6

/* Tolerance to compare doubles: */
#define DBL_TOL 1.0e-6

/* Check that error signals returned by functions are not ignored: */
#define WARN_UNUSED __attribute__((warn_unused_result))

/* Check for errors with a condition: */
#define PAMPA_CHECK(condition, message) { \
   if (condition) { \
      std::cout << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << "(): "; \
      std::cout << "error: " << message << "." << std::endl; \
      return 1; \
   } \
}

/* Check for errors with a condition and exit directly if there's an error: */
#define PAMPA_CHECK_EXIT(condition, message) { \
   if (condition) { \
      std::cout << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << "(): "; \
      std::cout << "error: " << message << "." << std::endl; \
      std::exit(1); \
   } \
}

/* Check for calls to virtual functions: */
#define PAMPA_CHECK_VIRTUAL { \
   PAMPA_CHECK(true, "virtual method called on the base class"); return 1; \
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
   enum Type {NONE, VACUUM, REFLECTIVE, ROBIN, DIRICHLET, ADIABATIC, CONVECTION};
}

/* The BoundaryCondition struct: */
struct BoundaryCondition {
   
   /* Boundary-condition type: */
   BC::Type type = BC::NONE;
   
   /* Boundary-condition parameters: */
   /*    - Robin boundary conditions: albedo factor. */
   /*    - Dirichlet boundary conditions: fixed value. */
   /*    - Convection boundary conditions: heat-transfer coefficient and free-stream temperature. */
   Array1D<Function> f;
   
};

/* The utils namespace: */
namespace utils {
   
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
   int WARN_UNUSED find(const std::string& name, const Array1D<T*>& v, int& index) {
      
      /* Find the object: */
      bool found = false;
      for (int i = 0; i < v.size(); i++) {
         if (v(i)->name == name) {
            PAMPA_CHECK(found, "duplicated name '" + name + "'");
            index = i;
            found = true;
         }
      }
      PAMPA_CHECK(!found, "undefined name '" + name + "'");
      
      return 0;
      
   }
   
   /* Find an object by its name: */
   template <typename T>
   int WARN_UNUSED find(const std::string& name, const Array1D<T*>& v, T** x) {
      
      /* Find the object: */
      bool found = false;
      for (int i = 0; i < v.size(); i++) {
         if (v(i)->name == name) {
            PAMPA_CHECK(found, "duplicated name '" + name + "'");
            *x = v(i);
            found = true;
         }
      }
      PAMPA_CHECK(!found, "undefined name '" + name + "'");
      
      return 0;
      
   }
   
   /* Find an object by its name: */
   template <typename T>
   int WARN_UNUSED find(const std::string& name, const Array1D<T>& v, T& x) {
      
      /* Find the object: */
      bool found = false;
      for (int i = 0; i < v.size(); i++) {
         if (v(i).name == name) {
            PAMPA_CHECK(found, "duplicated name '" + name + "'");
            x = v(i);
            found = true;
         }
      }
      PAMPA_CHECK(!found, "undefined name '" + name + "'");
      
      return 0;
      
   }
   
}
