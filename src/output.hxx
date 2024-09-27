#pragma once

#include "petsc.hxx"
#include "mpi.hxx"
#include "utils.hxx"

/* The output namespace: */
namespace output {
   
   /* Switch for verbose output: */
   extern bool verbose;
   
   /* Initialize: */
   int WARN_UNUSED initialize();
   
   /* Finalize: */
   int WARN_UNUSED finalize();
   
   /* Print a message to standard output from the main rank: */
   void print(const std::string& message, bool info = false);
   
   /* Print a variable to standard output from the main rank: */
   void print(const std::string& name, PetscScalar x, bool scientific, int precision, 
      bool info = false);
   
   /* Print a variable to standard output from the main rank: */
   void print(const std::string& name, PetscScalar min, PetscScalar max, bool scientific, 
      int precision, bool info = false);
   
}
