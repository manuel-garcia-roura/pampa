#pragma once

#include "petsc.hxx"
#include "mpi.hxx"
#include "utils.hxx"

/* The output namespace: */
namespace output {
   
   /* Switch for verbose output: */
   extern bool verbose;
   
   /* Switch for no output: */
   extern bool silent;
   
   /* Number of padding spaces for indented output: */
   extern int padding;
   
   /* Initialize: */
   int WARN_UNUSED initialize();
   
   /* Print a message to standard output from the main rank: */
   void print(const std::string& message, bool info = false);
   
   /* Print a variable to standard output from the main rank: */
   void print(const std::string& name, PetscScalar x, bool scientific, int precision, 
      bool info = false);
   
   /* Print a variable to standard output from the main rank: */
   void print(const std::string& name, PetscScalar min, PetscScalar max, bool scientific, 
      int precision, bool info = false);
   
   /* Indent: */
   void indent(bool info = false);
   
   /* Outdent: */
   void outdent(bool info = false);
   
}
