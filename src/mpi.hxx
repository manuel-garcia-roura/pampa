#pragma once

#include <string>

#include <mpi.h>
#include <petsc.h>

#include "utils.hxx"

/* The mpi namespace: */
namespace mpi {
   
   /* MPI size and rank: */
   extern int size, rank;
   
   /* Switch for verbose output: */
   extern bool verbose;
   
   /* Initialize MPI: */
   int WARN_UNUSED initialize(int argc, char* argv[]);
   
   /* Finalize MPI: */
   int WARN_UNUSED finalize();
   
   /* Get the path to the rank directory in parallel runs: */
   std::string get_path(const std::string& filename);
   
   /* Print a message to standard output from the main rank: */
   void print(const std::string& message, bool info = false);
   
   /* Print a variable to standard output from the main rank: */
   void print(const std::string& name, PetscScalar x);
   
}
