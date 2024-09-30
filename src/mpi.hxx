#pragma once

#include "utils.hxx"

/* The mpi namespace: */
namespace mpi {
   
   /* MPI size and rank: */
   extern int size, rank;
   
   /* Initialize: */
   int WARN_UNUSED initialize(int argc, char* argv[]);
   
   /* Finalize: */
   int WARN_UNUSED finalize();
   
   /* Get the path to the rank directory in parallel runs: */
   std::string get_path(const std::string& filename = std::string());
   
}
