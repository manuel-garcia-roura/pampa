#pragma once

#include <mpi.h>

#include "utils.hxx"

/* The mpi namespace: */
namespace mpi {
   
   /* MPI size and rank: */
   extern int size, rank;
   
   /* Initialize MPI: */
   int initialize(int argc, char* argv[]);
   
   /* Finalize MPI: */
   int finalize();
   
}
