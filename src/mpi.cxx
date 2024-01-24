#include "mpi.hxx"

/* MPI size and rank: */
namespace mpi {
   int size = 0, rank = 0;
}

/* Initialize MPI: */
int mpi::initialize(int argc, char* argv[]) {
   
   /* Initialize MPI: */
   MPI_CALL(MPI_Init(&argc, &argv));
   
   /* Get the MPI size and rank: */
   MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &size));
   MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
   
   return 0;
   
}

/* Finalize MPI: */
int mpi::finalize() {
   
   /* Finalize MPI: */
   MPI_CALL(MPI_Finalize());
   
   return 0;
   
}
