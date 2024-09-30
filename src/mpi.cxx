#include "mpi.hxx"

/* The mpi namespace: */
namespace mpi {
   
   /* MPI size and rank: */
   int size = 0, rank = 0;
   
}

/* Initialize: */
int mpi::initialize(int argc, char* argv[]) {
   
   /* Initialize MPI: */
   MPI_CALL(MPI_Init(&argc, &argv));
   
   /* Get the MPI size and rank: */
   MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &size));
   MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
   
   return 0;
   
}

/* Finalize: */
int mpi::finalize() {
   
   /* Finalize MPI: */
   MPI_CALL(MPI_Finalize());
   
   return 0;
   
}

/* Get the path to the rank directory in parallel runs: */
std::string mpi::get_path(const std::string& filename) {
   
   /* Get the path depending on the size and rank: */
   if (size > 1) {
      std::string dir = std::to_string(rank);
      utils::create_directory(dir);
      if (filename.empty())
         return dir;
      else
         return dir + "/" + filename;
   }
   else {
      if (filename.empty())
         return ".";
      else
         return filename;
   }
   
}
