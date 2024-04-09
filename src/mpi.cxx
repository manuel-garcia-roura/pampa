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

/* Get the path to the rank directory in parallel runs: */
std::string mpi::get_path(const std::string& filename) {
   
   /* Get the path depending on the size and rank: */
   if (size > 1) {
      std::string dir = std::to_string(rank);
      utils::create_directory(dir);
      return dir + "/" + filename;
   }
   else
      return filename;
   
}

/* Print a message to standard output from the main rank: */
void mpi::print(const std::string& message) {
   
   /* Print only from the main rank: */
   if (mpi::rank == 0)
      std::cout << message << std::endl;
   
}

/* Print a variable to standard output from the main rank: */
void mpi::print(const std::string& name, PetscScalar x) {
   
   /* Print only from the main rank: */
   if (mpi::rank == 0)
      std::cout << name << " = " << x << std::endl;
   
}
