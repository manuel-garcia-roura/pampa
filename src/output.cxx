#include "output.hxx"

/* The output namespace: */
namespace output {
   
   /* Switch for verbose output: */
   bool verbose = false;
   
}

/* Initialize: */
int output::initialize() {
   
   /* Get the switch for verbose output: */
   PAMPA_CHECK(petsc::get_switch("-verbose", verbose), "unable to get the 'verbose' switch");
   
   return 0;
   
}

/* Finalize: */
int output::finalize() {
   
   return 0;
   
}

/* Print a message to standard output from the main rank: */
void output::print(const std::string& message, bool info) {
   
   /* Print only from the main rank: */
   if (mpi::rank == 0 && (!info || verbose))
      std::cout << message << std::endl;
   
}

/* Print a variable to standard output from the main rank: */
void output::print(const std::string& name, PetscScalar x) {
   
   /* Print only from the main rank: */
   if (mpi::rank == 0)
      std::cout << name << " = " << x << std::endl;
   
}
