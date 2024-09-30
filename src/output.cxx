#include "output.hxx"

/* The output namespace: */
namespace output {
   
   /* Switch for verbose output: */
   bool verbose = false;
   
   /* Switch for no output: */
   bool silent = false;
   
   /* Number of padding spaces for indented output: */
   int padding = 0;
   
}

/* Initialize: */
int output::initialize() {
   
   /* Get the switch for verbose output: */
   PAMPA_CHECK(petsc::get("-verbose", verbose), "unable to get the 'verbose' switch");
   
   /* Get the switch for no output: */
   PAMPA_CHECK(petsc::get("-silent", silent), "unable to get the 'silent' switch");
   
   return 0;
   
}

/* Print a message to standard output from the main rank: */
void output::print(const std::string& message, bool info) {
   
   /* Print only from the main rank: */
   if (mpi::rank == 0 && ((!info || verbose) && !silent)) {
      std::cout << std::string(3*padding, ' ');
      std::cout << message << std::endl;
   }
   
}

/* Print a variable to standard output from the main rank: */
void output::print(const std::string& name, PetscScalar x, bool scientific, int precision, 
   bool info) {
   
   /* Print only from the main rank: */
   if (mpi::rank == 0 && ((!info || verbose) && !silent)) {
      if (scientific)
         std::cout << std::scientific;
      else
         std::cout << std::fixed;
      std::cout << std::setprecision(precision);
      std::cout << std::string(3*padding, ' ');
      std::cout << name << ": " << x << "." << std::endl;
   }
   
}

/* Print a variable to standard output from the main rank: */
void output::print(const std::string& name, PetscScalar min, PetscScalar max, bool scientific, 
   int precision, bool info) {
   
   /* Print only from the main rank: */
   if (mpi::rank == 0 && ((!info || verbose) && !silent)) {
      if (scientific)
         std::cout << std::scientific;
      else
         std::cout << std::fixed;
      std::cout << std::setprecision(precision);
      std::cout << std::string(3*padding, ' ');
      std::cout << name << ": " << min << " (min), " << max << " (max)." << std::endl;
   }
   
}

/* Indent: */
void output::indent(bool info) {
   
   /* Indent only from the main rank: */
   if (mpi::rank == 0 && ((!info || verbose) && !silent)) {
      padding++;
   }
   
}

/* Outdent: */
void output::outdent(bool info) {
   
   /* Outdent only from the main rank: */
   if (mpi::rank == 0 && ((!info || verbose) && !silent)) {
      padding--;
   }
   
}
