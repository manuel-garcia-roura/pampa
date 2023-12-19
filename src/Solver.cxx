#include "Solver.hxx"

/* The Solver constructor: */
Solver::Solver() {};

/* The Solver destructor: */
Solver::~Solver() {};

/* Initialize: */
int Solver::initialize(int argc, char* argv[]) {
   
   /* Initialize SLEPc: */
   static char help[] = "Solver for the generalized eigensystem R*phi = (1/k)*F*phi.\n";
   PETSC_CALL(SlepcInitialize(&argc, &argv, (char*)0, help), "unable to initialize SLEPc");
   
   /* Create the coefficient matrices: */
   PETSC_CALL(MatCreate(PETSC_COMM_WORLD, &R), "unable to create R");
   PETSC_CALL(MatCreate(PETSC_COMM_WORLD, &F), "unable to create F");
   PETSC_CALL(MatSetFromOptions(R), "unable to set up R");
   PETSC_CALL(MatSetFromOptions(F), "unable to set up F");
   
   /* Create the EPS context: */
   PETSC_CALL(EPSCreate(PETSC_COMM_WORLD, &eps), "unable to create EPS");
   
   return 0;
   
};

/* Finalize: */
int Solver::finalize() {
   
   /* Destroy the EPS context: */
   PETSC_CALL(EPSDestroy(&eps), "unable to destroy EPS");
   
   /* Destroy the coefficient matrices: */
   PETSC_CALL(MatDestroy(&R), "unable to create R");
   PETSC_CALL(MatDestroy(&F), "unable to create F");
   
   /* Finalize SLEPc: */
   PETSC_CALL(SlepcFinalize(), "unable to finalize SLEPc");
   
   return 0;
   
};
