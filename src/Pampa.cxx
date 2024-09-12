#include "Pampa.hxx"

/* Initialize the calculation: */
int Pampa::initialize(int argc, char* argv[], Array1D<double>& dt) {
   
   /* Get the input file name: */
   PAMPA_CHECK(argc < 2, 1, "missing input file");
   std::string filename(argv[1]);
   
   /* Initialize MPI: */
   PAMPA_CALL(mpi::initialize(argc, argv), "unable to initialize MPI");
   
   /* Initialize PETSc and SLEPc: */
   PAMPA_CALL(petsc::initialize(argc, argv), "unable to initialize PETSc and SLEPc");
   
   /* Read the main input file: */
   Parser parser;
   PAMPA_CALL(parser.read(filename, &mesh, &mesh_nodal, materials, solvers, dt), 
      "unable to parse " + filename);
   
   /* Remove boundary-condition materials from the mesh and swap the meshes: */
   bool remove_materials = false;
   for (int i = 0; i < materials.size(); i++)
      remove_materials |= materials(i)->isBC();
   if (remove_materials) {
      Mesh* mesh_new = nullptr;
      PAMPA_CALL(mesh->removeBCMats(materials, &mesh_new), 
         "unable to remove boundary-condition materials from the mesh");
      utils::free(&mesh);
      mesh = mesh_new;
   }
   
   /* Partition the mesh and swap the meshes: */
   if (mpi::size > 1 && !(mesh->isPartitioned())) {
      Mesh* submesh = nullptr;
      PAMPA_CALL(mesh->partition(&submesh), "unable to partition the mesh");
      utils::free(&mesh);
      mesh = submesh;
   }
   
   /* Set the mesh for all solvers: */
   for (int i = 0; i < solvers.size(); i++)
      solvers(i)->setMesh(mesh);
   
   /* Get the main solver: */
   PAMPA_CALL(getMainSolver(), "unable to get the main solver");
   
   /* Initialize the solver: */
   bool transient = !(dt.empty());
   PAMPA_CALL(solver->initialize(transient), "unable to initialize the solver");
   
   return 0;
   
}

/* Get the solution: */
int Pampa::solve(int n, double dt, double t) {
   
   /* Print the time-step number: */
   mpi::print("--------------------------------");
   mpi::print("n = " + std::to_string(n) + ":");
   
   /* Get the solution: */
   PAMPA_CALL(solver->solve(n, dt, t), "unable to get the solution");
   
   /* Output the solution: */
   std::string filename = "output_" + std::to_string(n) + ".vtk";
   PAMPA_CALL(solver->output(mpi::get_path(filename), n), "unable to output the solution");
   
   return 0;
   
}

/* Finalize the calculation: */
int Pampa::finalize() {
   
   /* Finalize the solver: */
   PAMPA_CALL(solver->finalize(), "unable to finalize the solver");
   
   /* Finalize PETSc and SLEPc: */
   PAMPA_CALL(petsc::finalize(), "unable to finalize PETSc and SLEPc");
   
   /* Finalize MPI: */
   PAMPA_CALL(mpi::finalize(), "unable to finalize MPI");
   
   /* Free the meshes, materials and solvers: */
   utils::free(&mesh);
   utils::free(&mesh_nodal);
   for (int i = 0; i < materials.size(); i++)
      utils::free(&materials(i));
   for (int i = 0; i < solvers.size(); i++)
      utils::free(&solvers(i));
   
   /* Print the last line: */
   mpi::print("--------------------------------");
   
   return 0;
   
}

/* Get the main solver: */
int Pampa::getMainSolver() {
   
   /* Check if there's at least one solver: */
   PAMPA_CHECK(solvers.empty(), 1, "no solvers defined");
   
   /* Look for the main solver: */
   if (solvers.size() == 1)
      solver = solvers(0);
   else {
      PAMPA_CALL(utils::find("main", solvers, &solver), "unable to find the main solver");
   }
   
   return 0;
   
}
