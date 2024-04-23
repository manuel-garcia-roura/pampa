#include "Parser.hxx"

/* Read a plain-text input file: */
int Parser::read(const std::string& filename, Mesh** mesh, Array1D<Material>& materials, 
   Array1D<Solver*>& solvers, Array1D<double>& dt) {
   
   /* Open the input file: */
   std::ifstream file(filename);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty()) break;
      
      /* Get the next keyword: */
      if (line[0] == "mesh") {
         
         /* Create the mesh depending on the mesh type: */
         std::string mesh_type = line[1];
         if (mesh_type == "cartesian")
            *mesh = new CartesianMesh();
         else if (mesh_type == "unstructured")
            *mesh = new UnstructuredExtrudedMesh();
         else if (mesh_type == "partitioned")
            *mesh = new PartitionedMesh();
         else {
            PAMPA_CHECK(true, 2, "wrong mesh type");
         }
         
         /* Read the mesh: */
         std::string mesh_filename = line[2];
         PAMPA_CALL((*mesh)->read(mesh_filename), "unable to read the mesh from " + mesh_filename);
         
         /* Build the mesh: */
         PAMPA_CALL((*mesh)->build(), "unable to build the mesh");
         
         /* Partition the mesh and swap the meshes: */
         if (mpi::size > 1 && !(*mesh)->isPartitioned()) {
            Mesh* submesh = NULL;
            PAMPA_CALL((*mesh)->partition(&submesh), "unable to partition the mesh");
            delete *mesh;
            *mesh = submesh;
         }
         
      }
      else if (line[0] == "material") {
         
         /* Get the material name: */
         std::string name = line[1];
         
         /* Create the material: */
         Material mat(name);
         
         /* Read the material: */
         if (line[2] == "{") {
            PAMPA_CALL(mat.read(file), "unable to read the material from " + filename);
         }
         else {
            std::string mat_filename = line[2];
            PAMPA_CALL(mat.read(mat_filename), "unable to read the material from " + mat_filename);
         }
         
         /* Keep the material definition: */
         materials.pushBack(mat);
         
      }
      else if (line[0] == "solver") {
         
         /* Create the solver depending on the solver type: */
         Solver* solver = NULL;
         unsigned int i = 1;
         std::string solver_type = line[i++];
         if (solver_type == "diffusion")
            solver = new DiffusionSolver(*mesh, materials);
         else if (solver_type == "sn")
            solver = new SNSolver(*mesh, materials);
         else if (solver_type == "conduction")
            solver = new HeatConductionSolver(*mesh, materials);
         else if (solver_type == "precursors")
            solver = new PrecursorSolver(*mesh, materials);
         else if (solver_type == "coupled") {
            std::string name = line[i++];
            solver = new CouplingSolver(name, *mesh);
         }
         else {
            PAMPA_CHECK(true, 3, "wrong solver type");
         }
         
         /* Read the solver: */
         if (i < line.size()) {
            std::string s = line[i++];
            if (s == "{") {
               PAMPA_CALL(solver->read(file, solvers), "unable to read the solver from " + filename);
            }
            else {
               std::string solver_filename = s;
               PAMPA_CALL(solver->read(solver_filename, solvers), 
                  "unable to read the solver from " + solver_filename);
            }
         }
         
         /* Keep the solver definition: */
         solvers.pushBack(solver);
         
      }
      else if (line[0] == "dt") {
         
         /* Get the dt values: */
         int nt;
         PAMPA_CALL(utils::read(nt, -INT_MAX, INT_MAX, line[1]), "wrong number of time steps");
         if (nt > 0) {
            PAMPA_CALL(utils::read(dt, nt, file), "wrong dt data");
         }
         else {
            PAMPA_CALL(utils::read(dt, 1, file), "wrong dt data");
            nt = -nt;
            dt.resize(nt, dt(0));
         }
         
      }
      else if (line[0] == "bc") {
         
         /* Get the solver name: */
         std::string name = line[1];
         
         /* Get the solver: */
         Solver* solver;
         PAMPA_CALL(utils::find(name, solvers, &solver), "unable to find solver");
         
         /* Get the boundary condition (1-based indexed): */
         int l, i = 2;
         BoundaryCondition bc;
         PAMPA_CALL(utils::read(l, 1, INT_MAX, line[i++]), "wrong boundary condition index");
         PAMPA_CALL(utils::read(bc, line, i), "wrong boundary condition");
         
         /* Add the boundary condition to the solver: */
         PAMPA_CALL(solver->addBoundaryCondition(bc, l), "unable to add the boundary condition");
         
      }
      else if (line[0] == "include") {
         
         /* Read an included input file: */
         std::string include_filename = line[1];
         PAMPA_CALL(read(include_filename, mesh, materials, solvers, dt), 
            "unable to parse " + include_filename);
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 5, "unrecognized keyword '" + line[0] + "'");
         
      }
      
   }
   
   return 0;
   
}
