#include "Parser.hxx"

/* Read a plain-text input file: */
int Parser::read(const std::string& filename, Mesh** mesh, Mesh** mesh_nodal, 
   Array1D<Material*>& materials, Array1D<Solver*>& solvers, Array1D<double>& dt) {
   
   /* Open the input file: */
   std::ifstream file(filename, std::ios_base::in);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty()) break;
      
      /* Get the next keyword: */
      unsigned int l = 0;
      if (line[l] == "mesh") {
         
         /* Create the mesh depending on the mesh type: */
         std::string mesh_type = line[++l];
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
         std::string mesh_filename = line[++l];
         PAMPA_CALL((*mesh)->read(mesh_filename), "unable to read the mesh from " + mesh_filename);
         
         /* Build the mesh: */
         PAMPA_CALL((*mesh)->build(), "unable to build the mesh");
         
      }
      else if (line[l] == "mesh-nodal") {
         
         /* Create the nodal mesh depending on the mesh type: */
         std::string mesh_nodal_type = line[++l];
         if (mesh_nodal_type == "cartesian")
            *mesh_nodal = new CartesianMesh();
         else if (mesh_nodal_type == "unstructured")
            *mesh_nodal = new UnstructuredExtrudedMesh();
         else {
            PAMPA_CHECK(true, 2, "wrong nodal mesh type");
         }
         
         /* Read the nodal mesh: */
         std::string mesh_nodal_filename = line[++l];
         PAMPA_CALL((*mesh_nodal)->read(mesh_nodal_filename), 
            "unable to read the nodal mesh from " + mesh_nodal_filename);
         
         /* Build the nodal mesh: */
         PAMPA_CALL((*mesh_nodal)->build(), "unable to build the nodal mesh");
         
      }
      else if (line[l] == "material") {
         
         /* Get the material name: */
         std::string name = line[++l];
         
         /* Check if this is a boundary-condition material: */
         bool bc = (line.size() > 2) && (line[2] == "bc");
         if (bc) (*mesh)->addBoundary(name);
         
         /* Create the material: */
         Material* mat = new Material(name, bc);
         
         /* Read the material: */
         if ((line.size() > 2) && (line[2] != "bc")) {
            std::string s = line[++l];
            if (s == "{") {
               PAMPA_CALL(mat->read(file), "unable to read the material from " + filename);
            }
            else {
               std::string mat_filename = s;
               PAMPA_CALL(mat->read(mat_filename), 
                  "unable to read the material from " + mat_filename);
            }
         }
         
         /* Keep the material definition: */
         materials.pushBack(mat);
         
      }
      else if (line[l] == "solver") {
         
         /* Create the solver depending on the solver type: */
         Solver* solver = nullptr;
         std::string solver_type = line[++l];
         if (solver_type == "diffusion")
            solver = new DiffusionSolver(materials);
         else if (solver_type == "sn")
            solver = new SNSolver(materials);
         else if (solver_type == "conduction")
            solver = new HeatConductionSolver(*mesh_nodal, materials);
         else if (solver_type == "precursors")
            solver = new PrecursorSolver(materials);
         else if (solver_type == "coupled") {
            std::string name = line[++l];
            solver = new CouplingSolver(name);
         }
         else {
            PAMPA_CHECK(true, 3, "wrong solver type");
         }
         
         /* Set the mesh: */
         solver->setMesh(*mesh);
         
         /* Read the solver: */
         if (l < line.size()+1) {
            std::string s = line[++l];
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
      else if (line[l] == "dt") {
         
         /* Get the dt values: */
         int nt;
         PAMPA_CALL(utils::read(nt, -INT_MAX, INT_MAX, line[++l]), "wrong number of time steps");
         if (nt > 0) {
            PAMPA_CALL(utils::read(dt, nt, file), "wrong dt data");
         }
         else {
            PAMPA_CALL(utils::read(dt, 1, file), "wrong dt data");
            nt = -nt;
            dt.resize(nt, dt(0));
         }
         
      }
      else if (line[l] == "include") {
         
         /* Read an included input file: */
         std::string include_filename = line[++l];
         PAMPA_CALL(read(include_filename, mesh, mesh_nodal, materials, solvers, dt), 
            "unable to parse " + include_filename);
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 5, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   return 0;
   
}
