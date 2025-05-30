#include "Parser.hxx"

/* Read a plain-text input file: */
int Parser::read(const std::string& filename, Mesh** mesh, Mesh** mesh_nodal, 
   Array1D<Material*>& materials, Array1D<Solver*>& solvers, Array1D<double>& dt) {
   
   /* Open the input file: */
   std::ifstream file(filename, std::ios_base::in);
   PAMPA_CHECK(!file.is_open(), "unable to open " + filename);
   
   /* Read the file line by line: */
   bool mesh_ready = false;
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = input::get_next_line(file);
      if (line.empty()) break;
      
      /* Get the next keyword: */
      unsigned int l = 0;
      if (line[l] == "mesh") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 3, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Create the mesh depending on the mesh type: */
         std::string mesh_type = line[++l];
         if (mesh_type == "cartesian")
            *mesh = new CartesianMesh();
         else if (mesh_type == "unstructured")
            *mesh = new UnstructuredExtrudedMesh();
         else if (mesh_type == "partitioned")
            *mesh = new PartitionedMesh();
         else {
            PAMPA_CHECK(true, "wrong mesh type");
         }
         
         /* Read the mesh: */
         std::string mesh_filename = line[++l];
         PAMPA_CHECK((*mesh)->read(mesh_filename), "unable to read the mesh from " + mesh_filename);
         
         /* Build the mesh: */
         PAMPA_CHECK((*mesh)->build(), "unable to build the mesh");
         
      }
      else if (line[l] == "mesh-nodal") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 3, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Create the nodal mesh depending on the mesh type: */
         std::string mesh_nodal_type = line[++l];
         if (mesh_nodal_type == "cartesian")
            *mesh_nodal = new CartesianMesh();
         else if (mesh_nodal_type == "unstructured")
            *mesh_nodal = new UnstructuredExtrudedMesh();
         else {
            PAMPA_CHECK(true, "wrong nodal mesh type");
         }
         
         /* Read the nodal mesh: */
         std::string mesh_nodal_filename = line[++l];
         PAMPA_CHECK((*mesh_nodal)->read(mesh_nodal_filename), 
            "unable to read the nodal mesh from " + mesh_nodal_filename);
         
         /* Build the nodal mesh: */
         PAMPA_CHECK((*mesh_nodal)->build(), "unable to build the nodal mesh");
         
      }
      else if (line[l] == "material") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 3, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the material name: */
         std::string name = line[++l];
         
         /* Create the material: */
         Material* mat = new Material(name);
         
         /* Read the material: */
         if (line.size() > 2) {
            std::string s = line[++l];
            if (s == "{") {
               PAMPA_CHECK(mat->read(file), "unable to read the material from " + filename);
            }
            else {
               std::string mat_filename = s;
               PAMPA_CHECK(mat->read(mat_filename), 
                  "unable to read the material from " + mat_filename);
            }
         }
         
         /* Check if this is a boundary-condition material: */
         if (mat->isBC()) (*mesh)->addBoundary(name);
         
         /* Keep the material definition: */
         materials.pushBack(mat);
         
      }
      else if (line[l] == "solver") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() < 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Process the mesh, if not done yet: */
         if (!mesh_ready) {
            
            /* Split materials in the mesh: */
            bool split_materials = false;
            for (int i = 0; i < materials.size(); i++)
               split_materials |= materials(i)->isSplit();
            if (split_materials) {
               PAMPA_CHECK((*mesh)->splitMaterials(materials), 
                  "unable to split materials in the mesh");
            }
            
            /* Remove boundary-condition materials from the mesh and swap the meshes: */
            bool remove_materials = false;
            for (int i = 0; i < materials.size(); i++)
               remove_materials |= materials(i)->isBC();
            if (remove_materials) {
               Mesh* mesh_new = nullptr;
               PAMPA_CHECK((*mesh)->removeBCMaterials(materials, &mesh_new), 
                  "unable to remove boundary-condition materials from the mesh");
               utils::free(mesh);
               *mesh = mesh_new;
            }
            
            /* Partition the mesh and swap the meshes: */
            if (mpi::size > 1 && !((*mesh)->isPartitioned())) {
               Mesh* submesh = nullptr;
               PAMPA_CHECK((*mesh)->partition(&submesh), "unable to partition the mesh");
               utils::free(mesh);
               *mesh = submesh;
            }
            
            /* Done processing the mesh: */
            mesh_ready = true;
#ifdef DEBUG
            PAMPA_CHECK((*mesh)->writeVTK("mesh_data", -1), "unable to write the mesh");
#endif
            
         }
         
         /* Create the solver depending on the solver type: */
         Solver* solver = nullptr;
         std::string solver_type = line[++l];
         if (solver_type == "diffusion")
            solver = new DiffusionSolver(*mesh, materials);
         else if (solver_type == "sn")
            solver = new SNSolver(*mesh, materials);
         else if (solver_type == "conduction")
            solver = new HeatConductionSolver(*mesh, *mesh_nodal, materials);
         else if (solver_type == "precursors")
            solver = new PrecursorSolver(*mesh, materials);
         else if (solver_type == "coupled") {
            /* Check the number of arguments: */
            PAMPA_CHECK(line.size() < 3, "wrong number of arguments for keyword '" + line[l] + "'");
            std::string name = line[++l];
            solver = new CouplingSolver(name, *mesh);
         }
         else {
            PAMPA_CHECK(true, "wrong solver type");
         }
         
         /* Read the solver: */
         if (l < line.size()+1) {
            std::string s = line[++l];
            if (s == "{") {
               PAMPA_CHECK(solver->read(file, solvers), 
                  "unable to read the solver from " + filename);
            }
            else {
               std::string solver_filename = s;
               PAMPA_CHECK(solver->read(solver_filename, solvers), 
                  "unable to read the solver from " + solver_filename);
            }
         }
         
         /* Keep the solver definition: */
         solvers.pushBack(solver);
         
      }
      else if (line[l] == "dt") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the dt values: */
         int nt;
         PAMPA_CHECK(input::read(nt, -INT_MAX, INT_MAX, line[++l]), "wrong number of time steps");
         if (nt > 0) {
            PAMPA_CHECK(input::read(dt, nt, 0.0, DBL_MAX, file), "wrong dt data");
         }
         else {
            PAMPA_CHECK(input::read(dt, 1, 0.0, DBL_MAX, file), "wrong dt data");
            nt = -nt;
            dt.resize(nt, dt(0));
         }
         
      }
      else if (line[l] == "vtk") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK((line.size() < 2) || (line.size() > 3), 
            "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the switch for .vtk output: */
         PAMPA_CHECK(input::read(vtk::on, line[++l]), "wrong switch for .vtk output");
         
         /* Get the output interval in time steps: */
         if (line.size() == 3) {
            PAMPA_CHECK(input::read(vtk::dn, 1, INT_MAX, line[++l]), "wrong .vtk output interval");
         }
         
      }
      else if (line[l] == "petsc") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 3, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the option name: */
         std::string name = line[++l];
         
         /* Set the PETSc/SLEPc option: */
         if (name == "verbose") {
            PAMPA_CHECK(input::read(petsc::verbose, line[++l]), 
               "wrong switch for verbose PETSc output");
         }
         else if (name == "dump") {
            PAMPA_CHECK(input::read(petsc::dump, line[++l]), 
               "wrong switch to write the solution in PETSc format");
         }
         else {
            std::string value = line[++l];
            PAMPA_CHECK(petsc::set(name, value), "wrong PETSc/SLEPc option");
         }
         
      }
      else if (line[l] == "include") {
         
         /* Read an included input file: */
         std::string include_filename = line[++l];
         PAMPA_CHECK(read(include_filename, mesh, mesh_nodal, materials, solvers, dt), 
            "unable to parse " + include_filename);
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   return 0;
   
}
