#include "Parser.hxx"

/* Read a plain-text input file: */
int Parser::read(const std::string& filename, Mesh** mesh, Array1D<Material>& materials, 
   Solver** solver) {
   
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
         
         /* Get the mesh data: */
         std::string mesh_type = line[1];
         std::string mesh_filename = line[2];
         
         /* Create the mesh: */
         if (mesh_type == "cartesian")
            *mesh = new CartesianMesh();
         else if (mesh_type == "unstructured")
            *mesh = new UnstructuredExtrudedMesh();
         else if (mesh_type == "partitioned")
            *mesh = new PartitionedMesh();
         else
            PAMPA_CHECK(true, 1, "wrong mesh type in " + filename);
         
         /* Read the mesh: */
         PAMPA_CALL((*mesh)->read(mesh_filename), "unable to read the mesh");
         
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
         
         /* Get the material data: */
         std::string material_filename = line[1];
         
         /* Create the new material: */
         Material material;
         
         /* Read the material: */
         PAMPA_CALL(material.read(material_filename), "unable to read the material");
         
         /* Keep the material definition: */
         materials.pushBack(material);
         
      }
      else if (line[0] == "solver") {
         
         /* Get the solver type: */
         int i = 1;
         std::string solver_type = line[i++];
         
         /* Create the solver depending on the solver type: */
         if (solver_type == "diffusion") {
            
            /* Get the number of energy groups: */
            int num_groups = std::stoi(line[i++]);
            
            /* Create the solver: */
            *solver = new DiffusionSolver(*mesh, materials, num_groups);
            
         }
         else if (solver_type == "sn") {
            
            /* Get the order (N) of the SN method: */
            int order = std::stoi(line[i++]);
            
            /* Get the number of energy groups: */
            int num_groups = std::stoi(line[i++]);
            
            /* Get the weight between upwind and linear schemes for face interpolation: */
            double face_interpolation_delta;
            PAMPA_CALL(utils::read(face_interpolation_delta, 0.0, 1.0, line[i++]), 
               "wrong weight between upwind and linear interpolation in " + filename);
            
            /* Get the switch to use the least-squares gradient for boundary interpolation: */
            bool boundary_interpolation_ls;
            PAMPA_CALL(utils::read(boundary_interpolation_ls, line[i++]), 
               "wrong switch for least-squares boundary interpolation in " + filename);
            
            /* Create the solver: */
            *solver = new SNSolver(*mesh, materials, num_groups, order, face_interpolation_delta, 
                             boundary_interpolation_ls);
            
         }
         else if (solver_type == "conduction") {
            
            /* Create the solver: */
            *solver = new HeatConductionSolver(*mesh, materials);
            
         }
         else
            PAMPA_CHECK(true, 1, "wrong solver in " + filename);
         
      }
      else if (line[0] == "include") {
         
         /* Read an included input file: */
         std::string include_filename = line[1];
         PAMPA_CALL(read(include_filename, mesh, materials, solver), 
            "unable to parse " + include_filename);
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 1, "unrecognized keyword '" + line[0] + "' in " + filename);
         
      }
      
   }
   
   return 0;
   
}
