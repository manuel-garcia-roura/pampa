#include "Parser.hxx"

/* Read a plain-text input file: */
int Parser::read(const std::string& filename, Mesh** mesh, Array1D<Material>& materials, 
   TransportMethod& method, GradientScheme& gradient) {
   
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
         
         /* Get the mesh info: */
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
         
      }
      else if (line[0] == "material") {
         
         /* Get the solver: */
         int i = 1;
         std::string solver = line[i++];
         
         /* Create the new material: */
         Material material;
         
         /* Get the material properties depending on the solver: */
         if (solver == "diffusion" || solver == "sn") {
            
            /* Get the number of energy groups: */
            int num_groups = std::stoi(line[i++]);
            material.num_groups = num_groups;
            
            /* Get the nuclear data: */
            PAMPA_CALL(utils::read(material.sigma_total, num_groups, file), 
               "wrong total cross sections in " + filename);
            PAMPA_CALL(utils::read(material.nu_sigma_fission, num_groups, file), 
               "wrong nu-fission cross sections in " + filename);
            PAMPA_CALL(utils::read(material.sigma_scattering, num_groups, num_groups, file), 
               "wrong scattering cross sections in " + filename);
            if (solver == "diffusion") {
               PAMPA_CALL(utils::read(material.diffusion_coefficient, num_groups, file), 
                  "wrong diffusion coefficients in " + filename);
            }
            PAMPA_CALL(utils::read(material.chi, num_groups, file), 
               "wrong fission spectrum in " + filename);
            
         }
         else if (solver == "conduction") {
            Array1D<double> thermal_properties;
            PAMPA_CALL(utils::read(thermal_properties, 3, file), 
               "wrong thermal properties in " + filename);
            material.rho = thermal_properties(0);
            material.cp = thermal_properties(1);
            material.k = thermal_properties(2);
         }
         else
            PAMPA_CHECK(true, 1, "wrong transport method in " + filename);
         
         /* Keep the material definition: */
         materials.pushBack(material);
         
      }
      else if (line[0] == "solver") {
         
         /* Get the transport method type: */
         int i = 1;
         std::string type = line[i++];
         if (type == "diffusion")
            method.type = TM::DIFFUSION;
         else if (type == "sn") {
            method.type = TM::SN;
            method.order = std::stoi(line[i++]);
         }
         else
            PAMPA_CHECK(true, 1, "wrong transport method in " + filename);
         
         /* Get the number of energy groups: */
         int num_groups = std::stoi(line[i++]);
         method.num_groups = num_groups;
         
      }
      else if (line[0] == "gradient") {
         
         /* Get the weight between upwind and linear interpolation: */
         gradient.delta = std::stod(line[1]);
         
         /* Get the boundary interpolation scheme: */
         std::string scheme = line[2];
         if (scheme == "upwind")
            gradient.boundary_interpolation = BI::UPWIND;
         else if (scheme == "ls") {
            gradient.boundary_interpolation = BI::LS;
         }
         else
            PAMPA_CHECK(true, 1, "wrong boundary interpolation scheme in " + filename);
         
      }
      else if (line[0] == "include") {
         
         /* Read an included input file: */
         std::string include_filename = line[1];
         PAMPA_CALL(read(include_filename, mesh, materials, method, gradient), 
            "unable to parse " + include_filename);
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 1, "unrecognized keyword '" + line[0] + "' in " + filename);
         
      }
      
   }
   
   return 0;
   
}
