#include "Parser.hxx"

/* Read a plain-text input file: */
int Parser::read(const std::string &filename, Model &model) {
   
   /* Open the input file: */
   std::ifstream file(filename);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty()) break;
      
      /* Get the next keyword: */
      if (line[0] == "model") {
         
         /* Create the transport method: */
         TransportMethod method;
         
         /* Get the transport method type: */
         std::string type = line[1];
         if (type == "diffusion")
            method.type = TM::DIFFUSION;
         else if (type == "sn")
            method.type = TM::SN;
         else
            PAMPA_CHECK(true, 1, "wrong transport method in " + filename);
         
         /* Get the number of energy groups: */
         int num_groups = std::stoi(line[2]);
         method.num_groups = num_groups;
         
         /* Keep the transport method: */
         model.setTransportMethod(method);
         
      }
      else if (line[0] == "mesh") {
         
         /* Get the mesh info: */
         std::string mesh_type = line[1];
         std::string mesh_filename = line[2];
         
         /* Create the mesh: */
         Mesh *mesh;
         if (mesh_type == "cartesian")
            mesh = new CartesianMesh();
         else if (mesh_type == "unstructured")
            mesh = new UnstructuredExtrudedMesh();
         else
            PAMPA_CHECK(true, 1, "wrong mesh type in " + filename);
         
         /* Read the mesh: */
         PAMPA_CALL(mesh->read(mesh_filename), "unable to read the mesh");
         
         /* Keep the mesh: */
         model.setMesh(mesh);
         
      }
      else if (line[0] == "material") {
         
         /* Create a new material: */
         Material material;
         
         /* Get the transport method type: */
         std::string type = line[1];
         if (type == "diffusion")
            material.method.type = TM::DIFFUSION;
         else if (type == "sn")
            material.method.type = TM::SN;
         else
            PAMPA_CHECK(true, 1, "wrong transport method in " + filename);
         
         /* Get the number of energy groups: */
         int num_groups = std::stoi(line[2]);
         material.method.num_groups = num_groups;
         
         /* Get the nuclear data: */
         PAMPA_CALL(utils::read(material.sigma_total, num_groups, file), 
            "wrong total cross-section data in " + filename);
         PAMPA_CALL(utils::read(material.nu_sigma_fission, num_groups, file), 
            "wrong nu-fission cross-section data in " + filename);
         PAMPA_CALL(utils::read(material.sigma_scattering, num_groups, num_groups, file), 
            "wrong scattering cross-section data in " + filename);
         switch (material.method.type) {
            case TM::DIFFUSION : {
               PAMPA_CALL(utils::read(material.diffusion_coefficient, num_groups, file), 
                  "wrong diffusion coefficient data in " + filename);
               break;
            }
            case TM::SN : {
               PAMPA_CALL(utils::read(material.sigma_scattering_s1, num_groups, num_groups, file), 
                  "wrong linear-anisotropic scattering cross-section data in " + filename);
               break;
            }
         }
         PAMPA_CALL(utils::read(material.chi, num_groups, file), 
            "wrong fission spectrum data in " + filename);
         
         /* Keep the material definition: */
         model.addMaterial(material);
         
      }
      else if (line[0] == "include") {
         
         /* Read an included input file: */
         std::string include_filename = line[1];
         PAMPA_CALL(read(include_filename, model), "unable to parse " + include_filename);
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 1, "unrecognized keyword '" + line[0] + "' in " + filename);
         
      }
      
   }
   
   return 0;
   
}
