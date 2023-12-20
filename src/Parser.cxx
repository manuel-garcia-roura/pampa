#include "Parser.hxx"

/* The Parser constructor: */
Parser::Parser() {};

/* The Parser destructor: */
Parser::~Parser() {};

/* Read a plain-text input file: */
int Parser::read(const std::string &filename, Model &model) {
   
   /* Open the input file: */
   std::ifstream file(filename);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line:*/
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty())
         break;
      
      /* Get the next keyword: */
      if (line[0] == "groups") {
         
         /* Get the number of energy groups: */
         int num_groups = std::stoi(line[1]);
         model.setNumEnergyGroups(num_groups);
         
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
         
         /* Get the number of energy groups: */
         int num_groups = model.getNumEnergyGroups();
         
         /* Create a new material: */
         Material material;
         
         /* Get the material name: */
         material.name = line[1];
         
         /* Get the nuclear data: */
         PAMPA_CALL(utils::read(material.sigma_absorption, num_groups, file), 
            "wrong absorption cross-section data in " + filename);
         PAMPA_CALL(utils::read(material.nu_sigma_fission, num_groups, file), 
            "wrong nu-fission cross-section data in " + filename);
         PAMPA_CALL(utils::read(material.sigma_scattering, num_groups, num_groups, file), 
            "wrong scattering cross-section data in " + filename);
         PAMPA_CALL(utils::read(material.diffusion_coefficient, num_groups, file), 
            "wrong diffusion coefficient data in " + filename);
         PAMPA_CALL(utils::read(material.chi, num_groups, file), 
            "wrong fission spectrum data in " + filename);
         
         /* Calculate the total cross sections: */
         material.sigma_total.resize(num_groups);
         material.sigma_removal.resize(num_groups);
         for (int g = 0; g < num_groups; g++) {
            material.sigma_total[g] = material.sigma_absorption[g] + material.nu_sigma_fission[g];
            for (int g2 = 0; g2 < num_groups; g2++)
               material.sigma_total[g] += material.sigma_scattering[g][g2];
            material.sigma_removal[g] = material.sigma_total[g] - material.sigma_scattering[g][g];
         }
         
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
   
};
