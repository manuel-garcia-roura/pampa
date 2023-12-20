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
         model.num_groups = std::stoi(line[1]);
         
      }
      else if (line[0] == "mesh") {
         
         /* Get the mesh info: */
         std::string mesh_type = line[1];
         std::string mesh_filename = line[2];
         if (mesh_type == "cartesian")
            model.mesh = new CartesianMesh();
         else if (mesh_type == "unstructured")
            model.mesh = new UnstructuredExtrudedMesh();
         else
            PAMPA_CHECK(true, 1, "wrong mesh type in " + filename);
         
         /* Read the mesh: */
         PAMPA_CALL((model.mesh)->read(mesh_filename), "unable to read the mesh");
         
      }
      else if (line[0] == "material") {
         
         /* Create a new material: */
         Material mat;
         
         /* Get the material name: */
         mat.name = line[1];
         
         /* Get the nuclear data: */
         PAMPA_CALL(utils::read(mat.sigma_absorption, model.num_groups, file), 
            "wrong absorption cross-section data in " + filename);
         PAMPA_CALL(utils::read(mat.nu_sigma_fission, model.num_groups, file), 
            "wrong nu-fission cross-section data in " + filename);
         PAMPA_CALL(utils::read(mat.sigma_scattering, model.num_groups, model.num_groups, file), 
            "wrong scattering cross-section data in " + filename);
         PAMPA_CALL(utils::read(mat.diffusion_coefficient, model.num_groups, file), 
            "wrong diffusion coefficient data in " + filename);
         PAMPA_CALL(utils::read(mat.chi, model.num_groups, file), 
            "wrong fission spectrum data in " + filename);
         
         /* Calculate the total cross sections: */
         mat.sigma_total.resize(model.num_groups);
         mat.sigma_removal.resize(model.num_groups);
         for (int g = 0; g < model.num_groups; g++) {
            mat.sigma_total[g] = mat.sigma_absorption[g] + mat.nu_sigma_fission[g];
            for (int g2 = 0; g2 < model.num_groups; g2++)
               mat.sigma_total[g] += mat.sigma_scattering[g][g2];
            mat.sigma_removal[g] = mat.sigma_total[g] - mat.sigma_scattering[g][g];
         }
         
         /* Keep the material definition: */
         model.materials.push_back(mat);
         
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
