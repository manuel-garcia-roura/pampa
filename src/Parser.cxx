#include "Parser.hxx"

/* The Parser constructor: */
Parser::Parser() {};

/* The Parser destructor: */
Parser::~Parser() {};

/* Read a plain-text input file: */
bool Parser::read(const std::string &filename, Model &model) {
   
   /* Open the input file: */
   std::ifstream file(filename);
   if (!file.is_open()) {
      std::cout << "Error: unable to open " << filename << "!\n";
      return false;
   }
   
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
         else {
            std::cout << "Error: wrong mesh type in " << filename << "!\n";
            return false;
         }
         
         /* Read the mesh: */
         if (!((model.mesh)->read(mesh_filename))) {
            std::cout << "Error: unable to read the mesh!" << std::endl;
            return false;
         }
         
      }
      else if (line[0] == "material") {
         
         /* Create a new material: */
         Material mat;
         
         /* Get the material name: */
         mat.name = line[1];
         
         /* Get the nuclear data: */
         if (!utils::read(mat.sigma_absorption, model.num_groups, file)) {
            std::cout << "Error: wrong absorption cross-section data in " << filename << "!\n";
            return false;
         }
         if (!utils::read(mat.nu_sigma_fission, model.num_groups, file)) {
            std::cout << "Error: wrong nu-fission cross-section data in " << filename << "!\n";
            return false;
         }
         if (!utils::read(mat.sigma_scattering, model.num_groups, model.num_groups, file)) {
            std::cout << "Error: wrong scattering cross-section data in " << filename << "!\n";
            return false;
         }
         if (!utils::read(mat.diffusion_coefficient, model.num_groups, file)) {
            std::cout << "Error: wrong diffusion coefficient data in " << filename << "!\n";
            return false;
         }
         if (!utils::read(mat.chi, model.num_groups, file)) {
            std::cout << "Error: wrong fission spectrum data in " << filename << "!\n";
            return false;
         }
         
         /* Keep the material definition: */
         model.materials.push_back(mat);
         
      }
      else if (line[0] == "include") {
         
         /* Read an included input file: */
         std::string include_filename = line[1];
         if (!read(include_filename, model)) {
            std::cout << "Error: unable to parse " << include_filename << "!\n";
            return false;
         }
         
      }
      else {
         
         /* Wrong keyword: */
         std::cout << "Error: wrong keyword in " << filename << "!\n";
         return false;
         
      }
      
   }
   
   return true;
   
};
