#include "Parser.hxx"

/* The Parser constructor: */
Parser::Parser() {};

/* The Parser destructor: */
Parser::~Parser() {};

/* Read a plain-text input file: */
bool Parser::read(const std::string &filename, Config &config, Mesh **mesh, 
   std::vector<Material> &materials) {
   
   /* Open the input file: */
   std::ifstream file(filename);
   if (!file.is_open()) {
      std::cout << "Error: unable to open " << filename << "!\n";
      return false;
   }
   
   /* Read the file line by line: */
   std::string line;
   while (std::getline(file, line)) {
      
      /* Skip empty lines and #-marked comments: */
      if (line.empty() || line.at(0) == '#')
         continue;
      
      /* Check for keywords: */
      utils::clean(line);
      std::istringstream iss(line);
      std::string s;
      std::getline(iss, s, ' ');
      if (s == "groups") {
         
         /* Get the number of energy groups: */
         std::getline(iss, s, ' ');
         config.num_groups = std::stoi(s);
         
      }
      else if (s == "mesh") {
         
         /* Get the mesh: */
         std::getline(iss, s, ' ');
         std::string mesh_type = s;
         if (mesh_type == "cartesian")
            *mesh = new CartesianMesh();
         else if (mesh_type == "unstructured")
            *mesh = new UnstructuredExtrudedMesh();
         else {
            std::cout << "Error: wrong mesh type!\n";
            return false;
         }
         std::getline(iss, s, ' ');
         std::string mesh_filename = s;
         if (!((*mesh)->read(mesh_filename))) std::cout << "Error reading the mesh!" << std::endl;
         
      }
      else if (s == "material") {
         
         /* Get the material name and the number of energy groups: */
         Material mat;
         std::getline(iss, s, ' ');
         mat.name = s;
         
         /* Get the nuclear data: */
         if (!utils::read(mat.sigma_absorption, config.num_groups, file)) {
            std::cout << "Error: wrong absorption cross-section data in " << filename << "!\n";
            return false;
         }
         if (!utils::read(mat.nu_sigma_fission, config.num_groups, file)) {
            std::cout << "Error: wrong nu-fission cross-section data in " << filename << "!\n";
            return false;
         }
         if (!utils::read(mat.sigma_scattering, config.num_groups, config.num_groups, file)) {
            std::cout << "Error: wrong scattering cross-section data in " << filename << "!\n";
            return false;
         }
         if (!utils::read(mat.diffusion_coefficient, config.num_groups, file)) {
            std::cout << "Error: wrong diffusion coefficient data in " << filename << "!\n";
            return false;
         }
         if (!utils::read(mat.chi, config.num_groups, file)) {
            std::cout << "Error: wrong fission spectrum data in " << filename << "!\n";
            return false;
         }
         
         /* Keep the material definition: */
         materials.push_back(mat);
         
      }
      else if (s == "include") {
         
         /* Read an included input file: */
         std::getline(iss, s, ' ');
         std::string filename_include = s;
         read(filename_include, config, mesh, materials);
         
      }
      else {
         std::cout << "Error: wrong keyword in " << filename << "!\n";
         return false;
      }
      
   }
   
   return true;
   
};
