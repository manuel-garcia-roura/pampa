#include "Material.hxx"

/* Read the material from a plain-text input file: */
int Material::read(const std::string& filename) {
   
   /* Open the input file: */
   std::ifstream file(filename);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty()) break;
      
      /* Get the next keyword: */
      if (line[0] == "energy-groups") {
         
         /* Get the number of energy groups: */
         line = utils::get_next_line(file);
         PAMPA_CALL(utils::read(num_groups, 1, INT_MAX, line[0]), "wrong number of energy groups");
         
      }
      else if (line[0] == "sigma-total") {
         
         /* Get the total cross sections: */
         PAMPA_CALL(utils::read(sigma_total, num_groups, file), "wrong total cross sections");
         
      }
      else if (line[0] == "nu-sigma-fission") {
         
         /* Get the nu-fission cross sections: */
         PAMPA_CALL(utils::read(nu_sigma_fission, num_groups, file), 
            "wrong nu-fission cross sections");
         
      }
      else if (line[0] == "sigma-scattering") {
         
         /* Get the scattering cross sections: */
         PAMPA_CALL(utils::read(sigma_scattering, num_groups, num_groups, file), 
            "wrong scattering cross sections");
         
      }
      else if (line[0] == "diffusion-coefficient") {
         
         /* Get the diffusion coefficients: */
         PAMPA_CALL(utils::read(diffusion_coefficient, num_groups, file), 
            "wrong diffusion coefficients");
         
      }
      else if (line[0] == "fission-spectrum") {
         
         /* Get the fission spectrum: */
         PAMPA_CALL(utils::read(chi, num_groups, file), "wrong fission spectrum");
         
      }
      else if (line[0] == "precursor-groups") {
         
         /* Get the number of delayed-neutron precursor groups: */
         line = utils::get_next_line(file);
         PAMPA_CALL(utils::read(num_precursor_groups, 1, INT_MAX, line[0]), 
            "wrong number of delayed-neutron precursor groups");
         
      }
      else if (line[0] == "precursor-lambda") {
         
         /* Get the precursor decay constants: */
         PAMPA_CALL(utils::read(lambda, num_precursor_groups, file), 
            "wrong precursor decay constants");
         
      }
      else if (line[0] == "precursor-beta") {
         
         /* Get the precursor fractions: */
         PAMPA_CALL(utils::read(beta, num_precursor_groups, file), "wrong precursor fractions");
         
      }
      else if (line[0] == "density") {
         
         /* Get the density: */
         line = utils::get_next_line(file);
         PAMPA_CALL(utils::read(rho, 0.0, DBL_MAX, line[0]), "wrong density");
         
      }
      else if (line[0] == "specific-heat-capacity") {
         
         /* Get the specific heat capacity: */
         line = utils::get_next_line(file);
         PAMPA_CALL(utils::read(cp, 0.0, DBL_MAX, line[0]), "wrong specific heat capacity");
         
      }
      else if (line[0] == "thermal-conductivity") {
         
         /* Get the thermal conductivity: */
         line = utils::get_next_line(file);
         PAMPA_CALL(utils::read(k, 0.0, DBL_MAX, line[0]), "wrong thermal conductivity");
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 1, "unrecognized keyword '" + line[0] + "'");
         
      }
      
   }
   
   return 0;
   
}
