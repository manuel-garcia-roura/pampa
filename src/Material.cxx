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
         PAMPA_CALL(utils::read(num_energy_groups, 1, INT_MAX, line[0]), 
            "wrong number of energy groups");
         
      }
      else if (line[0] == "sigma-total") {
         
         /* Get the total cross sections: */
         PAMPA_CALL(utils::read(sigma_total, num_energy_groups, file), 
            "wrong total cross sections");
         
      }
      else if (line[0] == "nu-sigma-fission") {
         
         /* Get the nu-fission cross sections: */
         PAMPA_CALL(utils::read(nu_sigma_fission, num_energy_groups, file), 
            "wrong nu-fission cross sections");
         
      }
      else if (line[0] == "kappa-sigma-fission") {
         
         /* Get the kappa-fission cross sections: */
         PAMPA_CALL(utils::read(kappa_sigma_fission, num_energy_groups, file), 
            "wrong kappa-fission cross sections");
         
      }
      else if (line[0] == "sigma-transport") {
         
         /* Get the transport cross sections: */
         PAMPA_CALL(utils::read(sigma_transport, num_energy_groups, file), 
            "wrong transport cross sections");
         
      }
      else if (line[0] == "sigma-scattering") {
         
         /* Get the scattering cross sections: */
         PAMPA_CALL(utils::read(sigma_scattering, num_energy_groups, num_energy_groups, file), 
            "wrong scattering cross sections");
         
      }
      else if (line[0] == "diffusion-coefficient") {
         
         /* Get the diffusion coefficients: */
         PAMPA_CALL(utils::read(diffusion_coefficient, num_energy_groups, file), 
            "wrong diffusion coefficients");
         
      }
      else if (line[0] == "fission-spectrum" || line[0] == "fission-spectrum-prompt") {
         
         /* Get the prompt fission spectrum: */
         PAMPA_CALL(utils::read(chi_prompt, num_energy_groups, file), 
            "wrong prompt fission spectrum");
         
      }
      else if (line[0] == "fission-spectrum-delayed") {
         
         /* Get the delayed fission spectrum: */
         PAMPA_CALL(utils::read(chi_delayed, num_energy_groups, file), 
            "wrong delayed fission spectrum");
         
      }
      else if (line[0] == "neutron-velocity") {
         
         /* Get the neutron velocity: */
         PAMPA_CALL(utils::read(velocity, num_energy_groups, file), "wrong neutron velocity");
         
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
         
         /* Get the total precursor fraction: */
         for (int g = 0; g < num_precursor_groups; g++)
            beta_total += beta(g);
         
      }
      else if (line[0] == "thermal-conductivity") {
         
         /* Get the thermal conductivity: */
         line = utils::get_next_line(file);
         PAMPA_CALL(utils::read(k, 0.0, DBL_MAX, line[0]), "wrong thermal conductivity");
         
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
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 2, "unrecognized keyword '" + line[0] + "'");
         
      }
      
   }
   
   /* Check non-fissile materials: */
   PAMPA_CHECK(!(chi_prompt.empty()) && nu_sigma_fission.empty(), 3, 
      "missing nu-fission cross sections");
   PAMPA_CHECK(!(nu_sigma_fission.empty()) && chi_prompt.empty(), 4, 
      "missing fission spectrum");
   
   /* Set zero nu-fission cross sections, if not given: */
   if (nu_sigma_fission.empty())
      nu_sigma_fission.resize(num_energy_groups, 0.0);
   
   /* Calculate the kappa-fission cross sections from the nu-fission ones, if not given: */
   if (kappa_sigma_fission.empty()) {
      kappa_sigma_fission.resize(num_energy_groups);
      for (int g = 0; g < num_energy_groups; g++)
         kappa_sigma_fission(g) = (kappa/nu) * nu_sigma_fission(g);
   }
   
   /* Set a zero fission spectrum, if not given: */
   if (chi_prompt.empty())
      chi_prompt.resize(num_energy_groups, 0.0);
   
   /* Use prompt fission spectrum for delayed neutrons, if not given: */
   if (chi_delayed.empty())
      chi_delayed = chi_prompt;
   
   /* Calculate the diffusion coefficients from the transport cross sections, if not given: */
   if (diffusion_coefficient.empty() && !(sigma_transport.empty())) {
      diffusion_coefficient.resize(num_energy_groups);
      for (int g = 0; g < num_energy_groups; g++)
         diffusion_coefficient(g) = 1.0 / (3.0*sigma_transport(g));
   }
   
   return 0;
   
}
