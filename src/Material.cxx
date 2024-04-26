#include "Material.hxx"

/* Read the material from a plain-text input file: */
int Material::read(const std::string& filename) {
   
   /* Open the input file: */
   std::ifstream file(filename);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Read the material: */
   PAMPA_CALL(read(file), "unable to read the material from " + filename);
   
   return 0;
   
}

/* Read the material from a plain-text input file: */
int Material::read(std::ifstream& file) {
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty() || line[0] == "}") break;
      
      /* Get the next keyword: */
      unsigned int l = 0;
      if (line[l] == "energy-groups") {
         
         /* Get the number of energy groups: */
         PAMPA_CALL(utils::read(num_energy_groups, 1, INT_MAX, line[++l]), 
            "wrong number of energy groups");
         
      }
      else if (line[l] == "sigma-total") {
         
         /* Get the total cross sections: */
         PAMPA_CALL(utils::read(sigma_total, num_energy_groups, file), 
            "wrong total cross sections");
         
      }
      else if (line[l] == "nu-sigma-fission") {
         
         /* Get the nu-fission cross sections: */
         PAMPA_CALL(utils::read(nu_sigma_fission, num_energy_groups, file), 
            "wrong nu-fission cross sections");
         
      }
      else if (line[l] == "kappa-sigma-fission") {
         
         /* Get the kappa-fission cross sections: */
         PAMPA_CALL(utils::read(kappa_sigma_fission, num_energy_groups, file), 
            "wrong kappa-fission cross sections");
         
      }
      else if (line[l] == "sigma-transport") {
         
         /* Get the transport cross sections: */
         PAMPA_CALL(utils::read(sigma_transport, num_energy_groups, file), 
            "wrong transport cross sections");
         
      }
      else if (line[l] == "sigma-scattering") {
         
         /* Get the scattering cross sections: */
         PAMPA_CALL(utils::read(sigma_scattering, num_energy_groups, num_energy_groups, file), 
            "wrong scattering cross sections");
         
      }
      else if (line[l] == "diffusion-coefficient") {
         
         /* Get the diffusion coefficients: */
         PAMPA_CALL(utils::read(diffusion_coefficient, num_energy_groups, file), 
            "wrong diffusion coefficients");
         
      }
      else if (line[l] == "fission-spectrum" || line[l] == "fission-spectrum-prompt") {
         
         /* Get the prompt fission spectrum: */
         PAMPA_CALL(utils::read(chi_prompt, num_energy_groups, file), 
            "wrong prompt fission spectrum");
         
      }
      else if (line[l] == "fission-spectrum-delayed") {
         
         /* Get the delayed fission spectrum: */
         PAMPA_CALL(utils::read(chi_delayed, num_energy_groups, file), 
            "wrong delayed fission spectrum");
         
      }
      else if (line[l] == "neutron-velocity") {
         
         /* Get the neutron velocity: */
         PAMPA_CALL(utils::read(velocity, num_energy_groups, file), "wrong neutron velocity");
         
      }
      else if (line[l] == "precursor-data") {
         
         /* Create the precursor data: */
         precursor_data = new PrecursorData();
         
         /* Read the precursor data: */
         PAMPA_CHECK(line[++l] != "{", 1, "missing opening '{' for precursor data");
         PAMPA_CALL(precursor_data->read(file), "unable to read the precursor data");
         
      }
      else if (line[l] == "thermal-properties") {
         
         /* Create the thermal properties depending on the type: */
         std::string thermal_properties_type = line[++l];
         if (thermal_properties_type == "constant") {
            double k0, rho0, cp0;
            PAMPA_CALL(utils::read(k0, 0.0, DBL_MAX, line[++l]), "wrong thermal conductivity");
            PAMPA_CALL(utils::read(rho0, 0.0, DBL_MAX, line[++l]), "wrong density");
            PAMPA_CALL(utils::read(cp0, 0.0, DBL_MAX, line[++l]), "wrong specific heat capacity");
            thermal_properties = new ConstantProperties(k0, rho0, cp0);
         }
         else if (thermal_properties_type == "graphite-h-451")
            thermal_properties = new GraphiteProperties();
         else if (thermal_properties_type == "graphite-matrix-a3-27")
            thermal_properties = new GraphiteMatrixProperties();
         else {
            PAMPA_CHECK(true, 1, "wrong mesh type");
         }
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 2, "unrecognized keyword '" + line[l] + "'");
         
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
   
   /* Calculate the effective fission spectrum for prompt and delayed neutrons: */
   if (hasPrecursorData()) {
      chi_eff.resize(num_energy_groups);
      for (int g = 0; g < num_energy_groups; g++)
         chi_eff(g) = (1.0-precursor_data->beta_total)*chi_prompt(g) + 
                         precursor_data->beta_total*chi_delayed(g);
   }
   else
      chi_eff = chi_prompt;
   
   /* Calculate the diffusion coefficients from the transport cross sections, if not given: */
   if (diffusion_coefficient.empty() && !(sigma_transport.empty())) {
      diffusion_coefficient.resize(num_energy_groups);
      for (int g = 0; g < num_energy_groups; g++)
         diffusion_coefficient(g) = 1.0 / (3.0*sigma_transport(g));
   }
   
   return 0;
   
}

/* Check the precursor data: */
int Material::checkPrecursorData(int num_precursor_groups) const {
   
   /* Check the precursor data: */
   PAMPA_CHECK(precursor_data->num_precursor_groups != num_precursor_groups, 1, 
      "wrong number of delayed-neutron precursor groups");
   PAMPA_CHECK((precursor_data->lambda).empty(), 2, "missing precursor decay constants");
   PAMPA_CHECK((precursor_data->beta).empty(), 3, "missing precursor fractions");
   
   return 0;
   
}
