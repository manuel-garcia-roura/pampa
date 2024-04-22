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
      else if (line[0] == "thermal-properties") {
         
         /* Get the thermal properties: */
         std::vector<std::string> names{"constant", "graphite-h-451", "graphite-matrix-a3-27"};
         PAMPA_CALL(utils::read<TH::Properties>(thermal_properties, names, line[1]), 
            "wrong thermal properties");
         
         /* Set reference properties: */
         k0 = k(1000.0);
         rho0 = rho(1000.0);
         cp0 = cp(1000.0);
         
      }
      else if (line[0] == "thermal-conductivity") {
         
         /* Get the thermal conductivity: */
         line = utils::get_next_line(file);
         PAMPA_CALL(utils::read(k0, 0.0, DBL_MAX, line[0]), "wrong thermal conductivity");
         
      }
      else if (line[0] == "density") {
         
         /* Get the density: */
         line = utils::get_next_line(file);
         PAMPA_CALL(utils::read(rho0, 0.0, DBL_MAX, line[0]), "wrong density");
         
      }
      else if (line[0] == "specific-heat-capacity") {
         
         /* Get the specific heat capacity: */
         line = utils::get_next_line(file);
         PAMPA_CALL(utils::read(cp0, 0.0, DBL_MAX, line[0]), "wrong specific heat capacity");
         
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
   
   /* Calculate the effective fission spectrum for prompt and delayed neutrons: */
   chi_eff.resize(num_energy_groups);
   for (int g = 0; g < num_energy_groups; g++)
      chi_eff(g) = (1.0-beta_total)*chi_prompt(g) + beta_total*chi_delayed(g);
   
   /* Calculate the diffusion coefficients from the transport cross sections, if not given: */
   if (diffusion_coefficient.empty() && !(sigma_transport.empty())) {
      diffusion_coefficient.resize(num_energy_groups);
      for (int g = 0; g < num_energy_groups; g++)
         diffusion_coefficient(g) = 1.0 / (3.0*sigma_transport(g));
   }
   
   return 0;
   
}

/* Get the thermal conductivity: */
double Material::k(double T) const {
   
   /* Return the thermal conductivity depending on the material: */
   double k_T = -1.0;
   switch (thermal_properties) {
      case TH::CONSTANT:
         k_T = k0;
         break;
      case TH::GRAPHITE_H_451:
         if (T < 500.0) T = 500.0;
         if (T > 1800.0) T = 1800.0;
         k_T = 3.28248e-5*std::pow(T, 2) - 1.24890e-1*T + 1.69214e2;
         k_T *= 1.0e-2;
         break;
      case TH::GRAPHITE_MATRIX_A3_27:
         k_T = 47.4 * (1.0-9.7556e-4*(T-100.0)*std::exp(-6.0360e-4*T));
         k_T *= 1.0e-2;
         break;
      default:
         break;
   }
   return k_T;
   
}

/* Get the density: */
double Material::rho(double T) const {
   
   /* Return the density depending on the material: */
   double rho_T = -1.0;
   switch (thermal_properties) {
      case TH::CONSTANT:
         rho_T = rho0;
         break;
      case TH::GRAPHITE_H_451:
         rho_T = 1.850e-3;
         break;
      case TH::GRAPHITE_MATRIX_A3_27:
         rho_T = 1.700e-3;
         break;
      default:
         break;
   }
   return rho_T;
   
}

/* Get the specific heat capacity: */
double Material::cp(double T) const {
   
   /* Return the specific heat capacity depending on the material: */
   double cp_T = -1.0;
   switch (thermal_properties) {
      case TH::CONSTANT:
         cp_T = cp0;
         break;
      case TH::GRAPHITE_H_451:
         if (T < 200.0) T = 200.0;
         if (T > 3500.0) T = 3500.0;
         cp_T = 0.54212 - 2.42667e-6*T - 90.2725*std::pow(T, -1) - 43449.3*std::pow(T, -3) + 
                   1.59309e7*std::pow(T, -3) - 1.43688e9*std::pow(T, -4);
         cp_T *= 4184.0;
         break;
      case TH::GRAPHITE_MATRIX_A3_27:
         if (T < 200.0) T = 200.0;
         if (T > 3500.0) T = 3500.0;
         cp_T = 0.54212 - 2.42667e-6*T - 90.2725*std::pow(T, -1) - 43449.3*std::pow(T, -3) + 
                   1.59309e7*std::pow(T, -3) - 1.43688e9*std::pow(T, -4);
         cp_T *= 4184.0;
         break;
      default:
         break;
   }
   return cp_T;
   
}
