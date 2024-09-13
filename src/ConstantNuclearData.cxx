#include "ConstantNuclearData.hxx"

/* Read the nuclear data from a plain-text input file: */
int ConstantNuclearData::read(std::ifstream& file) {
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty() || line[0] == "}") break;
      
      /* Get the next keyword: */
      unsigned int l = 0;
      if (line[l] == "energy-groups") {
         
         /* Get the number of energy groups: */
         PAMPA_CHECK(utils::read(num_energy_groups, 1, INT_MAX, line[++l]), 
            "wrong number of energy groups");
         
      }
      else if (line[l] == "sigma-total") {
         
         /* Get the total cross sections: */
         PAMPA_CHECK(utils::read(sigma_total, num_energy_groups, file), 
            "wrong total cross sections");
         
      }
      else if (line[l] == "nu-sigma-fission") {
         
         /* Get the nu-fission cross sections: */
         PAMPA_CHECK(utils::read(nu_sigma_fission, num_energy_groups, file), 
            "wrong nu-fission cross sections");
         
      }
      else if (line[l] == "kappa-sigma-fission") {
         
         /* Get the kappa-fission cross sections: */
         PAMPA_CHECK(utils::read(kappa_sigma_fission, num_energy_groups, file), 
            "wrong kappa-fission cross sections");
         
      }
      else if (line[l] == "sigma-transport") {
         
         /* Get the transport cross sections: */
         PAMPA_CHECK(utils::read(sigma_transport, num_energy_groups, file), 
            "wrong transport cross sections");
         
      }
      else if (line[l] == "sigma-scattering") {
         
         /* Get the scattering cross sections: */
         PAMPA_CHECK(utils::read(sigma_scattering, num_energy_groups, num_energy_groups, file), 
            "wrong scattering cross sections");
         
      }
      else if (line[l] == "diffusion-coefficient") {
         
         /* Get the diffusion coefficients: */
         PAMPA_CHECK(utils::read(diffusion_coefficient, num_energy_groups, file), 
            "wrong diffusion coefficients");
         
      }
      else if (line[l] == "fission-spectrum" || line[l] == "fission-spectrum-prompt") {
         
         /* Get the prompt fission spectrum: */
         PAMPA_CHECK(utils::read(chi_prompt, num_energy_groups, file), 
            "wrong prompt fission spectrum");
         
      }
      else if (line[l] == "fission-spectrum-delayed") {
         
         /* Get the delayed fission spectrum: */
         PAMPA_CHECK(utils::read(chi_delayed, num_energy_groups, file), 
            "wrong delayed fission spectrum");
         
      }
      else if (line[l] == "neutron-velocity") {
         
         /* Get the neutron velocity: */
         PAMPA_CHECK(utils::read(velocity, num_energy_groups, file), "wrong neutron velocity");
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   return 0;
   
}

/* Check the nuclear data after reading it: */
int ConstantNuclearData::check(double beta_total) {
   
   /* Check fissile materials: */
   if (!(nu_sigma_fission.empty()) || !(chi_prompt.empty())) {
      PAMPA_CHECK(nu_sigma_fission.empty(), "missing nu-fission cross sections");
      PAMPA_CHECK(chi_prompt.empty(), "missing fission spectrum");
   }
   
   /* Set zero nu-fission cross sections, if not given: */
   if (nu_sigma_fission.empty())
      nu_sigma_fission.resize(num_energy_groups, 0.0);
   
   /* Calculate the kappa-fission cross sections from the nu-fission ones, if not given: */
   if (kappa_sigma_fission.empty()) {
      kappa_sigma_fission = nu_sigma_fission;
      double f = kappa / nu;
      for (int g = 0; g < num_energy_groups; g++)
         kappa_sigma_fission(g) *= f;
   }
   
   /* Set a zero fission spectrum, if not given: */
   if (chi_prompt.empty())
      chi_prompt.resize(num_energy_groups, 0.0);
   
   /* Use prompt fission spectrum for delayed neutrons, if not given: */
   if (chi_delayed.empty())
      chi_delayed = chi_prompt;
   
   /* Calculate the effective fission spectrum for prompt and delayed neutrons: */
   if (beta_total > 0.0) {
      chi_effective.resize(num_energy_groups, 0.0);
      for (int g = 0; g < num_energy_groups; g++)
         chi_effective(g) = (1.0-beta_total)*chi_prompt(g) + beta_total*chi_delayed(g);
   }
   else
      chi_effective = chi_prompt;
   
   /* Calculate the diffusion coefficients from the transport cross sections, if not given: */
   if (!(sigma_transport.empty())) {
      PAMPA_CHECK(!(diffusion_coefficient.empty()), 
         "transport cross sections can only be defined if diffusion coefficients are not");
      diffusion_coefficient.resize(num_energy_groups, 0.0);
      for (int g = 0; g < num_energy_groups; g++)
         diffusion_coefficient(g) = 1.0 / (3.0*sigma_transport(g));
   }
   
   return 0;
   
}

/* Check the nuclear data to use it in a solver: */
int ConstantNuclearData::check(int num_energy_groups, bool diffusion, bool transient) const {
   
   /* Check the nuclear data: */
   PAMPA_CHECK(this->num_energy_groups != num_energy_groups, "wrong number of energy groups");
   PAMPA_CHECK(sigma_total.empty(), "missing total cross sections");
   PAMPA_CHECK(nu_sigma_fission.empty(), "missing nu-fission cross sections");
   PAMPA_CHECK(kappa_sigma_fission.empty(), "missing kappa-fission cross sections");
   PAMPA_CHECK(sigma_scattering.empty(), "missing scattering cross sections");
   PAMPA_CHECK(chi_effective.empty(), "missing effective fission spectrum");
   if (diffusion) {
      PAMPA_CHECK(diffusion_coefficient.empty(), "missing diffusion coefficients");
   }
   if (transient) {
      PAMPA_CHECK(chi_prompt.empty(), "missing prompt fission spectrum");
      PAMPA_CHECK(chi_delayed.empty(), "missing delayed fission spectrum");
      PAMPA_CHECK(velocity.empty(), "missing neutron velocities");
   }
   
   return 0;
   
}
