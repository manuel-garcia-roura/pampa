#include "NuclearData.hxx"

/* Read the nuclear data from a plain-text input file: */
int NuclearData::read(std::ifstream& file) {
   
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
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 2, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   return 0;
   
}
