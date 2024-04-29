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
      if (line[l] == "nuclear-data") {
         
         /* Create the nuclear data: */
         nuclear_data = new NuclearData();
         
         /* Read the nuclear data: */
         PAMPA_CHECK(line[++l] != "{", 1, "missing opening '{' for nuclear data");
         PAMPA_CALL(nuclear_data->read(file), "unable to read the nuclear data");
         
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
   PAMPA_CHECK(!((nuclear_data->chi_prompt).empty()) && (nuclear_data->nu_sigma_fission).empty(), 
      3, "missing nu-fission cross sections");
   PAMPA_CHECK(!((nuclear_data->nu_sigma_fission).empty()) && (nuclear_data->chi_prompt).empty(), 
      4, "missing fission spectrum");
   
   /* Set zero nu-fission cross sections, if not given: */
   if ((nuclear_data->nu_sigma_fission).empty())
      (nuclear_data->nu_sigma_fission).resize(nuclear_data->num_energy_groups, 0.0);
   
   /* Calculate the kappa-fission cross sections from the nu-fission ones, if not given: */
   if ((nuclear_data->kappa_sigma_fission).empty()) {
      (nuclear_data->kappa_sigma_fission).resize(nuclear_data->num_energy_groups);
      for (int g = 0; g < nuclear_data->num_energy_groups; g++)
         nuclear_data->kappa_sigma_fission(g) = (nuclear_data->kappa/nuclear_data->nu) * 
                                                   nuclear_data->nu_sigma_fission(g);
   }
   
   /* Set a zero fission spectrum, if not given: */
   if ((nuclear_data->chi_prompt).empty())
      (nuclear_data->chi_prompt).resize(nuclear_data->num_energy_groups, 0.0);
   
   /* Use prompt fission spectrum for delayed neutrons, if not given: */
   if ((nuclear_data->chi_delayed).empty())
      nuclear_data->chi_delayed = nuclear_data->chi_prompt;
   
   /* Calculate the effective fission spectrum for prompt and delayed neutrons: */
   if (hasPrecursorData()) {
      (nuclear_data->chi_eff).resize(nuclear_data->num_energy_groups);
      for (int g = 0; g < nuclear_data->num_energy_groups; g++)
         nuclear_data->chi_eff(g) = (1.0-precursor_data->beta_total)*nuclear_data->chi_prompt(g) + 
                                       precursor_data->beta_total*nuclear_data->chi_delayed(g);
   }
   else
      nuclear_data->chi_eff = nuclear_data->chi_prompt;
   
   /* Calculate the diffusion coefficients from the transport cross sections, if not given: */
   if ((nuclear_data->diffusion_coefficient).empty() && !(nuclear_data->sigma_transport).empty()) {
      (nuclear_data->diffusion_coefficient).resize(nuclear_data->num_energy_groups);
      for (int g = 0; g < nuclear_data->num_energy_groups; g++)
         nuclear_data->diffusion_coefficient(g) = 1.0 / (3.0*nuclear_data->sigma_transport(g));
   }
   
   return 0;
   
}

/* Check the nuclear data: */
int Material::checkNuclearData(int num_energy_groups, bool diffusion, bool transient) const {
   
   /* Check the nuclear data: */
   PAMPA_CHECK(nuclear_data->num_energy_groups != num_energy_groups, 1, 
      "wrong number of energy groups");
   PAMPA_CHECK((nuclear_data->sigma_total).empty(), 2, "missing total cross sections");
   PAMPA_CHECK((nuclear_data->nu_sigma_fission).empty(), 3, "missing nu-fission cross sections");
   PAMPA_CHECK((nuclear_data->kappa_sigma_fission).empty(), 4, 
      "missing kappa-fission cross sections");
   PAMPA_CHECK((nuclear_data->sigma_scattering).empty(), 5, "missing scattering cross sections");
   PAMPA_CHECK((nuclear_data->chi_prompt).empty(), 6, "missing prompt fission spectrum");
   if (diffusion) {
      PAMPA_CHECK((nuclear_data->diffusion_coefficient).empty(), 7, 
         "missing diffusion coefficients");
   }
   if (transient) {
      PAMPA_CHECK((nuclear_data->chi_delayed).empty(), 8, "missing delayed fission spectrum");
      PAMPA_CHECK((nuclear_data->velocity).empty(), 9, "missing neutron velocities");
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
