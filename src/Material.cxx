#include "Material.hxx"

/* Read the material from a plain-text input file: */
int Material::read(const std::string& filename) {
   
   /* Open the input file: */
   std::ifstream file(filename, std::ios_base::in);
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
         
         /* Create the constant nuclear data: */
         nuclear_data = new ConstantNuclearData();
         
         /* Read the nuclear data: */
         PAMPA_CHECK(line[++l] != "{", 1, "missing opening '{' for constant nuclear data");
         PAMPA_CALL(nuclear_data->read(file), "unable to read the constant nuclear data");
         
      }
      else if (line[l] == "nuclear-data-set") {
         
         /* Create the feedback nuclear data: */
         nuclear_data = new FeedbackNuclearData();
         
         /* Read the nuclear data: */
         PAMPA_CHECK(line[++l] != "{", 1, "missing opening '{' for feedback nuclear data");
         PAMPA_CALL(nuclear_data->read(file), "unable to read the feedback nuclear data");
         
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
      else if (line[l] == "fuel") {
         
         /* Get the switch for fuel materials: */
         PAMPA_CALL(utils::read(fuel, line[++l]), "wrong switch for fuel materials");
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 2, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   /* Check the nuclear data after reading it: */
   if (hasNuclearData()) {
      double beta_total = hasPrecursorData() ? precursor_data->beta_total : 0.0;
      PAMPA_CALL(nuclear_data->check(beta_total), "wrong nuclear data");
   }
   
   return 0;
   
}
