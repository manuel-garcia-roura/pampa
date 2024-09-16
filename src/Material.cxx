#include "Material.hxx"

/* Read the material from a plain-text input file: */
int Material::read(const std::string& filename) {
   
   /* Open the input file: */
   std::ifstream file(filename, std::ios_base::in);
   PAMPA_CHECK(!file.is_open(), "unable to open " + filename);
   
   /* Read the material: */
   PAMPA_CHECK(read(file), "unable to read the material from " + filename);
   
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
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Create the constant nuclear data: */
         nuclear_data = new ConstantNuclearData();
         
         /* Read the nuclear data: */
         PAMPA_CHECK(line[++l] != "{", "missing opening '{' for constant nuclear data");
         PAMPA_CHECK(nuclear_data->read(file), "unable to read the constant nuclear data");
         
      }
      else if (line[l] == "nuclear-data-set") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Create the feedback nuclear data: */
         nuclear_data = new FeedbackNuclearData();
         
         /* Read the nuclear data: */
         PAMPA_CHECK(line[++l] != "{", "missing opening '{' for feedback nuclear data");
         PAMPA_CHECK(nuclear_data->read(file), "unable to read the feedback nuclear data");
         
      }
      else if (line[l] == "precursor-data") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Create the precursor data: */
         precursor_data = new PrecursorData();
         
         /* Read the precursor data: */
         PAMPA_CHECK(line[++l] != "{", "missing opening '{' for precursor data");
         PAMPA_CHECK(precursor_data->read(file), "unable to read the precursor data");
         
      }
      else if (line[l] == "thermal-properties") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() < 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Create the thermal properties depending on the type: */
         std::string thermal_properties_type = line[++l];
         if (thermal_properties_type == "constant") {
            PAMPA_CHECK(line.size() != 5, "wrong number of arguments for keyword '" + line[l] + "'");
            double k0, rho0, cp0;
            PAMPA_CHECK(utils::read(k0, 0.0, DBL_MAX, line[++l]), "wrong thermal conductivity");
            PAMPA_CHECK(utils::read(rho0, 0.0, DBL_MAX, line[++l]), "wrong density");
            PAMPA_CHECK(utils::read(cp0, 0.0, DBL_MAX, line[++l]), "wrong specific heat capacity");
            thermal_properties = new ConstantProperties(k0, rho0, cp0);
         }
         else if (thermal_properties_type == "graphite-h-451")
            thermal_properties = new GraphiteProperties();
         else if (thermal_properties_type == "graphite-matrix-a3-27")
            thermal_properties = new GraphiteMatrixProperties();
         else {
            PAMPA_CHECK(true, "wrong thermal-property type");
         }
         
      }
      else if (line[l] == "fuel") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the switch for fuel materials: */
         PAMPA_CHECK(utils::read(fuel, line[++l]), "wrong switch for fuel materials");
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   /* Check the nuclear data after reading it: */
   if (hasNuclearData()) {
      double beta_total = hasPrecursorData() ? precursor_data->beta_total : 0.0;
      PAMPA_CHECK(nuclear_data->check(beta_total), "wrong nuclear data");
   }
   
   return 0;
   
}
