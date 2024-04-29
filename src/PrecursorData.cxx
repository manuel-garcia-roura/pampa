#include "PrecursorData.hxx"

/* Read the precursor data from a plain-text input file: */
int PrecursorData::read(std::ifstream& file) {
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty() || line[0] == "}") break;
      
      /* Get the next keyword: */
      unsigned int l = 0;
      if (line[l] == "precursor-groups") {
         
         /* Get the number of delayed-neutron precursor groups: */
         PAMPA_CALL(utils::read(num_precursor_groups, 1, INT_MAX, line[++l]), 
            "wrong number of delayed-neutron precursor groups");
         
      }
      else if (line[l] == "lambda") {
         
         /* Get the precursor decay constants: */
         PAMPA_CALL(utils::read(lambda, num_precursor_groups, file), 
            "wrong precursor decay constants");
         
      }
      else if (line[l] == "beta") {
         
         /* Get the precursor fractions: */
         PAMPA_CALL(utils::read(beta, num_precursor_groups, file), "wrong precursor fractions");
         
         /* Get the total precursor fraction: */
         for (int g = 0; g < num_precursor_groups; g++)
            beta_total += beta(g);
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 2, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   return 0;
   
}

/* Check the precursor data to use it in a solver: */
int PrecursorData::check(int num_precursor_groups) const {
   
   /* Check the precursor data: */
   PAMPA_CHECK(this->num_precursor_groups != num_precursor_groups, 1, 
      "wrong number of precursor groups");
   PAMPA_CHECK(lambda.empty(), 2, "missing precursor decay constants");
   PAMPA_CHECK(beta.empty(), 3, "missing precursor fractions");
   
   return 0;
   
}
