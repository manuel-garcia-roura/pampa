#include "PrecursorData.hxx"

/* Read the precursor data from a plain-text input file: */
int PrecursorData::read(std::ifstream& file) {
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = input::get_next_line(file);
      if (line.empty() || line[0] == "}") break;
      
      /* Get the next keyword: */
      unsigned int l = 0;
      if (line[l] == "precursor-groups") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the number of delayed-neutron precursor groups: */
         PAMPA_CHECK(input::read(num_precursor_groups, 1, INT_MAX, line[++l]), 
            "wrong number of delayed-neutron precursor groups");
         
      }
      else if (line[l] == "lambda") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 1, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the precursor decay constants: */
         PAMPA_CHECK(input::read(lambda, num_precursor_groups, 0.0, DBL_MAX, file), 
            "wrong precursor decay constants");
         
      }
      else if (line[l] == "beta") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 1, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the precursor fractions: */
         PAMPA_CHECK(input::read(beta, num_precursor_groups, 0.0, DBL_MAX, file), 
            "wrong precursor fractions");
         
         /* Get the total precursor fraction: */
         for (int g = 0; g < num_precursor_groups; g++)
            beta_total += beta(g);
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   return 0;
   
}

/* Check the precursor data to use it in a solver: */
int PrecursorData::check(int num_precursor_groups) const {
   
   /* Check the precursor data: */
   PAMPA_CHECK(this->num_precursor_groups != num_precursor_groups, 
      "wrong number of precursor groups");
   PAMPA_CHECK(lambda.empty(), "missing precursor decay constants");
   PAMPA_CHECK(beta.empty(), "missing precursor fractions");
   
   return 0;
   
}
