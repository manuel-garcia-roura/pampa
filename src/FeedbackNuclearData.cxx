#include "FeedbackNuclearData.hxx"

/* Read the nuclear data from a plain-text input file: */
int FeedbackNuclearData::read(std::ifstream& file) {
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = input::get_next_line(file);
      if (line.empty() || line[0] == "}") break;
      
      /* Get the next keyword: */
      unsigned int l = 0;
      if (line[l] == "temperature") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the reference temperature: */
         int n;
         PAMPA_CHECK(input::read(n, 0, INT_MAX, line[++l]), "wrong number of temperatures");
         PAMPA_CHECK(input::read(Tref, n, 0.0, DBL_MAX, file), "wrong temperature data");
         
      }
      else if (line[l] == "nuclear-data") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Create the nuclear data: */
         ConstantNuclearData* data = new ConstantNuclearData();
         
         /* Read the nuclear data: */
         PAMPA_CHECK(line[++l] != "{", "missing opening '{' for nuclear data");
         PAMPA_CHECK(data->read(file), "unable to read the nuclear data");
         
         /* Keep the nuclear data definition: */
         nuclear_data.pushBack(data);
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   return 0;
   
}
