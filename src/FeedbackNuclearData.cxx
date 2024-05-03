#include "FeedbackNuclearData.hxx"

/* Read the nuclear data from a plain-text input file: */
int FeedbackNuclearData::read(std::ifstream& file) {
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty() || line[0] == "}") break;
      
      /* Get the next keyword: */
      unsigned int l = 0;
      if (line[l] == "temperature") {
         
         /* Get the reference temperature: */
         int n;
         PAMPA_CALL(utils::read(n, 0, INT_MAX, line[++l]), "wrong number of temperatures");
         PAMPA_CALL(utils::read(Tref, n, file), "wrong temperature data");
         
      }
      else if (line[l] == "nuclear-data") {
         
         /* Create the nuclear data: */
         ConstantNuclearData* data = new ConstantNuclearData();
         
         /* Read the nuclear data: */
         PAMPA_CHECK(line[++l] != "{", 1, "missing opening '{' for nuclear data");
         PAMPA_CALL(data->read(file), "unable to read the nuclear data");
         
         /* Keep the nuclear data definition: */
         nuclear_data.pushBack(data);
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 2, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   return 0;
   
}
