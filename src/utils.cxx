#include "utils.hxx"

/* Trim and remove tabs and double spaces from a string: */
void utils::clean(std::string& s) {
   
   /* Trim: */
   int i1 = s.find_first_not_of(" ");
   int i2 = s.find_last_not_of(" ");
   s.erase(s.begin()+i2+1, s.end());
   s.erase(s.begin(), s.begin()+i1);
   
   /* Remove tabs and double spaces: */
   int i;
   while ((i = s.find("\t")) != std::string::npos)
      s.replace(i, 2, " ");
   while ((i = s.find("  ")) != std::string::npos)
      s.erase(i, 1);
   
}

/* Get the next line from a file stream: */
std::vector<std::string> utils::get_next_line(std::ifstream& file) {
   
   /* Read the file line by line: */
   std::string line;
   std::vector<std::string> words;
   while (std::getline(file, line)) {
      
      /* Skip empty lines and #-marked comments: */
      if (line.empty()) continue;
      if (line[0] == '#' || line[0] == ' ') continue;
      
      /* Split the new line: */
      clean(line);
      std::istringstream iss(line);
      std::string s;
      while (std::getline(iss, s, ' '))
         words.push_back(s);
      break;
      
   }
   
   return words;
   
}

/* Read a vector with n elements of type double from a file stream: */
int utils::read(std::vector<double>& v, int n, std::ifstream& file) {
   
   /* Read the elements: */
   v.reserve(n);
   while (v.size() < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), 1, "missing data");
      
      /* Read the elements in this line: */
      for (int i = 0; i < line.size(); i++) {
         PAMPA_CHECK(v.size() >= n, 2, "out-of-bounds data");
         v.push_back(std::stod(line[i]));
      }
      
   }
   
   return 0;
   
}

/* Read a vector with n elements of type int from a file stream: */
int utils::read(std::vector<int>& v, int n, std::ifstream& file) {
   
   /* Read the elements: */
   v.reserve(n);
   while (v.size() < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), 1, "missing data");
      
      /* Read the elements in this line: */
      for (int i = 0; i < line.size(); i++) {
         PAMPA_CHECK(v.size() >= n, 2, "out-of-bounds data");
         v.push_back(std::stoi(line[i]));
      }
      
   }
   
   return 0;
   
}

/* Read a vector with n elements of type std::vector<double> of size m from a file stream: */
int utils::read(std::vector<std::vector<double>>& v, int n, int m, std::ifstream& file) {
   
   /* Read the elements: */
   v.reserve(n);
   while (v.size() < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), 1, "missing data");
      
      /* Read the elements in this line: */
      PAMPA_CHECK(line.size() < m, 1, "missing data");
      PAMPA_CHECK(line.size() > m, 2, "out-of-bounds data");
      std::vector<double> p(m);
      for (int i = 0; i < line.size(); i++)
         p[i] = std::stod(line[i]);
      v.push_back(p);
      
   }
   
   return 0;
   
}

/* Read a vector with n elements of type std::vector<int> of unknown size from a file stream: */
int utils::read(std::vector<std::vector<int>>& v, int n, std::ifstream& file) {
   
   /* Read the elements: */
   v.reserve(n);
   while (v.size() < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), 1, "missing data");
      
      /* Read the elements in this line: */
      int m = line.size();
      std::vector<int> p(m);
      for (int i = 0; i < line.size(); i++)
         p[i] = std::stoi(line[i]);
      v.push_back(p);
      
   }
   
   return 0;
   
}
