#include "utils.hxx"

/* Trim and remove double spaces from a string: */
void utils::clean(std::string &s) {
   
   /* Trim: */
   int i1 = s.find_first_not_of(" ");
   int i2 = s.find_last_not_of(" ");
   s.erase(s.begin()+i2+1, s.end());
   s.erase(s.begin(), s.begin()+i1);
   
   /* Remove double spaces: */
   int i = s.find("  ");
   while (i != std::string::npos) {
      s.erase(i, 1);
      i = s.find("  ");
   }
   
};

/* Get the next line from a file stream: */
std::vector<std::string> utils::get_next_line(std::ifstream &file) {
   
   /* Read the file line by line: */
   std::string line;
   std::vector<std::string> words;
   while (std::getline(file, line)) {
      
      /* Skip empty lines and #-marked comments: */
      if (line.empty() || line.at(0) == '#')
         continue;
      
      /* Split the new line: */
      clean(line);
      std::istringstream iss(line);
      std::string s;
      while (std::getline(iss, s, ' '))
         words.push_back(s);
      break;
      
   }
   
   return words;
   
};

/* Read a vector with n elements of type double from a file stream: */
bool utils::read(std::vector<double> &v, int n, std::ifstream &file) {
   
   /* Read the elements: */
   v.reserve(n);
   while (v.size() < n) {
      
      /* Get the next line:*/
      std::vector<std::string> line = get_next_line(file);
      if (line.empty()) {
         std::cout << "Error: missing data!\n";
         return false;
      }
      
      /* Read the elements in this line: */
      for (int i = 0; i < line.size(); i++) {
         if (v.size() < n)
            v.push_back(std::stod(line[i]));
         else {
            std::cout << "Error: out-of-bounds data!\n";
            return false;
         }
      }
      
   }
   
   return true;
   
};

/* Read a vector with n elements of type int from a file stream: */
bool utils::read(std::vector<int> &v, int n, std::ifstream &file) {
   
   /* Read the elements: */
   v.reserve(n);
   while (v.size() < n) {
      
      /* Get the next line:*/
      std::vector<std::string> line = get_next_line(file);
      if (line.empty()) {
         std::cout << "Error: missing data!\n";
         return false;
      }
      
      /* Read the elements in this line: */
      for (int i = 0; i < line.size(); i++) {
         if (v.size() < n)
            v.push_back(std::stoi(line[i]));
         else {
            std::cout << "Error: out-of-bounds data!\n";
            return false;
         }
      }
      
   }
   
   return true;
   
};

/* Read a vector with n elements of type std::vector<double> of size m from a file stream: */
bool utils::read(std::vector<std::vector<double>> &v, int n, int m, std::ifstream &file) {
   
   /* Read the elements: */
   v.reserve(n);
   while (v.size() < n) {
      
      /* Get the next line:*/
      std::vector<std::string> line = get_next_line(file);
      if (line.empty()) {
         std::cout << "Error: missing data!\n";
         return false;
      }
      
      /* Read the elements in this line: */
      std::vector<double> p;
      p.reserve(m);
      for (int i = 0; i < line.size(); i++) {
         if (p.size() < m)
            p.push_back(std::stod(line[i]));
         else {
            std::cout << "Error: out-of-bounds data!\n";
            return false;
         }
      }
      v.push_back(p);
      
   }
   
   return true;
   
};

/* Read a vector with n elements of type std::vector<int> of unknown size from a file stream: */
bool utils::read(std::vector<std::vector<int>> &v, int n, std::ifstream &file) {
   
   /* Read the elements: */
   v.reserve(n);
   while (v.size() < n) {
      
      /* Get the next line:*/
      std::vector<std::string> line = get_next_line(file);
      if (line.empty()) {
         std::cout << "Error: missing data!\n";
         return false;
      }
      
      /* Read the elements in this line: */
      int m = std::stoi(line[0]);
      std::vector<int> p;
      p.reserve(m);
      for (int i = 1; i < line.size(); i++) {
         if (p.size() < m)
            p.push_back(std::stoi(line[i]));
         else {
            std::cout << "Error: out-of-bounds data!\n";
            return false;
         }
      }
      v.push_back(p);
      
   }
   
   return true;
   
};
