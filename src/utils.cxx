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

/* Get the next line from a file stream:*/
std::string utils::get_next_line() {
   
   std::string line;
   
   return line;
   
};

/* Read a vector with n elements from a file stream: */
bool utils::read(std::vector<double> &v, int n, std::ifstream &file) {
   
   /* Reserve the memory: */
   v.reserve(n);
   
   /* Read the elements: */
   while (v.size() < n) {
      
      /* Get a new line: */
      std::string line;
      if (!std::getline(file, line)) {
         std::cout << "Error: missing vector data!\n";
         return false;
      }
      
      /* Skip empty lines and #-marked comments: */
      if (line.empty() || line.at(0) == '#')
         continue;
      
      /* Read the elements in this line: */
      clean(line);
      std::istringstream iss(line);
      std::string s;
      while (std::getline(iss, s, ' ')) {
         if (v.size() < n)
            v.push_back(std::stod(s));
         else {
            std::cout << "Error: out-of-bounds vector data!\n";
            return false;
         }
      }
      
   }
   
   return true;
   
};

/* Read a vector with n elements of type std::vector<double> of size m from a file stream: */
bool utils::read(std::vector<std::vector<double>> &v, int n, int m, std::ifstream &file) {
   
   /* Reserve the memory: */
   v.reserve(n);
   
   /* Read the elements: */
   while (v.size() < n) {
      
      /* Get a new line: */
      std::string line;
      if (!std::getline(file, line)) {
         std::cout << "Error: missing vector data!\n";
         return false;
      }
      
      /* Skip empty lines and #-marked comments: */
      if (line.empty() || line.at(0) == '#')
         continue;
      
      /* Read the elements in this line: */
      clean(line);
      std::istringstream iss(line);
      std::string s;
      std::vector<double> p;
      p.reserve(m);
      while (std::getline(iss, s, ' ')) {
         if (p.size() < m)
            p.push_back(std::stod(s));
         else {
            std::cout << "Error: out-of-bounds vector data!\n";
            return false;
         }
      }
      v.push_back(p);
      
   }
   
   return true;
   
};

/* Read a vector with n elements of type std::vector<int> of unknown size from a file stream: */
bool utils::read(std::vector<std::vector<int>> &v, int n, std::ifstream &file) {
   
   /* Reserve the memory: */
   v.reserve(n);
   
   /* Read the elements: */
   while (v.size() < n) {
      
      /* Get a new line: */
      std::string line;
      if (!std::getline(file, line)) {
         std::cout << "Error: missing vector data!\n";
         return false;
      }
      
      /* Skip empty lines and #-marked comments: */
      if (line.empty() || line.at(0) == '#')
         continue;
      
      /* Read the elements in this line: */
      clean(line);
      std::istringstream iss(line);
      std::string s;
      std::getline(iss, s, ' ');
      int m = std::stoi(s);
      std::vector<int> p;
      p.reserve(m);
      while (std::getline(iss, s, ' ')) {
         if (p.size() < m)
            p.push_back(std::stoi(s));
         else {
            std::cout << "Error: out-of-bounds vector data!\n";
            return false;
         }
      }
      v.push_back(p);
      
   }
   
   return true;
   
};
