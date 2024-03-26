#include "utils.hxx"

/* Trim and remove tabs and double spaces from a string: */
void utils::clean(std::string& s) {
   
   /* Trim: */
   int i1 = s.find_first_not_of(" ");
   int i2 = s.find_last_not_of(" ");
   s.erase(s.begin()+i2+1, s.end());
   s.erase(s.begin(), s.begin()+i1);
   
   /* Remove tabs and double spaces: */
   size_t i;
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

/* Read an array with n elements of type double from a file stream: */
int utils::read(Array1D<double>& v, unsigned int n, std::ifstream& file) {
   
   /* Read the elements: */
   v.resize(n);
   unsigned int l = 0;
   while (l < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), 1, "missing data");
      
      /* Read the elements in this line: */
      for (unsigned int i = 0; i < line.size(); i++) {
         PAMPA_CHECK(l >= n, 2, "out-of-bounds data");
         v(l++) = std::stod(line[i]);
      }
      
   }
   
   return 0;
   
}

/* Read an array with n elements of type int from a file stream: */
int utils::read(Array1D<int>& v, unsigned int n, std::ifstream& file) {
   
   /* Read the elements: */
   v.resize(n);
   unsigned int l = 0;
   while (l < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), 1, "missing data");
      
      /* Read the elements in this line: */
      for (unsigned int i = 0; i < line.size(); i++) {
         PAMPA_CHECK(l >= n, 2, "out-of-bounds data");
         v(l++) = std::stoi(line[i]);
      }
      
   }
   
   return 0;
   
}

/* Read an array with (n, m) elements of type int from a file stream: */
int utils::read(Array2D<double>& v, unsigned int n, unsigned int m, std::ifstream& file) {
   
   /* Read the elements: */
   v.resize(n, m);
   unsigned int l = 0;
   while (l < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), 1, "missing data");
      
      /* Read the elements in this line: */
      PAMPA_CHECK(line.size() < m, 1, "missing data");
      PAMPA_CHECK(line.size() > m, 2, "out-of-bounds data");
      for (unsigned int i = 0; i < m; i++)
         v(l, i) = std::stod(line[i]);
      l++;
      
   }
   
   return 0;
   
}

/* Read a vector with n rows of type double and total size nt from a file stream: */
int utils::read(Vector2D<double>& v, unsigned int n, unsigned int nt, std::ifstream& file) {
   
   /* Read the elements: */
   v.reserve(nt);
   unsigned int l = 0;
   while (l < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), 1, "missing data");
      PAMPA_CHECK(v.size() + line.size() > nt, 2, "out-of-bounds data");
      
      /* Read the elements in this line: */
      unsigned int m = line.size();
      Array1D<double> p(m);
      for (unsigned int i = 0; i < line.size(); i++)
         p(i) = std::stod(line[i]);
      v.pushBack(p);
      l++;
      
   }
   
   return 0;
   
}

/* Read a vector with n rows of type int and total size nt from a file stream: */
int utils::read(Vector2D<int>& v, unsigned int n, unsigned int nt, std::ifstream& file) {
   
   /* Read the elements: */
   v.reserve(nt);
   unsigned int l = 0;
   while (l < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), 1, "missing data");
      PAMPA_CHECK(v.size() + line.size() > nt, 2, "out-of-bounds data");
      
      /* Read the elements in this line: */
      unsigned int m = line.size();
      Array1D<int> p(m);
      for (unsigned int i = 0; i < line.size(); i++)
         p(i) = std::stoi(line[i]);
      v.pushBack(p);
      l++;
      
   }
   
   return 0;
   
}

/* Read a vector with (n1, n2, n3) elements of type double from a file stream: */
int utils::read(Vector3D<double>& v, unsigned int n1, const Array1D<int>& n2, unsigned int n3, 
   std::ifstream& file) {
   
   /* Read the elements: */
   v.resize(n1, n2, n3);
   for (unsigned int i1 = 0; i1 < n1; i1++) {
      for (int i2 = 0; i2 < n2(i1); i2++) {
         
         /* Get the next line: */
         std::vector<std::string> line = get_next_line(file);
         PAMPA_CHECK(line.size() < n3, 1, "missing data");
         PAMPA_CHECK(line.size() > n3, 2, "out-of-bounds data");
         
         /* Read the elements in this line: */
         for (unsigned int i3 = 0; i3 < n3; i3++)
            v(i1, i2, i3) = std::stod(line[i3]);
         
      }
   }
   
   return 0;
   
}

/* Read a bool value from a string: */
int utils::read(bool& q, const std::string& s) {
   
   /* Get the value and check if it's either 0 (false) or 1 (true): */
   int x = std::stoi(s);
   PAMPA_CHECK(x != 0 && x != 1, 1, "wrong bool value");
   q = x;
   
   return 0;
   
}

/* Read an int value from a string: */
int utils::read(int& x, double x1, double x2, const std::string& s) {
   
   /* Get the value and check if it's within the limits: */
   x = std::stoi(s);
   PAMPA_CHECK(x < x1 || x > x2, 1, "out-of-bounds int value");
   
   return 0;
   
}

/* Read a double value from a string: */
int utils::read(double& x, double x1, double x2, const std::string& s) {
   
   /* Get the value and check if it's within the limits: */
   x = std::stod(s);
   PAMPA_CHECK(x < x1 || x > x2, 1, "out-of-bounds double value");
   
   return 0;
   
}

/* Read a boundary condition from a line: */
int utils::read(BoundaryCondition& bc, const std::vector<std::string>& line, int& i) {
   
   /* Get the boundary condition type: */
   bc.type = static_cast<BC::Type>(std::stoi(line[i++])-1);
   
   /* Get the albedo factor for Robin boundary conditions: */
   if (bc.type == BC::ROBIN) bc.a = std::stod(line[i++]);
   
   /* Get the fixed value for Dirichlet boundary conditions: */
   if (bc.type == BC::DIRICHLET) bc.x = std::stod(line[i++]);
   
   return 0;
   
}

/* Remove a directory: */
void utils::remove_directory(const std::string& dir) {
   
   /* Remove the directory if it exists: */
   struct stat st = {0};
   if (stat(dir.c_str(), &st) == 0)
      remove(dir.c_str());
   
}

/* Create a directory: */
void utils::create_directory(const std::string& dir) {
   
   /* Create the directory if it doesn't exist: */
   struct stat st = {0};
   if (stat(dir.c_str(), &st) == -1)
      mkdir(dir.c_str(), 0700);
   
}
