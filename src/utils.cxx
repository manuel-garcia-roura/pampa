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

/* Read an array with n elements of type double from a file stream: */
int utils::read(Array1D<double>& v, int n, std::ifstream& file) {
   
   /* Read the elements: */
   v.resize(n);
   int l = 0;
   while (l < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), 1, "missing data");
      
      /* Read the elements in this line: */
      for (int i = 0; i < line.size(); i++) {
         PAMPA_CHECK(l >= n, 2, "out-of-bounds data");
         v(l++) = std::stod(line[i]);
      }
      
   }
   
   return 0;
   
}

/* Read an array with n elements of type int from a file stream: */
int utils::read(Array1D<int>& v, int n, std::ifstream& file) {
   
   /* Read the elements: */
   v.resize(n);
   int l = 0;
   while (l < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), 1, "missing data");
      
      /* Read the elements in this line: */
      for (int i = 0; i < line.size(); i++) {
         PAMPA_CHECK(l >= n, 2, "out-of-bounds data");
         v(l++) = std::stoi(line[i]);
      }
      
   }
   
   return 0;
   
}

/* Read an array with (n, m) elements of type int from a file stream: */
int utils::read(Array2D<double>& v, int n, int m, std::ifstream& file) {
   
   /* Read the elements: */
   v.resize(n, m);
   int l = 0;
   while (l < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), 1, "missing data");
      
      /* Read the elements in this line: */
      PAMPA_CHECK(line.size() < m, 1, "missing data");
      PAMPA_CHECK(line.size() > m, 2, "out-of-bounds data");
      for (int i = 0; i < m; i++)
         v(l, i) = std::stod(line[i]);
      l++;
      
   }
   
   return 0;
   
}

/* Read a vector with n rows of type double and total size nt from a file stream: */
int utils::read(Vector2D<double>& v, int n, int nt, std::ifstream& file) {
   
   /* Read the elements: */
   v.reserve(nt);
   int l = 0;
   while (l < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), 1, "missing data");
      PAMPA_CHECK(v.size() + line.size() > nt, 2, "out-of-bounds data");
      
      /* Read the elements in this line: */
      int m = line.size();
      Array1D<double> p(m);
      for (int i = 0; i < line.size(); i++)
         p(i) = std::stod(line[i]);
      v.pushBack(p);
      l++;
      
   }
   
   return 0;
   
}

/* Read a vector with n rows of type int and total size nt from a file stream: */
int utils::read(Vector2D<int>& v, int n, int nt, std::ifstream& file) {
   
   /* Read the elements: */
   v.reserve(nt);
   int l = 0;
   while (l < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), 1, "missing data");
      PAMPA_CHECK(v.size() + line.size() > nt, 2, "out-of-bounds data");
      
      /* Read the elements in this line: */
      int m = line.size();
      Array1D<int> p(m);
      for (int i = 0; i < line.size(); i++)
         p(i) = std::stoi(line[i]);
      v.pushBack(p);
      l++;
      
   }
   
   return 0;
   
}

/* Read a vector with (n1, n2, n3) elements of type double from a file stream: */
int utils::read(Vector3D<double>& v, int n1, const Array1D<int>& n2, int n3, std::ifstream& file) {
   
   /* Read the elements: */
   v.resize(n1, n2, n3);
   for (int i1 = 0; i1 < n1; i1++) {
      for (int i2 = 0; i2 < n2(i1); i2++) {
         
         /* Get the next line: */
         std::vector<std::string> line = get_next_line(file);
         PAMPA_CHECK(line.size() < n3, 1, "missing data");
         PAMPA_CHECK(line.size() > n3, 2, "out-of-bounds data");
         
         /* Read the elements in this line: */
         for (int i3 = 0; i3 < n3; i3++)
            v(i1, i2, i3) = std::stod(line[i3]);
         
      }
   }
   
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
int utils::remove(const std::string& dir) {
   
   /* Remove the directory if it exists: */
   struct stat st = {0};
   if (stat(dir.c_str(), &st) == 0)
      remove(dir.c_str());
   
   return 0;
   
}

/* Create a directory: */
int utils::create(const std::string& dir) {
   
   /* Create the directory: */
   struct stat st = {0};
   if (stat(dir.c_str(), &st) == -1)
      mkdir(dir.c_str(), 0700);
   
   return 0;
   
}
