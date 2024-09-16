#include "input.hxx"

/* Trim and remove tabs and double spaces from a string: */
void input::clean(std::string& s) {
   
   /* Check for empty lines: */
   if (s.empty())
      return;
   
   /* Remove tabs and double spaces: */
   size_t i;
   while ((i = s.find("\t")) != std::string::npos)
      s.replace(i, 2, " ");
   while ((i = s.find("  ")) != std::string::npos)
      s.erase(i, 1);
   
   /* Check for empty lines: */
   if (s == " ") {
      s = std::string();
      return;
   }
   
   /* Trim: */
   int i1 = s.find_first_not_of(' ');
   int i2 = s.find_last_not_of(' ');
   s.erase(s.begin()+i2+1, s.end());
   s.erase(s.begin(), s.begin()+i1);
   
   /* Check for #-marked comments: */
   if (s[0] == '#') {
      s = std::string();
      return;
   }
   
}

/* Get the next line from a file stream: */
std::vector<std::string> input::get_next_line(std::ifstream& file) {
   
   /* Read the file line by line: */
   std::string line;
   std::vector<std::string> words;
   while (std::getline(file, line)) {
      
      /* Skip empty lines: */
      clean(line);
      if (line.empty()) continue;
      
      /* Split the new line: */
      std::istringstream iss(line);
      std::string s;
      while (std::getline(iss, s, ' '))
         words.push_back(s);
      break;
      
   }
   
   return words;
   
}

/* Read an array with n elements of type double from a file stream: */
int input::read(Array1D<double>& v, unsigned int n, double x1, double x2, std::ifstream& file) {
   
   /* Read the elements: */
   v.resize(n);
   unsigned int l = 0;
   while (l < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), "missing data");
      
      /* Read the elements in this line: */
      for (unsigned int i = 0; i < line.size(); i++) {
         PAMPA_CHECK(l >= n, "out-of-bounds data");
         PAMPA_CHECK(read(v(l++), x1, x2, line[i]), "wrong data");
      }
      
   }
   
   return 0;
   
}

/* Read an array with n elements of type int from a file stream: */
int input::read(Array1D<int>& v, unsigned int n, int x1, int x2, std::ifstream& file) {
   
   /* Read the elements: */
   v.resize(n);
   unsigned int l = 0;
   while (l < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), "missing data");
      
      /* Read the elements in this line: */
      for (unsigned int i = 0; i < line.size(); i++) {
         PAMPA_CHECK(l >= n, "out-of-bounds data");
         PAMPA_CHECK(read(v(l++), x1, x2, line[i]), "wrong data");
      }
      
   }
   
   return 0;
   
}

/* Read an array with (n, m) elements of type int from a file stream: */
int input::read(Array2D<double>& v, unsigned int n, unsigned int m, double x1, double x2, 
   std::ifstream& file) {
   
   /* Read the elements: */
   v.resize(n, m);
   unsigned int l = 0;
   while (l < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), "missing data");
      
      /* Read the elements in this line: */
      PAMPA_CHECK(line.size() < m, "missing data");
      PAMPA_CHECK(line.size() > m, "out-of-bounds data");
      for (unsigned int i = 0; i < m; i++) {
         PAMPA_CHECK(read(v(l, i), x1, x2, line[i]), "wrong data");
      }
      l++;
      
   }
   
   return 0;
   
}

/* Read a vector with n rows of type double and total size nt from a file stream: */
int input::read(Vector2D<double>& v, unsigned int n, unsigned int nt, double x1, double x2, 
   std::ifstream& file) {
   
   /* Read the elements: */
   v.reserve(nt);
   unsigned int l = 0;
   while (l < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), "missing data");
      PAMPA_CHECK(v.size() + line.size() > nt, "out-of-bounds data");
      
      /* Read the elements in this line: */
      unsigned int m = line.size();
      Array1D<double> p(m);
      for (unsigned int i = 0; i < line.size(); i++) {
         PAMPA_CHECK(read(p(i), x1, x2, line[i]), "wrong data");
      }
      v.pushBack(p);
      l++;
      
   }
   
   return 0;
   
}

/* Read a vector with n rows of type int and total size nt from a file stream: */
int input::read(Vector2D<int>& v, unsigned int n, unsigned int nt, int x1, int x2, 
   std::ifstream& file) {
   
   /* Read the elements: */
   v.reserve(nt);
   unsigned int l = 0;
   while (l < n) {
      
      /* Get the next line: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line.empty(), "missing data");
      PAMPA_CHECK(v.size() + line.size() > nt, "out-of-bounds data");
      
      /* Read the elements in this line: */
      unsigned int m = line.size();
      Array1D<int> p(m);
      for (unsigned int i = 0; i < line.size(); i++) {
         PAMPA_CHECK(read(p(i), x1, x2, line[i]), "wrong data");
      }
      v.pushBack(p);
      l++;
      
   }
   
   return 0;
   
}

/* Read a vector with (n1, n2, n3) elements of type double from a file stream: */
int input::read(Vector3D<double>& v, unsigned int n1, const Array1D<int>& n2, unsigned int n3, 
   double x1, double x2, std::ifstream& file) {
   
   /* Read the elements: */
   v.resize(n1, n2, n3);
   for (unsigned int i1 = 0; i1 < n1; i1++) {
      for (int i2 = 0; i2 < n2(i1); i2++) {
         
         /* Get the next line: */
         std::vector<std::string> line = get_next_line(file);
         PAMPA_CHECK(line.size() < n3, "missing data");
         PAMPA_CHECK(line.size() > n3, "out-of-bounds data");
         
         /* Read the elements in this line: */
         for (unsigned int i3 = 0; i3 < n3; i3++) {
            PAMPA_CHECK(read(v(i1, i2, i3), x1, x2, line[i3]), "wrong data");
         }
         
      }
   }
   
   return 0;
   
}

/* Read a bool value from a string: */
int input::read(bool& q, const std::string& s) {
   
   /* Get the value and check if it's either 0 (false) or 1 (true): */
   int x = std::stoi(s);
   PAMPA_CHECK(x != 0 && x != 1, "wrong bool value");
   q = x;
   
   return 0;
   
}

/* Read an int value from a string: */
int input::read(int& x, int x1, int x2, const std::string& s) {
   
   /* Get the value and check if it's within the limits: */
   x = std::stoi(s);
   PAMPA_CHECK(x < x1 || x > x2, "out-of-bounds int value");
   
   return 0;
   
}

/* Read a double value from a string: */
int input::read(double& x, double x1, double x2, const std::string& s) {
   
   /* Get the value and check if it's within the limits: */
   x = std::stod(s);
   PAMPA_CHECK(x < x1 || x > x2, "out-of-bounds double value");
   
   return 0;
   
}

/* Read a function from a file stream: */
int input::read(Function& f, double x1, double x2, const std::vector<std::string>& line, 
   unsigned int& i, std::ifstream& file) {
   
   /* Read the function depending on the number of arguments: */
   if (line.size() == i+1) {
      
      /* Read a single-value function: */
      double x;
      PAMPA_CHECK(read(x, x1, x2, line[i]), "wrong function data");
      f = Function(x);
      
   }
   else if (line.size() == i+2) {
      
      /* Get the number of time points: */
      int n;
      PAMPA_CHECK(read(n, 1, INT_MAX, line[i]), "wrong number of time points");
      
      /* Check the start of the function data: */
      PAMPA_CHECK(line[++i] != "{", "missing opening '{' for function data");
      
      /* Get the function data: */
      Array1D<double> t, x;
      PAMPA_CHECK(read(t, n, 0.0, DBL_MAX, file), "wrong time data");
      PAMPA_CHECK(read(x, n, x1, x2, file), "wrong function data");
      f = Function(t, x);
      
      /* Check the end of the function data: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line[0] != "}", "missing closing '}' for function data");
      
   }
   else {
      
      /* Wrong number of arguments: */
      PAMPA_CHECK(true, "wrong number of function arguments");
      
   }
   
   return 0;
   
}

/* Read multiple functions from a file stream: */
int input::read(Array1D<Function>& f, int m, double x1, double x2, 
   const std::vector<std::string>& line, unsigned int& i, std::ifstream& file) {
   
   /* Read the functions depending on the number of arguments: */
   if ((line.size() == i+m) && !((line.size() == i+2) && (line[i+1] == "{"))) {
      
      /* Read m single-value functions: */
      for (int j = 0; j < m; j++) {
         double x;
         PAMPA_CHECK(read(x, x1, x2, line[i++]), "wrong function data");
         f(j) = Function(x);
      }
      
   }
   else if (line.size() == i+2) {
      
      /* Get the number of time points: */
      int n;
      PAMPA_CHECK(read(n, 1, INT_MAX, line[i]), "wrong number of time points");
      
      /* Check the start of the function data: */
      PAMPA_CHECK(line[++i] != "{", "missing opening '{' for function data");
      
      /* Get the function data: */
      Array1D<double> t;
      PAMPA_CHECK(read(t, n, 0.0, DBL_MAX, file), "wrong time data");
      for (int j = 0; j < m; j++) {
         Array1D<double> x;
         PAMPA_CHECK(read(x, n, x1, x2, file), "wrong function data");
         f(j) = Function(t, x);
      }
      
      /* Check the end of the function data: */
      std::vector<std::string> line = get_next_line(file);
      PAMPA_CHECK(line[0] != "}", "missing closing '}' for function data");
      
   }
   else {
      
      /* Wrong number of arguments: */
      PAMPA_CHECK(true, "wrong number of function arguments");
      
   }
   
   return 0;
   
}

/* Read a boundary condition from a line: */
int input::read(BoundaryCondition& bc, const std::vector<std::string>& line, unsigned int& i, 
   std::ifstream& file) {
   
   /* Get the boundary-condition type: */
   std::string bc_type = line[i++];
   if (bc_type == "vacuum")
      bc.type = BC::VACUUM;
   else if (bc_type == "reflective")
      bc.type = BC::REFLECTIVE;
   else if (bc_type == "robin")
      bc.type = BC::ROBIN;
   else if (bc_type == "dirichlet")
      bc.type = BC::DIRICHLET;
   else if (bc_type == "adiabatic")
      bc.type = BC::ADIABATIC;
   else if (bc_type == "convection")
      bc.type = BC::CONVECTION;
   else {
      PAMPA_CHECK(true, "wrong boundary-condition type");
   }
   
   /* Get a single parameter for Robin or Dirichlet boundary conditions: */
   if (bc.type == BC::ROBIN || bc.type == BC::DIRICHLET) {
      bc.f.resize(1);
      PAMPA_CHECK(read(bc.f(0), -DBL_MAX, DBL_MAX, line, i, file), 
         "wrong boundary-condition parameter");
   }
   
   /* Get two parameters for convection boundary conditions: */
   if (bc.type == BC::CONVECTION) {
      bc.f.resize(2);
      PAMPA_CHECK(read(bc.f, 2, -DBL_MAX, DBL_MAX, line, i, file), 
         "wrong boundary-condition parameter");
   }
   
   return 0;
   
}
