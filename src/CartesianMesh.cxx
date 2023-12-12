#include "CartesianMesh.hxx"

/* The CartesianMesh constructor: */
CartesianMesh::CartesianMesh() {};

/* The CartesianMesh destructor: */
CartesianMesh::~CartesianMesh() {};

/* Read the mesh from a plain-text file: */
bool CartesianMesh::read(const std::string &filename) {
   
   /* Open the input file: */
   std::ifstream file(filename);
   
   /* Read the file line by line: */
   std::string line;
   std::vector<double> *d;
   int i, n = -1;
   while (std::getline(file, line)) {
      
      /* Skip empty lines and #-marked comments: */
      if (line.empty() || line.at(0) == '#')
         continue;
      
      /* Read the line word by word: */
      std::istringstream iss(line);
      std::string s;
      while (std::getline(iss, s, ' ')) {
         
         /* Check for dx, dy and dz to manage the data: */
         if (s == "dx" || s == "dy" || s == "dz") {
            if (n > 0 && i != n) {
               std::cout << "Error: wrong discretization!\n";
               return false;
            }
            i = 0;
            if (s == "dx") {
               std::getline(iss, s, ' ');
               nx = std::stoi(s);
               dx.reserve(nx);
               d = &dx;
               n = nx;
            }
            else if (s == "dy") {
               std::getline(iss, s, ' ');
               ny = std::stoi(s);
               dy.reserve(ny);
               d = &dy;
               n = ny;
            }
            else if (s == "dz") {
               std::getline(iss, s, ' ');
               nz = std::stoi(s);
               dz.reserve(nz);
               d = &dz;
               n = nz;
            }
         }
         
         /* Get the dx, dy or dz values: */
         else {
            if (i < n)
               (*d)[i++] = std::stod(s);
            else {
               std::cout << "Error: wrong discretization!\n";
               return false;
            }
         }
         
      }
      
   }
   
   return true;
   
};

/* Build the mesh: */
bool CartesianMesh::build() {
   
   /* Build the mesh points: */
   std::vector<double> x(nx+1), y(ny+1), z(nz+1);
   x[0] = 0.0;
   for (int i = 0; i < nx; i++)
      x[i+1] = x[i] + dx[i];
   y[0] = 0.0;
   for (int j = 0; j < ny; j++)
      y[j+1] = y[j] + dy[j];
   z[0] = 0.0;
   for (int k = 0; k < nz; k++)
      z[k+1] = z[k] + dz[k];
   points.reserve((nx+1)*(ny+1)*(nz+1));
   for (int k = 0; k < nz+1; k++) {
      for (int j = 0; j < ny+1; j++) {
         for (int i = 0; i < nx+1; i++) {
            points.push_back(std::array<double, 3>{x[i], y[j], z[k]});
         }
      }
   }
   
   /* Build the mesh cells: */
   cells.reserve(nx*ny*nz);
   for (int k = 0; k < nz; k++) {
      for (int j = 0; j < ny; j++) {
         for (int i = 0; i < nx; i++) {
            int p1 = i + j*(nx+1) + k*(nx+1)*(ny+1);
            int p2 = (i+1) + j*(nx+1) + k*(nx+1)*(ny+1);
            int p3 = (i+1) + (j+1)*(nx+1) + k*(nx+1)*(ny+1);
            int p4 = i + (j+1)*(nx+1) + k*(nx+1)*(ny+1);
            int p5 = i + j*(nx+1) + (k+1)*(nx+1)*(ny+1);
            int p6 = (i+1) + j*(nx+1) + (k+1)*(nx+1)*(ny+1);
            int p7 = (i+1) + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
            int p8 = i + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
            cells.push_back(std::vector<int>{p1, p2, p3, p4, p5, p6, p7, p8});
         }
      }
   }
   
   return true;
   
};
