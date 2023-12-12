#include "CartesianMesh.hxx"

/* The CartesianMesh constructor: */
CartesianMesh::CartesianMesh() {};

/* The CartesianMesh destructor: */
CartesianMesh::~CartesianMesh() {};

/* Read the mesh from a plain-text file: */
bool CartesianMesh::read(const std::string &filename) {
   
   /* Open the input file: */
   std::ifstream file(filename);
   if (!file.is_open()) {
      std::cout << "Error: unable to open " << filename << "!\n";
      return false;
   }
   
   /* Read the file line by line: */
   std::string line;
   while (std::getline(file, line)) {
      
      /* Skip empty lines and #-marked comments: */
      if (line.empty() || line.at(0) == '#')
         continue;
      
      /* Check for dx, dy and dz: */
      std::istringstream iss(line);
      std::string s;
      std::getline(iss, s, ' ');
      if (s == "dx") {
         
         /* Get the dx values: */
         std::getline(iss, s, ' ');
         nx = std::stoi(s);
         if (!utils::read(dx, nx, file)) {
            std::cout << "Error: wrong dx data in " << filename << "!\n";
            return false;
         }
         
      }
      else if (s == "dy") {
         
         /* Get the dy values: */
         std::getline(iss, s, ' ');
         ny = std::stoi(s);
         if (!utils::read(dy, ny, file)) {
            std::cout << "Error: wrong dy data in " << filename << "!\n";
            return false;
         }
         
      }
      else if (s == "dz") {
         
         /* Get the dz values: */
         std::getline(iss, s, ' ');
         nz = std::stoi(s);
         if (!utils::read(dz, nz, file)) {
            std::cout << "Error: wrong dz data in " << filename << "!\n";
            return false;
         }
         
      }
      else {
         std::cout << "Error: wrong keyword in " << filename << "!\n";
         return false;
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
   num_points = (nx+1) * (ny+1) * (nz+1);
   points.reserve(num_points);
   for (int k = 0; k < nz+1; k++) {
      for (int j = 0; j < ny+1; j++) {
         for (int i = 0; i < nx+1; i++) {
            points.push_back(std::vector<double>{x[i], y[j], z[k]});
         }
      }
   }
   
   /* Build the mesh cells: */
   num_cells = nx * ny * nz;
   cells.reserve(num_cells);
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
