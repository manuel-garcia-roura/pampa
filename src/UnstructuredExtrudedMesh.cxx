#include "UnstructuredExtrudedMesh.hxx"

/* The UnstructuredExtrudedMesh constructor: */
UnstructuredExtrudedMesh::UnstructuredExtrudedMesh() {};

/* The UnstructuredExtrudedMesh destructor: */
UnstructuredExtrudedMesh::~UnstructuredExtrudedMesh() {};

/* Read the mesh from a plain-text file: */
bool UnstructuredExtrudedMesh::read(const std::string &filename) {
   
   /* Open the input file: */
   std::ifstream file(filename);
   
   /* Read the file line by line: */
   std::string line;
   while (std::getline(file, line)) {
      
      /* Skip empty lines and #-marked comments: */
      if (line.empty() || line.at(0) == '#')
         continue;
      
      /* Check for points, cells and dz to manage the data: */
      std::istringstream iss(line);
      std::string s;
      std::getline(iss, s, ' ');
      if (s == "points") {
         
         /* Get the coordinates for each point: */
         std::getline(iss, s, ' ');
         num_xy_points = std::stoi(s);
         xy_points.reserve(num_xy_points);
         while (xy_points.size() < num_xy_points) {
            
            /* Skip empty lines and #-marked comments: */
            if (line.empty() || line.at(0) == '#')
               continue;
            
            /* Get the coordinates for this point: */
            if (!std::getline(file, line))
               return false;
            std::istringstream iss2(line);
            double x, y;
            if (!(iss2 >> x >> y))
               return false;
            xy_points.push_back(std::array<double, 2>{x, y});
            
         }
         
      }
      else if (s == "cells") {
         
         /* Get the indices for each cell: */
         std::getline(iss, s, ' ');
         num_xy_cells = std::stoi(s);
         xy_cells.reserve(num_xy_cells);
         while (xy_cells.size() < num_xy_cells) {
            
            /* Skip empty lines and #-marked comments: */
            if (line.empty() || line.at(0) == '#')
               continue;
            
            /* Get the indices for this cell: */
            if (!std::getline(file, line))
               return false;
            std::istringstream iss2(line);
            std::getline(iss2, s, ' ');
            int n = std::stoi(s);
            std::vector<int> cell(n);
            for (int i = 0; i < n; i++)
               if (!(iss2 >> cell[i]))
                  return false;
            xy_cells.push_back(cell);
            
         }
         
      }
      else if (s == "dz") {
         
         /* Get the dz values: */
         std::getline(iss, s, ' ');
         nz = std::stoi(s);
         dz.reserve(nz);
         while (dz.size() < nz) {
            
            /* Skip empty lines and #-marked comments: */
            if (line.empty() || line.at(0) == '#')
               continue;
            
            /* Get the dz values for this line: */
            if (!std::getline(file, line))
               return false;
            std::istringstream iss2(line);
            while (std::getline(iss2, s, ' ')) {
               if (dz.size() < nz)
                  dz.push_back(std::stod(s));
               else {
                  std::cout << "Error: wrong discretization!\n";
                  return false;
               }
            }
            
         }
         
      }
      else {
         return false;
      }
      
   }
   
   return true;
   
};

/* Build the mesh: */
bool UnstructuredExtrudedMesh::build() {
   
   /* Build the mesh points: */
   std::vector<double> z(nz+1);
   z[0] = 0.0;
   for (int k = 0; k < nz; k++)
      z[k+1] = z[k] + dz[k];
   points.reserve(num_xy_points*(nz+1));
   for (int k = 0; k < nz+1; k++) {
      for (int i = 0; i < num_xy_points; i++) {
         points.push_back(std::array<double, 3>{xy_points[i][0], xy_points[i][1], z[k]});
      }
   }
   
   /* Build the mesh cells: */
   cells.reserve(num_xy_cells*nz);
   for (int k = 0; k < nz; k++) {
      for (int i = 0; i < num_xy_cells; i++) {
         int n = xy_cells[i].size();
         std::vector<int> cell;
         cell.reserve(2*n);
         for (int dk = 0; dk < 2; dk++) {
            for (int l = 0; l < n; l++) {
               cell.push_back(xy_cells[i][l]+(k+dk)*num_xy_points);
            }
         }
         cells.push_back(cell);
      }
   }
   
   return true;
   
};
