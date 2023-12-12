#include "UnstructuredExtrudedMesh.hxx"

/* The UnstructuredExtrudedMesh constructor: */
UnstructuredExtrudedMesh::UnstructuredExtrudedMesh() {};

/* The UnstructuredExtrudedMesh destructor: */
UnstructuredExtrudedMesh::~UnstructuredExtrudedMesh() {};

/* Read the mesh from a plain-text file: */
bool UnstructuredExtrudedMesh::read(const std::string &filename) {
   
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
      
      /* Check for points, cells and dz: */
      std::istringstream iss(line);
      std::string s;
      std::getline(iss, s, ' ');
      if (s == "points") {
         
         /* Get the point coordinates: */
         std::getline(iss, s, ' ');
         num_xy_points = std::stoi(s);
         if (!utils::read(xy_points, num_xy_points, 2, file)) {
            std::cout << "Error: wrong point data in " << filename << "!\n";
            return false;
         }
         
      }
      else if (s == "cells") {
         
         /* Get the cell indices: */
         std::getline(iss, s, ' ');
         num_xy_cells = std::stoi(s);
         if (!utils::read(xy_cells, num_xy_cells, file)) {
            std::cout << "Error: wrong cell data in " << filename << "!\n";
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
bool UnstructuredExtrudedMesh::build() {
   
   /* Build the mesh points: */
   std::vector<double> z(nz+1);
   z[0] = 0.0;
   for (int k = 0; k < nz; k++)
      z[k+1] = z[k] + dz[k];
   num_points = num_xy_points * (nz+1);
   points.reserve(num_points);
   for (int k = 0; k < nz+1; k++) {
      for (int i = 0; i < num_xy_points; i++) {
         points.push_back(std::vector<double>{xy_points[i][0], xy_points[i][1], z[k]});
      }
   }
   
   /* Build the mesh cells: */
   num_cells = num_xy_cells * nz;
   cells.reserve(num_cells);
   for (int k = 0; k < nz; k++) {
      for (int i = 0; i < num_xy_cells; i++) {
         int num_xy_indices = xy_cells[i].size();
         std::vector<int> cell;
         cell.reserve(2*num_xy_indices);
         for (int dk = 0; dk < 2; dk++) {
            for (int l = 0; l < num_xy_indices; l++) {
               cell.push_back(xy_cells[i][l]+(k+dk)*num_xy_points);
            }
         }
         cells.push_back(cell);
      }
   }
   
   return true;
   
};
