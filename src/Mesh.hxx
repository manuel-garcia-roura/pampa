#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>

/* The Mesh class: */
class Mesh {
   
   protected:
      
      /* Mesh dimensions: */
      int num_points, num_cells;
      
      /* Mesh points: */
      std::vector<std::vector<double>> points;
      
      /* Mesh cells: */
      std::vector<std::vector<int>> cells;
   
   public:
      
      /* The Mesh constructor: */
      Mesh();
      
      /* The Mesh destructor: */
      ~Mesh();
      
      /* Read the mesh from a plain-text file: */
      virtual bool read(const std::string &filename);
      
      /* Build the mesh: */
      virtual bool build();
      
      /* Write the mesh to a plain-text file in .vtk format: */
      bool write(const std::string &filename);
   
};
