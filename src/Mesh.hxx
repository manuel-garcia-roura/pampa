#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "utils.hxx"

/* The bc namespace: */
namespace bc {
  const int vacuum = 1;
  const int reflective = 2;
  const int robin = 3; 
}

/* The Cells struct: */
struct Cells {
   
   /* Cell points: */
   std::vector<std::vector<int>> points;
   
   /* Cell volumes: */
   std::vector<double> volumes;
   
   /* Cell centroids: */
   std::vector<std::vector<double>> centroids;
   
   /* Cell materials: */
   std::vector<int> materials;
   
};

/* The Faces struct: */
struct Faces {
   
   /* Face points: */
   std::vector<std::vector<std::vector<int>>> points;
   
   /* Face areas: */
   std::vector<std::vector<double>> areas;
   
   /* Face centroids: */
   std::vector<std::vector<std::vector<double>>> centroids;
   
   /* Face normals: */
   std::vector<std::vector<std::vector<double>>> normals;
   
   /* Face neighboring cell: */
   std::vector<std::vector<int>> neighbour;
   
};

/* The Mesh class: */
class Mesh {
   
   protected:
      
      /* Mesh dimensions: */
      int num_points, num_cells;
      
      /* Mesh points: */
      std::vector<std::vector<double>> points;
      
      /* Mesh cells: */
      Cells cells;
      
      /* Mesh faces: */
      Faces faces;
   
   public:
      
      /* The Mesh constructor: */
      Mesh();
      
      /* The Mesh destructor: */
      ~Mesh();
      
      /* Read the mesh from a plain-text file: */
      virtual int read(const std::string &filename);
      
      /* Build the mesh: */
      virtual int build();
      
      /* Write the mesh to a plain-text file in .vtk format: */
      int write(const std::string &filename);
   
};
