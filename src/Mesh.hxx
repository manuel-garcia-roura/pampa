#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "utils.hxx"

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
   
   /* Face neighboring cells (non-negative) or boundary conditions (negative, 1-based): */
   std::vector<std::vector<int>> neighbours;
   
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
      
      /* Boundary conditions (1-based indexed): */
      std::vector<BoundaryCondition> bcs;
   
   public:
      
      /* The Mesh constructor: */
      Mesh() {}
      
      /* The Mesh destructor: */
      ~Mesh() {}
      
      /* Get the number of cells: */
      int getNumCells() const {return num_cells;}
      
      /* Get the mesh cells: */
      const Cells& getCells() const {return cells;}
      
      /* Get the mesh faces: */
      const Faces& getFaces() const {return faces;}
      
      /* Get the mesh boundary conditions: */
      const std::vector<BoundaryCondition>& getBoundaryConditions() const {return bcs;}
      
      /* Read the mesh from a plain-text input file: */
      virtual int read(const std::string& filename);
      
      /* Build the mesh: */
      virtual int build();
      
      /* Write the mesh to a plain-text file in .vtk format: */
      int writeVTK(const std::string& filename) const;
      
      /* Write all the mesh data to a plain-text file: */
      int writeData(const std::string& filename) const;
   
};
