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
   Vector2D<int> points;
   
   /* Cell volumes: */
   Array1D<double> volumes;
   
   /* Cell centroids: */
   Array2D<double> centroids;
   
   /* Cell materials: */
   Array1D<int> materials;
   
};

/* The Faces struct: */
struct Faces {
   
   /* Number of cell faces: */
   Array1D<int> num_faces;
   
   /* Face areas: */
   Array2D<double> areas;
   
   /* Face centroids: */
   Array3D<double> centroids;
   
   /* Face normals: */
   Array3D<double> normals;
   
   /* Face neighboring cells (non-negative) or boundary conditions (negative, 1-based): */
   Array2D<int> neighbours;
   
};

/* The Mesh class: */
class Mesh {
   
   protected:
      
      /* Mesh dimensions: */
      int num_points, num_cells;
      
      /* Mesh points: */
      Array2D<double> points;
      
      /* Mesh cells: */
      Cells cells;
      
      /* Mesh faces: */
      Faces faces;
      
      /* Boundary conditions (1-based indexed): */
      Array1D<BoundaryCondition> bcs;
   
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
      const Array1D<BoundaryCondition>& getBoundaryConditions() const {return bcs;}
      
      /* Read the mesh from a plain-text input file: */
      virtual int read(const std::string& filename) 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Build the mesh: */
      virtual int build() 
         {PAMPA_CHECK(true, 1, "virtual method called on the base class"); return 1;}
      
      /* Write the mesh to a plain-text file in .vtk format: */
      int writeVTK(const std::string& filename) const;
      
      /* Write all the mesh data to a plain-text file: */
      int writeData(const std::string& filename) const;
   
};
