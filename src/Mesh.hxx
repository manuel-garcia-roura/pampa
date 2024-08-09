#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#ifdef WITH_METIS
#include <metis.h>
#endif

#include "vtk.hxx"
#include "mpi.hxx"
#include "math.hxx"
#include "utils.hxx"

#ifdef WITH_METIS
#define METIS_RECURSIVE 1
#define METIS_KWAY 2
#define METIS_PARTGRAPH METIS_KWAY
#endif

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
   
   /* Cell indices in the global mesh: */
   Array1D<int> indices;
   
};

/* The Faces struct: */
struct Faces {
   
   /* Number of cell faces: */
   Array1D<int> num_faces;
   
   /* Face areas: */
   Vector2D<double> areas;
   
   /* Face centroids: */
   Vector3D<double> centroids;
   
   /* Face normals: */
   Vector3D<double> normals;
   
   /* Face neighboring cells (non-negative) or boundary conditions (negative, 1-based): */
   Vector2D<int> neighbours;
   
};

/* The Mesh class: */
class Mesh {
   
   private:
      
      /* Get the domain indices to partition the mesh: */
      int WARN_UNUSED getDomainIndices(Array1D<int>& part, Array1D<int>& size);
   
   protected:
      
      /* Flag to know if the mesh has been partitioned or not: */
      bool partitioned = false;
      
      /* Mesh dimensions: */
      int num_dims = 0, num_points = 0, num_cells = 0, num_ghost_cells = 0, num_faces_max = 0, 
         num_boundaries = 0;
      
      /* Global mesh dimensions: */
      int num_cells_global = 0;
      
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
      virtual ~Mesh() {}
      
      /* Check if the mesh has been partitioned: */
      bool isPartitioned() const {return partitioned;}
      
      /* Get the number of dimensions: */
      int getNumDimensions() const {return num_dims;}
      
      /* Get the number of cells: */
      int getNumCells() const {return num_cells;}
      
      /* Get the maximum number of faces: */
      int getNumFacesMax() const {return num_faces_max;}
      
      /* Get the number of boundaries: */
      int getNumBoundaries() const {return num_boundaries;}
      
      /* Get the number of cells in the global mesh: */
      int getNumCellsGlobal() const {return num_cells_global;}
      
      /* Get the mesh cells: */
      const Cells& getCells() const {return cells;}
      
      /* Get the mesh faces: */
      const Faces& getFaces() const {return faces;}
      
      /* Get the mesh boundary conditions: */
      const Array1D<BoundaryCondition>& getBoundaryConditions() const {return bcs;}
      
      /* Read the mesh from a plain-text input file: */
      virtual int WARN_UNUSED read(const std::string& filename) {PAMPA_CHECK_VIRTUAL}
      
      /* Build the mesh: */
      virtual int WARN_UNUSED build() {PAMPA_CHECK_VIRTUAL}
      
      /* Partition the mesh: */
      int WARN_UNUSED partition(Mesh** submesh);
      
      /* Write the mesh to a plain-text file in .vtk format: */
      int WARN_UNUSED writeVTK(const std::string& filename) const;
      
      /* Write all the mesh data to a plain-text file: */
      int WARN_UNUSED writeData(const std::string& filename) const;
   
};
