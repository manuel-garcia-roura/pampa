#pragma once

#include "Material.hxx"
#include "vtk.hxx"
#include "input.hxx"
#include "output.hxx"
#include "mpi.hxx"
#include "math.hxx"
#include "utils.hxx"

/* The Cells struct: */
struct Cells {
   
   /* Points: */
   Vector2D<int> points;
   
   /* Volumes: */
   Array1D<double> volumes;
   
   /* Centroids: */
   Array2D<double> centroids;
   
   /* Materials: */
   Array1D<int> materials;
   
   /* Indices in the nodal mesh: */
   Array1D<int> nodal_indices;
   
   /* Indices in the global mesh: */
   Array1D<int> global_indices;
   
};

/* The Faces struct: */
struct Faces {
   
   /* Number of cell faces: */
   Array1D<int> num_faces;
   
   /* Areas: */
   Vector2D<double> areas;
   
   /* Centroids: */
   Vector3D<double> centroids;
   
   /* Normals: */
   Vector3D<double> normals;
   
   /* Neighboring cells (non-negative) or boundary conditions (negative, 1-based): */
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
      
      /* Dimensions: */
      int num_dims = 0, num_points = 0, num_cells = 0, num_ghost_cells = 0, num_faces_max = 0;
      
      /* Global dimensions: */
      int num_cells_global = 0;
      
      /* Points: */
      Array2D<double> points;
      
      /* Cells: */
      Cells cells;
      
      /* Faces: */
      Faces faces;
      
      /* Boundary names: */
      Array1D<std::string> boundaries;
      
      /* Boundary conditions (1-based indexed): */
      Array1D<BoundaryCondition> bcs;
   
   public:
      
      /* The Mesh constructor: */
      Mesh() {}
      
      /* The Mesh destructor: */
      virtual ~Mesh() {}
      
      /* Add a mesh boundary: */
      void addBoundary(const std::string& name) {boundaries.pushBack(name);}
      
      /* Check if the mesh has been partitioned: */
      bool isPartitioned() const {return partitioned;}
      
      /* Get the number of dimensions: */
      int getNumDimensions() const {return num_dims;}
      
      /* Get the number of cells: */
      int getNumCells() const {return num_cells;}
      
      /* Get the maximum number of faces: */
      int getNumFacesMax() const {return num_faces_max;}
      
      /* Get the number of cells in the global mesh: */
      int getNumCellsGlobal() const {return num_cells_global;}
      
      /* Get the mesh cells: */
      const Cells& getCells() const {return cells;}
      
      /* Get the mesh faces: */
      const Faces& getFaces() const {return faces;}
      
      /* Get the mesh boundaries: */
      const Array1D<std::string>& getBoundaries() const {return boundaries;}
      
      /* Get the mesh boundary conditions: */
      const Array1D<BoundaryCondition>& getBoundaryConditions() const {return bcs;}
      
      /* Read the mesh from a plain-text input file: */
      virtual int WARN_UNUSED read(const std::string& filename) {PAMPA_CHECK_VIRTUAL}
      
      /* Build the mesh: */
      virtual int WARN_UNUSED build() {PAMPA_CHECK_VIRTUAL}
      
      /* Remove boundary-condition materials from the mesh: */
      int WARN_UNUSED removeBCMats(const Array1D<Material*>& materials, Mesh** mesh);
      
      /* Partition the mesh: */
      int WARN_UNUSED partition(Mesh** submesh);
      
      /* Write the mesh to a plain-text file in .vtk format: */
      int WARN_UNUSED writeVTK(const std::string& filename) const;
      
      /* Write all the mesh data to a plain-text file: */
      int WARN_UNUSED writeData(const std::string& filename) const;
   
};
