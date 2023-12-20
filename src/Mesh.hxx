#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "utils.hxx"
#include "Material.hxx"

/* The bc namespace: */
namespace bc {
  const int vacuum = -1;
  const int reflective = -2;
  const int robin = -3; 
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
      
      /* Model materials: */
      const std::vector<Material> *materials;
   
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
      
      /* Set the model materials: */
      void setModelMaterials(const std::vector<Material> *materials);
      
      /* Get the number of cells: */
      int getNumCells() const;
      
      /* Get the material for cell i: */
      int getMaterial(int i) const;
      
      /* Get the integral of the total cross section for cell i and group g: */
      double getVolumeIntegralSigmaTotal(int i, int g) const;
      
      /* Get the integral of the removal cross section for cell i and group g: */
      double getVolumeIntegralSigmaRemoval(int i, int g) const;
      
      /* Get the integral of the absorption cross section for cell i and group g: */
      double getVolumeIntegralSigmaAbsorption(int i, int g) const;
      
      /* Get the integral of the nu-fission cross section for cell i and group g: */
      double getVolumeIntegralSigmaNuFission(int i, int g) const;
      
      /* Get the integral of the scattering cross section for cell i and group g: */
      double getVolumeIntegralSigmaScattering(int i, int g, int g2) const;
      
      /* Get the integral of the diffusion coeff. for face f of cell i and group g: */
      double getSurfaceIntegralDiffusionCoefficient(int i, int f, int g) const;
      
      /* Get the integral of the diffusion coeff. of cell i2 for face f of cell i and group g: */
      double getSurfaceIntegralDiffusionCoefficient(int i2, int i, int f, int g) const;
      
      /* Get the fission spectrum for cell i and group g: */
      double getChi(int i, int g) const;
      
      /* Get the centroid of cell i: */
      const std::vector<double>& getCellCentroid(int i) const;
      
      /* Get the centroid of face j of cell i: */
      const std::vector<double>& getFaceCentroid(int i, int j) const;
      
      /* Get the normal of face j of cell i: */
      const std::vector<double>& getFaceNormal(int i, int j) const;
      
      /* Get the neighboring cells of all faces of cell i: */
      const std::vector<int>& getFaceNeighbours(int i) const;
   
};
