#include "Mesh.hxx"

/* The Mesh constructor: */
Mesh::Mesh() : materials(NULL) {};

/* The Mesh destructor: */
Mesh::~Mesh() {};

/* Read the mesh from a plain-text file: */
int Mesh::read(const std::string &filename) {
   
   /* This method should never be called: */
   PAMPA_CHECK(true, 1, "Mesh.read() method called on the base class");
   
   return 0;
   
};

/* Build the mesh: */
int Mesh::build() {
   
   /* This method should never be called: */
   PAMPA_CHECK(true, 1, "Mesh.build() method called on the base class");
   
   return 0;
   
};

/* Set the model materials: */
void Mesh::setModelMaterials(const std::vector<Material> *materials) {
   
   /* Keep a pointer to the materials: */
   this->materials = materials;
   
};

/* Get the number of cells: */
int Mesh::getNumCells() const {
   
   return num_cells;
   
};

/* Get the boundary conditions: */
const std::vector<BoundaryCondition>& Mesh::getBoundaryConditions() const {
   
   return bcs;
   
};

/* Get the volume of cell i: */
double Mesh::getCellVolume(int i) const {
   
   return cells.volumes[i];
   
};

/* Get the centroid of cell i: */
const std::vector<double>& Mesh::getCellCentroid(int i) const {
   
   return cells.centroids[i];
   
};

/* Get the area of face j of cell i: */
double Mesh::getFaceArea(int i, int j) const {
   
   return faces.areas[i][j];
   
};

/* Get the centroid of face j of cell i: */
const std::vector<double>& Mesh::getFaceCentroid(int i, int j) const {
   
   return faces.centroids[i][j];
   
};

/* Get the normal of face j of cell i: */
const std::vector<double>& Mesh::getFaceNormal(int i, int j) const {
   
   return faces.normals[i][j];
   
};

/* Get the neighboring cells for all faces of cell i: */
const std::vector<int>& Mesh::getFaceNeighbours(int i) const {
   
   return faces.neighbours[i];
   
};

/* Get the material for cell i: */
int Mesh::getMaterial(int i) const {
   
   return cells.materials[i];
   
};

/* Get the integral of the total cross section for cell i and group g: */
double Mesh::getVolIntSigmaTotal(int i, int g) const {
   
   /* Calculate the integral from the mean value and the volume: */
   return (*materials)[cells.materials[i]].sigma_total[g] * cells.volumes[i];
   
};

/* Get the integral of the nu-fission cross section for cell i and group g: */
double Mesh::getVolIntNuSigmaFission(int i, int g) const {
   
   /* Calculate the integral from the mean value and the volume: */
   return (*materials)[cells.materials[i]].nu_sigma_fission[g] * cells.volumes[i];
   
};

/* Get the integral of the scattering cross section for cell i and group g: */
double Mesh::getVolIntSigmaScattering(int i, int g, int g2) const {
   
   /* Calculate the integral from the mean value and the volume: */
   return (*materials)[cells.materials[i]].sigma_scattering[g][g2] * cells.volumes[i];
   
};

/* Get the integral of the removal cross section for cell i and group g: */
double Mesh::getVolIntSigmaRemoval(int i, int g) const {
   
   /* Calculate the integral from the mean value and the volume: */
   return (*materials)[cells.materials[i]].sigma_removal[g] * cells.volumes[i];
   
};

/* Get the integral of the diffusion coeff. for face f of cell i and group g: */
double Mesh::getSurfIntDiffusionCoefficient(int i, int f, int g) const {
   
   /* Calculate the integral from the mean value and the area: */
   return (*materials)[cells.materials[i]].diffusion_coefficient[g] * faces.areas[i][f];
   
};

/* Get the integral of the diffusion coeff. of cell i2 for face f of cell i and group g: */
double Mesh::getSurfIntDiffusionCoefficient(int i2, int i, int f, int g) const {
   
   /* Calculate the integral from the mean value and the area: */
   return (*materials)[cells.materials[i2]].diffusion_coefficient[g] * faces.areas[i][f];
   
};

/* Get the fission spectrum for cell i and group g: */
double Mesh::getChi(int i, int g) const {
   
   return (*materials)[cells.materials[i]].chi[g];
   
};

/* Write the mesh to a plain-text file in .vtk format: */
int Mesh::write(const std::string &filename) {
   
   /* Open the output file: */
   std::ofstream file(filename);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Set the precision for the point coordinates: */
   file << std::fixed;
   file << std::setprecision(3);
   
   /* Write the header: */
   file << "# vtk DataFile Version 3.0" << std::endl;
   file << "FVM mesh" << std::endl;
   file << "ASCII" << std::endl;
   file << "DATASET UNSTRUCTURED_GRID" << std::endl;
   file << std::endl;
   
   /* Write the point coordinates: */
   file << "POINTS " << num_points << " double" << std::endl;
   for (int i = 0; i < num_points; i++)
      file << points[i][0] << " " << points[i][1] << " " << points[i][2] << std::endl;
   file << std::endl;
   
   /* Write the cell indices: */
   int num_indices = num_cells;
   for (int i = 0; i < num_cells; i++)
      num_indices += cells.points[i].size();
   file << "CELLS " << num_cells << " " << num_indices << std::endl;
   for (int i = 0; i < num_cells; i++) {
      num_indices = cells.points[i].size();
      file << num_indices;
      for (int j = 0; j < num_indices; j++) {
         file << " " << cells.points[i][j];
         if (j == num_indices-1) file << std::endl;
      }
   }
   file << std::endl;
   
   /* Write the cell types (TODO: write this for cells other than hexahedrons!): */
   file << "CELL_TYPES " << num_cells << std::endl;
   for (int i = 0; i < num_cells; i++) {
      num_indices = cells.points[i].size();
      if (num_indices == 8)
         file << "12" << std::endl;
      else
         PAMPA_CHECK(true, 1, "wrong cell type");
   }
   file << std::endl;
   
   /* Write the cell data: */
   file << "CELL_DATA " << num_cells << std::endl << std::endl;
   
   /* Write the cell materials: */
   file << "SCALARS materials double 1" << std::endl;
   file << "LOOKUP_TABLE default" << std::endl;
   for (int i = 0; i < num_cells; i++)
      file << cells.materials[i]+1 << std::endl;
   file << std::endl;
   
   return 0;
   
};
