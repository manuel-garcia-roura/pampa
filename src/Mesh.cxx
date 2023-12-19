#include "Mesh.hxx"

/* The Mesh constructor: */
Mesh::Mesh() {};

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
   for (int i = 0; i < num_cells; i++)
      file << "12" << std::endl;
   file << std::endl;
   
   /* Write the cell materials: */
   file << "CELL_DATA " << num_cells << std::endl;
   file << "SCALARS materials double 1" << std::endl;
   file << "LOOKUP_TABLE default" << std::endl;
   for (int i = 0; i < num_cells; i++)
      file << cells.materials[i] << std::endl;
   
   return 0;
   
};
