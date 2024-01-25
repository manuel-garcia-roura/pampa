#include "Mesh.hxx"

/* Write the mesh to a plain-text file in .vtk format: */
int Mesh::writeVTK(const std::string& filename) const {
   
   /* Open the output file: */
   std::ofstream file(filename);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Set the precision: */
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
   for (int i = 0; i < num_points; i++) {
      file << points[i][0] << " ";
      file << points[i][1] << " ";
      file << points[i][2] << std::endl;
   }
   file << std::endl;
   
   /* Write the cell points: */
   int num_cell_points = num_cells;
   for (int i = 0; i < num_cells; i++)
      num_cell_points += cells.points[i].size();
   file << "CELLS " << num_cells << " " << num_cell_points << std::endl;
   for (int i = 0; i < num_cells; i++) {
      num_cell_points = cells.points[i].size();
      file << num_cell_points;
      for (int j = 0; j < num_cell_points; j++) {
         file << " " << cells.points[i][j];
         if (j == num_cell_points-1) file << std::endl;
      }
   }
   file << std::endl;
   
   /* Write the cell types: */
   file << "CELL_TYPES " << num_cells << std::endl;
   for (int i = 0; i < num_cells; i++) {
      num_cell_points = cells.points[i].size();
      if (num_cell_points == 2)
         file << "3" << std::endl;
      else if (num_cell_points == 4)
         file << "9" << std::endl;
      else if (num_cell_points == 8)
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
   
}

/* Write all the mesh data to a plain-text file: */
int Mesh::writeData(const std::string& filename) const {
   
   /* Open the output file: */
   std::ofstream file(filename);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Set the precision: */
   file << std::fixed;
   file << std::setprecision(3);
   
   /* Write the point coordinates: */
   file << "POINTS " << num_points << " double" << std::endl;
   for (int i = 0; i < num_points; i++) {
      file << points[i][0] << " ";
      file << points[i][1] << " ";
      file << points[i][2] << std::endl;
   }
   file << std::endl;
   
   /* Write the cell points: */
   int num_cell_points = num_cells;
   for (int i = 0; i < num_cells; i++)
      num_cell_points += cells.points[i].size();
   file << "CELLS " << num_cells << " " << num_cell_points << std::endl;
   for (int i = 0; i < num_cells; i++) {
      num_cell_points = cells.points[i].size();
      file << num_cell_points;
      for (int j = 0; j < num_cell_points; j++) {
         file << " " << cells.points[i][j];
         if (j == num_cell_points-1) file << std::endl;
      }
   }
   file << std::endl;
   
   /* Write the cell volumes: */
   file << "CELL_VOLUMES " << num_cells << std::endl;
   for (int i = 0; i < num_cells; i++)
      file << cells.volumes[i] << std::endl;
   file << std::endl;
   
   /* Write the cell centroids: */
   file << "CELL_CENTROIDS " << num_cells << std::endl;
   for (int i = 0; i < num_cells; i++) {
      file << cells.centroids[i][0] << " ";
      file << cells.centroids[i][1] << " ";
      file << cells.centroids[i][2] << std::endl;
   }
   file << std::endl;
   
   /* Write the cell materials: */
   file << "CELL_MATERIALS " << num_cells << std::endl;
   for (int i = 0; i < num_cells; i++)
      file << cells.materials[i] << std::endl;
   file << std::endl;
   
   /* Write the face points: */
   int num_face_points = 0;
   for (int i = 0; i < num_cells; i++)
      for (int f = 0; f < faces.points[i].size(); f++)
         num_face_points += 3 + faces.points[i][f].size();
   file << "FACES " << num_cells << " " << num_face_points << std::endl;
   for (int i = 0; i < num_cells; i++) {
      for (int f = 0; f < faces.points[i].size(); f++) {
         num_face_points = faces.points[i][f].size();
         file << i << " " << f << " " << num_face_points;
         for (int j = 0; j < num_face_points; j++) {
            file << " " << faces.points[i][f][j];
            if (j == num_face_points-1) file << std::endl;
         }
      }
   }
   file << std::endl;
   
   /* Write the face areas: */
   int num_faces = 0;
   for (int i = 0; i < num_cells; i++)
      num_faces += faces.areas[i].size();
   file << "FACE_AREAS " << num_cells << " " << num_faces << std::endl;
   for (int i = 0; i < num_cells; i++)
      for (int f = 0; f < faces.areas[i].size(); f++)
         file << i << " " << f << " " << faces.areas[i][f] << std::endl;
   file << std::endl;
   
   /* Write the face centroids: */
   file << "FACE_CENTROIDS " << num_cells << " " << num_faces << std::endl;
   for (int i = 0; i < num_cells; i++) {
      for (int f = 0; f < faces.centroids[i].size(); f++) {
         file << i << " " << f << " ";
         file << faces.centroids[i][f][0] << " ";
         file << faces.centroids[i][f][1] << " ";
         file << faces.centroids[i][f][2] << std::endl;
      }
   }
   file << std::endl;
   
   /* Write the face normals: */
   file << "FACE_NORMALS " << num_cells << " " << num_faces << std::endl;
   for (int i = 0; i < num_cells; i++) {
      for (int f = 0; f < faces.normals[i].size(); f++) {
         file << i << " " << f << " ";
         file << faces.normals[i][f][0] << " ";
         file << faces.normals[i][f][1] << " ";
         file << faces.normals[i][f][2] << std::endl;
      }
   }
   file << std::endl;
   
   /* Write the face neighbours: */
   file << "FACE_NEIGHBOURS " << num_cells << " " << num_faces << std::endl;
   for (int i = 0; i < num_cells; i++) {
      for (int f = 0; f < faces.neighbours[i].size(); f++) {
         int i2 = faces.neighbours[i][f];
         if (i2 < 0)
            file << i << " " << f << " " << bcs[-i2].type << std::endl;
         else
            file << i << " " << f << " " << i2 << std::endl;
      }
   }
   file << std::endl;
   
   return 0;
   
}
