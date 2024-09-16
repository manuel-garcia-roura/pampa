#include "PartitionedMesh.hxx"

/* Read the mesh from a plain-text input file: */
int PartitionedMesh::read(const std::string& filename) {
   
   /* Open the input file: */
   std::string path = mpi::get_path(filename);
   std::ifstream file(path, std::ios_base::in);
   PAMPA_CHECK(!file.is_open(), "unable to open " + path);
   
   /* Read the file line by line: */
   int num_cell_faces = 0;
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty()) break;
      
      /* Get the next keyword: */
      unsigned int l = 0;
      if (line[l] == "points") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the point coordinates: */
         PAMPA_CHECK(utils::read(num_points, 1, INT_MAX, line[++l]), "wrong number of points");
         PAMPA_CHECK(utils::read(points, num_points, 3, file), "wrong point data");
         
      }
      else if (line[l] == "cells") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 4, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the number of cells: */
         PAMPA_CHECK(utils::read(num_cells, 1, INT_MAX, line[++l]), "wrong number of cells");
         PAMPA_CHECK(utils::read(num_ghost_cells, 1, INT_MAX, line[++l]), 
            "wrong number of ghost cells");
         PAMPA_CHECK(utils::read(num_cells_global, 1, INT_MAX, line[++l]), 
            "wrong global number of cells");
         
      }
      else if (line[l] == "cell-points") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 3, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the cell points: */
         int num_rows, num_cell_points;
         PAMPA_CHECK(utils::read(num_rows, num_cells, num_cells, line[++l]), 
            "wrong number of cells");
         PAMPA_CHECK(utils::read(num_cell_points, 1, INT_MAX, line[++l]), 
            "wrong number of cell points");
         PAMPA_CHECK(utils::read(cells.points, num_cells, num_cell_points, file), 
            "wrong cell-point data");
         
      }
      else if (line[l] == "cell-volumes") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the cell volumes: */
         int num_rows;
         PAMPA_CHECK(utils::read(num_rows, num_cells, num_cells, line[++l]), 
            "wrong number of cell volumes");
         PAMPA_CHECK(utils::read(cells.volumes, num_cells, file), "wrong cell-volume data");
         
      }
      else if (line[l] == "cell-centroids") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the cell centroids: */
         int num_rows, num_cells_total = num_cells + num_ghost_cells;
         PAMPA_CHECK(utils::read(num_rows, num_cells_total, num_cells_total, line[++l]), 
            "wrong number of cell centroids");
         PAMPA_CHECK(utils::read(cells.centroids, num_cells_total, 3, file), 
            "wrong cell-centroid data");
         
         /* Get the number of dimensions: */
         num_dims = 1;
         for (int i = 1; i < num_cells; i++) {
            if (cells.centroids(i, 1) != cells.centroids(0, 1)) {
               num_dims++;
               break;
            }
         }
         for (int i = 1; i < num_cells; i++) {
            if (cells.centroids(i, 2) != cells.centroids(0, 2)) {
               num_dims++;
               break;
            }
         }
         
      }
      else if (line[l] == "cell-materials") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the cell materials (1-based indexed): */
         int num_rows, num_cells_total = num_cells + num_ghost_cells;
         PAMPA_CHECK(utils::read(num_rows, num_cells_total, num_cells_total, line[++l]), 
            "wrong number of cell materials");
         PAMPA_CHECK(utils::read(cells.materials, num_cells_total, file), 
            "wrong cell-material data");
         
         /* Switch the material index from 1-based to 0-based: */
         for (int i = 0; i < num_cells+num_ghost_cells; i++)
            cells.materials(i)--;
         
      }
      else if (line[l] == "cell-nodal-indices") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the cell nodal indices: */
         int num_rows;
         PAMPA_CHECK(utils::read(num_rows, num_cells, num_cells, line[++l]), 
            "wrong number of cell nodal indices");
         PAMPA_CHECK(utils::read(cells.nodal_indices, num_cells, file), "wrong nodal-index data");
         
      }
      else if (line[l] == "cell-global-indices") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the cell global indices: */
         int num_rows, num_cells_total = num_cells + num_ghost_cells;
         PAMPA_CHECK(utils::read(num_rows, num_cells_total, num_cells_total, line[++l]), 
            "wrong number of cell global indices");
         PAMPA_CHECK(utils::read(cells.global_indices, num_cells_total, file), 
            "wrong global-index data");
         
      }
      else if (line[l] == "faces") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the number of cell faces: */
         int num_rows;
         PAMPA_CHECK(utils::read(num_rows, num_cells, num_cells, line[++l]), 
            "wrong number of cell faces");
         PAMPA_CHECK(utils::read(faces.num_faces, num_cells, file), "wrong cell-face data");
         
         /* Get the maximum and the total number of cell faces: */
         num_cell_faces = 0;
         for (int i = 0; i < num_cells; i++) {
            int n = faces.num_faces(i);
            num_cell_faces += n;
            num_faces_max = std::max(n, num_faces_max);
         }
         
      }
      else if (line[l] == "face-areas") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 3, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the face areas: */
         int num_rows, num_elements;
         PAMPA_CHECK(utils::read(num_rows, num_cells, num_cells, line[++l]), 
            "wrong number of cell faces");
         PAMPA_CHECK(utils::read(num_elements, num_cell_faces, num_cell_faces, line[++l]), 
            "wrong number of face areas");
         PAMPA_CHECK(utils::read(faces.areas, num_cells, num_cell_faces, file), 
            "wrong face-area data");
         
      }
      else if (line[l] == "face-centroids") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 3, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the face centroids: */
         int num_rows, num_elements;
         PAMPA_CHECK(utils::read(num_rows, num_cells, num_cells, line[++l]), 
            "wrong number of cell faces");
         PAMPA_CHECK(utils::read(num_elements, 3*num_cell_faces, 3*num_cell_faces, line[++l]), 
            "wrong number of face-centroid coordinates");
         PAMPA_CHECK(utils::read(faces.centroids, num_cells, faces.num_faces, 3, file), 
            "wrong face-centroid data");
         
      }
      else if (line[l] == "face-normals") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 3, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the face normals: */
         int num_rows, num_elements;
         PAMPA_CHECK(utils::read(num_rows, num_cells, num_cells, line[++l]), 
            "wrong number of cell faces");
         PAMPA_CHECK(utils::read(num_elements, 3*num_cell_faces, 3*num_cell_faces, line[++l]), 
            "wrong number of face-normal coordinates");
         PAMPA_CHECK(utils::read(faces.normals, num_cells, faces.num_faces, 3, file), 
            "wrong face-normal data");
         
      }
      else if (line[l] == "face-neighbours") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 3, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the face neighbours: */
         int num_rows, num_elements;
         PAMPA_CHECK(utils::read(num_rows, num_cells, num_cells, line[++l]), 
            "wrong number of cell faces");
         PAMPA_CHECK(utils::read(num_elements, num_cell_faces, num_cell_faces, line[++l]), 
            "wrong number of face neighbours");
         PAMPA_CHECK(utils::read(faces.neighbours, num_cells, num_cell_faces, file), 
            "wrong face-neighbour data");
         
      }
      else if (line[l] == "boundary") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the boundary name: */
         std::string name = line[++l];
         boundaries.pushBack(name);
         
      }
      else if (line[l] == "bc") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() < 3, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Initialize the boundary-condition array, if not done yet: */
         if (bcs.empty()) bcs.resize(1+boundaries.size());
         
         /* Get the boundary name and index: */
         std::string name = line[++l];
         int i = boundaries.find(name);
         PAMPA_CHECK(i < 0, "wrong boundary name");
         
         /* Get the boundary condition (1-based indexed): */
         PAMPA_CHECK(utils::read(bcs(i+1), line, ++l, file), "wrong boundary condition");
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   return 0;
   
}
