#include "PartitionedMesh.hxx"

/* Read the mesh from a plain-text input file: */
int PartitionedMesh::read(const std::string& filename) {
   
   /* Open the input file: */
   std::string path = mpi::get_path(filename);
   std::ifstream file(path);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + path);
   
   /* Read the file line by line: */
   int num_cell_faces = 0;
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty()) break;
      
      /* Get the next keyword: */
      if (line[0] == "points") {
         
         /* Get the point coordinates: */
         PAMPA_CALL(utils::read(num_points, 1, INT_MAX, line[1]), "wrong number of points");
         PAMPA_CALL(utils::read(points, num_points, 3, file), "wrong point data");
         
      }
      else if (line[0] == "cells") {
         
         /* Get the number of cells: */
         PAMPA_CALL(utils::read(num_cells, 1, INT_MAX, line[1]), "wrong number of cells");
         PAMPA_CALL(utils::read(num_ghost_cells, 1, INT_MAX, line[2]), 
            "wrong number of ghost cells");
         PAMPA_CALL(utils::read(num_cells_global, 1, INT_MAX, line[3]), 
            "wrong global number of cells");
         
      }
      else if (line[0] == "cell-points") {
         
         /* Get the cell points: */
         int num_rows, num_cell_points;
         PAMPA_CALL(utils::read(num_rows, num_cells, num_cells, line[1]), "wrong number of cells");
         PAMPA_CALL(utils::read(num_cell_points, 1, INT_MAX, line[2]), 
            "wrong number of cell points");
         PAMPA_CALL(utils::read(cells.points, num_cells, num_cell_points, file), 
            "wrong cell-point data");
         
      }
      else if (line[0] == "cell-volumes") {
         
         /* Get the cell volumes: */
         int num_rows;
         PAMPA_CALL(utils::read(num_rows, num_cells, num_cells, line[1]), 
            "wrong number of cell volumes");
         PAMPA_CALL(utils::read(cells.volumes, num_cells, file), "wrong cell-volume data");
         
      }
      else if (line[0] == "cell-centroids") {
         
         /* Get the cell centroids: */
         int num_rows, num_cells_total = num_cells + num_ghost_cells;
         PAMPA_CALL(utils::read(num_rows, num_cells_total, num_cells_total, line[1]), 
            "wrong number of cell centroids");
         PAMPA_CALL(utils::read(cells.centroids, num_cells_total, 3, file), 
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
      else if (line[0] == "cell-materials") {
         
         /* Get the cell materials (1-based indexed): */
         int num_rows, num_cells_total = num_cells + num_ghost_cells;
         PAMPA_CALL(utils::read(num_rows, num_cells_total, num_cells_total, line[1]), 
            "wrong number of cell materials");
         PAMPA_CALL(utils::read(cells.materials, num_cells_total, file), 
            "wrong cell-material data");
         
         /* Switch the material index from 1-based to 0-based: */
         for (int i = 0; i < num_cells+num_ghost_cells; i++)
            cells.materials(i)--;
         
      }
      else if (line[0] == "cell-indices") {
         
         /* Get the cell indices: */
         int num_rows, num_cells_total = num_cells + num_ghost_cells;
         PAMPA_CALL(utils::read(num_rows, num_cells_total, num_cells_total, line[1]), 
            "wrong number of cell indices");
         PAMPA_CALL(utils::read(cells.indices, num_cells_total, file), "wrong cell-index data");
         
      }
      else if (line[0] == "faces") {
         
         /* Get the number of cell faces: */
         int num_rows;
         PAMPA_CALL(utils::read(num_rows, num_cells, num_cells, line[1]), 
            "wrong number of cell faces");
         PAMPA_CALL(utils::read(faces.num_faces, num_cells, file), "wrong cell-face data");
         
         /* Get the maximum and the total number of cell faces: */
         num_cell_faces = 0;
         for (int i = 0; i < num_cells; i++) {
            int n = faces.num_faces(i);
            num_cell_faces += n;
            num_faces_max = std::max(n, num_faces_max);
         }
         
      }
      else if (line[0] == "face-areas") {
         
         /* Get the face areas: */
         int num_rows, num_elements;
         PAMPA_CALL(utils::read(num_rows, num_cells, num_cells, line[1]), 
            "wrong number of cell faces");
         PAMPA_CALL(utils::read(num_elements, num_cell_faces, num_cell_faces, line[2]), 
            "wrong number of face areas");
         PAMPA_CALL(utils::read(faces.areas, num_cells, num_cell_faces, file), 
            "wrong face-area data");
         
      }
      else if (line[0] == "face-centroids") {
         
         /* Get the face centroids: */
         int num_rows, num_elements;
         PAMPA_CALL(utils::read(num_rows, num_cells, num_cells, line[1]), 
            "wrong number of cell faces");
         PAMPA_CALL(utils::read(num_elements, 3*num_cell_faces, 3*num_cell_faces, line[2]), 
            "wrong number of face-centroid coordinates");
         PAMPA_CALL(utils::read(faces.centroids, num_cells, faces.num_faces, 3, file), 
            "wrong face-centroid data");
         
      }
      else if (line[0] == "face-normals") {
         
         /* Get the face normals: */
         int num_rows, num_elements;
         PAMPA_CALL(utils::read(num_rows, num_cells, num_cells, line[1]), 
            "wrong number of cell faces");
         PAMPA_CALL(utils::read(num_elements, 3*num_cell_faces, 3*num_cell_faces, line[2]), 
            "wrong number of face-normal coordinates");
         PAMPA_CALL(utils::read(faces.normals, num_cells, faces.num_faces, 3, file), 
            "wrong face-normal data");
         
      }
      else if (line[0] == "face-neighbours") {
         
         /* Get the face neighbours: */
         int num_rows, num_elements;
         PAMPA_CALL(utils::read(num_rows, num_cells, num_cells, line[1]), 
            "wrong number of cell faces");
         PAMPA_CALL(utils::read(num_elements, num_cell_faces, num_cell_faces, line[2]), 
            "wrong number of face neighbours");
         PAMPA_CALL(utils::read(faces.neighbours, num_cells, num_cell_faces, file), 
            "wrong face-neighbour data");
         
      }
      else if (line[0] == "bc") {
         
         /* Initialize the boundary-condition array if not done yet: */
         if (bcs.empty()) bcs.resize(1);
         
         /* Get the boundary condition (1-based indexed): */
         int l, i = 1;
         PAMPA_CALL(utils::read(l, bcs.size(), bcs.size(), line[i++]), 
            "wrong boundary condition index");
         BoundaryCondition bc;
         PAMPA_CALL(utils::read(bc, line, i), "wrong boundary condition");
         bcs.pushBack(bc);
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 1, "unrecognized keyword '" + line[0] + "'");
         
      }
      
   }
   
   return 0;
   
}
