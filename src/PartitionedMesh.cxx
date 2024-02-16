#include "PartitionedMesh.hxx"

/* Read the mesh from a plain-text input file: */
int PartitionedMesh::read(const std::string& filename) {
   
   /* Read from a rank directory in parallel runs: */
   std::string path;
   if (mpi::size > 1) {
      std::string dir = std::to_string(mpi::rank);
      path = dir + "/" + filename;
   }
   else
      path = filename;
   
   /* Open the input file: */
   std::ifstream file(path);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + path);
   
   /* Read the file line by line: */
   int num_cell_faces;
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty()) break;
      
      /* Get the next keyword: */
      if (line[0] == "points") {
         
         /* Get the point coordinates: */
         num_points = std::stoi(line[1]);
         PAMPA_CALL(utils::read(points, num_points, 3, file), "wrong point data in " + filename);
         
      }
      else if (line[0] == "cells") {
         
         /* Get the number of cells: */
         num_cells = std::stoi(line[1]);
         num_ghost_cells = std::stoi(line[2]);
         num_cells_global = std::stoi(line[3]);
         
      }
      else if (line[0] == "cell-points") {
         
         /* Get the cell points: */
         int num_rows = std::stoi(line[1]);
         int num_cell_points = std::stoi(line[2]);
         PAMPA_CHECK(num_rows != num_cells, 1, "wrong number of cell points in " + filename);
         PAMPA_CALL(utils::read(cells.points, num_cells, num_cell_points, file), 
            "wrong cell-point data in " + filename);
         
      }
      else if (line[0] == "cell-volumes") {
         
         /* Get the cell volumes: */
         int num_volumes = std::stoi(line[1]);
         PAMPA_CHECK(num_volumes != num_cells, 1, "wrong number of cell volumes in " + filename);
         PAMPA_CALL(utils::read(cells.volumes, num_cells, file), 
            "wrong cell-volume data in " + filename);
         
      }
      else if (line[0] == "cell-centroids") {
         
         /* Get the cell centroids: */
         int num_centroids = std::stoi(line[1]);
         PAMPA_CHECK(num_centroids != num_cells+num_ghost_cells, 1, 
            "wrong number of cell centroids in " + filename);
         PAMPA_CALL(utils::read(cells.centroids, num_cells+num_ghost_cells, 3, file), 
            "wrong cell-centroid data in " + filename);
         
      }
      else if (line[0] == "cell-materials") {
         
         /* Get the cell materials (1-based indexed): */
         int num_materials = std::stoi(line[1]);
         PAMPA_CHECK(num_materials != num_cells+num_ghost_cells, 1, 
            "wrong number of cell materials in " + filename);
         PAMPA_CALL(utils::read(cells.materials, num_cells+num_ghost_cells, file), 
            "wrong cell-material data in " + filename);
         
         /* Switch the material index from 1-based to 0-based: */
         for (int i = 0; i < num_cells+num_ghost_cells; i++)
            cells.materials(i)--;
         
      }
      else if (line[0] == "cell-indices") {
         
         /* Get the cell indices: */
         int num_indices = std::stoi(line[1]);
         PAMPA_CHECK(num_indices != num_cells+num_ghost_cells, 1, 
            "wrong number of cell indices in " + filename);
         PAMPA_CALL(utils::read(cells.indices, num_cells+num_ghost_cells, file), 
            "wrong cell-index data in " + filename);
         
      }
      else if (line[0] == "faces") {
         
         /* Get the number of cell faces: */
         int num_faces = std::stoi(line[1]);
         PAMPA_CHECK(num_faces != num_cells, 1, "wrong number of cell faces in " + filename);
         PAMPA_CALL(utils::read(faces.num_faces, num_cells, file), 
            "wrong cell-face data in " + filename);
         
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
         int num_rows = std::stoi(line[1]);
         int num_elements = std::stoi(line[2]);
         PAMPA_CHECK(num_rows != num_cells, 1, "wrong number of face areas in " + filename);
         PAMPA_CHECK(num_elements != num_cell_faces, 1, 
            "wrong number of face areas in " + filename);
         PAMPA_CALL(utils::read(faces.areas, num_cells, num_cell_faces, file), 
            "wrong face-area data in " + filename);
         
      }
      else if (line[0] == "face-centroids") {
         
         /* Get the face centroids: */
         int num_rows = std::stoi(line[1]);
         int num_elements = std::stoi(line[2]);
         PAMPA_CHECK(num_rows != num_cells, 1, "wrong number of face centroids in " + filename);
         PAMPA_CHECK(num_elements != 3*num_cell_faces, 1, 
            "wrong number of face-centroid coordinates in " + filename);
         PAMPA_CALL(utils::read(faces.centroids, num_cells, faces.num_faces, 3, file), 
            "wrong face-centroid data in " + filename);
         
      }
      else if (line[0] == "face-normals") {
         
         /* Get the face normals: */
         int num_rows = std::stoi(line[1]);
         int num_elements = std::stoi(line[2]);
         PAMPA_CHECK(num_rows != num_cells, 1, "wrong number of face normals in " + filename);
         PAMPA_CHECK(num_elements != 3*num_cell_faces, 1, 
            "wrong number of face-normal coordinates in " + filename);
         PAMPA_CALL(utils::read(faces.normals, num_cells, faces.num_faces, 3, file), 
            "wrong face-normal data in " + filename);
         
      }
      else if (line[0] == "face-neighbours") {
         
         /* Get the face neighbours: */
         int num_rows = std::stoi(line[1]);
         int num_elements = std::stoi(line[2]);
         PAMPA_CHECK(num_rows != num_cells, 1, "wrong number of face neighbours in " + filename);
         PAMPA_CHECK(num_elements != num_cell_faces, 1, 
            "wrong number of face neighbours in " + filename);
         PAMPA_CALL(utils::read(faces.neighbours, num_cells, num_cell_faces, file), 
            "wrong face-neighbour data in " + filename);
         
      }
      else if (line[0] == "bc") {
         
         /* Initialize the boundary-condition array if not done yet: */
         if (bcs.empty()) bcs.resize(1);
         
         /* Get the boundary condition (1-based indexed): */
         int i = 1, l = std::stoi(line[i++]);
         PAMPA_CHECK(l != bcs.size(), 1, "wrong boundary condition in " + filename);
         BoundaryCondition bc;
         PAMPA_CALL(utils::read(bc, line, i), "wrong boundary condition in " + filename);
         bcs.pushBack(bc);
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 1, "unrecognized keyword '" + line[0] + "' in " + filename);
         
      }
      
   }
   
   return 0;
   
}
