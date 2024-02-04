#include "UnstructuredExtrudedMesh.hxx"

/* Read the mesh from a plain-text input file: */
int UnstructuredExtrudedMesh::read(const std::string& filename) {
   
   /* Open the input file: */
   std::ifstream file(filename);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty()) break;
      
      /* Get the next keyword: */
      if (line[0] == "points") {
         
         /* Get the point coordinates: */
         num_xy_points = std::stoi(line[1]);
         PAMPA_CALL(utils::read(xy_points, num_xy_points, 2, file), 
            "wrong point data in " + filename);
         
      }
      else if (line[0] == "cells") {
         
         /* Get the cell indices: */
         num_xy_cells = std::stoi(line[1]);
         PAMPA_CALL(utils::read(xy_cells, num_xy_cells, file), "wrong cell data in " + filename);
         
      }
      else if (line[0] == "dz") {
         
         /* Get the dz values: */
         nz = std::stoi(line[1]);
         if (nz > 0) {
            PAMPA_CALL(utils::read(dz, nz, file), "wrong dz data in " + filename);
         }
         else {
            PAMPA_CALL(utils::read(dz, 1, file), "wrong dz data in " + filename);
            nz = -nz;
            dz = Array1D<double>(nz, dz(0));
         }
         
      }
      else if (line[0] == "boundary") {
         
         /* Get the cell indices: */
         Array1D<int> xy_boundary;
         int num_xy_boundary_points = std::stoi(line[1]);
         PAMPA_CALL(utils::read(xy_boundary, num_xy_boundary_points, file), 
            "wrong boundary data in " + filename);
         xy_boundaries.push_back(xy_boundary);
         num_xy_boundaries++;
         
      }
      else if (line[0] == "bc") {
         
         /* Get the boundary conditions (1-based indexed): */
         /* Note: [1, nxy] = xy-plane, nxy+1 = -z, nxy+2 = +z (nxy = num_xy_boundaries). */
         bcs.resize(1+num_xy_boundaries+2);
         int i = 1;
         std::string dir = line[i++];
         if (dir == "z") {
            int l = num_xy_boundaries+1;
            bcs[l].type = static_cast<BC::Type>(std::stoi(line[i++])-1);
            if (bcs[l].type == BC::ROBIN) bcs[l].a = std::stod(line[i++]);
            l++;
            bcs[l].type = static_cast<BC::Type>(std::stoi(line[i++])-1);
            if (bcs[l].type == BC::ROBIN) bcs[l].a = std::stod(line[i++]);
         }
         else {
            int l = std::stoi(dir);
            PAMPA_CHECK(l > num_xy_boundaries, 1, "wrong boundary condition in " + filename);
            bcs[l].type = static_cast<BC::Type>(std::stoi(line[i++])-1);
            if (bcs[l].type == BC::ROBIN) bcs[l].a = std::stod(line[i++]);
         }
         
      }
      else if (line[0] == "materials") {
         
         /* Get the material distribution: */
         int num_materials = std::stoi(line[1]);
         int num_cells = num_xy_cells * std::max(nz, 1);
         PAMPA_CHECK(num_materials != num_cells, 1, "wrong number of materials in " + filename);
         PAMPA_CALL(utils::read(cells.materials, num_materials, file), 
            "wrong material data in " + filename);
         
         /* Switch the material index from 1-based to 0-based: */
         for (int i = 0; i < num_materials; i++)
            cells.materials(i)--;
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 1, "unrecognized keyword '" + line[0] + "' in " + filename);
         
      }
      
   }
   
   return 0;
   
}

/* Build the mesh: */
int UnstructuredExtrudedMesh::build() {
   
   /* Build the mesh points: */
   std::vector<double> z(nz+1);
   for (int k = 0; k < nz; k++)
      z[k+1] = z[k] + dz(k);
   num_points = num_xy_points * (nz+1);
   int ip = 0;
   points = Array2D<double>(num_points, 3);
   for (int k = 0; k < nz+1; k++) {
      for (int i = 0; i < num_xy_points; i++) {
         points(ip, 0) = xy_points(i, 0);
         points(ip, 1) = xy_points(i, 1);
         points(ip, 2) = z[k];
         ip++;
      }
   }
   
   /* Build the mesh cells: */
   /* Note: the cell points are ordered according to the gmsh convention. */
   num_cells = num_xy_cells * std::max(nz, 1);
   int ic = 0;
   int num_xy_cell_points_max = 0;
   std::vector<int> num_cell_points(num_cells);
   for (int k = 0; k < std::max(nz, 1); k++)
      for (int i = 0; i < num_xy_cells; i++)
         num_cell_points[ic++] = (nz > 0) ? 2 * xy_cells[i].size() : xy_cells[i].size();
   ic = 0;
   cells.points = Vector2D<int>(num_cells, num_cell_points);
   cells.volumes = Array1D<double>(num_cells);
   cells.centroids = Array2D<double>(num_cells, 3);
   for (int k = 0; k < std::max(nz, 1); k++) {
      for (int i = 0; i < num_xy_cells; i++) {
         
         /* Get the cell points: */
         int n = xy_cells[i].size();
         for (int l = 0; l < n; l++)
            cells.points(ic, l) = xy_cells[i][l];
         if (nz > 0)
            for (int l = 0; l < n; l++)
               cells.points(ic, n+l) = xy_cells[i][l] + (k+1)*num_xy_points;
         num_xy_cell_points_max = std::max(int(xy_cells[i].size()), num_xy_cell_points_max);
         
         /* Get the cell volume: */
         double a = math::get_area(xy_points, xy_cells[i]);
         double v = (nz > 0) ? a * dz(k) : a;
         cells.volumes(ic) = v;
         
         /* Get the cell centroid: */
         std::vector<double> p0 = math::get_centroid(xy_points, xy_cells[i], a);
         cells.centroids(ic, 0) = p0[0];
         cells.centroids(ic, 1) = p0[1];
         cells.centroids(ic, 2) = z[k] + 0.5*dz(k);
         
         /* Move to the next cell: */
         ic++;
         
      }
   }
   
   /* Get the cells for each point in the xy-plane: */
   std::vector<std::vector<int>> xy_points_to_cells(num_xy_points);
   for (int i = 0; i < num_xy_cells; i++)
      for (int j = 0; j < xy_cells[i].size(); j++)
         xy_points_to_cells[xy_cells[i][j]].push_back(i);
   
   /* Set the boundary conditions (1-based indexed): */
   for (int i = 0; i < num_xy_boundaries; i++)
      for (int j = 0; j < xy_boundaries[i].size(); j++)
         xy_points_to_cells[xy_boundaries[i](j)].push_back(-i-1);
   
   /* Get the neighboring cell for each cell face in the xy-plane: */
   std::vector<std::vector<int>> xy_neighbours(num_xy_cells);
   bool found;
   for (int i = 0; i < num_xy_cells; i++) {
      
      /* Get the neighbouring cell for each face of this cell: */
      int num_xy_cell_points = xy_cells[i].size();
      xy_neighbours[i].resize(num_xy_cell_points);
      for (int f = 0; f < num_xy_cell_points; f++) {
         
         /* Get the cells for both points in this face: */
         const std::vector<int>& cls1 = xy_points_to_cells[xy_cells[i][f]];
         const std::vector<int>& cls2 = xy_points_to_cells[xy_cells[i][(f+1)%num_xy_cell_points]];
         
         /* Find the cell connected to both points in this face: */
         found = false;
         for (int l1 = 0; l1 < cls1.size(); l1++) {
            for (int l2 = 0; l2 < cls2.size(); l2++) {
               int i1 = cls1[l1];
               int i2 = cls2[l2];
               if (i1 == i2 && i1 != i) {
                  PAMPA_CHECK(xy_neighbours[i][f] > 0, 1, "wrong mesh connectivity");
                  xy_neighbours[i][f] = i1;
                  found = true;
               }
            }
         }
         PAMPA_CHECK(!found, 1, "wrong mesh connectivity");
         
      }
   }
   
   /* Build the mesh faces: */
   /* Note: the face points are ordered counterclockwise so that the normal points outward.*/
   int num_faces_max = (nz > 0) ? num_xy_cell_points_max + 2 : num_xy_cell_points_max;
   faces.num_faces = Array1D<int>(num_cells);
   faces.areas = Array2D<double>(num_cells, num_faces_max);
   faces.centroids = Array3D<double>(num_cells, num_faces_max, 3);
   faces.normals = Array3D<double>(num_cells, num_faces_max, 3);
   faces.neighbours = Array2D<int>(num_cells, num_faces_max);
   ic = 0;
   for (int k = 0; k < std::max(nz, 1); k++) {
      for (int i = 0; i < num_xy_cells; i++) {
         
         /* Initialize the face data for this cell: */
         int nxy = xy_cells[i].size();
         faces.num_faces(ic) = (nz > 0) ? nxy + 2 : nxy;
         
         /* xy-plane faces: */
         for (int f = 0; f < nxy; f++) {
		      int f2 = (f+1)%nxy;
            faces.areas(ic, f) = math::get_distance(xy_points, xy_cells[i][f], xy_cells[i][f2], 2);
            if (nz > 0) faces.areas(ic, f) *= dz(k);
            std::vector<double> p0 = math::get_midpoint(xy_points, xy_cells[i][f], xy_cells[i][f2], 2);
            faces.centroids(ic, f, 0) = p0[0];
            faces.centroids(ic, f, 1) = p0[1];
            faces.centroids(ic, f, 2) = z[k]+0.5*dz(k);
            std::vector<double> n = math::get_normal(xy_points, xy_cells[i][f], xy_cells[i][f2]);
            faces.normals(ic, f, 0) = n[0];
            faces.normals(ic, f, 1) = n[1];
            faces.normals(ic, f, 2) = 0.0;
            faces.neighbours(ic, f) = xy_neighbours[i][f];
            if (faces.neighbours(ic, f) >= 0) faces.neighbours(ic, f) += k * num_xy_cells;
         }
         
         /* -z face: */
         if (nz > 0) {
            faces.areas(ic, nxy) = math::get_area(xy_points, xy_cells[i]);
            std::vector<double> p0 = math::get_centroid(xy_points, xy_cells[i], faces.areas(ic, nxy));
            faces.centroids(ic, nxy, 0) = p0[0];
            faces.centroids(ic, nxy, 1) = p0[1];
            faces.centroids(ic, nxy, 2) = z[k];
            faces.normals(ic, nxy, 0) = 0.0;
            faces.normals(ic, nxy, 1) = 0.0;
            faces.normals(ic, nxy, 2) = -1.0;
            faces.neighbours(ic, nxy) = (k == 0) ? -num_xy_boundaries - 1 : ic - num_xy_cells;
         }
         
         /* +z face: */
         if (nz > 0) {
            faces.areas(ic, nxy+1) = math::get_area(xy_points, xy_cells[i]);
            std::vector<double> p0 = math::get_centroid(xy_points, xy_cells[i], faces.areas(ic, nxy+1));
            faces.centroids(ic, nxy+1, 0) = p0[0];
            faces.centroids(ic, nxy+1, 1) = p0[1];
            faces.centroids(ic, nxy+1, 2) = z[k]+dz(k);
            faces.normals(ic, nxy+1, 0) = 0.0;
            faces.normals(ic, nxy+1, 1) = 0.0;
            faces.normals(ic, nxy+1, 2) = 1.0;
            faces.neighbours(ic, nxy+1) = (k == nz-1) ? -num_xy_boundaries - 2 : ic + num_xy_cells;
         }
         
         /* Move to the next cell: */
         ic++;
         
      }
   }
   
   /* Write all the mesh data for debugging: */
#ifdef DEBUG
   writeData("mesh_data.pmp");
#endif
   
   return 0;
   
}
