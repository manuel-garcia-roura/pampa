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
            dz.resize(nz, dz(0));
         }
         
      }
      else if (line[0] == "boundary") {
         
         /* Get the cell indices: */
         Array1D<int> xy_boundary;
         int num_xy_boundary_points = std::stoi(line[1]);
         PAMPA_CALL(utils::read(xy_boundary, num_xy_boundary_points, file), 
            "wrong boundary data in " + filename);
         xy_boundaries.pushBack(xy_boundary);
         num_xy_boundaries++;
         
      }
      else if (line[0] == "bc") {
         
         /* Get the boundary conditions (1-based indexed): */
         /* Note: [1, n] = xy-plane, n+1 = -z, n+2 = +z (n = num_xy_boundaries). */
         if (bcs.empty()) bcs.resize(1+num_xy_boundaries+2);
         int i = 1;
         std::string dir = line[i++];
         if (dir == "z") {
            int l = num_xy_boundaries+1;
            PAMPA_CALL(utils::read(bcs(l++), line, i), "wrong boundary condition in " + filename);
            PAMPA_CALL(utils::read(bcs(l), line, i), "wrong boundary condition in " + filename);
         }
         else {
            int l = std::stoi(dir);
            PAMPA_CHECK(l > num_xy_boundaries, 1, "wrong boundary condition in " + filename);
            PAMPA_CALL(utils::read(bcs(l), line, i), "wrong boundary condition in " + filename);
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
   
   /* Build the axis coordinates: */
   Array1D<double> z(nz+1);
   for (int k = 0; k < nz; k++)
      z(k+1) = z(k) + dz(k);
   
   /* Build the mesh points: */
   num_points = num_xy_points * (nz+1);
   points.resize(num_points, 3);
   for (int ip = 0, k = 0; k < nz+1; k++) {
      for (int i = 0; i < num_xy_points; i++) {
         points(ip, 0) = xy_points(i, 0);
         points(ip, 1) = xy_points(i, 1);
         points(ip, 2) = z(k);
         ip++;
      }
   }
   
   /* Get the number of points for each cell in the xy-plane: */
   Array1D<int> num_xy_cell_points(num_xy_cells);
   int num_xy_cell_points_max = 0;
   for (int i = 0; i < num_xy_cells; i++) {
      int n = xy_cells.size(i);
      num_xy_cell_points(i) = n;
      num_xy_cell_points_max = std::max(n, num_xy_cell_points_max);
   }
   
   /* Get the number of points for each cell: */
   num_cells = num_xy_cells * std::max(nz, 1);
   Array1D<int> num_cell_points(num_cells);
   for (int ic = 0, k = 0; k < std::max(nz, 1); k++) {
      for (int i = 0; i < num_xy_cells; i++) {
         int n = xy_cells.size(i);
         num_cell_points(ic++) = (nz > 0) ? 2*n : n;
      }
   }
   
   /* Build the mesh cells: */
   /* Note: the cell points are ordered according to the gmsh convention. */
   cells.points.resize(num_cells, num_cell_points);
   cells.volumes.resize(num_cells);
   cells.centroids.resize(num_cells, 3);
   cells.indices.resize(num_cells);
   for (int ic = 0, k = 0; k < std::max(nz, 1); k++) {
      for (int i = 0; i < num_xy_cells; i++) {
         
         /* Get the cell points: */
         int n = xy_cells.size(i);
         for (int l = 0; l < n; l++)
            cells.points(ic, l) = xy_cells(i, l);
         if (nz > 0)
            for (int l = 0; l < n; l++)
               cells.points(ic, n+l) = xy_cells(i, l) + (k+1)*num_xy_points;
         
         /* Get the cell volume: */
         double a = math::area(xy_points, xy_cells(i), n);
         double v = (nz > 0) ? a*dz(k) : a;
         cells.volumes(ic) = v;
         
         /* Get the cell centroid: */
         math::centroid(cells.centroids(ic), xy_points, xy_cells(i), n, a);
         cells.centroids(ic, 2) = z(k) + 0.5*dz(k);
         
         /* Get the cell index in the global mesh: */
         cells.indices(ic) = ic;
         
         ic++;
         
      }
   }
   
   /* Get the number of cells and boundary conditions for each point in the xy-plane: */
   Array1D<int> num_xy_point_cells(num_xy_points);
   for (int i = 0; i < num_xy_cells; i++)
      for (int j = 0; j < xy_cells.size(i); j++)
         num_xy_point_cells(xy_cells(i, j))++;
   for (int i = 0; i < num_xy_boundaries; i++)
      for (int j = 0; j < xy_boundaries.size(i); j++)
         num_xy_point_cells(xy_boundaries(i, j))++;
   
   /* Get the cells and boundary conditions (1-based indexed) for each point in the xy-plane: */
   Vector2D<int> xy_point_cells(num_xy_points, num_xy_point_cells);
   Array1D<int> ipc(num_xy_points);
   for (int i = 0; i < num_xy_cells; i++)
      for (int j = 0; j < xy_cells.size(i); j++)
         xy_point_cells(xy_cells(i, j), ipc(xy_cells(i, j))++) = i;
   for (int i = 0; i < num_xy_boundaries; i++)
      for (int j = 0; j < xy_boundaries.size(i); j++)
         xy_point_cells(xy_boundaries(i, j), ipc(xy_boundaries(i, j))++) = -i-1;
   
   /* Get the neighboring cell for each face in the xy-plane: */
   Vector2D<int> xy_neighbours(num_xy_cells, num_xy_cell_points);
   for (int i = 0; i < num_xy_cells; i++) {
      for (int f = 0; f < num_xy_cell_points(i); f++) {
         
         /* Get the cell indices for both points in this face: */
         int ic1 = xy_cells(i, f);
         int ic2 = xy_cells(i, (f+1)%num_xy_cell_points(i));
         
         /* Find the cell connected to both points in this face: */
         bool found = false;
         for (int l1 = 0; l1 < xy_point_cells.size(ic1); l1++) {
            for (int l2 = 0; l2 < xy_point_cells.size(ic2); l2++) {
               int i1 = xy_point_cells(ic1, l1);
               int i2 = xy_point_cells(ic2, l2);
               if (i1 == i2 && i1 != i) {
                  PAMPA_CHECK(xy_neighbours(i, f) > 0, 1, "wrong mesh connectivity");
                  xy_neighbours(i, f) = i1;
                  found = true;
               }
            }
         }
         PAMPA_CHECK(!found, 1, "wrong mesh connectivity");
         
      }
   }
   
   /* Build the mesh faces: */
   /* Note: the face points are ordered counterclockwise so that the normal points outward.*/
   num_faces_max = (nz > 0) ? num_xy_cell_points_max+2 : num_xy_cell_points_max;
   faces.num_faces.resize(num_cells);
   faces.areas.resize(num_cells, num_faces_max);
   faces.centroids.resize(num_cells, num_faces_max, 3);
   faces.normals.resize(num_cells, num_faces_max, 3);
   faces.neighbours.resize(num_cells, num_faces_max);
   for (int ic = 0, k = 0; k < std::max(nz, 1); k++) {
      for (int i = 0; i < num_xy_cells; i++) {
         
         /* Get the number of faces for this cell: */
         int n = xy_cells.size(i);
         faces.num_faces(ic) = (nz > 0) ? n+2 : n;
         
         /* xy-plane faces: */
         int f;
         for (f = 0; f < n; f++) {
		      int f2 = (f+1) % n;
            faces.areas(ic, f) = math::distance(xy_points, xy_cells(i, f), xy_cells(i, f2), 2);
            if (nz > 0) faces.areas(ic, f) *= dz(k);
            math::midpoint(faces.centroids(ic, f), xy_points, xy_cells(i, f), xy_cells(i, f2), 2);
            faces.centroids(ic, f, 2) = z(k)+0.5*dz(k);
            math::normal(faces.normals(ic, f), xy_points, xy_cells(i, f), xy_cells(i, f2));
            faces.normals(ic, f, 2) = 0.0;
            faces.neighbours(ic, f) = xy_neighbours(i, f);
            if (faces.neighbours(ic, f) >= 0) faces.neighbours(ic, f) += k * num_xy_cells;
         }
         
         /* -z face: */
         if (nz > 0) {
            faces.areas(ic, f) = math::area(xy_points, xy_cells(i), n);
            math::centroid(faces.centroids(ic, f), xy_points, xy_cells(i), n, faces.areas(ic, f));
            faces.centroids(ic, f, 2) = z(k);
            faces.normals(ic, f, 0) = 0.0;
            faces.normals(ic, f, 1) = 0.0;
            faces.normals(ic, f, 2) = -1.0;
            faces.neighbours(ic, f) = (k == 0) ? -num_xy_boundaries-1 : ic-num_xy_cells;
            f++;
         }
         
         /* +z face: */
         if (nz > 0) {
            faces.areas(ic, f) = math::area(xy_points, xy_cells(i), n);
            math::centroid(faces.centroids(ic, f), xy_points, xy_cells(i), n, faces.areas(ic, f));
            faces.centroids(ic, f, 2) = z(k)+dz(k);
            faces.normals(ic, f, 0) = 0.0;
            faces.normals(ic, f, 1) = 0.0;
            faces.normals(ic, f, 2) = 1.0;
            faces.neighbours(ic, f) = (k == nz-1) ? -num_xy_boundaries-2 : ic+num_xy_cells;
            f++;
         }
         
         ic++;
         
      }
   }
   
   /* Write all the mesh data for debugging: */
#ifdef DEBUG
   writeData("mesh_data.pmp");
#endif
   
   return 0;
   
}
