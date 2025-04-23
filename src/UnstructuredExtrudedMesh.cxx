#include "UnstructuredExtrudedMesh.hxx"

/* Read the mesh from a plain-text input file: */
int UnstructuredExtrudedMesh::read(const std::string& filename) {
   
   /* Open the input file: */
   std::ifstream file(filename, std::ios_base::in);
   PAMPA_CHECK(!file.is_open(), "unable to open " + filename);
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = input::get_next_line(file);
      if (line.empty()) break;
      
      /* Get the next keyword: */
      unsigned int l = 0;
      if (line[l] == "points") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the point coordinates: */
         PAMPA_CHECK(input::read(num_xy_points, 1, INT_MAX, line[++l]), "wrong number of points");
         PAMPA_CHECK(input::read(xy_points, num_xy_points, 2, -DBL_MAX, DBL_MAX, file), 
            "wrong point data");
         
      }
      else if (line[l] == "cells") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 3, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the cell points: */
         int num_xy_cell_points;
         PAMPA_CHECK(input::read(num_xy_cells, 1, INT_MAX, line[++l]), "wrong number of cells");
         PAMPA_CHECK(input::read(num_xy_cell_points, 1, INT_MAX, line[++l]), 
            "wrong number of cell points");
         PAMPA_CHECK(input::read(xy_cells, num_xy_cells, num_xy_cell_points, 0, INT_MAX, file), 
            "wrong cell data");
         num_dims += 2;
         
      }
      else if (line[l] == "dz") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the dz values: */
         PAMPA_CHECK(input::read(nz, -INT_MAX, INT_MAX, line[++l]), "wrong nz value");
         if (nz > 0) {
            PAMPA_CHECK(input::read(dz, nz, 0, DBL_MAX, file), "wrong dz data");
         }
         else {
            PAMPA_CHECK(input::read(dz, 1, 0, DBL_MAX, file), "wrong dz data");
            nz = -nz;
            dz.resize(nz, dz(0));
         }
         if (nz > 1) num_dims++;
         
         /* Add the z boundaries: */
         boundaries.pushBack("-z");
         boundaries.pushBack("+z");
         
      }
      else if (line[l] == "boundary") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 3, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the boundary name: */
         std::string name = line[++l];
         boundaries.pushBack(name);
         xy_boundary_names.pushBack(name);
         
         /* Get the cell indices: */
         Array1D<int> xy_boundary;
         int num_xy_boundary_points;
         PAMPA_CHECK(input::read(num_xy_boundary_points, 0, INT_MAX, line[++l]), 
            "wrong number of boundary points");
         if (num_xy_boundary_points > 0) {
            PAMPA_CHECK(input::read(xy_boundary, num_xy_boundary_points, 0, INT_MAX, file), 
               "wrong boundary data");
         }
         else
            xy_default_boundary = num_xy_boundaries;
         xy_boundary_points.pushBack(xy_boundary);
         num_xy_boundaries++;
         
      }
      else if (line[l] == "bc") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() < 3, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Initialize the boundary-condition array, if not done yet: */
         if (bcs.empty()) bcs.resize(1+boundaries.size());
         
         /* Get the boundary name and index: */
         std::string name = line[++l];
         int ibc = boundaries.find(name);
         PAMPA_CHECK(ibc < 0, "wrong boundary name");
         
         /* Get the boundary condition (1-based indexed): */
         PAMPA_CHECK(input::read(bcs(ibc+1), line, ++l, file), "wrong boundary condition");
         
      }
      else if (line[l] == "materials") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the material distribution (1-based indexed): */
         int num_materials, num_cells = num_xy_cells * std::max(nz, 1);
         PAMPA_CHECK(input::read(num_materials, num_cells, num_cells, line[++l]), 
            "wrong number of materials");
         PAMPA_CHECK(input::read(cells.materials, num_cells, 0, INT_MAX, file), 
            "wrong material data");
         
         /* Switch the material index from 1-based to 0-based: */
         for (int i = 0; i < num_cells; i++)
            cells.materials(i)--;
         
      }
      else if (line[l] == "nodal-indices") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the nodal indices: */
         int num_indices, num_cells = num_xy_cells * std::max(nz, 1);
         PAMPA_CHECK(input::read(num_indices, num_cells, num_cells, line[++l]), 
            "wrong number of nodal indices");
         PAMPA_CHECK(input::read(cells.nodal_indices, num_cells, -INT_MAX, INT_MAX, file), 
            "wrong nodal-index data");
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, "unrecognized keyword '" + line[l] + "'");
         
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
   num_cells_global = num_cells;
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
   cells.global_indices.resize(num_cells);
   for (int ic = 0, k = 0; k < std::max(nz, 1); k++) {
      for (int i = 0; i < num_xy_cells; i++) {
         
         /* Get the cell points: */
         int n = xy_cells.size(i);
         for (int l = 0; l < n; l++)
            cells.points(ic, l) = xy_cells(i, l) + k*num_xy_points;
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
         cells.global_indices(ic) = ic;
         
         ic++;
         
      }
   }
   
   /* Get the number of cells and boundary conditions for each point in the xy-plane: */
   Array1D<int> num_xy_point_cells(num_xy_points);
   for (int i = 0; i < num_xy_cells; i++)
      for (int j = 0; j < xy_cells.size(i); j++)
         num_xy_point_cells(xy_cells(i, j))++;
   for (int i = 0; i < num_xy_boundaries; i++)
      for (int j = 0; j < xy_boundary_points.size(i); j++)
         num_xy_point_cells(xy_boundary_points(i, j))++;
   
   /* Get the cells and boundary conditions (1-based indexed) for each point in the xy-plane: */
   Vector2D<int> xy_point_cells(num_xy_points, num_xy_point_cells);
   Array1D<int> ipc(num_xy_points);
   for (int i = 0; i < num_xy_cells; i++)
      for (int j = 0; j < xy_cells.size(i); j++)
         xy_point_cells(xy_cells(i, j), ipc(xy_cells(i, j))++) = i;
   for (int i = 0; i < num_xy_boundaries; i++) {
      int boundary_index = boundaries.find(xy_boundary_names(i));
      PAMPA_CHECK(boundary_index < 0, "wrong boundary name");
      for (int j = 0; j < xy_boundary_points.size(i); j++)
         xy_point_cells(xy_boundary_points(i, j), ipc(xy_boundary_points(i, j))++) = 
            -boundary_index - 1;
   }
   
   /* Get the neighboring cell for each face in the xy-plane: */
   Vector2D<int> xy_neighbors(num_xy_cells, num_xy_cell_points);
   for (int i = 0; i < num_xy_cells; i++) {
      for (int f = 0; f < num_xy_cell_points(i); f++) {
         
         /* Get the cell indices for both points in this face: */
         int ip1 = xy_cells(i, f);
         int ip2 = xy_cells(i, (f+1)%num_xy_cell_points(i));
         
         /* Find the cell connected to both points in this face: */
         bool found = false;
         for (int l1 = 0; l1 < xy_point_cells.size(ip1); l1++) {
            for (int l2 = 0; l2 < xy_point_cells.size(ip2); l2++) {
               int i1 = xy_point_cells(ip1, l1);
               int i2 = xy_point_cells(ip2, l2);
               if (i1 == i2 && i1 != i) {
                  PAMPA_CHECK(xy_neighbors(i, f) > 0, "wrong mesh connectivity");
                  xy_neighbors(i, f) = i1;
                  found = true;
               }
            }
         }
         if (!found && xy_default_boundary >= 0) {
            xy_neighbors(i, f) = -xy_default_boundary - 1;
            found = true;
         }
         PAMPA_CHECK(!found, "wrong mesh connectivity");
         
      }
   }
   
   /* Get the number of cell faces: */
   num_faces_max = (nz > 0) ? num_xy_cell_points_max+2 : num_xy_cell_points_max;
   faces.num_faces.resize(num_cells);
   for (int ic = 0, k = 0; k < std::max(nz, 1); k++) {
      for (int i = 0; i < num_xy_cells; i++) {
         int n = xy_cells.size(i);
         faces.num_faces(ic++) = (nz > 0) ? n+2 : n;
      }
   }
   
   /* Get the boundaries for the -z and +z directions: */
   Array1D<std::string> directions;
   if (nz > 0) {
      directions.pushBack("-z");
      directions.pushBack("+z");
   }
   Array1D<int> dir_indices(directions.size());
   for (int i = 0; i < directions.size(); i++) {
      dir_indices(i) = boundaries.find(directions(i));
      PAMPA_CHECK(dir_indices(i) < 0, "wrong boundary name");
   }
   
   /* Build the mesh faces: */
   /* Note: the face points are ordered counterclockwise so that the normal points outward.*/
   faces.areas.resize(num_cells, faces.num_faces);
   faces.centroids.resize(num_cells, faces.num_faces, 3);
   faces.normals.resize(num_cells, faces.num_faces, 3);
   faces.neighbors.resize(num_cells, faces.num_faces);
   for (int ic = 0, k = 0; k < std::max(nz, 1); k++) {
      for (int i = 0; i < num_xy_cells; i++) {
         
         /* Get the number of faces for this cell: */
         int n = xy_cells.size(i);
         
         /* xy-plane faces: */
         int f;
         for (f = 0; f < n; f++) {
		      int f2 = (f+1) % n;
            faces.areas(ic, f) = math::distance(xy_points, xy_cells(i, f), xy_cells(i, f2), 2);
            if (nz > 0) faces.areas(ic, f) *= dz(k);
            math::midpoint(faces.centroids(ic, f), xy_points, xy_cells(i, f), xy_cells(i, f2), 2);
            faces.centroids(ic, f, 2) = z(k) + 0.5*dz(k);
            math::normal(faces.normals(ic, f), xy_points, xy_cells(i, f), xy_cells(i, f2));
            faces.normals(ic, f, 2) = 0.0;
            faces.neighbors(ic, f) = xy_neighbors(i, f);
            if (faces.neighbors(ic, f) >= 0) faces.neighbors(ic, f) += k * num_xy_cells;
         }
         
         /* -z face: */
         if (nz > 0) {
            faces.areas(ic, f) = math::area(xy_points, xy_cells(i), n);
            math::centroid(faces.centroids(ic, f), xy_points, xy_cells(i), n, faces.areas(ic, f));
            faces.centroids(ic, f, 2) = z(k);
            faces.normals(ic, f, 0) = 0.0;
            faces.normals(ic, f, 1) = 0.0;
            faces.normals(ic, f, 2) = -1.0;
            faces.neighbors(ic, f) = (k == 0) ? -dir_indices(0)-1 : ic - num_xy_cells;
            f++;
         }
         
         /* +z face: */
         if (nz > 0) {
            faces.areas(ic, f) = math::area(xy_points, xy_cells(i), n);
            math::centroid(faces.centroids(ic, f), xy_points, xy_cells(i), n, faces.areas(ic, f));
            faces.centroids(ic, f, 2) = z(k) + dz(k);
            faces.normals(ic, f, 0) = 0.0;
            faces.normals(ic, f, 1) = 0.0;
            faces.normals(ic, f, 2) = 1.0;
            faces.neighbors(ic, f) = (k == nz-1) ? -dir_indices(1)-1 : ic + num_xy_cells;
            f++;
         }
         
         ic++;
         
      }
   }
   
   /* Write all the mesh data for debugging: */
#ifdef DEBUG
   PAMPA_CHECK(writeData("mesh_data.pmp"), "unable to write the mesh data");
#endif
   
   return 0;
   
}
