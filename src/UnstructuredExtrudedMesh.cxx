#include "UnstructuredExtrudedMesh.hxx"

/* The UnstructuredExtrudedMesh constructor: */
UnstructuredExtrudedMesh::UnstructuredExtrudedMesh() {};

/* The UnstructuredExtrudedMesh destructor: */
UnstructuredExtrudedMesh::~UnstructuredExtrudedMesh() {};

/* Read the mesh from a plain-text input file: */
int UnstructuredExtrudedMesh::read(const std::string &filename) {
   
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
         dz.clear();
         if (nz > 0) {
            PAMPA_CALL(utils::read(dz, nz, file), "wrong dz data in " + filename);
         }
         else {
            PAMPA_CALL(utils::read(dz, 1, file), "wrong dz data in " + filename);
            nz = -nz;
            dz.reserve(nz);
            for (int k = 1; k < nz; k++)
               dz.push_back(dz[0]);
         }
         
      }
      else if (line[0] == "boundary") {
         
         /* Get the cell indices: */
         std::vector<int> xy_boundary;
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
            cells.materials[i]--;
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 1, "unrecognized keyword '" + line[0] + "' in " + filename);
         
      }
      
   }
   
   return 0;
   
};

/* Build the mesh: */
int UnstructuredExtrudedMesh::build() {
   
   /* Build the mesh points: */
   std::vector<double> z(nz+1);
   for (int k = 0; k < nz; k++)
      z[k+1] = z[k] + dz[k];
   num_points = num_xy_points * (nz+1);
   points.reserve(num_points);
   for (int k = 0; k < nz+1; k++)
      for (int i = 0; i < num_xy_points; i++)
         points.push_back(std::vector<double>{xy_points[i][0], xy_points[i][1], z[k]});
   
   /* Build the mesh cells: */
   /* Note: the cell points are ordered according to the gmsh convention. */
   num_cells = num_xy_cells * std::max(nz, 1);
   cells.points.reserve(num_cells);
   cells.volumes.reserve(num_cells);
   cells.centroids.reserve(num_cells);
   for (int k = 0; k < std::max(nz, 1); k++)
      for (int i = 0; i < num_xy_cells; i++) {
         
         /* Get the cell points: */
         if (nz > 0) {
            int n = xy_cells[i].size();
            std::vector<int> pts;
            pts.reserve(2*n);
            for (int dk = 0; dk < 2; dk++)
               for (int l = 0; l < n; l++)
                  pts.push_back(xy_cells[i][l]+(k+dk)*num_xy_points);
            cells.points.push_back(pts);
         }
         else
            cells.points.push_back(xy_cells[i]);
         
         /* Get the cell volume: */
         double a = math::get_area(points, xy_cells[i]);
         double v = (nz > 0) ? a * dz[k] : a;
         cells.volumes.push_back(v);
         
         /* Get the cell centroid: */
         std::vector<double> p0 = math::get_centroid(points, xy_cells[i], a);
         p0.push_back(z[k]+0.5*dz[k]);
         cells.centroids.push_back(p0);
         
      }
   
   /* Get the cells for each point in the xy-plane: */
   std::vector<std::vector<int>> xy_points_to_cells(num_xy_points);
   for (int i = 0; i < num_xy_cells; i++)
      for (int j = 0; j < xy_cells[i].size(); j++)
         xy_points_to_cells[xy_cells[i][j]].push_back(i);
   
   /* Set the boundary conditions (1-based indexed): */
   for (int i = 0; i < num_xy_boundaries; i++)
      for (int j = 0; j < xy_boundaries[i].size(); j++)
         xy_points_to_cells[xy_boundaries[i][j]].push_back(-i-1);
   
   /* Get the neighboring cell for each cell face in the xy-plane: */
   std::vector<std::vector<int>> xy_neighbours(num_xy_cells);
   bool found;
   for (int i = 0; i < num_xy_cells; i++) {
      
      /* Get the neighbouring cell for each face of this cell: */
      int num_xy_cell_points = xy_cells[i].size();
      xy_neighbours[i].resize(num_xy_cell_points);
      for (int f = 0; f < num_xy_cell_points; f++) {
         
         /* Get the cells for both points in this face: */
         const std::vector<int> &cls1 = xy_points_to_cells[xy_cells[i][f]];
         const std::vector<int> &cls2 = xy_points_to_cells[xy_cells[i][(f+1)%num_xy_cell_points]];
         
         /* Find the cell connected to both points in this face: */
         found = false;
         for (int l1 = 0; l1 < cls1.size(); l1++)
            for (int l2 = 0; l2 < cls2.size(); l2++) {
               int i1 = cls1[l1];
               int i2 = cls2[l2];
               if (i1 == i2 && i1 != i) {
                  PAMPA_CHECK(xy_neighbours[i][f] > 0, 1, "wrong mesh connectivity");
                  xy_neighbours[i][f] = i1;
                  found = true;
               }
            }
         PAMPA_CHECK(!found, 1, "wrong mesh connectivity");
         
      }
   }
   
   /* Build the mesh faces: */
   /* Note: the face points are ordered counterclockwise so that the normal points outward.*/
   faces.points.reserve(num_cells);
   faces.areas.reserve(num_cells);
   faces.centroids.reserve(num_cells);
   faces.normals.reserve(num_cells);
   faces.neighbours.reserve(num_cells);
   int ic = 0;
   for (int k = 0; k < std::max(nz, 1); k++)
      for (int i = 0; i < num_xy_cells; i++) {
         
         /* Initialize the face data for this cell: */
         int nxy = xy_cells[i].size();
         int num_faces = (nz > 0) ? nxy + 2 : nxy;
         std::vector<std::vector<int>> pts(num_faces);
         std::vector<double> a(num_faces);
         std::vector<std::vector<double>> p0(num_faces);
         std::vector<std::vector<double>> n(num_faces);
         std::vector<int> ic2(num_faces);
         
         /* xy-plane faces: */
         for (int f = 0; f < nxy; f++) {
            if (nz > 0)
               pts[f] = math::extrude_edge(cells.points[ic], f, nxy);
            else
               pts[f] = std::vector<int>{cells.points[ic][f], cells.points[ic][(f+1)%nxy]};
            a[f] = math::get_distance(points, xy_cells[i][f], xy_cells[i][(f+1)%nxy], 2);
            if (nz > 0) a[f] *= dz[k];
            p0[f] = math::get_midpoint(points, xy_cells[i][f], xy_cells[i][(f+1)%nxy], 2);
            p0[f].push_back(z[k]+0.5*dz[k]);
            n[f] = math::get_normal(points, xy_cells[i][f], xy_cells[i][(f+1)%nxy]);
            n[f].push_back(0.0);
            ic2[f] = xy_neighbours[i][f];
            if (ic2[f] >= 0) ic2[f] += k * num_xy_cells;
         }
         
         /* -z face: */
         if (nz > 0) {
            pts[nxy] = std::vector<int>(cells.points[ic].begin(), cells.points[ic].begin()+nxy);
            std::reverse(pts[nxy].begin(), pts[nxy].end());
            a[nxy] = math::get_area(points, xy_cells[i]);
            p0[nxy] = math::get_centroid(points, xy_cells[i], a[nxy]);
            p0[nxy].push_back(z[k]);
            n[nxy] = std::vector<double>{0.0, 0.0, -1.0};
            ic2[nxy] = (k == 0) ? -num_xy_boundaries - 1 : ic - num_xy_cells;
         }
         
         /* +z face: */
         if (nz > 0) {
            pts[nxy+1] = std::vector<int>(cells.points[ic].begin()+nxy, cells.points[ic].end());
            a[nxy+1] = math::get_area(points, xy_cells[i]);
            p0[nxy+1] = math::get_centroid(points, xy_cells[i], a[nxy+1]);
            p0[nxy+1].push_back(z[k]+dz[k]);
            n[nxy+1] = std::vector<double>{0.0, 0.0, 1.0};
            ic2[nxy+1] = (k == nz-1) ? -num_xy_boundaries - 2 : ic + num_xy_cells;
         }
         
         /* Keep the data for this cell: */
         faces.points.push_back(pts);
         faces.areas.push_back(a);
         faces.centroids.push_back(p0);
         faces.normals.push_back(n);
         faces.neighbours.push_back(ic2);
         ic++;
         
      }
   
   /* Write all the mesh data for debugging: */
#ifdef DEBUG
   writeData("mesh_data.pmp");
#endif
   
   return 0;
   
};
