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
      if (line.empty())
         break;
      
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
         PAMPA_CALL(utils::read(dz, nz, file), "wrong dz data in " + filename);
         
      }
      else if (line[0] == "boundary") {
         
         /* Get the cell indices: */
         std::vector<int> boundary;
         int num_xy_bc_points = std::stoi(line[1]);
         PAMPA_CALL(utils::read(boundary, num_xy_bc_points, file), 
            "wrong boundary data in " + filename);
         boundaries.push_back(boundary);
         bcs.resize(boundaries.size());
         
      }
      else if (line[0] == "bc") {
         
         /* Get the boundary condition: */
         std::string dir = line[1];
         if (dir == "z") {
            bc_z[0] = std::stoi(line[2]);
            bc_z[1] = std::stoi(line[3]);
         }
         else {
            int i = std::stoi(dir);
            PAMPA_CHECK(i >= boundaries.size(), 1, 
               "wrong boundary condition direction in " + filename)
            bcs[i] = std::stoi(line[2]);
         }
         
      }
      else if (line[0] == "materials") {
         
         /* Get the material distribution: */
         int num_materials = std::stoi(line[1]);
         PAMPA_CHECK(num_materials != num_xy_cells*nz, 1, 
            "wrong number of materials in " + filename);
         PAMPA_CALL(utils::read(cells.materials, num_materials, file), 
            "wrong material data in " + filename);
         
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
   for (int k = 0; k < nz+1; k++) {
      for (int i = 0; i < num_xy_points; i++) {
         points.push_back(std::vector<double>{xy_points[i][0], xy_points[i][1], z[k]});
      }
   }
   
   /* Build the mesh cells: */
   /* Note: the cell points are ordered according to the gmsh convention. */
   num_cells = num_xy_cells * nz;
   cells.points.reserve(num_cells);
   cells.volumes.reserve(num_cells);
   cells.centroids.reserve(num_cells);
   for (int k = 0; k < nz; k++) {
      for (int i = 0; i < num_xy_cells; i++) {
         int n = xy_cells[i].size();
         std::vector<int> pts;
         pts.reserve(2*n);
         for (int dk = 0; dk < 2; dk++) {
            for (int l = 0; l < n; l++) {
               pts.push_back(xy_cells[i][l]+(k+dk)*num_xy_points);
            }
         }
         double a = math::get_area(points, xy_cells[i]);
         double v = a * dz[k];
         std::vector<double> p0 = math::get_centroid(points, xy_cells[i], a);
         p0.push_back(z[k]+0.5*dz[k]);
         cells.points.push_back(pts);
         cells.volumes.push_back(v);
         cells.centroids.push_back(p0);
      }
   }
   
   /* Get the cells for each point in the xy-plane: */
   std::vector<std::vector<int>> xy_points_to_cells(num_xy_points);
   for (int i = 0; i < num_xy_cells; i++) {
      for (int j = 0; j < xy_cells[i].size(); j++) {
         xy_points_to_cells[xy_cells[i][j]].push_back(i);
      }
   }
   
   /* Set the boundary conditions: */
   for (int i = 0; i < boundaries.size(); i++) {
      for (int j = 0; j < boundaries[i].size(); j++) {
         xy_points_to_cells[boundaries[i][j]].push_back(-bcs[i]);
      }
   }
   
   /* Get the neighbour for each cell face in the xy-plane: */
   std::vector<std::vector<int>> xy_neighbours(num_xy_cells);
   bool found;
   for (int i = 0; i < num_xy_cells; i++) {
      int nxy = xy_cells[i].size();
      xy_neighbours[i].resize(nxy);
      for (int f = 0; f < nxy; f++) {
         found = false;
         for (int l1 = 0; l1 < xy_points_to_cells[xy_cells[i][f]].size(); l1++) {
            for (int l2 = 0; l2 < xy_points_to_cells[xy_cells[i][(f+1)%nxy]].size(); l2++) {
               int i1 = xy_points_to_cells[xy_cells[i][f]][l1];
               int i2 = xy_points_to_cells[xy_cells[i][(f+1)%nxy]][l2];
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
   faces.points.reserve(num_cells);
   faces.areas.reserve(num_cells);
   faces.centroids.reserve(num_cells);
   faces.normals.reserve(num_cells);
   faces.neighbours.reserve(num_cells);
   int l = 0;
   for (int k = 0; k < nz; k++) {
      for (int i = 0; i < num_xy_cells; i++) {
         
         /* Initialize the face data for this cell: */
         int nxy = xy_cells[i].size();
         std::vector<std::vector<int>> pts(nxy+2);
         std::vector<double> a(nxy+2);
         std::vector<std::vector<double>> p0(nxy+2);
         std::vector<std::vector<double>> n(nxy+2);
         std::vector<int> l2(nxy+2);
         
         /* xy-plane faces: */
         for (int f = 0; f < nxy; f++) {
            pts[f] = math::extrude_edge(cells.points[l], f, nxy);
            a[f] = math::get_distance(points, xy_cells[i][f], xy_cells[i][(f+1)%nxy], 2);
            a[f] *= dz[k];
            p0[f] = math::get_midpoint(points, xy_cells[i][f], xy_cells[i][(f+1)%nxy], 2);
            p0[f].push_back(z[k]+0.5*dz[k]);
            n[f] = math::get_normal(points, xy_cells[i][f], xy_cells[i][(f+1)%nxy]);
            n[f].push_back(0.0);
            l2[f] = xy_neighbours[i][f];
            if (l2[f] > 0)
               l2[f] += k * num_xy_cells;
         }
         
         /* -z face: */
         pts[nxy] = std::vector<int>(cells.points[l].rbegin()+nxy, cells.points[l].rend()+2*nxy);
         a[nxy] = math::get_area(points, xy_cells[i]);
         p0[nxy] = math::get_centroid(points, xy_cells[i], a[nxy]);
         p0[nxy].push_back(z[k]);
         n[nxy] = std::vector<double>{0.0, 0.0, -1.0};
         l2[nxy] = (k == 0) ? -bc_z[0] : l - num_xy_cells;
         
         /* +z face: */
         pts[nxy+1] = std::vector<int>(cells.points[l].begin()+nxy, cells.points[l].end());
         a[nxy+1] = math::get_area(points, xy_cells[i]);
         p0[nxy+1] = math::get_centroid(points, xy_cells[i], a[nxy+1]);
         p0[nxy+1].push_back(z[k]+dz[k]);
         n[nxy+1] = std::vector<double>{0.0, 0.0, 1.0};
         l2[nxy+1] = (k == nz-1) ? -bc_z[1] : l + num_xy_cells;
         
         /* Keep the data for this cell: */
         faces.points.push_back(pts);
         faces.areas.push_back(a);
         faces.centroids.push_back(p0);
         faces.normals.push_back(n);
         faces.neighbours.push_back(l2);
         l++;
         
      }
   }
   
   return 0;
   
};
