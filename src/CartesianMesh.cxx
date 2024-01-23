#include "CartesianMesh.hxx"

/* Read the mesh from a plain-text input file: */
int CartesianMesh::read(const std::string& filename) {
   
   /* Open the input file: */
   std::ifstream file(filename);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty()) break;
      
      /* Get the next keyword: */
      if (line[0] == "dx") {
         
         /* Get the dx values: */
         nx = std::stoi(line[1]);
         dx.clear();
         if (nx > 0) {
            PAMPA_CALL(utils::read(dx, nx, file), "wrong dx data in " + filename);
         }
         else {
            PAMPA_CALL(utils::read(dx, 1, file), "wrong dx data in " + filename);
            nx = -nx;
            dx.reserve(nx);
            for (int i = 1; i < nx; i++)
               dx.push_back(dx[0]);
         }
         
      }
      else if (line[0] == "dy") {
         
         /* Get the dy values: */
         ny = std::stoi(line[1]);
         dy.clear();
         if (ny > 0) {
            PAMPA_CALL(utils::read(dy, ny, file), "wrong dy data in " + filename);
         }
         else {
            PAMPA_CALL(utils::read(dy, 1, file), "wrong dy data in " + filename);
            ny = -ny;
            dy.reserve(ny);
            for (int j = 1; j < ny; j++)
               dy.push_back(dy[0]);
         }
         
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
      else if (line[0] == "bc") {
         
         /* Get the boundary conditions (1-based indexed): */
         /* Note: 1 = -x, 2 = +x, 3 = -y, 4 = +y, 5 = -z, 6 = +z. */
         bcs.resize(7);
         int i = 1;
         std::string dir = line[i++];
         if (dir == "x") {
            bcs[1].type = static_cast<BC::Type>(std::stoi(line[i++])-1);
            if (bcs[1].type == BC::ROBIN) bcs[1].a = std::stod(line[i++]);
            bcs[2].type = static_cast<BC::Type>(std::stoi(line[i++])-1);
            if (bcs[2].type == BC::ROBIN) bcs[2].a = std::stod(line[i++]);
         }
         else if (dir == "y") {
            bcs[3].type = static_cast<BC::Type>(std::stoi(line[i++])-1);
            if (bcs[3].type == BC::ROBIN) bcs[3].a = std::stod(line[i++]);
            bcs[4].type = static_cast<BC::Type>(std::stoi(line[i++])-1);
            if (bcs[4].type == BC::ROBIN) bcs[4].a = std::stod(line[i++]);
         }
         else if (dir == "z") {
            bcs[5].type = static_cast<BC::Type>(std::stoi(line[i++])-1);
            if (bcs[5].type == BC::ROBIN) bcs[5].a = std::stod(line[i++]);
            bcs[6].type = static_cast<BC::Type>(std::stoi(line[i++])-1);
            if (bcs[6].type == BC::ROBIN) bcs[6].a = std::stod(line[i++]);
         }
         else
            PAMPA_CHECK(true, 1, "wrong boundary condition in " + filename);
         
      }
      else if (line[0] == "materials") {
         
         /* Get the material distribution: */
         int num_materials = std::stoi(line[1]);
         int num_cells = nx * std::max(ny, 1) * std::max(nz, 1);
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
   
}

/* Build the mesh: */
int CartesianMesh::build() {
   
   /* Build the mesh points: */
   std::vector<double> x(nx+1), y(ny+1), z(nz+1);
   for (int i = 0; i < nx; i++)
      x[i+1] = x[i] + dx[i];
   for (int j = 0; j < ny; j++)
      y[j+1] = y[j] + dy[j];
   for (int k = 0; k < nz; k++)
      z[k+1] = z[k] + dz[k];
   num_points = (nx+1) * (ny+1) * (nz+1);
   points.reserve(num_points);
   for (int k = 0; k < nz+1; k++)
      for (int j = 0; j < ny+1; j++)
         for (int i = 0; i < nx+1; i++)
            points.push_back(std::vector<double>{x[i], y[j], z[k]});
   
   /* Get the number of physical cells and check the material definition: */
   int num_xy_cells = 0;
   std::vector<std::array<int, 2>> num_x_void_cells(std::max(ny, 1));
   int im = 0;
   for (int k = 0; k < std::max(nz, 1); k++) {
      for (int j = 0; j < std::max(ny, 1); j++) {
         int l = 0;
         for (int i = 0; i < nx; i++) {
            int mat = cells.materials[im];
            if (k == 0) {
               if (mat != -1) {
                  if (l == 0) l++;
                  num_xy_cells++;
               }
               else
                  num_x_void_cells[j][l]++;
            }
            else {
               int mat0 = cells.materials[im-nx*std::max(ny, 1)];
               PAMPA_CHECK((mat != -1) && (mat0 == -1) || (mat == -1) && (mat0 != -1), 1, 
                  "wrong material definition");
            }
            im++;
         }
      }
   }
   num_cells = num_xy_cells * std::max(nz, 1);
   
   /* Build the mesh cells: */
   /* Note: the cell points are ordered according to the gmsh convention. */
   cells.points.reserve(num_cells);
   cells.volumes.reserve(num_cells);
   cells.centroids.reserve(num_cells);
   im = 0;
   for (int k = 0; k < std::max(nz, 1); k++) {
      for (int j = 0; j < std::max(ny, 1); j++) {
         for (int i = 0; i < nx; i++) {
            
            /* Build only physical cells: */
            if (cells.materials[im] != -1) {
               
               /* Get the cell points: */
               int p1 = i + j*(nx+1) + k*(nx+1)*(ny+1);
               int p2 = (i+1) + j*(nx+1) + k*(nx+1)*(ny+1);
               if (ny > 0) {
                  int p3 = (i+1) + (j+1)*(nx+1) + k*(nx+1)*(ny+1);
                  int p4 = i + (j+1)*(nx+1) + k*(nx+1)*(ny+1);
                  if (nz > 0) {
                     int p5 = i + j*(nx+1) + (k+1)*(nx+1)*(ny+1);
                     int p6 = (i+1) + j*(nx+1) + (k+1)*(nx+1)*(ny+1);
                     int p7 = (i+1) + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
                     int p8 = i + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
                     cells.points.push_back(std::vector<int>{p1, p2, p3, p4, p5, p6, p7, p8});
                  }
                  else
                     cells.points.push_back(std::vector<int>{p1, p2, p3, p4});
               }
               else
                  cells.points.push_back(std::vector<int>{p1, p2});
               
               /* Get the cell volume: */
               double v = (nz > 0) ? dx[i] * dy[j] * dz[k] : (ny > 0) ? dx[i] * dy[j] : dx[i];
               cells.volumes.push_back(v);
               
               /* Get the cell centroid: */
               double x0 = x[i] + 0.5*dx[i];
               double y0 = y[j] + 0.5*dy[j];
               double z0 = z[k] + 0.5*dz[k];
               cells.centroids.push_back(std::vector<double>{x0, y0, z0});
               
            }
            
            /* Move to the next material: */
            im++;
            
         }
      }
   }
   
   /* Build the mesh faces: */
   /* Note: the face points are ordered counterclockwise so that the normal points outward.*/
   faces.points.reserve(num_cells);
   faces.areas.reserve(num_cells);
   faces.centroids.reserve(num_cells);
   faces.normals.reserve(num_cells);
   faces.neighbours.reserve(num_cells);
   im = 0;
   int ic = 0;
   for (int k = 0; k < std::max(nz, 1); k++) {
      for (int j = 0; j < std::max(ny, 1); j++) {
         for (int i = 0; i < nx; i++) {
            
            /* Build only physical cells: */
            if (cells.materials[im] != -1) {
               
               /* Initialize the face data for this cell: */
               int num_faces = (nz > 0) ? 6 : (ny > 0) ? 4 : 2, f = 0;
               std::vector<std::vector<int>> pts(num_faces);
               std::vector<double> a(num_faces);
               std::vector<std::vector<double>> p0(num_faces);
               std::vector<std::vector<double>> n(num_faces);
               std::vector<int> ic2(num_faces);
               
               /* -y face: */
               if (ny > 0) {
                  if (nz > 0)
                     pts[f] = math::extrude_edge(cells.points[ic], 0, 4);
                  else
                     pts[f] = std::vector<int>{cells.points[ic][0], cells.points[ic][1]};
                  a[f] = (nz > 0) ? dx[i] * dz[k] : dx[i];
                  p0[f] = std::vector<double>{x[i]+0.5*dx[i], y[j], z[k]+0.5*dz[k]};
                  n[f] = std::vector<double>{0.0, -1.0, 0.0};
                  if (j == 0)
                     ic2[f] = -3;
                  else {
                     if (cells.materials[im-nx] == -1)
                        ic2[f] = -3;
                     else
                        ic2[f] = ic - nx + num_x_void_cells[j][0] + num_x_void_cells[j-1][1];
                  }
                  f++;
               }
               
               /* +x face: */
               if (nz > 0)
                  pts[f] = math::extrude_edge(cells.points[ic], 1, 4);
               else {
                  if (ny > 0)
                     pts[f] = std::vector<int>{cells.points[ic][1], cells.points[ic][2]};
                  else
                     pts[f] = std::vector<int>{cells.points[ic][1]};
               }
               a[f] = (nz > 0) ? dy[j] * dz[k] : (ny > 0) ? dy[j] : 1.0;
               p0[f] = std::vector<double>{x[i]+dx[i], y[j]+0.5*dy[j], z[k]+0.5*dz[k]};
               n[f] = std::vector<double>{1.0, 0.0, 0.0};
               if (i == nx-1)
                  ic2[f] = -2;
               else {
                  if (cells.materials[im+1] == -1)
                     ic2[f] = -2;
                  else
                     ic2[f] = ic + 1;
               }
               f++;
               
               /* +y face: */
               if (ny > 0) {
                  if (nz > 0)
                     pts[f] = math::extrude_edge(cells.points[ic], 2, 4);
                  else
                     pts[f] = std::vector<int>{cells.points[ic][2], cells.points[ic][3]};
                  a[f] = (nz > 0) ? dx[i] * dz[k] : dx[i];
                  p0[f] = std::vector<double>{x[i]+0.5*dx[i], y[j]+dy[j], z[k]+0.5*dz[k]};
                  n[f] = std::vector<double>{0.0, 1.0, 0.0};
                  if (j == ny-1)
                     ic2[f] = -4;
                  else {
                     if (cells.materials[im+nx] == -1)
                        ic2[f] = -4;
                     else
                        ic2[f] = ic + nx - num_x_void_cells[j][1] - num_x_void_cells[j+1][0];
                  }
                  f++;
               }
               
               /* -x face: */
               if (nz > 0)
                  pts[f] = math::extrude_edge(cells.points[ic], 3, 4);
               else {
                  if (ny > 0)
                     pts[f] = std::vector<int>{cells.points[ic][3], cells.points[ic][0]};
                  else
                     pts[f] = std::vector<int>{cells.points[ic][0]};
               }
               a[f] = (nz > 0) ? dy[j] * dz[k] : (ny > 0) ? dy[j] : 1.0;
               p0[f] = std::vector<double>{x[i], y[j]+0.5*dy[j], z[k]+0.5*dz[k]};
               n[f] = std::vector<double>{-1.0, 0.0, 0.0};
               if (i == 0)
                  ic2[f] = -1;
               else {
                  if (cells.materials[im-1] == -1)
                     ic2[f] = -1;
                  else
                     ic2[f] = ic - 1;
               }
               f++;
               
               /* -z face: */
               if (nz > 0) {
                  pts[f] = std::vector<int>(cells.points[ic].begin(), cells.points[ic].begin()+4);
                  std::reverse(pts[f].begin(), pts[f].end());
                  a[f] = dx[i] * dy[j];
                  p0[f] = std::vector<double>{x[i]+0.5*dx[i], y[j]+0.5*dy[j], z[k]};
                  n[f] = std::vector<double>{0.0, 0.0, -1.0};
                  ic2[f] = (k == 0) ? -5 : ic - num_xy_cells;
                  f++;
               }
               
               /* +z face: */
               if (nz > 0) {
                  pts[f] = std::vector<int>(cells.points[ic].begin()+4, cells.points[ic].end());
                  a[f] = dx[i] * dy[j];
                  p0[f] = std::vector<double>{x[i]+0.5*dx[i], y[j]+0.5*dy[j], z[k]+dz[k]};
                  n[f] = std::vector<double>{0.0, 0.0, 1.0};
                  ic2[f] = (k == nz-1) ? -6 : ic + num_xy_cells;
                  f++;
               }
               
               /* Keep the data for this cell: */
               faces.points.push_back(pts);
               faces.areas.push_back(a);
               faces.centroids.push_back(p0);
               faces.normals.push_back(n);
               faces.neighbours.push_back(ic2);
               ic++;
               
            }
            
            /* Move to the next material: */
            im++;
            
         }
      }
   }
   
   /* Remove the unused materials to get the indexing right: */
   cells.materials.erase(std::remove(cells.materials.begin(), cells.materials.end(), -1), 
      cells.materials.end());
   
   /* Write all the mesh data for debugging: */
#ifdef DEBUG
   writeData("mesh_data.pmp");
#endif
   
   return 0;
   
}
