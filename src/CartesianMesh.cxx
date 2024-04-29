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
      unsigned int l = 0;
      if (line[l] == "dx") {
         
         /* Get the dx values: */
         PAMPA_CALL(utils::read(nx, -INT_MAX, INT_MAX, line[++l]), "wrong nx value");
         if (nx > 0) {
            PAMPA_CALL(utils::read(dx, nx, file), "wrong dx data");
         }
         else {
            PAMPA_CALL(utils::read(dx, 1, file), "wrong dx data");
            nx = -nx;
            dx.resize(nx, dx(0));
         }
         num_dims++;
         
      }
      else if (line[l] == "dy") {
         
         /* Get the dy values: */
         PAMPA_CALL(utils::read(ny, -INT_MAX, INT_MAX, line[++l]), "wrong ny value");
         if (ny > 0) {
            PAMPA_CALL(utils::read(dy, ny, file), "wrong dy data");
         }
         else {
            PAMPA_CALL(utils::read(dy, 1, file), "wrong dy data");
            ny = -ny;
            dy.resize(ny, dy(0));
         }
         if (ny > 1) num_dims++;
         
      }
      else if (line[l] == "dz") {
         
         /* Get the dz values: */
         PAMPA_CALL(utils::read(nz, -INT_MAX, INT_MAX, line[++l]), "wrong nz value");
         if (nz > 0) {
            PAMPA_CALL(utils::read(dz, nz, file), "wrong dz data");
         }
         else {
            PAMPA_CALL(utils::read(dz, 1, file), "wrong dz data");
            nz = -nz;
            dz.resize(nz, dz(0));
         }
         if (nz > 1) num_dims++;
         
      }
      else if (line[l] == "bc") {
         
         /* Initialize the boundary-condition array, if not done yet: */
         if (bcs.empty()) bcs.resize(7);
         
         /* Get the boundary conditions (1-based indexed): */
         /* Note: 1 = -x, 2 = +x, 3 = -y, 4 = +y, 5 = -z, 6 = +z. */
         std::string dir = line[++l];
         if (dir == "x") {
            PAMPA_CALL(utils::read(bcs(1), line, ++l, file), "wrong boundary condition");
            PAMPA_CALL(utils::read(bcs(2), line, l, file), "wrong boundary condition");
         }
         else if (dir == "y") {
            PAMPA_CALL(utils::read(bcs(3), line, ++l, file), "wrong boundary condition");
            PAMPA_CALL(utils::read(bcs(4), line, l, file), "wrong boundary condition");
         }
         else if (dir == "z") {
            PAMPA_CALL(utils::read(bcs(5), line, ++l, file), "wrong boundary condition");
            PAMPA_CALL(utils::read(bcs(6), line, l, file), "wrong boundary condition");
         }
         else {
            PAMPA_CHECK(true, 2, "wrong boundary condition");
         }
         
      }
      else if (line[l] == "materials") {
         
         /* Get the material distribution (1-based indexed): */
         int num_materials, num_cells = nx * std::max(ny, 1) * std::max(nz, 1);
         PAMPA_CALL(utils::read(num_materials, num_cells, num_cells, line[++l]), 
            "wrong number of materials");
         PAMPA_CALL(utils::read(cells.materials, num_cells, file), "wrong material data");
         
         /* Switch the material index from 1-based to 0-based: */
         for (int i = 0; i < num_cells; i++)
            cells.materials(i)--;
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 3, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   return 0;
   
}

/* Build the mesh: */
int CartesianMesh::build() {
   
   /* Build the axis coordinates: */
   Array1D<double> x(nx+1), y(ny+1), z(nz+1);
   for (int i = 0; i < nx; i++)
      x(i+1) = x(i) + dx(i);
   for (int j = 0; j < ny; j++)
      y(j+1) = y(j) + dy(j);
   for (int k = 0; k < nz; k++)
      z(k+1) = z(k) + dz(k);
   
   /* Build the mesh points: */
   num_points = (nx+1) * (ny+1) * (nz+1);
   points.resize(num_points, 3);
   for (int ip = 0, k = 0; k < nz+1; k++) {
      for (int j = 0; j < ny+1; j++) {
         for (int i = 0; i < nx+1; i++) {
            points(ip, 0) = x(i);
            points(ip, 1) = y(j);
            points(ip, 2) = z(k);
            ip++;
         }
      }
   }
   
   /* Get the number of physical cells and check the material definition: */
   Array2D<int> num_x_void_cells(std::max(ny, 1), 2);
   int num_xy_cells = 0;
   for (int im = 0, k = 0; k < std::max(nz, 1); k++) {
      for (int j = 0; j < std::max(ny, 1); j++) {
         int l = 0;
         for (int i = 0; i < nx; i++) {
            
            /* Check the material: */
            int mat = cells.materials(im);
            if (k == 0) {
               if (mat != -1) {
                  if (l == 0) l++;
                  num_xy_cells++;
               }
               else
                  num_x_void_cells(j, l)++;
            }
            else {
               int mat0 = cells.materials(im-nx*std::max(ny, 1));
               PAMPA_CHECK(((mat != -1) && (mat0 == -1)) || ((mat == -1) && (mat0 != -1)), 1, 
                  "wrong material definition");
            }
            
            im++;
            
         }
      }
   }
   num_cells = num_xy_cells * std::max(nz, 1);
   num_cells_global = num_cells;
   
   /* Build the mesh cells: */
   /* Note: the cell points are ordered according to the gmsh convention. */
   int num_cell_points = (nz > 0) ? 8 : (ny > 0) ? 4 : 2;
   cells.points.resize(num_cells, Array1D<int>(num_cells, num_cell_points));
   cells.volumes.resize(num_cells);
   cells.centroids.resize(num_cells, 3);
   cells.indices.resize(num_cells);
   for (int ic = 0, im = 0, k = 0; k < std::max(nz, 1); k++) {
      for (int j = 0; j < std::max(ny, 1); j++) {
         for (int i = 0; i < nx; i++) {
            
            /* Build only physical cells: */
            if (cells.materials(im) != -1) {
               
               /* Get the cell points: */
               cells.points(ic, 0) = i + j*(nx+1) + k*(nx+1)*(ny+1);
               cells.points(ic, 1) = (i+1) + j*(nx+1) + k*(nx+1)*(ny+1);
               if (ny > 0) {
                  cells.points(ic, 2) = (i+1) + (j+1)*(nx+1) + k*(nx+1)*(ny+1);
                  cells.points(ic, 3) = i + (j+1)*(nx+1) + k*(nx+1)*(ny+1);
                  if (nz > 0) {
                     cells.points(ic, 4) = i + j*(nx+1) + (k+1)*(nx+1)*(ny+1);
                     cells.points(ic, 5) = (i+1) + j*(nx+1) + (k+1)*(nx+1)*(ny+1);
                     cells.points(ic, 6) = (i+1) + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
                     cells.points(ic, 7) = i + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
                  }
               }
               
               /* Get the cell volume: */
               cells.volumes(ic) = (nz > 0) ? dx(i)*dy(j)*dz(k) : (ny > 0) ? dx(i)*dy(j) : dx(i);
               
               /* Get the cell centroid: */
               cells.centroids(ic, 0) = x(i) + 0.5*dx(i);
               cells.centroids(ic, 1) = y(j) + 0.5*dy(j);
               cells.centroids(ic, 2) = z(k) + 0.5*dz(k);
               
               /* Get the cell index in the global mesh: */
               cells.indices(ic) = ic;
               
               ic++;
               
            }
            
            im++;
            
         }
      }
   }
   
   /* Get the number of cell faces: */
   num_faces_max = (nz > 0) ? 6 : (ny > 0) ? 4 : 2;
   faces.num_faces.resize(num_cells, num_faces_max);
   
   /* Build the mesh faces: */
   /* Note: the face points are ordered counterclockwise so that the normal points outward.*/
   faces.areas.resize(num_cells, faces.num_faces);
   faces.centroids.resize(num_cells, faces.num_faces, 3);
   faces.normals.resize(num_cells, faces.num_faces, 3);
   faces.neighbours.resize(num_cells, faces.num_faces);
   for (int ic = 0, im = 0, k = 0; k < std::max(nz, 1); k++) {
      for (int j = 0; j < std::max(ny, 1); j++) {
         for (int i = 0; i < nx; i++) {
            
            /* Build only physical cells: */
            int f = 0;
            if (cells.materials(im) != -1) {
               
               /* -y face: */
               if (ny > 0) {
                  faces.areas(ic, f) = (nz > 0) ? dx(i)*dz(k) : dx(i);
                  faces.centroids(ic, f, 0) = x(i) + 0.5*dx(i);
                  faces.centroids(ic, f, 1) = y(j);
                  faces.centroids(ic, f, 2) = z(k) + 0.5*dz(k);
                  faces.normals(ic, f, 0) = 0.0;
                  faces.normals(ic, f, 1) = -1.0;
                  faces.normals(ic, f, 2) = 0.0;
                  if (j == 0)
                     faces.neighbours(ic, f) = -3;
                  else {
                     if (cells.materials(im-nx) == -1)
                        faces.neighbours(ic, f) = -3;
                     else
                        faces.neighbours(ic, f) = 
                           ic - nx + num_x_void_cells(j, 0) + num_x_void_cells(j-1, 1);
                  }
                  f++;
               }
               
               /* +x face: */
               faces.areas(ic, f) = (nz > 0) ? dy(j)*dz(k) : (ny > 0) ? dy(j) : 1.0;
               faces.centroids(ic, f, 0) = x(i) + dx(i);
               faces.centroids(ic, f, 1) = y(j) + 0.5*dy(j);
               faces.centroids(ic, f, 2) = z(k) + 0.5*dz(k);
               faces.normals(ic, f, 0) = 1.0;
               faces.normals(ic, f, 1) = 0.0;
               faces.normals(ic, f, 2) = 0.0;
               if (i == nx-1)
                  faces.neighbours(ic, f) = -2;
               else {
                  if (cells.materials(im+1) == -1)
                     faces.neighbours(ic, f) = -2;
                  else
                     faces.neighbours(ic, f) = ic + 1;
               }
               f++;
               
               /* +y face: */
               if (ny > 0) {
                  faces.areas(ic, f) = (nz > 0) ? dx(i)*dz(k) : dx(i);
                  faces.centroids(ic, f, 0) = x(i) + 0.5*dx(i);
                  faces.centroids(ic, f, 1) = y(j) + dy(j);
                  faces.centroids(ic, f, 2) = z(k) + 0.5*dz(k);
                  faces.normals(ic, f, 0) = 0.0;
                  faces.normals(ic, f, 1) = 1.0;
                  faces.normals(ic, f, 2) = 0.0;
                  if (j == ny-1)
                     faces.neighbours(ic, f) = -4;
                  else {
                     if (cells.materials(im+nx) == -1)
                        faces.neighbours(ic, f) = -4;
                     else
                        faces.neighbours(ic, f) = 
                           ic + nx - num_x_void_cells(j, 1) - num_x_void_cells(j+1, 0);
                  }
                  f++;
               }
               
               /* -x face: */
               faces.areas(ic, f) = (nz > 0) ? dy(j)*dz(k) : (ny > 0) ? dy(j) : 1.0;
               faces.centroids(ic, f, 0) = x(i);
               faces.centroids(ic, f, 1) = y(j) + 0.5*dy(j);
               faces.centroids(ic, f, 2) = z(k) + 0.5*dz(k);
               faces.normals(ic, f, 0) = -1.0;
               faces.normals(ic, f, 1) = 0.0;
               faces.normals(ic, f, 2) = 0.0;
               if (i == 0)
                  faces.neighbours(ic, f) = -1;
               else {
                  if (cells.materials(im-1) == -1)
                     faces.neighbours(ic, f) = -1;
                  else
                     faces.neighbours(ic, f) = ic - 1;
               }
               f++;
               
               /* -z face: */
               if (nz > 0) {
                  faces.areas(ic, f) = dx(i) * dy(j);
                  faces.centroids(ic, f, 0) = x(i) + 0.5*dx(i);
                  faces.centroids(ic, f, 1) = y(j) + 0.5*dy(j);
                  faces.centroids(ic, f, 2) = z(k);
                  faces.normals(ic, f, 0) = 0.0;
                  faces.normals(ic, f, 1) = 0.0;
                  faces.normals(ic, f, 2) = -1.0;
                  faces.neighbours(ic, f) = (k == 0) ? -5 : ic-num_xy_cells;
                  f++;
               }
               
               /* +z face: */
               if (nz > 0) {
                  faces.areas(ic, f) = dx(i) * dy(j);
                  faces.centroids(ic, f, 0) = x(i) + 0.5*dx(i);
                  faces.centroids(ic, f, 1) = y(j) + 0.5*dy(j);
                  faces.centroids(ic, f, 2) = z(k) + dz(k);
                  faces.normals(ic, f, 0) = 0.0;
                  faces.normals(ic, f, 1) = 0.0;
                  faces.normals(ic, f, 2) = 1.0;
                  faces.neighbours(ic, f) = (k == nz-1) ? -6 : ic+num_xy_cells;
                  f++;
               }
               
               ic++;
               
            }
            
            im++;
            
         }
      }
   }
   
   /* Remove the unused materials to get the indexing right: */
   cells.materials.remove(-1);
   
   /* Write all the mesh data for debugging: */
#ifdef DEBUG
   writeData("mesh_data.pmp");
#endif
   
   return 0;
   
}
