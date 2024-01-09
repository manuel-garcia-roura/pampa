#include "CartesianMesh.hxx"

/* The CartesianMesh constructor: */
CartesianMesh::CartesianMesh() {};

/* The CartesianMesh destructor: */
CartesianMesh::~CartesianMesh() {};

/* Read the mesh from a plain-text input file: */
int CartesianMesh::read(const std::string &filename) {
   
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
         int num_cells = (nz > 0) ? nx * ny * nz : nx * ny;
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
   for (int k = 0; k < nz+1; k++) {
      for (int j = 0; j < ny+1; j++) {
         for (int i = 0; i < nx+1; i++) {
            points.push_back(std::vector<double>{x[i], y[j], z[k]});
         }
      }
   }
   
   /* Build the mesh cells: */
   /* Note: the cell points are ordered according to the gmsh convention. */
   num_cells = (nz > 0) ? nx * ny * nz : nx * ny;
   cells.points.reserve(num_cells);
   cells.volumes.reserve(num_cells);
   cells.centroids.reserve(num_cells);
   for (int k = 0; k < num_cells/(nx*ny); k++) {
      for (int j = 0; j < ny; j++) {
         for (int i = 0; i < nx; i++) {
            
            /* Get the cell points: */
            int p1 = i + j*(nx+1) + k*(nx+1)*(ny+1);
            int p2 = (i+1) + j*(nx+1) + k*(nx+1)*(ny+1);
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
            
            /* Get the cell volume: */
            double v = (nz > 0) ? dx[i] * dy[j] * dz[k] : dx[i] * dy[j];
            cells.volumes.push_back(v);
            
            /* Get the cell centroid: */
            double x0 = x[i] + 0.5*dx[i];
            double y0 = y[j] + 0.5*dy[j];
            double z0 = z[k] + 0.5*dz[k];
            cells.centroids.push_back(std::vector<double>{x0, y0, z0});
            
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
   int l = 0;
   for (int k = 0; k < num_cells/(nx*ny); k++) {
      for (int j = 0; j < ny; j++) {
         for (int i = 0; i < nx; i++) {
            
            /* Initialize the face data for this cell: */
            int num_faces = (nz > 0) ? 6 : 4;
            std::vector<std::vector<int>> pts(num_faces);
            std::vector<double> a(num_faces);
            std::vector<std::vector<double>> p0(num_faces);
            std::vector<std::vector<double>> n(num_faces);
            std::vector<int> l2(num_faces);
            
            /* -y face: */
            pts[0] = (nz > 0) ? math::extrude_edge(cells.points[l], 0, 4) : 
                        std::vector<int>{cells.points[l][0], cells.points[l][1]};
            a[0] = (nz > 0) ? dx[i] * dz[k] : dx[i];
            p0[0] = std::vector<double>{x[i]+0.5*dx[i], y[j], z[k]+0.5*dz[k]};
            n[0] = std::vector<double>{0.0, -1.0, 0.0};
            l2[0] = (j == 0) ? -3 : l - nx;
            
            /* +x face: */
            pts[1] = (nz > 0) ? math::extrude_edge(cells.points[l], 1, 4) : 
                        std::vector<int>{cells.points[l][1], cells.points[l][2]};
            a[1] = (nz > 0) ? dy[j] * dz[k] : dy[j];
            p0[1] = std::vector<double>{x[i]+dx[i], y[j]+0.5*dy[j], z[k]+0.5*dz[k]};
            n[1] = std::vector<double>{1.0, 0.0, 0.0};
            l2[1] = (i == nx-1) ? -2 : l + 1;
            
            /* +y face: */
            pts[2] = (nz > 0) ? math::extrude_edge(cells.points[l], 2, 4) : 
                        std::vector<int>{cells.points[l][2], cells.points[l][3]};
            a[2] = (nz > 0) ? dx[i] * dz[k] : dx[i];
            p0[2] = std::vector<double>{x[i]+0.5*dx[i], y[j]+dy[j], z[k]+0.5*dz[k]};
            n[2] = std::vector<double>{0.0, 1.0, 0.0};
            l2[2] = (j == ny-1) ? -4 : l + nx;
            
            /* -x face: */
            pts[3] = (nz > 0) ? math::extrude_edge(cells.points[l], 3, 4) : 
                        std::vector<int>{cells.points[l][3], cells.points[l][0]};
            a[3] = (nz > 0) ? dy[j] * dz[k] : dy[j];
            p0[3] = std::vector<double>{x[i], y[j]+0.5*dy[j], z[k]+0.5*dz[k]};
            n[3] = std::vector<double>{-1.0, 0.0, 0.0};
            l2[3] = (i == 0) ? -1 : l - 1;
            
            /* -z face: */
            if (nz > 0) {
               pts[4] = std::vector<int>(cells.points[l].rbegin()+4, cells.points[l].rend()+8);
               a[4] = dx[i] * dy[j];
               p0[4] = std::vector<double>{x[i]+0.5*dx[i], y[j]+0.5*dy[j], z[k]};
               n[4] = std::vector<double>{0.0, 0.0, -1.0};
               l2[4] = (k == 0) ? -5 : l - nx*ny;
            }
            
            /* +z face: */
            if (nz > 0) {
               pts[5] = std::vector<int>(cells.points[l].begin()+4, cells.points[l].end());
               a[5] = dx[i] * dy[j];
               p0[5] = std::vector<double>{x[i]+0.5*dx[i], y[j]+0.5*dy[j], z[k]+dz[k]};
               n[5] = std::vector<double>{0.0, 0.0, 1.0};
               l2[5] = (k == nz-1) ? -6 : l + nx*ny;
            }
            
            /* Keep the data for this cell: */
            faces.points.push_back(pts);
            faces.areas.push_back(a);
            faces.centroids.push_back(p0);
            faces.normals.push_back(n);
            faces.neighbours.push_back(l2);
            l++;
            
         }
      }
   }
   
   return 0;
   
};
