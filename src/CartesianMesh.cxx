#include "CartesianMesh.hxx"

/* The CartesianMesh constructor: */
CartesianMesh::CartesianMesh() {};

/* The CartesianMesh destructor: */
CartesianMesh::~CartesianMesh() {};

/* Read the mesh from a plain-text file: */
bool CartesianMesh::read(const std::string &filename) {
   
   /* Open the input file: */
   std::ifstream file(filename);
   if (!file.is_open()) {
      std::cout << "Error: unable to open " << filename << "!\n";
      return false;
   }
   
   /* Read the file line by line: */
   std::string line;
   while (std::getline(file, line)) {
      
      /* Skip empty lines and #-marked comments: */
      if (line.empty() || line.at(0) == '#')
         continue;
      
      /* Check for dx, dy and dz: */
      std::istringstream iss(line);
      std::string s;
      std::getline(iss, s, ' ');
      if (s == "dx") {
         
         /* Get the dx values: */
         std::getline(iss, s, ' ');
         nx = std::stoi(s);
         if (!utils::read(dx, nx, file)) {
            std::cout << "Error: wrong dx data in " << filename << "!\n";
            return false;
         }
         
      }
      else if (s == "dy") {
         
         /* Get the dy values: */
         std::getline(iss, s, ' ');
         ny = std::stoi(s);
         if (!utils::read(dy, ny, file)) {
            std::cout << "Error: wrong dy data in " << filename << "!\n";
            return false;
         }
         
      }
      else if (s == "dz") {
         
         /* Get the dz values: */
         std::getline(iss, s, ' ');
         nz = std::stoi(s);
         if (!utils::read(dz, nz, file)) {
            std::cout << "Error: wrong dz data in " << filename << "!\n";
            return false;
         }
         
      }
      else {
         std::cout << "Error: wrong keyword in " << filename << "!\n";
         return false;
      }
      
   }
   
   return true;
   
};

/* Build the mesh: */
bool CartesianMesh::build() {
   
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
   num_cells = nx * ny * nz;
   cells.points.reserve(num_cells);
   cells.volumes.reserve(num_cells);
   cells.centroids.reserve(num_cells);
   for (int k = 0; k < nz; k++) {
      for (int j = 0; j < ny; j++) {
         for (int i = 0; i < nx; i++) {
            int p1 = i + j*(nx+1) + k*(nx+1)*(ny+1);
            int p2 = (i+1) + j*(nx+1) + k*(nx+1)*(ny+1);
            int p3 = (i+1) + (j+1)*(nx+1) + k*(nx+1)*(ny+1);
            int p4 = i + (j+1)*(nx+1) + k*(nx+1)*(ny+1);
            int p5 = i + j*(nx+1) + (k+1)*(nx+1)*(ny+1);
            int p6 = (i+1) + j*(nx+1) + (k+1)*(nx+1)*(ny+1);
            int p7 = (i+1) + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
            int p8 = i + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
            double v = dx[i] * dy[j] * dz[k];
            double x0 = x[i] + 0.5*dx[i];
            double y0 = y[j] + 0.5*dy[j];
            double z0 = z[k] + 0.5*dz[k];
            cells.points.push_back(std::vector<int>{p1, p2, p3, p4, p5, p6, p7, p8});
            cells.volumes.push_back(v);
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
   faces.neighbour.reserve(num_cells);
   int l = 0;
   for (int k = 0; k < nz; k++) {
      for (int j = 0; j < ny; j++) {
         for (int i = 0; i < nx; i++) {
            
            /* Initialize the face data for this cell: */
            std::vector<std::vector<int>> pts(6);
            std::vector<double> a(6);
            std::vector<std::vector<double>> p0(6);
            std::vector<std::vector<double>> n(6);
            std::vector<int> l2(6);
            
            /* -y face: */
            pts[0] = math::extrude_edge(cells.points[l], 0, 4);
            a[0] = dx[i] * dz[k];
            p0[0] = std::vector<double>{x[i]+0.5*dx[i], y[j], z[k]+0.5*dz[k]};
            n[0] = std::vector<double>{0.0, -1.0, 0.0};
            l2[0] = (j == 0) ? -1 : j-1;
            
            /* +x face: */
            pts[1] = math::extrude_edge(cells.points[l], 1, 4);
            a[1] = dy[j] * dz[k];
            p0[1] = std::vector<double>{x[i]+dx[i], y[j]+0.5*dy[j], z[k]+0.5*dz[k]};
            n[1] = std::vector<double>{1.0, 0.0, 0.0};
            l2[1] = (i == nx-1) ? -1 : i+1;
            
            /* +y face: */
            pts[2] = math::extrude_edge(cells.points[l], 2, 4);
            a[2] = dx[i] * dz[k];
            p0[2] = std::vector<double>{x[i]+0.5*dx[i], y[j]+dy[j], z[k]+0.5*dz[k]};
            n[2] = std::vector<double>{0.0, 1.0, 0.0};
            l2[2] = (j == ny-1) ? -1 : j+1;
            
            /* -x face: */
            pts[3] = math::extrude_edge(cells.points[l], 3, 4);;
            a[3] = dy[j] * dz[k];
            p0[3] = std::vector<double>{x[i], y[j]+0.5*dy[j], z[k]+0.5*dz[k]};
            n[3] = std::vector<double>{-1.0, 0.0, 0.0};
            l2[3] = (i == 0) ? -1 : i-1;
            
            /* -z face: */
            pts[4] = std::vector<int>(cells.points[l].rbegin()+4, cells.points[l].rend()+8);
            a[4] = dx[i] * dy[j];
            p0[4] = std::vector<double>{x[i]+0.5*dx[i], y[j]+0.5*dy[j], z[k]};
            n[4] = std::vector<double>{0.0, 0.0, -1.0};
            l2[4] = (k == 0) ? -1 : k-1;
            
            /* +z face: */
            pts[5] = std::vector<int>(cells.points[l].begin()+4, cells.points[l].end());
            a[5] = dx[i] * dy[j];
            p0[5] = std::vector<double>{x[i]+0.5*dx[i], y[j]+0.5*dy[j], z[k]+dz[k]};
            n[5] = std::vector<double>{0.0, 0.0, 1.0};
            l2[5] = (k == nz-1) ? -1 : k+1;
            
            /* Keep the data for this cell: */
            faces.points.push_back(pts);
            faces.areas.push_back(a);
            faces.centroids.push_back(p0);
            faces.normals.push_back(n);
            faces.neighbour.push_back(l2);
            l++;
            
            if (l == num_cells) {
               for (int h = 0; h < 6; h++) { 
                  std::cout << pts[h][0] << " " << pts[h][1] << " " << pts[h][2] << " " << pts[h][3] << std::endl;
                  std::cout << a[h] << std::endl;
                  std::cout << p0[h][0] << " " << p0[h][1] << " " << p0[h][2] << std::endl;
                  std::cout << n[h][0] << " " << n[h][1] << " " << n[h][2] << std::endl;
                  std::cout << l2[h] << std::endl;
               }
            }
            
         }
      }
   }
   
   return true;
   
};
