#include "CartesianMesh.hxx"

/* The CartesianMesh constructor: */
CartesianMesh::CartesianMesh() {};

/* The CartesianMesh destructor: */
CartesianMesh::~CartesianMesh() {};

/* Read the mesh from a plain-text input file: */
bool CartesianMesh::read(const std::string &filename) {
   
   /* Open the input file: */
   std::ifstream file(filename);
   if (!file.is_open()) {
      std::cout << "Error: unable to open " << filename << "!\n";
      return false;
   }
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line:*/
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty())
         break;
      
      /* Get the next keyword: */
      if (line[0] == "dx") {
         
         /* Get the dx values: */
         nx = std::stoi(line[1]);
         if (!utils::read(dx, nx, file)) {
            std::cout << "Error: wrong dx data in " << filename << "!\n";
            return false;
         }
         
      }
      else if (line[0] == "dy") {
         
         /* Get the dy values: */
         ny = std::stoi(line[1]);
         if (!utils::read(dy, ny, file)) {
            std::cout << "Error: wrong dy data in " << filename << "!\n";
            return false;
         }
         
      }
      else if (line[0] == "dz") {
         
         /* Get the dz values: */
         nz = std::stoi(line[1]);
         if (!utils::read(dz, nz, file)) {
            std::cout << "Error: wrong dz data in " << filename << "!\n";
            return false;
         }
         
      }
      else {
         
         /* Wrong keyword: */
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
            l2[0] = (j == 0) ? -1 : l - nx;
            
            /* +x face: */
            pts[1] = math::extrude_edge(cells.points[l], 1, 4);
            a[1] = dy[j] * dz[k];
            p0[1] = std::vector<double>{x[i]+dx[i], y[j]+0.5*dy[j], z[k]+0.5*dz[k]};
            n[1] = std::vector<double>{1.0, 0.0, 0.0};
            l2[1] = (i == nx-1) ? -1 : l + 1;
            
            /* +y face: */
            pts[2] = math::extrude_edge(cells.points[l], 2, 4);
            a[2] = dx[i] * dz[k];
            p0[2] = std::vector<double>{x[i]+0.5*dx[i], y[j]+dy[j], z[k]+0.5*dz[k]};
            n[2] = std::vector<double>{0.0, 1.0, 0.0};
            l2[2] = (j == ny-1) ? -1 : l + nx;
            
            /* -x face: */
            pts[3] = math::extrude_edge(cells.points[l], 3, 4);;
            a[3] = dy[j] * dz[k];
            p0[3] = std::vector<double>{x[i], y[j]+0.5*dy[j], z[k]+0.5*dz[k]};
            n[3] = std::vector<double>{-1.0, 0.0, 0.0};
            l2[3] = (i == 0) ? -1 : l - 1;
            
            /* -z face: */
            pts[4] = std::vector<int>(cells.points[l].rbegin()+4, cells.points[l].rend()+8);
            a[4] = dx[i] * dy[j];
            p0[4] = std::vector<double>{x[i]+0.5*dx[i], y[j]+0.5*dy[j], z[k]};
            n[4] = std::vector<double>{0.0, 0.0, -1.0};
            l2[4] = (k == 0) ? -1 : l - nx*ny;
            
            /* +z face: */
            pts[5] = std::vector<int>(cells.points[l].begin()+4, cells.points[l].end());
            a[5] = dx[i] * dy[j];
            p0[5] = std::vector<double>{x[i]+0.5*dx[i], y[j]+0.5*dy[j], z[k]+dz[k]};
            n[5] = std::vector<double>{0.0, 0.0, 1.0};
            l2[5] = (k == nz-1) ? -1 : l + nx*ny;
            
            /* Keep the data for this cell: */
            faces.points.push_back(pts);
            faces.areas.push_back(a);
            faces.centroids.push_back(p0);
            faces.normals.push_back(n);
            faces.neighbour.push_back(l2);
            l++;
            
         }
      }
   }
   
   return true;
   
};
