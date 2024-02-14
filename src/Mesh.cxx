#include "Mesh.hxx"

/* Partition the mesh: */
int Mesh::partition(Mesh** submesh) {
   
   /* Get the domain indices: */
   Array1D<int> part;
   PAMPA_CALL(getDomainIndices(part), "unable to get the domain indices");
   
   /* Create the submesh: */
   *submesh = new Mesh();
   
   /* Get the number of physical and ghost cells and their cell indices in the global mesh: */
   Array1D<int> mesh_cell_indices, mesh_ghost_cell_indices;
   (*submesh)->num_cells = 0;
   (*submesh)->num_ghost_cells = 0;
   for (int im = 0; im < num_cells; im++) {
      if (part(im) == mpi::rank) {
         mesh_cell_indices.pushBack(im);
         (*submesh)->num_cells++;
         for (int f = 0; f < faces.num_faces(im); f++) {
            int j = faces.neighbours(im, f);
            if (j >= 0) {
               if (part(j) != mpi::rank && !(mesh_ghost_cell_indices.find(j))) {
                  mesh_ghost_cell_indices.pushBack(j);
                  (*submesh)->num_ghost_cells++;
               }
            }
         }
      }
   }
   
   /* Add the global indices for ghost cells to the ones for physical cells: */
   int num_cells_total = (*submesh)->num_cells + (*submesh)->num_ghost_cells;
   mesh_cell_indices.resize(num_cells_total);
   for (int ism = (*submesh)->num_cells, i = 0; i < (*submesh)->num_ghost_cells; i++)
      mesh_cell_indices(ism++) = mesh_ghost_cell_indices(i);
   
   /* Get mesh points that are actually used in the submesh: */
   Array1D<int> used(num_points, 0);
   for (int ism = 0; ism < num_cells_total; ism++) {
      int im = mesh_cell_indices(ism);
      for (int j = 0; j < cells.points.size(im); j++)
         used(cells.points(im, j)) = 1;
   }
   
   /* Get the number of points in the submesh: */
   (*submesh)->num_points = 0;
   for (int i = 0; i < num_points; i++)
      if (used(i))
         (*submesh)->num_points++;
   
   /* Build the submesh points and get the mapping from the mesh to the submesh: */
   Array1D<int> submesh_point_indices(num_points, -1);
   (*submesh)->points.resize((*submesh)->num_points, 3);
   for (int ism = 0, im = 0; im < num_points; im++) {
      if (used(im)) {
         (*submesh)->points(ism, 0) = points(im, 0);
         (*submesh)->points(ism, 1) = points(im, 1);
         (*submesh)->points(ism, 2) = points(im, 2);
         submesh_point_indices(im) = ism;
         ism++;
      }
   }
   
   /* Get the number of points for each cell in the submesh: */
   Array1D<int> num_cell_points(num_cells_total);
   for (int ism = 0; ism < num_cells_total; ism++) {
      int im = mesh_cell_indices(ism);
      num_cell_points(ism) = cells.points.size(im);
   }
   
   /* Build the submesh cells and get the mapping from the mesh to the submesh: */
   Array1D<int> submesh_cell_indices(num_cells, -1);
   ((*submesh)->cells).points.resize(num_cells_total, num_cell_points);
   ((*submesh)->cells).volumes.resize(num_cells_total);
   ((*submesh)->cells).centroids.resize(num_cells_total, 3);
   ((*submesh)->cells).materials.resize(num_cells_total);
   for (int ism = 0; ism < num_cells_total; ism++) {
      int im = mesh_cell_indices(ism);
      for (int j = 0; j < cells.points.size(im); j++)
         ((*submesh)->cells).points(ism, j) = submesh_point_indices(cells.points(im, j));
      ((*submesh)->cells).volumes(ism) = cells.volumes(im);
      ((*submesh)->cells).centroids(ism, 0) = cells.centroids(im, 0);
      ((*submesh)->cells).centroids(ism, 1) = cells.centroids(im, 1);
      ((*submesh)->cells).centroids(ism, 2) = cells.centroids(im, 2);
      ((*submesh)->cells).materials(ism) = cells.materials(im);
      submesh_cell_indices(im) = ism;
   }
   
   /* Get the maximum number of faces: */
   (*submesh)->num_faces_max = 0;
   for (int ism = 0; ism < (*submesh)->num_cells; ism++) {
      int im = mesh_cell_indices(ism);
      (*submesh)->num_faces_max = std::max((*submesh)->num_faces_max, faces.num_faces(im));
   }
   
   /* Build the submesh faces: */
   /* Note: only the faces for the physical cells are needed. */
   ((*submesh)->faces).num_faces.resize((*submesh)->num_cells);
   ((*submesh)->faces).areas.resize((*submesh)->num_cells, (*submesh)->num_faces_max);
   ((*submesh)->faces).centroids.resize((*submesh)->num_cells, (*submesh)->num_faces_max, 3);
   ((*submesh)->faces).normals.resize((*submesh)->num_cells, (*submesh)->num_faces_max, 3);
   ((*submesh)->faces).neighbours.resize((*submesh)->num_cells, (*submesh)->num_faces_max);
   for (int ism = 0; ism < (*submesh)->num_cells; ism++) {
      int im = mesh_cell_indices(ism);
      ((*submesh)->faces).num_faces(ism) = faces.num_faces(im);
      for (int f = 0; f < faces.num_faces(im); f++) {
         ((*submesh)->faces).areas(ism, f) = faces.areas(im, f);
         ((*submesh)->faces).centroids(ism, f, 0) = faces.centroids(im, f, 0);
         ((*submesh)->faces).centroids(ism, f, 1) = faces.centroids(im, f, 1);
         ((*submesh)->faces).centroids(ism, f, 2) = faces.centroids(im, f, 2);
         ((*submesh)->faces).normals(ism, f, 0) = faces.normals(im, f, 0);
         ((*submesh)->faces).normals(ism, f, 1) = faces.normals(im, f, 1);
         ((*submesh)->faces).normals(ism, f, 2) = faces.normals(im, f, 2);
         ((*submesh)->faces).neighbours(ism, f) = submesh_cell_indices(faces.neighbours(im, f));
      }
   }
   
   /* Copy the Boundary conditions: */
   (*submesh)->bcs = bcs;
   
   std::string filename = "mesh" + std::to_string(mpi::rank) + ".vtk";
   (*submesh)->writeVTK(filename);
   
   return 0;
   
}

/* Write the mesh to a plain-text file in .vtk format: */
int Mesh::writeVTK(const std::string& filename) const {
   
   /* Open the output file: */
   std::ofstream file(filename);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Set the precision: */
   file << std::fixed;
   file << std::setprecision(3);
   
   /* Write the header: */
   file << "# vtk DataFile Version 3.0" << std::endl;
   file << "FVM mesh" << std::endl;
   file << "ASCII" << std::endl;
   file << "DATASET UNSTRUCTURED_GRID" << std::endl;
   file << std::endl;
   
   /* Write the point coordinates: */
   file << "POINTS " << num_points << " double" << std::endl;
   for (int i = 0; i < num_points; i++) {
      file << points(i, 0) << " ";
      file << points(i, 1) << " ";
      file << points(i, 2) << std::endl;
   }
   file << std::endl;
   
   /* Write the cell points: */
   int num_cell_points = num_cells;
   for (int i = 0; i < num_cells; i++)
      num_cell_points += cells.points.size(i);
   file << "CELLS " << num_cells << " " << num_cell_points << std::endl;
   for (int i = 0; i < num_cells; i++) {
      num_cell_points = cells.points.size(i);
      file << num_cell_points;
      for (int j = 0; j < num_cell_points; j++) {
         file << " " << cells.points(i, j);
         if (j == num_cell_points-1) file << std::endl;
      }
   }
   file << std::endl;
   
   /* Write the cell types: */
   file << "CELL_TYPES " << num_cells << std::endl;
   for (int i = 0; i < num_cells; i++) {
      num_cell_points = cells.points.size(i);
      if (num_cell_points == 2)
         file << "3" << std::endl;
      else if (num_cell_points == 4)
         file << "9" << std::endl;
      else if (num_cell_points == 8)
         file << "12" << std::endl;
      else
         PAMPA_CHECK(true, 1, "wrong cell type");
   }
   file << std::endl;
   
   /* Write the cell data: */
   file << "CELL_DATA " << num_cells << std::endl << std::endl;
   
   /* Write the cell materials: */
   file << "SCALARS materials double 1" << std::endl;
   file << "LOOKUP_TABLE default" << std::endl;
   for (int i = 0; i < num_cells; i++)
      file << cells.materials(i)+1 << std::endl;
   file << std::endl;
   
   return 0;
   
}

/* Write all the mesh data to a plain-text file: */
int Mesh::writeData(const std::string& filename) const {
   
   /* Open the output file: */
   std::ofstream file(filename);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Set the precision: */
   file << std::fixed;
   file << std::setprecision(3);
   
   /* Write the point coordinates: */
   file << "POINTS " << num_points << " double" << std::endl;
   for (int i = 0; i < num_points; i++) {
      file << points(i, 0) << " ";
      file << points(i, 1) << " ";
      file << points(i, 2) << std::endl;
   }
   file << std::endl;
   
   /* Write the cell points: */
   int num_cell_points = num_cells;
   for (int i = 0; i < num_cells; i++)
      num_cell_points += cells.points.size(i);
   file << "CELLS " << num_cells << " " << num_cell_points << std::endl;
   for (int i = 0; i < num_cells; i++) {
      num_cell_points = cells.points.size(i);
      file << num_cell_points;
      for (int j = 0; j < num_cell_points; j++) {
         file << " " << cells.points(i, j);
         if (j == num_cell_points-1) file << std::endl;
      }
   }
   file << std::endl;
   
   /* Write the cell volumes: */
   file << "CELL_VOLUMES " << num_cells << std::endl;
   for (int i = 0; i < num_cells; i++)
      file << cells.volumes(i) << std::endl;
   file << std::endl;
   
   /* Write the cell centroids: */
   file << "CELL_CENTROIDS " << num_cells << std::endl;
   for (int i = 0; i < num_cells; i++) {
      file << cells.centroids(i, 0) << " ";
      file << cells.centroids(i, 1) << " ";
      file << cells.centroids(i, 2) << std::endl;
   }
   file << std::endl;
   
   /* Write the cell materials: */
   file << "CELL_MATERIALS " << num_cells << std::endl;
   for (int i = 0; i < num_cells; i++)
      file << cells.materials(i)+1 << std::endl;
   file << std::endl;
   
   /* Write the face areas: */
   int num_faces = 0;
   for (int i = 0; i < num_cells; i++)
      num_faces += faces.num_faces(i);
   file << "FACE_AREAS " << num_cells << " " << num_faces << std::endl;
   for (int i = 0; i < num_cells; i++)
      for (int f = 0; f < faces.num_faces(i); f++)
         file << i << " " << f << " " << faces.areas(i, f) << std::endl;
   file << std::endl;
   
   /* Write the face centroids: */
   file << "FACE_CENTROIDS " << num_cells << " " << num_faces << std::endl;
   for (int i = 0; i < num_cells; i++) {
      for (int f = 0; f < faces.num_faces(i); f++) {
         file << i << " " << f << " ";
         file << faces.centroids(i, f, 0) << " ";
         file << faces.centroids(i, f, 1) << " ";
         file << faces.centroids(i, f, 2) << std::endl;
      }
   }
   file << std::endl;
   
   /* Write the face normals: */
   file << "FACE_NORMALS " << num_cells << " " << num_faces << std::endl;
   for (int i = 0; i < num_cells; i++) {
      for (int f = 0; f < faces.num_faces(i); f++) {
         file << i << " " << f << " ";
         file << faces.normals(i, f, 0) << " ";
         file << faces.normals(i, f, 1) << " ";
         file << faces.normals(i, f, 2) << std::endl;
      }
   }
   file << std::endl;
   
   /* Write the face neighbours: */
   file << "FACE_NEIGHBOURS " << num_cells << " " << num_faces << std::endl;
   for (int i = 0; i < num_cells; i++) {
      for (int f = 0; f < faces.num_faces(i); f++) {
         int i2 = faces.neighbours(i, f);
         if (i2 < 0)
            file << i << " " << f << " " << bcs(-i2).type << std::endl;
         else
            file << i << " " << f << " " << i2 << std::endl;
      }
   }
   file << std::endl;
   
   return 0;
   
}

/* Get the domain indices to partition the mesh: */
#ifdef WITH_METIS
int Mesh::getDomainIndices(Array1D<int>& part) {
   
   /* Check the number of MPI ranks: */
   if (mpi::size == 1)
      return 0;
   
   /* Get the number of graph vertices, i.e. cells, and the total number of edges, i.e. faces: */
   /* Note: only faces with physical neighbouring cells need to be considered here. */
   idx_t nvtxs = num_cells, ncon = 0;
   for (int i = 0; i < num_cells; i++)
      for (int f = 0; f < faces.num_faces(i); f++)
         if (faces.neighbours(i, f) >= 0)
            ncon++;
   
   /* Get the adjacency structure of the graph in CSR format: */
   Array1D<idx_t> xadj(nvtxs+1), adjncy(ncon);
   for (int icon = 0, i = 0; i < num_cells; i++) {
      xadj(i) = icon;
      for (int f = 0; f < faces.num_faces(i); f++) {
         if (faces.neighbours(i, f) >= 0) {
            adjncy(icon) = faces.neighbours(i, f);
            icon++;
         }
      }
   }
   xadj(nvtxs) = ncon;
   
   /* Set the METIS options: */
   Array1D<idx_t> options(METIS_NOPTIONS);
   METIS_SetDefaultOptions(&options(0));
   options(METIS_OPTION_CTYPE) = METIS_CTYPE_RM;
   options(METIS_OPTION_IPTYPE) = METIS_IPTYPE_GROW;
   options(METIS_OPTION_RTYPE) = METIS_RTYPE_FM;
   options(METIS_OPTION_NO2HOP) = 0;
   options(METIS_OPTION_DBGLVL) = 0;
   #if METIS_PARTGRAPH == METIS_KWAY
   options(METIS_OPTION_OBJTYPE) = METIS_OBJTYPE_CUT;
   options(METIS_OPTION_MINCONN) = 1;
   #endif
   
   /* Initialize the graph partition vector: */
   idx_t nparts = mpi::size, objval;
   part.resize(nvtxs);
   
   /* Perform the partitioning: */
   #if METIS_PARTGRAPH == METIS_RECURSIVE
   METIS_CALL(METIS_PartGraphRecursive(&nvtxs, &ncon, &xadj(0), &adjncy(0), NULL, NULL, NULL, 
      &nparts, NULL, NULL, &options(0), &objval, &part(0)));
   #elif METIS_PARTGRAPH == METIS_KWAY
   METIS_CALL(METIS_PartGraphKway(&nvtxs, &ncon, &xadj(0), &adjncy(0), NULL, NULL, NULL, 
      &nparts, NULL, NULL, &options(0), &objval, &part(0)));
   #else
      #error "Wrong METIS graph partitioning routine."
   #endif
   
   return 0;
   
}
#else
int Mesh::getDomainIndices(Array1D<int>& part) {
   
   /* Check the number of MPI ranks: */
   if (mpi::size == 1)
      return 0;
   
   /* Get the number of domains and the number of cells for each domain: */
   int num_domains = mpi::size;
   Array1D<int> num_domain_cells(num_domains, num_cells/num_domains);
   for (int i = 0; i < num_domains; i++)
      if (i < num_cells%num_domains)
         num_domain_cells(i)++;
   
   /* Get the domain indices: */
   part.resize(num_cells);
   for (int ic = 0, i = 0; i < num_domains; i++)
      for (int j = 0; j < num_domain_cells(i); j++)
         part(ic++) = i;
   
   return 0;
   
}
#endif
