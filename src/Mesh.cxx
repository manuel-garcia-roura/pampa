#include "Mesh.hxx"

/* Partition the mesh: */
int Mesh::partition(Mesh** submesh) {
   
   /* Get the domain indices: */
   Array1D<int> domain_indices, num_cells_local;
   PAMPA_CALL(getDomainIndices(domain_indices, num_cells_local), 
      "unable to get the domain indices");
   
   /* Create the submesh: */
   *submesh = new Mesh();
   
   /* Get the number of dimensions: */
   (*submesh)->num_dims = num_dims;
   
   /* Get the number of local and global physical cells in the partitioned mesh: */
   (*submesh)->num_cells = num_cells_local(mpi::rank);
   (*submesh)->num_cells_global = num_cells;
   
   /* Get the cell mapping from each submesh in the partitioned mesh to the original mesh: */
   /* Note: ism_to_im(ip, ism) is the index in the original mesh for cell ism in partition ip. */
   Array1D<int> ism(mpi::size);
   Vector2D<int> ism_to_im(mpi::size, num_cells_local);
   for (int im = 0; im < num_cells; im++)
      ism_to_im(domain_indices(im), ism(domain_indices(im))++) = im;
   
   /* Get the number of ghost cells in this submesh and their mapping to the original mesh: */
   /* Note: ism_to_im_ghost(ism) is the index in the original mesh for ghost cell ism. */
   Array1D<int> ism_to_im_ghost;
   (*submesh)->num_ghost_cells = 0;
   for (int im = 0; im < num_cells; im++) {
      if (domain_indices(im) == mpi::rank) {
         for (int f = 0; f < faces.num_faces(im); f++) {
            int im2 = faces.neighbours(im, f);
            if (im2 >= 0) {
               if (domain_indices(im2) != mpi::rank && !(ism_to_im_ghost.find(im2))) {
                  ism_to_im_ghost.pushBack(im2);
                  (*submesh)->num_ghost_cells++;
               }
            }
         }
      }
   }
   
   /* Get the cell mapping from the original mesh to each submesh in the partitioned mesh: */
   /* Note: im_to_ism(im) is the index in its submesh for cell im in the original mesh: */
   Array1D<int> im_to_ism(num_cells);
   for (int ip = 0; ip < mpi::size; ip++) {
      for (int ism = 0; ism < num_cells_local(ip); ism++) {
         int im = ism_to_im(ip, ism);
         im_to_ism(im) = ism;
      }
   }
   
   /* Get the total number of physical and ghost cells in this submesh: */
   int num_cells_total = (*submesh)->num_cells + (*submesh)->num_ghost_cells;
   
   /* Get the cell mapping from this submesh to the original mesh for physical and ghost cells: */
   /* Note: ism_to_im_total(ism) is the index in the original mesh for cell ism. */
   Array1D<int> ism_to_im_total(num_cells_total);
   for (int ism = 0; ism < (*submesh)->num_cells; ism++)
      ism_to_im_total(ism) = ism_to_im(mpi::rank, ism);
   for (int ism = 0; ism < (*submesh)->num_ghost_cells; ism++)
      ism_to_im_total((*submesh)->num_cells+ism) = ism_to_im_ghost(ism);
   
   /* Get the cell mapping from the original mesh to this submesh for physical and ghost cells: */
   /* Note: im_to_ism_total(im) is the index in this submesh for cell ism in the original mesh. */
   Array1D<int> im_to_ism_total(num_cells, -1);
   for (int ism = 0; ism < num_cells_total; ism++) {
      int im = ism_to_im_total(ism);
      im_to_ism_total(im) = ism;
   }
   
   /* Get the mesh points that are actually used in the submesh: */
   Array1D<int> used(num_points);
   for (int ism = 0; ism < (*submesh)->num_cells; ism++) {
      int im = ism_to_im(mpi::rank, ism);
      for (int j = 0; j < cells.points.size(im); j++)
         used(cells.points(im, j)) = 1;
   }
   
   /* Get the number of points in the submesh: */
   (*submesh)->num_points = 0;
   for (int i = 0; i < num_points; i++)
      if (used(i))
         (*submesh)->num_points++;
   
   /* Build the submesh points and get the mapping from the mesh to the submesh: */
   Array1D<int> im_to_ism_points(num_points, -1);
   (*submesh)->points.resize((*submesh)->num_points, 3);
   for (int ism = 0, im = 0; im < num_points; im++) {
      if (used(im)) {
         (*submesh)->points(ism, 0) = points(im, 0);
         (*submesh)->points(ism, 1) = points(im, 1);
         (*submesh)->points(ism, 2) = points(im, 2);
         im_to_ism_points(im) = ism;
         ism++;
      }
   }
   
   /* Get the number of points for each physical cell in the submesh: */
   Array1D<int> num_cell_points((*submesh)->num_cells);
   for (int ism = 0; ism < (*submesh)->num_cells; ism++) {
      int im = ism_to_im_total(ism);
      num_cell_points(ism) = cells.points.size(im);
   }
   
   /* Build the submesh cells: */
   /* Note: for ghost cells only the centroids and materials are needed. */
   ((*submesh)->cells).points.resize((*submesh)->num_cells, num_cell_points);
   ((*submesh)->cells).volumes.resize((*submesh)->num_cells);
   ((*submesh)->cells).centroids.resize(num_cells_total, 3);
   ((*submesh)->cells).materials.resize(num_cells_total);
   for (int ism = 0; ism < num_cells_total; ism++) {
      int im = ism_to_im_total(ism);
      if (ism < (*submesh)->num_cells) {
         for (int j = 0; j < cells.points.size(im); j++)
            ((*submesh)->cells).points(ism, j) = im_to_ism_points(cells.points(im, j));
         ((*submesh)->cells).volumes(ism) = cells.volumes(im);
      }
      ((*submesh)->cells).centroids(ism, 0) = cells.centroids(im, 0);
      ((*submesh)->cells).centroids(ism, 1) = cells.centroids(im, 1);
      ((*submesh)->cells).centroids(ism, 2) = cells.centroids(im, 2);
      ((*submesh)->cells).materials(ism) = cells.materials(im);
   }
   
   /* Build the submesh cell indices in the global mesh for physical and ghost cells: */
   Array1D<int> ism0(mpi::size);
   for (int i = 1; i < mpi::size; i++)
      ism0(i) = ism0(i-1) + num_cells_local(i-1);
   ((*submesh)->cells).indices.resize(num_cells_total);
   for (int ism = 0; ism < (*submesh)->num_cells; ism++)
      ((*submesh)->cells).indices(ism) = ism0(mpi::rank) + ism;
   for (int ism = 0; ism < (*submesh)->num_ghost_cells; ism++) {
      int im = ism_to_im_ghost(ism);
      int ip = domain_indices(im);
      ((*submesh)->cells).indices((*submesh)->num_cells+ism) = ism0(ip) + im_to_ism(im);
   }
   
   /* Get the number of cell faces: */
   (*submesh)->num_faces_max = num_faces_max;
   ((*submesh)->faces).num_faces.resize((*submesh)->num_cells);
   for (int ism = 0; ism < (*submesh)->num_cells; ism++) {
      int im = ism_to_im_total(ism);
      ((*submesh)->faces).num_faces(ism) = faces.num_faces(im);
   }
   
   /* Build the submesh faces: */
   /* Note: only the faces for the physical cells are needed. */
   ((*submesh)->faces).areas.resize((*submesh)->num_cells, ((*submesh)->faces).num_faces);
   ((*submesh)->faces).centroids.resize((*submesh)->num_cells, ((*submesh)->faces).num_faces, 3);
   ((*submesh)->faces).normals.resize((*submesh)->num_cells, ((*submesh)->faces).num_faces, 3);
   ((*submesh)->faces).neighbours.resize((*submesh)->num_cells, ((*submesh)->faces).num_faces);
   for (int ism = 0; ism < (*submesh)->num_cells; ism++) {
      int im = ism_to_im_total(ism);
      for (int f = 0; f < faces.num_faces(im); f++) {
         ((*submesh)->faces).areas(ism, f) = faces.areas(im, f);
         ((*submesh)->faces).centroids(ism, f, 0) = faces.centroids(im, f, 0);
         ((*submesh)->faces).centroids(ism, f, 1) = faces.centroids(im, f, 1);
         ((*submesh)->faces).centroids(ism, f, 2) = faces.centroids(im, f, 2);
         ((*submesh)->faces).normals(ism, f, 0) = faces.normals(im, f, 0);
         ((*submesh)->faces).normals(ism, f, 1) = faces.normals(im, f, 1);
         ((*submesh)->faces).normals(ism, f, 2) = faces.normals(im, f, 2);
         if (faces.neighbours(im, f) >= 0)
            ((*submesh)->faces).neighbours(ism, f) = im_to_ism_total(faces.neighbours(im, f));
         else
            ((*submesh)->faces).neighbours(ism, f) = faces.neighbours(im, f);
      }
   }
   
   /* Copy the Boundary conditions: */
   (*submesh)->bcs = bcs;
   
   /* Write all the mesh data to a plain-text file: */
   PAMPA_CALL(writeData("mesh.pmp"), "unable to write the mesh data");
   
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
   
   /* Get the total number of cell points: */
   int num_cell_points = num_cells;
   for (int i = 0; i < num_cells; i++)
      num_cell_points += cells.points.size(i);
   
   /* Write the cell points: */
   file << "CELLS " << num_cells << " " << num_cell_points << std::endl;
   for (int i = 0; i < num_cells; i++) {
      num_cell_points = cells.points.size(i);
      file << num_cell_points;
      for (int j = 0; j < num_cell_points; j++)
         file << " " << cells.points(i, j);
      file << std::endl;
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
   
   /* Write to a rank directory in parallel runs: */
   std::string path;
   if (mpi::size > 1) {
      std::string dir = std::to_string(mpi::rank);
      PAMPA_CALL(utils::create(dir), "unable to create the output directory");
      path = dir + "/" + filename;
   }
   else
      path = filename;
   
   /* Open the output file: */
   std::ofstream file(path);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + path);
   
   /* Set the precision: */
   file << std::fixed;
   file << std::setprecision(3);
   
   /* Write the point coordinates: */
   file << "# points:" << std::endl;
   file << "points " << num_points << std::endl;
   for (int i = 0; i < num_points; i++) {
      file << points(i, 0) << " ";
      file << points(i, 1) << " ";
      file << points(i, 2) << std::endl;
   }
   file << std::endl;
   
   /* Write the number of cells: */
   file << "# number of cells:" << std::endl;
   file << "cells " << num_cells << " " << num_ghost_cells << " " << num_cells_global << std::endl;
   
   /* Get the total number of cell points: */
   int num_cell_points = 0;
   for (int i = 0; i < num_cells; i++)
      num_cell_points += cells.points.size(i);
   
   /* Write the cell points: */
   file << "# cell points:" << std::endl;
   file << "cell-points " << num_cells << " " << num_cell_points << std::endl;
   for (int i = 0; i < num_cells; i++) {
      for (int j = 0; j < cells.points.size(i); j++) {
         if (j > 0) file << " ";
         file << cells.points(i, j);
      }
      file << std::endl;
   }
   file << std::endl;
   
   /* Write the cell volumes: */
   file << "# cell volumes:" << std::endl;
   file << "cell-volumes " << num_cells << std::endl;
   for (int i = 0; i < num_cells; i++)
      file << cells.volumes(i) << std::endl;
   file << std::endl;
   
   /* Write the cell centroids: */
   file << "# cell centroids:" << std::endl;
   file << "cell-centroids " << num_cells+num_ghost_cells << std::endl;
   for (int i = 0; i < num_cells+num_ghost_cells; i++) {
      file << cells.centroids(i, 0) << " ";
      file << cells.centroids(i, 1) << " ";
      file << cells.centroids(i, 2) << std::endl;
   }
   file << std::endl;
   
   /* Write the cell materials (1-based indexed): */
   file << "# cell materials:" << std::endl;
   file << "cell-materials " << num_cells+num_ghost_cells << std::endl;
   for (int i = 0; i < num_cells+num_ghost_cells; i++)
      file << cells.materials(i)+1 << std::endl;
   file << std::endl;
   
   /* Write the cell indices in the global mesh: */
   file << "# cell indices:" << std::endl;
   file << "cell-indices " << num_cells+num_ghost_cells << std::endl;
   for (int i = 0; i < num_cells+num_ghost_cells; i++)
      file << cells.indices(i) << std::endl;
   file << std::endl;
   
   /* Write the number of cell faces: */
   file << "# number of cell faces:" << std::endl;
   file << "faces " << num_cells << std::endl;
   for (int i = 0; i < num_cells; i++)
      file << faces.num_faces(i) << std::endl;
   file << std::endl;
   
   /* Get the total number of cell faces: */
   int num_cell_faces = 0;
   for (int i = 0; i < num_cells; i++)
      num_cell_faces += faces.num_faces(i);
   
   /* Write the face areas: */
   file << "# face areas:" << std::endl;
   file << "face-areas " << num_cells << " " << num_cell_faces << std::endl;
   for (int i = 0; i < num_cells; i++) {
      for (int j = 0; j < faces.num_faces(i); j++) {
         if (j > 0) file << " ";
         file << faces.areas(i, j);
      }
      file << std::endl;
   }
   file << std::endl;
   
   /* Write the face centroids: */
   file << "# face centroids:" << std::endl;
   file << "face-centroids " << num_cells << " " << 3*num_cell_faces << std::endl;
   for (int i = 0; i < num_cells; i++) {
      for (int f = 0; f < faces.num_faces(i); f++) {
         file << faces.centroids(i, f, 0) << " ";
         file << faces.centroids(i, f, 1) << " ";
         file << faces.centroids(i, f, 2) << std::endl;
      }
   }
   file << std::endl;
   
   /* Write the face normals: */
   file << "# face normals:" << std::endl;
   file << "face-normals " << num_cells << " " << 3*num_cell_faces << std::endl;
   for (int i = 0; i < num_cells; i++) {
      for (int f = 0; f < faces.num_faces(i); f++) {
         file << faces.normals(i, f, 0) << " ";
         file << faces.normals(i, f, 1) << " ";
         file << faces.normals(i, f, 2) << std::endl;
      }
   }
   file << std::endl;
   
   /* Write the face neighbours: */
   file << "# face neighbours:" << std::endl;
   file << "face-neighbours " << num_cells << " " << num_cell_faces << std::endl;
   for (int i = 0; i < num_cells; i++) {
      for (int j = 0; j < faces.num_faces(i); j++) {
         if (j > 0) file << " ";
         file << faces.neighbours(i, j);
      }
      file << std::endl;
   }
   file << std::endl;
   
   /* Write the boundary conditions (1-based indexed): */
   file << "# boundary conditions:" << std::endl;
   for (int i = 1; i < bcs.size(); i++) {
      file << "bc " << i << " " << bcs(i).type+1;
      if (bcs(i).type == BC::ROBIN) file << " " << bcs(i).a;
      file << std::endl;
   }
   file << std::endl;
   
   return 0;
   
}

/* Get the domain indices to partition the mesh: */
#ifdef WITH_METIS
int Mesh::getDomainIndices(Array1D<int>& part, Array1D<int>& size) {
   
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
   
   /* Get the number of cells for each domain: */
   size.resize(mpi::size);
   for (int i = 0; i < num_cells; i++)
      size(part(i))++;
   
   return 0;
   
}
#else
int Mesh::getDomainIndices(Array1D<int>& part, Array1D<int>& size) {
   
   /* Get the number of cells for each domain: */
   size.resize(mpi::size, num_cells/mpi::size);
   for (int i = 0; i < mpi::size; i++)
      if (i < num_cells%mpi::size)
         size(i)++;
   
   /* Get the domain indices: */
   part.resize(num_cells);
   for (int ic = 0, i = 0; i < mpi::size; i++)
      for (int j = 0; j < size(i); j++)
         part(ic++) = i;
   
   return 0;
   
}
#endif
