#include "Mesh.hxx"

#ifdef WITH_METIS
#include <metis.h>
#define METIS_RECURSIVE 1
#define METIS_KWAY 2
#define METIS_PARTGRAPH METIS_KWAY
#endif

/* Remove boundary-condition materials from the mesh: */
int Mesh::removeBCMats(const Array1D<Material*>& materials, Mesh** mesh) {
   
   /* Check if the mesh has already been partitioned: */
   PAMPA_CHECK(partitioned, "unable to remove materials from a partitioned mesh");
   
   /* Create the new mesh: */
   *mesh = new Mesh();
   
   /* Get the number of dimensions: */
   (*mesh)->num_dims = num_dims;
   
   /* Get the number of cells after removing the materials: */
   (*mesh)->num_cells = 0;
   for (int i = 0; i < num_cells; i++)
      if (!(materials(cells.materials(i))->isBC()))
         (*mesh)->num_cells++;
   (*mesh)->num_cells_global = (*mesh)->num_cells;
   
   /* Get the number of faces after removing the materials: */
   (*mesh)->num_faces_max = 0;
   ((*mesh)->faces).num_faces.resize((*mesh)->num_cells);
   for (int ic = 0, i = 0; i < num_cells; i++) {
      if (!(materials(cells.materials(i))->isBC())) {
         (*mesh)->num_faces_max = std::max((*mesh)->num_faces_max, faces.num_faces(i));
         ((*mesh)->faces).num_faces(ic++) = faces.num_faces(i);
      }
   }
   
   /* Get the mesh points that are actually used after removing the materials: */
   Array1D<int> used(num_points);
   for (int i = 0; i < num_cells; i++)
      if (!(materials(cells.materials(i))->isBC()))
         for (int j = 0; j < cells.points.size(i); j++)
            used(cells.points(i, j)) = 1;
   
   /* Get the number of points after removing the materials: */
   (*mesh)->num_points = 0;
   for (int i = 0; i < num_points; i++)
      if (used(i))
         (*mesh)->num_points++;
   
   /* Build the new mesh points and get the mapping from the original to the new mesh: */
   Array1D<int> i_to_ip(num_points, -1);
   (*mesh)->points.resize((*mesh)->num_points, 3);
   for (int ip = 0, i = 0; i < num_points; i++) {
      if (used(i)) {
         (*mesh)->points(ip, 0) = points(i, 0);
         (*mesh)->points(ip, 1) = points(i, 1);
         (*mesh)->points(ip, 2) = points(i, 2);
         i_to_ip(i) = ip;
         ip++;
      }
   }
   
   /* Get the number of points for each cell after removing the materials: */
   Array1D<int> num_cell_points((*mesh)->num_cells);
   for (int ic = 0, i = 0; i < num_cells; i++)
      if (!(materials(cells.materials(i))->isBC()))
         num_cell_points(ic++) = cells.points.size(i);
   
   /* Copy the boundaries: */
   (*mesh)->boundaries = boundaries;
   
   /* Create the new boundaries: */
   Array1D<int> bcmat_indices(materials.size(), -1);
   for (int i = 0; i < materials.size(); i++) {
      if (materials(i)->isBC()) {
         bcmat_indices(i) = ((*mesh)->boundaries).find(materials(i)->name);
         PAMPA_CHECK(bcmat_indices(i) < 0, "wrong boundary name");
      }
   }
   
   /* Get the cell mapping from the original mesh to the new mesh: */
   Array1D<int> i_to_ic(num_cells, -1);
   for (int ic = 0, i = 0; i < num_cells; i++) {
      if (!(materials(cells.materials(i))->isBC())) {
         i_to_ic(i) = ic;
         ic++;
      }
   }
   
   /* Get all the cells with physical materials: */
   ((*mesh)->cells).points.resize((*mesh)->num_cells, num_cell_points);
   ((*mesh)->cells).volumes.resize((*mesh)->num_cells);
   ((*mesh)->cells).centroids.resize((*mesh)->num_cells, 3);
   ((*mesh)->cells).materials.resize((*mesh)->num_cells);
   if (!(cells.nodal_indices.empty()))
      ((*mesh)->cells).nodal_indices.resize((*mesh)->num_cells);
   ((*mesh)->cells).global_indices.resize((*mesh)->num_cells);
   ((*mesh)->faces).areas.resize((*mesh)->num_cells, ((*mesh)->faces).num_faces);
   ((*mesh)->faces).centroids.resize((*mesh)->num_cells, ((*mesh)->faces).num_faces, 3);
   ((*mesh)->faces).normals.resize((*mesh)->num_cells, ((*mesh)->faces).num_faces, 3);
   ((*mesh)->faces).neighbours.resize((*mesh)->num_cells, ((*mesh)->faces).num_faces);
   for (int ic = 0, i = 0; i < num_cells; i++) {
      if (!(materials(cells.materials(i))->isBC())) {
         
         /* Get the cell data: */
         for (int j = 0; j < cells.points.size(i); j++)
            ((*mesh)->cells).points(ic, j) = i_to_ip(cells.points(i, j));
         ((*mesh)->cells).volumes(ic) = cells.volumes(i);
         ((*mesh)->cells).centroids(ic, 0) = cells.centroids(i, 0);
         ((*mesh)->cells).centroids(ic, 1) = cells.centroids(i, 1);
         ((*mesh)->cells).centroids(ic, 2) = cells.centroids(i, 2);
         ((*mesh)->cells).materials(ic) = cells.materials(i);
         if (!(((*mesh)->cells).nodal_indices.empty()))
            ((*mesh)->cells).nodal_indices(ic) = cells.nodal_indices(i);
         ((*mesh)->cells).global_indices(ic) = ic;
         
         /* Get the face data and replace boundary-condition materials with boundaries: */
         for (int f = 0; f < faces.num_faces(i); f++) {
            ((*mesh)->faces).areas(ic, f) = faces.areas(i, f);
            ((*mesh)->faces).centroids(ic, f, 0) = faces.centroids(i, f, 0);
            ((*mesh)->faces).centroids(ic, f, 1) = faces.centroids(i, f, 1);
            ((*mesh)->faces).centroids(ic, f, 2) = faces.centroids(i, f, 2);
            ((*mesh)->faces).normals(ic, f, 0) = faces.normals(i, f, 0);
            ((*mesh)->faces).normals(ic, f, 1) = faces.normals(i, f, 1);
            ((*mesh)->faces).normals(ic, f, 2) = faces.normals(i, f, 2);
            int i2 = faces.neighbours(i, f);
            if (i2 >= 0) {
               if (materials(cells.materials(i2))->isBC())
                  ((*mesh)->faces).neighbours(ic, f) = -bcmat_indices(cells.materials(i2)) - 1;
               else
                  ((*mesh)->faces).neighbours(ic, f) = i_to_ic(i2);
            }
            else
               ((*mesh)->faces).neighbours(ic, f) = i2;
         }
         
         ic++;
         
      }
   }
   
   /* Copy the boundary conditions: */
   (*mesh)->bcs = bcs;
   
   return 0;
   
}

/* Partition the mesh: */
int Mesh::partition(Mesh** submesh) {
   
   /* Get the domain indices: */
   Array1D<int> domain_indices, num_cells_local;
   PAMPA_CHECK(getDomainIndices(domain_indices, num_cells_local), 
      "unable to get the domain indices");
   
   /* Create the submesh: */
   *submesh = new Mesh();
   
   /* Mark the submesh as partitioned: */
   (*submesh)->partitioned = true;
   
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
            if (im2 >= 0 && domain_indices(im2) != mpi::rank && ism_to_im_ghost.find(im2) < 0) {
               ism_to_im_ghost.pushBack(im2);
               (*submesh)->num_ghost_cells++;
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
   if (!(cells.nodal_indices.empty()))
      ((*submesh)->cells).nodal_indices.resize((*submesh)->num_cells);
   for (int ism = 0; ism < num_cells_total; ism++) {
      int im = ism_to_im_total(ism);
      if (ism < (*submesh)->num_cells) {
         for (int j = 0; j < cells.points.size(im); j++)
            ((*submesh)->cells).points(ism, j) = im_to_ism_points(cells.points(im, j));
         ((*submesh)->cells).volumes(ism) = cells.volumes(im);
         if (!(cells.nodal_indices.empty()))
            ((*submesh)->cells).nodal_indices(ism) = cells.nodal_indices(im);
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
   ((*submesh)->cells).global_indices.resize(num_cells_total);
   for (int ism = 0; ism < (*submesh)->num_cells; ism++)
      ((*submesh)->cells).global_indices(ism) = ism0(mpi::rank) + ism;
   for (int ism = 0; ism < (*submesh)->num_ghost_cells; ism++) {
      int im = ism_to_im_ghost(ism);
      int ip = domain_indices(im);
      ((*submesh)->cells).global_indices((*submesh)->num_cells+ism) = ism0(ip) + im_to_ism(im);
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
   
   /* Copy the boundaries: */
   (*submesh)->boundaries = boundaries;
   
   /* Copy the boundary conditions: */
   (*submesh)->bcs = bcs;
   
   /* Write all the submesh data to a plain-text file: */
   PAMPA_CHECK((*submesh)->writeData("mesh.pmp"), "unable to write the submesh data");
   
   return 0;
   
}

/* Write the mesh to a plain-text file in .vtk format: */
int Mesh::writeVTK(const std::string& filename) const {
   
   /* Write the mesh in .vtk format: */
   PAMPA_CHECK(vtk::write(filename, points, num_points, cells.points, num_cells, cells.materials, 
      cells.nodal_indices), "unable to write the mesh");
   
   return 0;
   
}

/* Write all the mesh data to a plain-text file: */
int Mesh::writeData(const std::string& filename) const {
   
   /* Open the output file: */
   std::string path = mpi::get_path(filename);
   std::ofstream file(path, std::ios_base::out);
   PAMPA_CHECK(!file.is_open(), "unable to open " + path);
   
   /* Set the precision: */
   file << std::fixed;
   file << std::setprecision(3);
   
   /* Write the point coordinates: */
   file << "points " << num_points << std::endl;
   for (int i = 0; i < num_points; i++) {
      file << points(i, 0) << " ";
      file << points(i, 1) << " ";
      file << points(i, 2) << std::endl;
   }
   file << std::endl;
   
   /* Write the number of cells: */
   file << "cells " << num_cells << " " << num_ghost_cells << " " << num_cells_global << std::endl;
   
   /* Get the total number of cell points: */
   int num_cell_points = 0;
   for (int i = 0; i < num_cells; i++)
      num_cell_points += cells.points.size(i);
   
   /* Write the cell points: */
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
   file << "cell-volumes " << num_cells << std::endl;
   for (int i = 0; i < num_cells; i++)
      file << cells.volumes(i) << std::endl;
   file << std::endl;
   
   /* Write the cell centroids: */
   file << "cell-centroids " << num_cells+num_ghost_cells << std::endl;
   for (int i = 0; i < num_cells+num_ghost_cells; i++) {
      file << cells.centroids(i, 0) << " ";
      file << cells.centroids(i, 1) << " ";
      file << cells.centroids(i, 2) << std::endl;
   }
   file << std::endl;
   
   /* Write the cell materials (1-based indexed): */
   file << "cell-materials " << num_cells+num_ghost_cells << std::endl;
   for (int i = 0; i < num_cells+num_ghost_cells; i++)
      file << cells.materials(i)+1 << std::endl;
   file << std::endl;
   
   /* Write the cell indices in the nodal mesh: */
   if (!(cells.nodal_indices.empty())) {
      file << "cell-nodal-indices " << num_cells << std::endl;
      for (int i = 0; i < num_cells; i++)
         file << cells.nodal_indices(i) << std::endl;
      file << std::endl;
   }
   
   /* Write the cell indices in the global mesh: */
   file << "cell-global-indices " << num_cells+num_ghost_cells << std::endl;
   for (int i = 0; i < num_cells+num_ghost_cells; i++)
      file << cells.global_indices(i) << std::endl;
   file << std::endl;
   
   /* Write the number of cell faces: */
   file << "faces " << num_cells << std::endl;
   for (int i = 0; i < num_cells; i++)
      file << faces.num_faces(i) << std::endl;
   file << std::endl;
   
   /* Get the total number of cell faces: */
   int num_cell_faces = 0;
   for (int i = 0; i < num_cells; i++)
      num_cell_faces += faces.num_faces(i);
   
   /* Write the face areas: */
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
   file << "face-neighbours " << num_cells << " " << num_cell_faces << std::endl;
   for (int i = 0; i < num_cells; i++) {
      for (int j = 0; j < faces.num_faces(i); j++) {
         if (j > 0) file << " ";
         file << faces.neighbours(i, j);
      }
      file << std::endl;
   }
   file << std::endl;
   
   /* Write the boundary names: */
   for (int i = 0; i < boundaries.size(); i++)
      file << "boundary " << boundaries(i) << std::endl;
   
   /* Write the boundary conditions (1-based indexed): */
   for (int i = 1; i < bcs.size(); i++) {
      file << "bc " << boundaries(i-1);
      if (bcs(i).type == BC::VACUUM)
         file << " vacuum";
      else if (bcs(i).type == BC::REFLECTIVE)
         file << " reflective";
      else if (bcs(i).type == BC::ROBIN)
         file << " robin";
      else if (bcs(i).type == BC::DIRICHLET)
         file << " dirichlet";
      else if (bcs(i).type == BC::ADIABATIC)
         file << " adiabatic";
      else if (bcs(i).type == BC::CONVECTION)
         file << " convection";
      else {
         PAMPA_CHECK(true, "wrong boundary-condition type");
      }
      if (bcs(i).type == BC::ROBIN || bcs(i).type == BC::DIRICHLET || bcs(i).type == BC::CONVECTION)
         file << " " << bcs(i).f(0)(0.0);
      if (bcs(i).type == BC::CONVECTION)
         file << " " << bcs(i).f(1)(0.0);
      file << std::endl;
   }
   file << std::endl;
   
   return 0;
   
}

/* Get the domain indices to partition the mesh: */
#ifdef WITH_METIS
int Mesh::getDomainIndices(Array1D<int>& part, Array1D<int>& size) {
   
   /* Get the number of graph vertices, i.e. cells, and the total number of edges, i.e. faces: */
   /* Note: only faces with physical neighboring cells need to be considered here. */
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
   METIS_CALL(METIS_PartGraphRecursive(&nvtxs, &ncon, &xadj(0), &adjncy(0), nullptr, nullptr, 
      nullptr, &nparts, nullptr, nullptr, &options(0), &objval, &part(0)));
   #elif METIS_PARTGRAPH == METIS_KWAY
   METIS_CALL(METIS_PartGraphKway(&nvtxs, &ncon, &xadj(0), &adjncy(0), nullptr, nullptr, nullptr, 
      &nparts, nullptr, nullptr, &options(0), &objval, &part(0)));
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
