#include "DiffusionSolver.hxx"

/* Build the coefficient matrices and solution vectors: */
int DiffusionSolver::build() {
   
   /* Build the coefficient matrices: */
   PAMPA_CALL(buildMatrices(), "unable to build the coefficient matrices");
   
   /* Build the solution vectors: */
   PAMPA_CALL(buildVectors(), "unable to build the solution vectors");
   
   return 0;
   
}

/* Build the coefficient matrices: */
int DiffusionSolver::buildMatrices() {
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   int num_cells_global = mesh->getNumCellsGlobal();
   int num_faces_max = mesh->getNumFacesMax();
   const Cells& cells = mesh->getCells();
   const Faces& faces = mesh->getFaces();
   const Array1D<BoundaryCondition>& bcs = mesh->getBoundaryConditions();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Create, preallocate and set up the coefficient matrices: */
   int size_local = num_cells * num_groups;
   int size_global = num_cells_global * num_groups;
   int num_r_nonzero_max = num_faces_max + num_groups;
   int num_f_nonzero = num_groups;
   petsc::create_matrix(R, size_local, size_global, num_r_nonzero_max);
   petsc::create_matrix(F, size_local, size_global, num_f_nonzero);
   
   /* Get the local ownership range: */
   int lmin, lmax, f_lmin, f_lmax;
   MatGetOwnershipRange(R, &lmin, &lmax);
   MatGetOwnershipRange(F, &f_lmin, &f_lmax);
   PAMPA_CHECK((f_lmin != lmin) || (f_lmax != lmax), 1, "wrong local ownership range");
   
   /* Initialize the matrix rows for R and F: */
   PetscInt r_l2[num_r_nonzero_max];
   PetscInt f_l2[num_f_nonzero];
   PetscScalar r_l_l2[num_r_nonzero_max];
   PetscScalar f_l_l2[num_f_nonzero];
   
   /* Calculate the coefficients for each cell i: */
   for (int i = 0; i < num_cells; i++) {
      
      /* Get the material for cell i: */
      const Material& mat = materials(cells.materials(i));
      
      /* Calculate the coefficients for each group g: */
      for (int g = 0; g < num_groups; g++) {
         
         /* Get the matrix index for cell i and group g: */
         PetscInt l = cells.indices(i)*num_groups + g;
         int r_i = 1, f_i = 0;
         
         /* Set the total-reaction term: */
         r_l2[0] = l;
         r_l_l2[0] = mat.sigma_total(g) * cells.volumes(i);
         
         /* Set the group-to-group coupling terms: */
         for (int g2 = 0; g2 < num_groups; g2++) {
            
            /* Get the matrix index for cell i and group g2: */
            PetscInt l2 = cells.indices(i)*num_groups + g2;
            
            /* Set the (g2 -> g) scattering term: */
            if (l2 == l)
               r_l_l2[0] += -mat.sigma_scattering(g2, g) * cells.volumes(i);
            else {
               r_l2[r_i] = l2;
               r_l_l2[r_i++] = -mat.sigma_scattering(g2, g) * cells.volumes(i);
            }
            
            /* Set the (g2 -> g) fission term: */
            f_l2[f_i] = l2;
            f_l_l2[f_i++] = mat.chi(g) * mat.nu_sigma_fission(g2) * cells.volumes(i);
            
         }
         
         /* Set the cell-to-cell coupling terms: */
         double r;
         for (int f = 0; f < faces.num_faces(i); f++) {
            
            /* Get the index for cell i2 (actual cell or boundary condition): */
            /* Note: boundary conditions have negative, 1-based indexes: */
            int i2 = faces.neighbours(i, f);
            
            /* Set the boundary conditions: */
            if (i2 < 0) {
               
               /* Check the boundary-condition type: */
               switch (bcs(-i2).type) {
                  
                  /* Set vacuum (zero-flux) boundary conditions: */
                  case BC::VACUUM : {
                     
                     /* Get the surface leakage factor: */
                     double w = math::surface_leakage_factor(cells.centroids(i), 
                                   faces.centroids(i, f), faces.normals(i, f));
                     
                     /* Set the leakage term for cell i: */
                     r_l_l2[0] += w * mat.diffusion_coefficient(g) * faces.areas(i, f);
                     
                     break;
                     
                  }
                  
                  /* Set reflective (zero-current) boundary conditions (nothing to be done): */
                  case BC::REFLECTIVE : {
                     
                     break;
                     
                  }
                  
                  /* Set Robin boundary conditions: */
                  case BC::ROBIN : {
                     
                     /* Set the leakage term for cell i: */
                     r_l_l2[0] -= bcs(-i2).a * faces.areas(i, f);
                     
                     break;
                     
                  }
                  
               }
               
            }
            
            /* Set the cell-to-cell coupling terms depending on the neighbour material: */
            else {
               
               /* Get the matrix index for cell i2 and group g: */
               PetscInt l2 = cells.indices(i2)*num_groups + g;
               
               /* Get the material for cell i2: */
               const Material& mat2 = materials(cells.materials(i2));
               
               /* Set the terms for cells with the same materials: */
               if (&mat2 == &mat) {
                  
                  /* Get the surface leakage factor: */
                  double w = math::surface_leakage_factor(cells.centroids(i), cells.centroids(i2), 
                                faces.normals(i, f));
                  
                  /* Get the leakage term for cell i2: */
                  r = -w * mat.diffusion_coefficient(g) * faces.areas(i, f);
                  
               }
               
               /* Set the terms for cells with different materials: */
               else {
                  
                  /* Get the surface leakage factor and the weight for cell i: */
                  double w_i_i2 = math::surface_leakage_factor(cells.centroids(i), 
                                     faces.centroids(i, f), faces.normals(i, f));
                  w_i_i2 *= mat.diffusion_coefficient(g) * faces.areas(i, f);
                  
                  /* Get the surface leakage factor and the weight for cell i2: */
                  double w_i2_i = math::surface_leakage_factor(cells.centroids(i2), 
                                     faces.centroids(i, f), faces.normals(i, f));
                  w_i2_i *= -mat2.diffusion_coefficient(g) * faces.areas(i, f);
                  
                  /* Get the leakage term for cell i2: */
                  r = -(w_i_i2*w_i2_i) / (w_i_i2+w_i2_i);
                  
               }
               
               /* Set the leakage term for cell i: */
               r_l_l2[0] -= r;
               
               /* Set the leakage term for cell i2: */
               r_l2[r_i] = l2;
               r_l_l2[r_i++] = r;
               
            }
            
         }
         
         /* Set the matrix rows for R and F: */
         PETSC_CALL(MatSetValues(R, 1, &l, r_i, r_l2, r_l_l2, INSERT_VALUES));
         PETSC_CALL(MatSetValues(F, 1, &l, f_i, f_l2, f_l_l2, INSERT_VALUES));
         
      }
      
   }
   
   /* Assembly the coefficient matrices: */
   PETSC_CALL(MatAssemblyBegin(R, MAT_FINAL_ASSEMBLY));
   PETSC_CALL(MatAssemblyBegin(F, MAT_FINAL_ASSEMBLY));
   PETSC_CALL(MatAssemblyEnd(R, MAT_FINAL_ASSEMBLY));
   PETSC_CALL(MatAssemblyEnd(F, MAT_FINAL_ASSEMBLY));
   
   return 0;
   
}

/* Build the solution vectors: */
int DiffusionSolver::buildVectors() {
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   int num_cells_global = mesh->getNumCellsGlobal();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Create the scalar-flux MPI vector: */
   PETSC_CALL(MatCreateVecs(R, NULL, &phi_mpi));
   
   /* Create the scalar-flux sequential vector: */
   PETSC_CALL(VecCreateSeq(MPI_COMM_SELF, num_cells_global*num_groups, &phi_seq));
   PETSC_CALL(VecZeroEntries(phi_seq));
   
   return 0;
   
}

/* Get the solution after solving the eigensystem: */
int DiffusionSolver::getSolution() {
   
   /* Get the scalar flux from the EPS context: */
   double lambda;
   PETSC_CALL(EPSGetEigenpair(eps, 0, &lambda, NULL, phi_mpi, NULL));
   keff = 1.0 / lambda;
   
   /* Gather the scalar flux from all ranks: */
   VecScatter context;
   PETSC_CALL(VecScatterCreateToAll(phi_mpi, &context, &phi_seq));
   PETSC_CALL(VecScatterBegin(context, phi_mpi, phi_seq, INSERT_VALUES, SCATTER_FORWARD));
   PETSC_CALL(VecScatterEnd(context, phi_mpi, phi_seq, INSERT_VALUES, SCATTER_FORWARD));
   PETSC_CALL(VecScatterDestroy(&context));
   
   /* Normalize the scalar flux: */
   PAMPA_CALL(normalizeScalarFlux(), "unable to normalize the scalar flux");
   
   return 0;
   
}

/* Write the solution to a plain-text file in .vtk format: */
int DiffusionSolver::writeVTK(const std::string& filename) const {
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   const Cells& cells = mesh->getCells();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Get the array for the scalar flux: */
   PetscScalar* data_phi;
   PETSC_CALL(VecGetArray(phi_seq, &data_phi));
   
   /* Write the mesh: */
   PAMPA_CALL(mesh->writeVTK(filename), "unable to write the mesh");
   
   /* Open the output file: */
   std::ofstream file(filename, std::ios_base::app);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Write the scalar flux: */
   for (int g = 0; g < num_groups; g++) {
      file << "SCALARS flux_" << (g+1) << " double 1" << std::endl;
      file << "LOOKUP_TABLE default" << std::endl;
      for (int i = 0; i < num_cells; i++)
         file << data_phi[cells.indices(i)*num_groups+g] << std::endl;
      file << std::endl;
   }
   
   /* Restore the array for the scalar flux: */
   PETSC_CALL(VecRestoreArray(phi_seq, &data_phi));
   
   return 0;
   
}

/* Write the solution to a binary file in PETSc format: */
int DiffusionSolver::writePETSc(const std::string& filename) const {
   
   /* Write the solution to a binary file: */
   PetscViewer viewer;
   PETSC_CALL(PetscViewerBinaryOpen(MPI_COMM_WORLD, "flux.ptc", FILE_MODE_WRITE, &viewer));
   PETSC_CALL(VecView(phi_mpi, viewer));
   PETSC_CALL(PetscViewerDestroy(&viewer));
   
   return 0;
   
}

/* Destroy the solution vectors: */
int DiffusionSolver::destroyVectors() {
   
   /* Destroy the solution vectors: */
   PETSC_CALL(VecDestroy(&phi_mpi));
   PETSC_CALL(VecDestroy(&phi_seq));
   
   return 0;
   
}
