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
   const Cells& cells = mesh->getCells();
   const Faces& faces = mesh->getFaces();
   const std::vector<BoundaryCondition>& bcs = mesh->getBoundaryConditions();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Create, preallocate and set up the coefficient matrices: */
   petsc::create_matrix(R, num_cells*num_groups, 6+num_groups);
   petsc::create_matrix(F, num_cells*num_groups, num_groups);
   
   /* Get the local ownership range: */
   int l1, l2, f_l1, f_l2;
   MatGetOwnershipRange(R, &l1, &l2);
   MatGetOwnershipRange(F, &f_l1, &f_l2);
   PAMPA_CHECK((f_l1 != l1) || (f_l2 != l2), 1, "wrong local ownership range");
   
   /* Calculate the coefficients for each cell i: */
   for (int i = 0; i < num_cells; i++) {
      
      /* Get the material for cell i: */
      const Material& mat = materials[cells.materials[i]];
      
      /* Calculate the coefficients for each group g: */
      for (int g = 0; g < num_groups; g++) {
         
         /* Get the matrix index for cell i and group g: */
         int l = i*num_groups + g;
         
         /* Check if the matrix index is local: */
         if ((l >= l2) || (l < l1))
            continue;
         
         /* Set the total-reaction term: */
         double r_l_l = mat.sigma_total[g] * cells.volumes[i];
         
         /* Set the group-to-group coupling terms: */
         for (int g2 = 0; g2 < num_groups; g2++) {
            
            /* Get the matrix index for cell i and group g2: */
            int l2 = i*num_groups + g2;
            
            /* Set the (g2 -> g) scattering term: */
            double r_l_l2 = -mat.sigma_scattering[g2][g] * cells.volumes[i];
            if (l2 == l)
               r_l_l += r_l_l2;
            else {
               PETSC_CALL(MatSetValues(R, 1, &l, 1, &l2, &r_l_l2, INSERT_VALUES));
            }
            
            /* Set the (g2 -> g) fission term: */
            double f_l_l2 = mat.chi[g] * mat.nu_sigma_fission[g2] * cells.volumes[i];
            PETSC_CALL(MatSetValues(F, 1, &l, 1, &l2, &f_l_l2, INSERT_VALUES));
            
         }
         
         /* Set the cell-to-cell coupling terms: */
         double r_l_l2;
         for (int f = 0; f < faces.neighbours[i].size(); f++) {
            
            /* Get the index for cell i2 (actual cell or boundary condition): */
            /* Note: boundary conditions have negative, 1-based indexes: */
            int i2 = faces.neighbours[i][f];
            
            /* Set the boundary conditions: */
            if (i2 < 0) {
               
               /* Check the boundary-condition type: */
               switch (bcs[-i2].type) {
                  
                  /* Set vacuum (zero-flux) boundary conditions: */
                  case BC::VACUUM : {
                     
                     /* Get the geometrical data: */
                     const std::vector<double>& p_i = cells.centroids[i];
                     const std::vector<double>& p_f = faces.centroids[i][f];
                     const std::vector<double>& n_i_f = faces.normals[i][f];
                     
                     /* Get the surface leakage factor: */
                     double w = math::surface_leakage_factor(p_i, p_f, n_i_f);
                     
                     /* Set the leakage term for cell i: */
                     r_l_l += w * mat.diffusion_coefficient[g] * faces.areas[i][f];
                     
                     break;
                     
                  }
                  
                  /* Set reflective (zero-current) boundary conditions (nothing to be done): */
                  case BC::REFLECTIVE : {
                     
                     break;
                     
                  }
                  
                  /* Set Robin boundary conditions: */
                  case BC::ROBIN : {
                     
                     /* Set the leakage term for cell i: */
                     r_l_l -= bcs[-i2].a * faces.areas[i][f];
                     
                     break;
                     
                  }
                  
               }
               
            }
            
            /* Set the cell-to-cell coupling terms depending on the neighbour material: */
            else {
               
               /* Get the matrix index for cell i2 and group g: */
               int l2 = i2*num_groups + g;
               
               /* Get the material for cell i2: */
               const Material& mat2 = materials[cells.materials[i2]];
               
               /* Set the terms for cells with the same materials: */
               if (&mat2 == &mat) {
                  
                  /* Get the geometrical data: */
                  const std::vector<double>& p_i = cells.centroids[i];
                  const std::vector<double>& p_i2 = cells.centroids[i2];
                  const std::vector<double>& n_i_f = faces.normals[i][f];
                  
                  /* Get the surface leakage factor: */
                  double w = math::surface_leakage_factor(p_i, p_i2, n_i_f);
                  
                  /* Get the leakage term for cell i2: */
                  r_l_l2 = -w * mat.diffusion_coefficient[g] * faces.areas[i][f];
                  
               }
               
               /* Set the terms for cells with different materials: */
               else {
                  
                  /* Get the geometrical data: */
                  const std::vector<double>& p_i = cells.centroids[i];
                  const std::vector<double>& p_i2 = cells.centroids[i2];
                  const std::vector<double>& p_f = faces.centroids[i][f];
                  const std::vector<double>& n_i_f = faces.normals[i][f];
                  
                  /* Get the surface leakage factor and the weight for cell i: */
                  double w_i_i2 = math::surface_leakage_factor(p_i, p_f, n_i_f);
                  w_i_i2 *= mat.diffusion_coefficient[g] * faces.areas[i][f];
                  
                  /* Get the surface leakage factor and the weight for cell i2: */
                  double w_i2_i = math::surface_leakage_factor(p_i2, p_f, n_i_f);
                  w_i2_i *= -mat2.diffusion_coefficient[g] * faces.areas[i][f];
                  
                  /* Get the leakage term for cell i2: */
                  r_l_l2 = -(w_i_i2*w_i2_i) / (w_i_i2+w_i2_i);
                  
               }
               
               /* Set the leakage term for cell i: */
               r_l_l -= r_l_l2;
               
               /* Set the leakage term for cell i2: */
               PETSC_CALL(MatSetValues(R, 1, &l, 1, &l2, &r_l_l2, INSERT_VALUES));
               
            }
            
         }
         
         /* Set the diagonal coefficient: */
         PETSC_CALL(MatSetValues(R, 1, &l, 1, &l, &r_l_l, INSERT_VALUES));
         
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
   
   /* Get the number of cells: */
   int num_cells = mesh->getNumCells();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Create the scalar-flux MPI vector: */
   PETSC_CALL(MatCreateVecs(R, NULL, &phi_mpi));
   
   /* Create the scalar-flux sequential vector: */
   PETSC_CALL(VecCreateSeq(MPI_COMM_SELF, num_cells*num_groups, &phi_seq));
   PETSC_CALL(VecZeroEntries(phi_seq));
   
   return 0;
   
}

/* Get the solution after solving the eigensystem: */
int DiffusionSolver::getSolution() {
   
   /* Get the number of cells: */
   int num_cells = mesh->getNumCells();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
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
   
   /* Get the number of cells: */
   int num_cells = mesh->getNumCells();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Get the array for the scalar flux: */
   PetscScalar* data_phi;
   PETSC_CALL(VecGetArray(phi_seq, &data_phi));
   
   /* Check the MPI rank: */
   if (mpi::rank == 0) {
      
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
            file << data_phi[i*num_groups+g] << std::endl;
         file << std::endl;
      }
      
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
