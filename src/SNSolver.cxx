#include "SNSolver.hxx"

/* Build the coefficient matrices and solution vectors: */
int SNSolver::build() {
   
   /* Build the angular quadrature set: */
   quadrature = AngularQuadratureSet(method.order);
   PAMPA_CALL(quadrature.build(), "unable to build the angular quadrature set");
   
   /* Build the coefficient matrices: */
   PAMPA_CALL(buildMatrices(), "unable to build the coefficient matrices");
   
   /* Build the solution vectors: */
   PAMPA_CALL(buildVectors(), "unable to build the solution vectors");
   
   return 0;
   
}

/* Build the coefficient matrices: */
int SNSolver::buildMatrices() {
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   const Cells& cells = mesh->getCells();
   const Faces& faces = mesh->getFaces();
   const std::vector<BoundaryCondition>& bcs = mesh->getBoundaryConditions();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Get the angular quadrature data: */
   int num_directions = quadrature.getNumDirections();
   const std::vector<std::vector<double>>& directions = quadrature.getDirections();
   const std::vector<double>& weights = quadrature.getWeights();
   
   /* Create, preallocate and set up the coefficient matrices: */
   petsc::create_matrix(R, num_cells*num_groups*num_directions, 6+num_groups*num_directions);
   petsc::create_matrix(F, num_cells*num_groups*num_directions, num_groups*num_directions);
   
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
         
         /* Calculate the coefficients for each direction m: */
         for (int m = 0; m < num_directions; m++) {
            
            /* Get the matrix index for cell i, group g and direction m: */
            int l = i*num_directions*num_groups + g*num_directions + m;
            
            /* Check if the matrix index is local: */
            if ((l >= l2) || (l < l1))
               continue;
            
            /* Set the total-reaction term: */
            double r_l_l = mat.sigma_total[g] * cells.volumes[i];
            
            /* Set the group-to-group coupling terms: */
            for (int g2 = 0; g2 < num_groups; g2++) {
               
               /* Set the direction-to-direction coupling terms: */
               for (int m2 = 0; m2 < num_directions; m2++) {
                  
                  /* Get the matrix index for cell i, group g2 and direction m2: */
                  int l2 = i*num_directions*num_groups + g2*num_directions + m2;
                  
                  /* Set the (g2 -> g, m2 -> m) isotropic scattering term: */
                  double r_l_l2 = -mat.sigma_scattering[g2][g] * weights[m2] * cells.volumes[i];
                  if (l2 == l)
                     r_l_l += r_l_l2;
                  else {
                     PETSC_CALL(MatSetValues(R, 1, &l, 1, &l2, &r_l_l2, INSERT_VALUES));
                  }
                  
                  /* Set the (g2 -> g, m2 -> m) fission term: */
                  double f_l_l2 = mat.chi[g] * mat.nu_sigma_fission[g2] * weights[m2] * 
                                     cells.volumes[i];
                  PETSC_CALL(MatSetValues(F, 1, &l, 1, &l2, &f_l_l2, INSERT_VALUES));
                  
               }
               
            }
            
            /* Set the cell-to-cell coupling terms: */
            for (int f = 0; f < faces.neighbours[i].size(); f++) {
               
               /* Get the index for cell i2 (actual cell or boundary condition): */
               /* Note: boundary conditions have negative, 1-based indexes: */
               int i2 = faces.neighbours[i][f];
               
               /* Set the boundary conditions: */
               if (i2 < 0) {
                  
                  /* Check the boundary-condition type: */
                  switch (bcs[-i2].type) {
                     
                     /* Set vacuum (zero-incomming-current) boundary conditions: */
                     case BC::VACUUM : {
                        
                        /* Get the geometrical data: */
                        const std::vector<double>& n_i_f = faces.normals[i][f];
                        
                        /* Get the dot product between the direction and the face normal: */
                        double w = math::dot_product(directions[m], n_i_f, 3);
                        
                        /* Set the leakage term for cell i for outgoing directions: */
                        if (w > 0.0) r_l_l += w * faces.areas[i][f];
                        
                        break;
                        
                     }
                     
                     /* Set reflective (zero-current) boundary conditions (not implemented): */
                     case BC::REFLECTIVE : {
                        
                        /* Not implemented: */
                        PAMPA_CHECK(true, 1, "reflective boundary conditions not implemented");
                        
                        break;
                        
                     }
                     
                     /* Set Robin boundary conditions (not implemented): */
                     case BC::ROBIN : {
                        
                        /* Not implemented: */
                        PAMPA_CHECK(true, 1, "Robin boundary conditions not implemented");
                        
                        break;
                        
                     }
                     
                  }
                  
               }
               
               /* Set the cell-to-cell coupling terms: */
               else {
                  
                  /* Get the matrix index for cell i2, group g and direction m: */
                  int l2 = i2*num_directions*num_groups + g*num_directions + m;
                  
                  /* Get the geometrical data: */
                  const std::vector<double>& p_i = cells.centroids[i];
                  const std::vector<double>& p_i2 = cells.centroids[i2];
                  const std::vector<double>& p_f = faces.centroids[i][f];
                  const std::vector<double>& n_i_f = faces.normals[i][f];
                  
                  /* Get the distances between the cell centers and the face: */
                  double r_i = math::l2_norm(math::subtract(p_f, p_i, 3), 3);
                  double r_i2 = math::l2_norm(math::subtract(p_f, p_i2, 3), 3);
                  
                  /* Get the dot product between the direction and the face normal: */
                  double w = math::dot_product(directions[m], n_i_f, 3);
                  
                  /* Get the flux weights for the face flux: */
                  double delta = 0.01;
                  double w_i, w_i2;
                  if (w > 0.0) {
                     w_i = (r_i2+delta*r_i) / (r_i+r_i2);
                     w_i2 = (1.0-delta)*r_i / (r_i+r_i2);
                  }
                  else {
                     w_i = (1.0-delta)*r_i2 / (r_i+r_i2);
                     w_i2 = (r_i+delta*r_i2) / (r_i+r_i2);
                  }
                  
                  /* Set the leakage term for cell i: */
                  r_l_l += w_i * w * faces.areas[i][f];
                  
                  /* Set the leakage term for cell i2: */
                  double r_l_l2 = w_i2 * w * faces.areas[i][f];
                  PETSC_CALL(MatSetValues(R, 1, &l, 1, &l2, &r_l_l2, INSERT_VALUES));
                  
               }
               
            }
            
            /* Set the diagonal coefficient: */
            PETSC_CALL(MatSetValues(R, 1, &l, 1, &l, &r_l_l, INSERT_VALUES));
            
         }
         
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
int SNSolver::buildVectors() {
   
   /* Get the number of cells: */
   int num_cells = mesh->getNumCells();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Get the number of directions: */
   int num_directions = quadrature.getNumDirections();
   
   /* Create the scalar-flux MPI vector: */
   PETSC_CALL(VecCreate(MPI_COMM_WORLD, &phi_mpi));
   PETSC_CALL(VecSetSizes(phi_mpi, PETSC_DECIDE, num_cells*num_groups));
   PETSC_CALL(VecSetFromOptions(phi_mpi));
   
   /* Create the angular-flux MPI vector: */
   PETSC_CALL(MatCreateVecs(R, NULL, &psi_mpi));
   
   /* Create the scalar-flux sequential vector: */
   PETSC_CALL(VecCreateSeq(MPI_COMM_SELF, num_cells*num_groups, &phi_seq));
   PETSC_CALL(VecZeroEntries(phi_seq));
   
   /* Create the angular-flux sequential vector: */
   PETSC_CALL(VecCreateSeq(MPI_COMM_SELF, num_cells*num_groups*num_directions, &psi_seq));
   PETSC_CALL(VecZeroEntries(psi_seq));
   
   return 0;
   
}

/* Get the solution after solving the eigensystem: */
int SNSolver::getSolution() {
   
   /* Get the number of cells: */
   int num_cells = mesh->getNumCells();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Get the number of directions: */
   int num_directions = quadrature.getNumDirections();
   
   /* Get the angular flux from the EPS context: */
   double lambda;
   PETSC_CALL(EPSGetEigenpair(eps, 0, &lambda, NULL, psi_mpi, NULL));
   keff = 1.0 / lambda;
   
   /* Gather the angular flux from all ranks: */
   VecScatter context;
   PETSC_CALL(VecScatterCreateToAll(psi_mpi, &context, &psi_seq));
   PETSC_CALL(VecScatterBegin(context, psi_mpi, psi_seq, INSERT_VALUES, SCATTER_FORWARD));
   PETSC_CALL(VecScatterEnd(context, psi_mpi, psi_seq, INSERT_VALUES, SCATTER_FORWARD));
   PETSC_CALL(VecScatterDestroy(&context));
   
   /* Calculate the scalar flux: */
   PAMPA_CALL(calculateScalarFlux(), "unable to calculate the scalar flux");
   
   /* Normalize the scalar flux: */
   PAMPA_CALL(normalizeScalarFlux(), "unable to normalize the scalar flux");
   
   /* Normalize the angular flux: */
   PAMPA_CALL(normalizeAngularFlux(), "unable to normalize the angular flux");
   
   return 0;
   
}

/* Write the solution to a plain-text file in .vtk format: */
int SNSolver::writeVTK(const std::string& filename) const {
   
   /* Get the number of cells: */
   int num_cells = mesh->getNumCells();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Get the number of directions: */
   int num_directions = quadrature.getNumDirections();
   
   /* Get the arrays for the scalar and angular fluxes: */
   PetscScalar *data_phi, *data_psi;
   PETSC_CALL(VecGetArray(phi_seq, &data_phi));
   PETSC_CALL(VecGetArray(psi_seq, &data_psi));
   
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
      
      /* Write the angular flux: */
      for (int g = 0; g < num_groups; g++) {
         for (int m = 0; m < num_directions; m++) {
            file << "SCALARS flux_" << (g+1) << "_" << (m+1) << " double 1" << std::endl;
            file << "LOOKUP_TABLE default" << std::endl;
            for (int i = 0; i < num_cells; i++)
               file << data_psi[i*num_directions*num_groups+g*num_directions+m] << std::endl;
            file << std::endl;
         }
      }
      
   }
   
   /* Restore the arrays for the scalar and angular fluxes: */
   PETSC_CALL(VecRestoreArray(phi_seq, &data_phi));
   PETSC_CALL(VecRestoreArray(psi_seq, &data_psi));
   
   return 0;
   
}

/* Write the solution to a binary file in PETSc format: */
int SNSolver::writePETSc(const std::string& filename) const {
   
   /* Write the solution to a binary file: */
   PetscViewer viewer;
   PETSC_CALL(PetscViewerBinaryOpen(MPI_COMM_WORLD, "flux.ptc", FILE_MODE_WRITE, &viewer));
   PETSC_CALL(VecView(psi_mpi, viewer));
   PETSC_CALL(PetscViewerDestroy(&viewer));
   
   return 0;
   
}

/* Destroy the solution vectors: */
int SNSolver::destroyVectors() {
   
   /* Destroy the solution vectors: */
   PETSC_CALL(VecDestroy(&phi_mpi));
   PETSC_CALL(VecDestroy(&phi_seq));
   PETSC_CALL(VecDestroy(&psi_mpi));
   PETSC_CALL(VecDestroy(&psi_seq));
   
   return 0;
   
}

/* Calculate the scalar flux: */
int SNSolver::calculateScalarFlux() {
   
   /* Get the number of cells: */
   int num_cells = mesh->getNumCells();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Get the angular quadrature data: */
   int num_directions = quadrature.getNumDirections();
   const std::vector<double>& weights = quadrature.getWeights();
   
   /* Get the arrays for the scalar and angular fluxes: */
   PetscScalar *data_phi, *data_psi;
   PETSC_CALL(VecGetArray(phi_seq, &data_phi));
   PETSC_CALL(VecGetArray(psi_seq, &data_psi));
   
   /* Integrate the angular flux over all directions: */
   for (int i = 0; i < num_cells; i++) {
      for (int g = 0; g < num_groups; g++) {
         data_phi[i*num_groups+g] = 0.0;
         for (int m = 0; m < num_directions; m++)
            data_phi[i*num_groups+g] += weights[m] * 
                                           data_psi[i*num_directions*num_groups+g*num_directions+m];
         data_phi[i*num_groups+g] *= 4.0 * M_PI;
      }
   }
   
   /* Restore the arrays for the scalar and angular fluxes: */
   PETSC_CALL(VecRestoreArray(phi_seq, &data_phi));
   PETSC_CALL(VecRestoreArray(psi_seq, &data_psi));
   
   return 0;
   
}

/* Normalize the angular flux: */
int SNSolver::normalizeAngularFlux() {
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   const Cells& cells = mesh->getCells();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Get the angular quadrature data: */
   int num_directions = quadrature.getNumDirections();
   const std::vector<double>& weights = quadrature.getWeights();
   
   /* Get the array for the angular flux: */
   PetscScalar* data_psi;
   PETSC_CALL(VecGetArray(psi_seq, &data_psi));
   
   /* Normalize the angular flux (TODO: normalize correctly with the power!): */
   double vol = 0.0;
   for (int i = 0; i < num_cells; i++)
      if (materials[cells.materials[i]].nu_sigma_fission[1] > 0.0)
         vol += cells.volumes[i];
   double sum = 0.0;
   for (int g = 0; g < num_groups; g++)
      for (int i = 0; i < num_cells; i++)
         for (int m = 0; m < num_directions; m++)
            sum += weights[m] * data_psi[i*num_directions*num_groups+g*num_directions+m] * 
                      materials[cells.materials[i]].nu_sigma_fission[g] * cells.volumes[i];
   double f = vol / sum;
   for (int g = 0; g < num_groups; g++)
      for (int i = 0; i < num_cells; i++)
         for (int m = 0; m < num_directions; m++)
            data_psi[i*num_directions*num_groups+g*num_directions+m] *= f;
   
   /* Restore the array for the angular flux: */
   PETSC_CALL(VecRestoreArray(psi_seq, &data_psi));
   
   return 0;
   
}
