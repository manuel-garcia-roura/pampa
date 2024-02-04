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
   const Array2D<double>& directions = quadrature.getDirections();
   const Array1D<double>& weights = quadrature.getWeights();
   const Array2D<double>& axes = quadrature.getAxes();
   const Array2D<int>& reflected_directions = quadrature.getReflectedDirections();
   
   /* Create, preallocate and set up the coefficient matrices: */
   petsc::create_matrix(R, num_cells*num_groups*num_directions, 6+num_groups*num_directions);
   petsc::create_matrix(F, num_cells*num_groups*num_directions, num_groups*num_directions);
   
   /* Get the local ownership range: */
   int lmin, lmax, f_lmin, f_lmax;
   MatGetOwnershipRange(R, &lmin, &lmax);
   MatGetOwnershipRange(F, &f_lmin, &f_lmax);
   PAMPA_CHECK((f_lmin != lmin) || (f_lmax != lmax), 1, "wrong local ownership range");
   
   /* Initialize the matrix rows for R and F: */
   PetscInt r_l2[6+num_groups*num_directions];
   PetscInt f_l2[num_groups*num_directions];
   PetscScalar r_l_l2[6+num_groups*num_directions];
   PetscScalar f_l_l2[num_groups*num_directions];
   
   /* Calculate the coefficients for each cell i: */
   for (int i = 0; i < num_cells; i++) {
      
      /* Get the material for cell i: */
      const Material& mat = materials[cells.materials[i]];
      
      /* Calculate the coefficients for each group g: */
      for (int g = 0; g < num_groups; g++) {
         
         /* Calculate the coefficients for each direction m: */
         for (int m = 0; m < num_directions; m++) {
            
            /* Get the matrix index for cell i, group g and direction m: */
            PetscInt l = i*num_directions*num_groups + g*num_directions + m;
            int r_i = 1, f_i = 0;
            
            /* Check if the matrix index is local: */
            if ((l >= lmax) || (l < lmin))
               continue;
            
            /* Set the total-reaction term: */
            r_l2[0] = l;
            r_l_l2[0] = mat.sigma_total(g) * cells.volumes(i);
            
            /* Set the group-to-group coupling terms: */
            for (int g2 = 0; g2 < num_groups; g2++) {
               
               /* Set the direction-to-direction coupling terms: */
               for (int m2 = 0; m2 < num_directions; m2++) {
                  
                  /* Get the matrix index for cell i, group g2 and direction m2: */
                  PetscInt l2 = i*num_directions*num_groups + g2*num_directions + m2;
                  
                  /* Set the (g2 -> g, m2 -> m) isotropic scattering term: */
                  if (l2 == l)
                     r_l_l2[0] += -mat.sigma_scattering[g2][g] * weights(m2) * cells.volumes(i);
                  else {
                     r_l2[r_i] = l2;
                     r_l_l2[r_i++] = -mat.sigma_scattering[g2][g] * weights(m2) * cells.volumes(i);
                  }
                  
                  /* Set the (g2 -> g, m2 -> m) fission term: */
                  f_l2[f_i] = l2;
                  f_l_l2[f_i++] = mat.chi(g) * mat.nu_sigma_fission(g2) * weights(m2) * 
                                     cells.volumes(i);
                  
               }
               
            }
            
            /* Set the cell-to-cell coupling terms: */
            for (int f = 0; f < faces.num_faces(i); f++) {
               
               /* Get the index for cell i2 (actual cell or boundary condition): */
               /* Note: boundary conditions have negative, 1-based indexes: */
               int i2 = faces.neighbours(i, f);
               
               /* Set the boundary conditions: */
               if (i2 < 0) {
                  
                  /* Check the boundary-condition type: */
                  switch (bcs[-i2].type) {
                     
                     /* Set vacuum (zero-incomming-current) boundary conditions: */
                     case BC::VACUUM : {
                        
                        /* Get the dot product between the direction and the face normal: */
                        double w = math::dot_product(directions(m), faces.normals(i, f), 3);
                        
                        /* Set the leakage term for cell i for outgoing directions: */
                        if (w > 0.0) r_l_l2[0] += w * faces.areas(i, f);
                        
                        break;
                        
                     }
                     
                     /* Set reflective (zero-current) boundary conditions: */
                     case BC::REFLECTIVE : {
                        
                        /* Get the dot product between the direction and the face normal: */
                        double w = math::dot_product(directions(m), faces.normals(i, f), 3);
                        
                        /* Set the leakage term for cell i depending on the direction: */
                        if (w > 0.0)
                           
                           /* Set the leakage term for cell i for outgoing directions: */
                           r_l_l2[0] += w * faces.areas(i, f);
                           
                        else {
                           
                           /* Get the reflected outgoing direction for incoming directions: */
                           int m2 = -1;
                           for (int i = 0; i < 3; i++) {
                              double d = math::dot_product(axes(i), faces.normals(i, f), 3);
                              if (fabs(d) > 1.0-TOL)
                                 m2 = reflected_directions(m, i);
                           }
                           PAMPA_CHECK(m2 == -1, 1, "reflected direction not found");
                           
                           /* Get the matrix index for cell i, group g and direction m2: */
                           PetscInt l2 = i*num_directions*num_groups + g*num_directions + m2;
                           
                           /* Set the leakage term for cell i for incoming directions: */
                           r_l2[r_i] = l2;
                           r_l_l2[r_i++] = w * faces.areas(i, f);
                           
                        }
                        
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
                  PetscInt l2 = i2*num_directions*num_groups + g*num_directions + m;
                  
                  /* Get the distances between the cell centers and the face: */
                  double r_i_f = math::get_distance(faces.centroids(i, f), cells.centroids(i), 3);
                  double r_i2_f = math::get_distance(faces.centroids(i, f), cells.centroids(i2), 3);
                  double r_i_i2 = math::get_distance(cells.centroids(i), cells.centroids(i2), 3);
                  
                  /* Get the dot product between the direction and the face normal: */
                  double w = math::dot_product(directions(m), faces.normals(i, f), 3);
                  
                  /* Get the flux weights for the face flux: */
                  double delta = 0.01;
                  double w_i, w_i2;
                  if (w > 0.0) {
                     w_i = (r_i2_f+delta*r_i_f) / r_i_i2;
                     w_i2 = (1.0-delta)*r_i_f / r_i_i2;
                  }
                  else {
                     w_i = (1.0-delta)*r_i2_f / r_i_i2;
                     w_i2 = (r_i_f+delta*r_i2_f) / r_i_i2;
                  }
                  
                  /* Set the leakage term for cell i: */
                  r_l_l2[0] += w_i * w * faces.areas(i, f);
                  
                  /* Set the leakage term for cell i2: */
                  r_l2[r_i] = l2;
                  r_l_l2[r_i++] = w_i2 * w * faces.areas(i, f);
                  
               }
               
            }
            
            /* Set the matrix rows for R and F: */
            PETSC_CALL(MatSetValues(R, 1, &l, r_i, r_l2, r_l_l2, INSERT_VALUES));
            PETSC_CALL(MatSetValues(F, 1, &l, f_i, f_l2, f_l_l2, INSERT_VALUES));
            
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
   const Array1D<double>& weights = quadrature.getWeights();
   
   /* Get the arrays for the scalar and angular fluxes: */
   PetscScalar *data_phi, *data_psi;
   PETSC_CALL(VecGetArray(phi_seq, &data_phi));
   PETSC_CALL(VecGetArray(psi_seq, &data_psi));
   
   /* Integrate the angular flux over all directions: */
   for (int i = 0; i < num_cells; i++) {
      for (int g = 0; g < num_groups; g++) {
         data_phi[i*num_groups+g] = 0.0;
         for (int m = 0; m < num_directions; m++)
            data_phi[i*num_groups+g] += weights(m) * 
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
   const Array1D<double>& weights = quadrature.getWeights();
   
   /* Get the array for the angular flux: */
   PetscScalar* data_psi;
   PETSC_CALL(VecGetArray(psi_seq, &data_psi));
   
   /* Normalize the angular flux (TODO: normalize correctly with the power!): */
   double vol = 0.0;
   for (int i = 0; i < num_cells; i++)
      if (materials[cells.materials[i]].nu_sigma_fission(1) > 0.0)
         vol += cells.volumes(i);
   double sum = 0.0;
   for (int g = 0; g < num_groups; g++)
      for (int i = 0; i < num_cells; i++)
         for (int m = 0; m < num_directions; m++)
            sum += weights(m) * data_psi[i*num_directions*num_groups+g*num_directions+m] * 
                      materials[cells.materials[i]].nu_sigma_fission(g) * cells.volumes(i);
   double f = vol / sum;
   for (int g = 0; g < num_groups; g++)
      for (int i = 0; i < num_cells; i++)
         for (int m = 0; m < num_directions; m++)
            data_psi[i*num_directions*num_groups+g*num_directions+m] *= f;
   
   /* Restore the array for the angular flux: */
   PETSC_CALL(VecRestoreArray(psi_seq, &data_psi));
   
   return 0;
   
}
