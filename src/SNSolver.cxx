#include "SNSolver.hxx"

/* Check the material data: */
int SNSolver::checkMaterials() {
   
   /* Check the materials: */
   for (int i = 0; i < materials.size(); i++) {
      PAMPA_CHECK(materials(i).sigma_total.empty(), 1, "missing total cross sections");
      PAMPA_CHECK(materials(i).nu_sigma_fission.empty(), 1, "missing nu-fission cross sections");
      PAMPA_CHECK(materials(i).sigma_scattering.empty(), 1, "missing scattering cross sections");
      PAMPA_CHECK(materials(i).chi.empty(), 1, "missing fission spectrum");
   }
   
   return 0;
   
}

/* Build the coefficient matrices and solution vectors: */
int SNSolver::build() {
   
   /* Build the angular quadrature set: */
   quadrature = AngularQuadratureSet(order);
   PAMPA_CALL(quadrature.build(), "unable to build the angular quadrature set");
   
   /* Build the cell-to-cell coupling coefficients for the gradient-discretization scheme: */
   PAMPA_CALL(buildGaussGradientScheme(grad_coefs, false), 
      "unable to build the Gauss gradient discretization");
   if (gradient.boundary_interpolation == BI::LS) {
      PAMPA_CALL(buildLSGradientScheme(grad_coefs_bc, true), 
         "unable to build the least-squares gradient discretization for boundary cells");
   }
   
   /* Build the coefficient matrices: */
   PAMPA_CALL(buildMatrices(), "unable to build the coefficient matrices");
   
   /* Build the solution vectors: */
   PAMPA_CALL(buildVectors(), "unable to build the solution vectors");
   
   return 0;
   
}

/* Get the mapping and the number of faces for boundary cells: */
int SNSolver::getBoundaryCells(Array1D<int>& num_faces_bc) {
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   const Faces& faces = mesh->getFaces();
   
   /* Get the number of boundary cells: */
   int num_cells_bc = 0;
   for (int i = 0; i < num_cells; i++) {
      for (int f = 0; f < faces.num_faces(i); f++) {
         if (faces.neighbours(i, f) < 0) {
            num_cells_bc++;
            break;
         }
      }
   }
   
   /* Get the boundary cells: */
   ic_to_ibc.resize(num_cells, -1);
   num_faces_bc.resize(num_cells_bc);
   for (int ibc = 0, i = 0; i < num_cells; i++) {
      int num_faces = faces.num_faces(i);
      for (int f = 0; f < num_faces; f++) {
         if (faces.neighbours(i, f) < 0) {
            ic_to_ibc(i) = ibc;
            num_faces_bc(ibc++) = num_faces;
            break;
         }
      }
   }
   
   return 0;
   
}

/* Build the coefficients for the Gauss gradient-discretization scheme: */
int SNSolver::buildGaussGradientScheme(Vector3D<double>& coefs, bool bc) {
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   const Cells& cells = mesh->getCells();
   const Faces& faces = mesh->getFaces();
   
   /* Get the weight between upwind and linear interpolation: */
   double delta = gradient.delta;
   
   /* Initialize the coefficients: */
   if (bc) {
      Array1D<int> num_faces_bc;
      PAMPA_CALL(getBoundaryCells(num_faces_bc), "unable to get the boundary cells");
      coefs.resize(num_faces_bc.size(), num_faces_bc, 4);
   }
   else
      coefs.resize(num_cells, faces.num_faces, 4);
   
   /* Build the coefficients: */
   for (int i = 0; i < num_cells; i++) {
      
      /* Check if the coefficients for this cell are needed: */
      if (bc && ic_to_ibc(i) < 0) continue;
      int ic = (bc) ? ic_to_ibc(i) : i;
      
      /* Get the cell-to-cell coupling terms: */
      for (int f = 0; f < faces.num_faces(i); f++) {
         
         /* Get the index for cell i2 (actual cell or boundary condition): */
         /* Note: boundary conditions have negative, 1-based indexes: */
         int i2 = faces.neighbours(i, f);
         
         /* Get the coefficients (only needed for physical cells): */
         if (i2 >= 0) {
            
            /* Get the distances between the cell centers and the face: */
            double r_i_f = math::distance(faces.centroids(i, f), cells.centroids(i), 3);
            double r_i2_f = math::distance(faces.centroids(i, f), cells.centroids(i2), 3);
            double r_i_i2 = math::distance(cells.centroids(i), cells.centroids(i2), 3);
            
            /* Get the flux weights for the face flux for outgoing directions: */
            coefs(ic, f, 0) = (r_i2_f+delta*r_i_f) / r_i_i2;
            coefs(ic, f, 1) = (1.0-delta)*r_i_f / r_i_i2;
            
            /* Get the flux weights for the face flux for incoming directions: */
            coefs(ic, f, 2) = (1.0-delta)*r_i2_f / r_i_i2;
            coefs(ic, f, 3) = (r_i_f+delta*r_i2_f) / r_i_i2;
            
         }
         
      }
      
   }
   
   return 0;
   
}

/* Build the coefficients for the least-squares gradient-discretization scheme: */
int SNSolver::buildLSGradientScheme(Vector3D<double>& coefs, bool bc) {
   #include "DiffusionSolver.hxx"
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   int num_dims = mesh->getNumDimensions();
   const Cells& cells = mesh->getCells();
   const Faces& faces = mesh->getFaces();
   
   /* Initialize the coefficients: */
   if (bc) {
      Array1D<int> num_faces_bc;
      PAMPA_CALL(getBoundaryCells(num_faces_bc), "unable to get the boundary cells");
      coefs.resize(num_faces_bc.size(), num_faces_bc, 3);
   }
   else
      coefs.resize(num_cells, faces.num_faces, 3);
   
   /* Build the coefficients: */
   for (int i = 0; i < num_cells; i++) {
      
      /* Check if the coefficients for this cell are needed: */
      if (bc && ic_to_ibc(i) < 0) continue;
      int ic = (bc) ? ic_to_ibc(i) : i;
      
      /* Get the centroid and the number of faces for this cell: */
      const double* c_i = cells.centroids(i);
      int num_faces = faces.num_faces(i);
      
      /* Get the d matrix with the cell-to-cell distances for all neighbours: */
      /* Note: for boundary faces the face centroid is used. */
      Array2D<double> d(num_faces, num_dims);
      for (int f = 0; f < num_faces; f++) {
         int i2 = faces.neighbours(i, f);
         const double* c_i2 = (i2 < 0) ? faces.centroids(i, f) : cells.centroids(i2);
         for (int id = 0; id < num_dims; id++)
            d(f, id) = c_i2[id] - c_i[id];
      }
      
      /* Get the G = d^T * d matrix: */
      Eigen::MatrixXd G(num_dims, num_dims);
      for (int jg = 0; jg < num_dims; jg++)
         for (int ig = 0; ig < num_dims; ig++)
            for (int f = 0; f < num_faces; f++)
               G(jg, ig) += d(jg, f) * d(f, ig);
      
      /* Invert the G matrix: */
      Eigen::MatrixXd Ginv(num_dims, num_dims);
      Ginv = G.inverse();
      
      /* Get the M = Ginv * d^T coupling matrix: */
      for (int f = 0; f < num_faces; f++)
         for (int id = 0; id < num_dims; id++)
            for (int jd = 0; jd < num_dims; jd++)
               coefs(ic, f, id) += Ginv(id, jd) * d(jd, f);
      
   }
   
   return 0;
   
}

/* Build the coefficient matrices: */
int SNSolver::buildMatrices() {
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   int num_cells_global = mesh->getNumCellsGlobal();
   int num_faces_max = mesh->getNumFacesMax();
   const Cells& cells = mesh->getCells();
   const Faces& faces = mesh->getFaces();
   const Array1D<BoundaryCondition>& bcs = mesh->getBoundaryConditions();
   
   /* Get the angular quadrature data: */
   int num_directions = quadrature.getNumDirections();
   const Array2D<double>& directions = quadrature.getDirections();
   const Array1D<double>& weights = quadrature.getWeights();
   const Array2D<double>& axes = quadrature.getAxes();
   const Array2D<int>& reflected_directions = quadrature.getReflectedDirections();
   
   /* Create, preallocate and set up the coefficient matrices: */
   int size_local = num_cells * num_groups * num_directions;
   int size_global = num_cells_global * num_groups * num_directions;
   int num_r_nonzero_max = num_groups*num_directions + num_faces_max;
   int num_f_nonzero = num_groups * num_directions;
   petsc::create_matrix(R, size_local, size_global, num_r_nonzero_max);
   petsc::create_matrix(F, size_local, size_global, num_f_nonzero);
   
   /* Initialize the matrix rows for R and F: */
   PetscInt r_l2[num_groups*num_directions];
   PetscInt f_l2[num_groups*num_directions];
   PetscScalar r_l_l2[num_groups*num_directions];
   PetscScalar f_l_l2[num_groups*num_directions];
   
   /* Calculate the coefficients for each cell i: */
   for (int i = 0; i < num_cells; i++) {
      
      /* Get the material for cell i: */
      const Material& mat = materials(cells.materials(i));
      
      /* Calculate the coefficients for each group g: */
      for (int g = 0; g < num_groups; g++) {
         
         /* Calculate the coefficients for each direction m: */
         for (int m = 0; m < num_directions; m++) {
            
            /* Get the matrix index for cell i, group g and direction m: */
            PetscInt l = index(cells.indices(i), g, m, num_groups, num_directions);
            int r_i = 1, f_i = 0;
            
            /* Set the total-reaction term: */
            r_l2[0] = l;
            r_l_l2[0] = mat.sigma_total(g) * cells.volumes(i);
            
            /* Set the group-to-group coupling terms: */
            for (int g2 = 0; g2 < num_groups; g2++) {
               
               /* Set the direction-to-direction coupling terms: */
               for (int m2 = 0; m2 < num_directions; m2++) {
                  
                  /* Get the matrix index for cell i, group g2 and direction m2: */
                  PetscInt l2 = index(cells.indices(i), g2, m2, num_groups, num_directions);
                  
                  /* Set the (g2 -> g, m2 -> m) isotropic scattering term: */
                  if (l2 == l)
                     r_l_l2[0] += -mat.sigma_scattering(g2, g) * weights(m2) * cells.volumes(i);
                  else {
                     r_l2[r_i] = l2;
                     r_l_l2[r_i++] = -mat.sigma_scattering(g2, g) * weights(m2) * cells.volumes(i);
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
               
               /* Get the dot product between the direction and the face normal: */
               double w = math::dot_product(directions(m), faces.normals(i, f), 3);
               
               /* Set the boundary conditions: */
               if (i2 < 0) {
                  
                  /* Check the boundary-condition type: */
                  switch (bcs(-i2).type) {
                     
                     /* Set vacuum (zero-incomming-current) boundary conditions: */
                     case BC::VACUUM : {
                        
                        /* Set the leakage term for outgoing directions: */
                        if (w > 0.0) {
                           
                           /* Set the upwind contribution for cell i: */
                           r_l_l2[0] += w * faces.areas(i, f);
                           
                           /* Set the correction of the face flux using the LS gradient: */
                           if (gradient.boundary_interpolation == BI::LS) {
							  
							  /* Get the vector difference between the face and cell centroids: */
                              double dp[3];
                              math::subtract(dp, faces.centroids(i, f), cells.centroids(i), 3);
                              
                              /* Set the coupling terms for neighboring cells: */
                              for (int f2 = 0; f2 < faces.num_faces(i); f2++) {
                                 
                                 /* Get the index for cell i3: */
                                 int i3 = faces.neighbours(i, f2);
                                 
                                 /* Get the contribution from neighboring physical cells: */
                                 if (i3 >= 0) {
                                    
                                    /* Get the matrix index for cell i3, group g and direction m: */
                                    PetscInt l3 = index(cells.indices(i3), g, m, num_groups, 
                                       num_directions);
                                    
                                    /* Get the boundary-cell index: */
                                    int ibc = ic_to_ibc(i);
                                    
                                    /* Set the LS contribution for cell i: */
                                    double w_i3 = math::dot_product(dp, grad_coefs_bc(ibc, f2), 3);
                                    r_l_l2[0] += -w_i3 * w * faces.areas(i, f);
                                    
                                    /* Set the LS contribution for cell i3: */
                                    double r = w_i3 * w * faces.areas(i, f);
                                    PETSC_CALL(MatSetValues(R, 1, &l, 1, &l3, &r, ADD_VALUES));
                                    
                                 }
                                 
                              }
                              
						   }
                           
                        }
                        
                        break;
                        
                     }
                     
                     /* Set reflective (zero-current) boundary conditions: */
                     /* TODO: the LS correction should be implemented here as well. */
                     case BC::REFLECTIVE : {
                        
                        /* Set the leakage term for outgoing directions: */
                        if (w > 0.0)
                           r_l_l2[0] += w * faces.areas(i, f);
                        
                        /* Set the leakage term for incoming directions: */
                        else {
                           
                           /* Get the reflected outgoing direction: */
                           int m2 = -1;
                           for (int i = 0; i < 3; i++) {
                              double d = math::dot_product(axes(i), faces.normals(i, f), 3);
                              if (fabs(d) > 1.0-TOL)
                                 m2 = reflected_directions(m, i);
                           }
                           PAMPA_CHECK(m2 == -1, 1, "reflected direction not found");
                           
                           /* Get the matrix index for cell i, group g and direction m2: */
                           PetscInt l2 = index(cells.indices(i), g, m2, num_groups, num_directions);
                           
                           /* Set the leakage term for cell i: */
                           double r = w * faces.areas(i, f);
                           PETSC_CALL(MatSetValues(R, 1, &l, 1, &l2, &r, ADD_VALUES));
                           
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
                  PetscInt l2 = index(cells.indices(i2), g, m, num_groups, num_directions);
                  
                  /* Get the flux weights for the face flux: */
                  double w_i, w_i2;
                  if (w > 0.0) {
                     w_i = grad_coefs(i, f, 0);
                     w_i2 = grad_coefs(i, f, 1);
                  }
                  else {
                     w_i = grad_coefs(i, f, 2);
                     w_i2 = grad_coefs(i, f, 3);
                  }
                  
                  /* Set the leakage term for cell i: */
                  r_l_l2[0] += w_i * w * faces.areas(i, f);
                  
                  /* Set the leakage term for cell i2: */
                  double r = w_i2 * w * faces.areas(i, f);
                  PETSC_CALL(MatSetValues(R, 1, &l, 1, &l2, &r, ADD_VALUES));
                  
               }
               
            }
            
            /* Set the matrix rows for R and F: */
            PETSC_CALL(MatSetValues(R, 1, &l, r_i, r_l2, r_l_l2, ADD_VALUES));
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
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   int num_cells_global = mesh->getNumCellsGlobal();
   
   /* Create the scalar-flux vector: */
   PETSC_CALL(VecCreate(MPI_COMM_WORLD, &phi));
   PETSC_CALL(VecSetSizes(phi, num_cells*num_groups, num_cells_global*num_groups));
   PETSC_CALL(VecSetFromOptions(phi));
   
   /* Create the angular-flux vector: */
   PETSC_CALL(MatCreateVecs(R, NULL, &psi));
   
   return 0;
   
}

/* Get the solution after solving the eigensystem: */
int SNSolver::getSolution() {
   
   /* Get the angular flux from the EPS context: */
   double lambda;
   PETSC_CALL(EPSGetEigenpair(eps, 0, &lambda, NULL, psi, NULL));
   keff = 1.0 / lambda;
   
   /* Calculate the scalar flux: */
   PAMPA_CALL(calculateScalarFlux(), "unable to calculate the scalar flux");
   
   /* Normalize the scalar flux: */
   PAMPA_CALL(normalizeScalarFlux(), "unable to normalize the scalar flux");
   
   /* Normalize the angular flux: */
   PAMPA_CALL(normalizeAngularFlux(), "unable to normalize the angular flux");
   
   return 0;
   
}

/* Calculate the scalar flux: */
int SNSolver::calculateScalarFlux() {
   
   /* Get the number of cells: */
   int num_cells = mesh->getNumCells();
   
   /* Get the angular quadrature data: */
   int num_directions = quadrature.getNumDirections();
   const Array1D<double>& weights = quadrature.getWeights();
   
   /* Get the arrays for the scalar and angular fluxes: */
   PetscScalar *data_phi, *data_psi;
   PETSC_CALL(VecGetArray(phi, &data_phi));
   PETSC_CALL(VecGetArray(psi, &data_psi));
   
   /* Integrate the angular flux over all directions: */
   for (int iphi = 0, ipsi = 0, i = 0; i < num_cells; i++) {
      for (int g = 0; g < num_groups; g++) {
         data_phi[iphi] = 0.0;
         for (int m = 0; m < num_directions; m++)
            data_phi[iphi] += weights(m) * data_psi[ipsi++];
         data_phi[iphi] *= 4.0 * M_PI;
         iphi++;
      }
   }
   
   /* Restore the arrays for the scalar and angular fluxes: */
   PETSC_CALL(VecRestoreArray(phi, &data_phi));
   PETSC_CALL(VecRestoreArray(psi, &data_psi));
   
   return 0;
   
}

/* Normalize the angular flux: */
int SNSolver::normalizeAngularFlux() {
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   const Cells& cells = mesh->getCells();
   
   /* Get the angular quadrature data: */
   int num_directions = quadrature.getNumDirections();
   const Array1D<double>& weights = quadrature.getWeights();
   
   /* Get the array for the angular flux: */
   PetscScalar* data_psi;
   PETSC_CALL(VecGetArray(psi, &data_psi));
   
   /* Normalize the angular flux: */
   /* TODO: normalize correctly with the power. */
   double vol = 0.0;
   for (int i = 0; i < num_cells; i++)
      if (materials(cells.materials(i)).nu_sigma_fission(1) > 0.0)
         vol += cells.volumes(i);
   MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
   double sum = 0.0;
   for (int ipsi = 0, i = 0; i < num_cells; i++)
      for (int g = 0; g < num_groups; g++)
         for (int m = 0; m < num_directions; m++)
            sum += weights(m) * data_psi[ipsi++] * 
                      materials(cells.materials(i)).nu_sigma_fission(g) * cells.volumes(i);
   MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
   double f = vol / sum;
   for (int ipsi = 0, i = 0; i < num_cells; i++)
      for (int g = 0; g < num_groups; g++)
         for (int m = 0; m < num_directions; m++)
            data_psi[ipsi++] *= f;
   
   /* Check for negative fluxes: */
   for (int ipsi = 0, i = 0; i < num_cells; i++) {
      for (int g = 0; g < num_groups; g++) {
         for (int m = 0; m < num_directions; m++) {
            PAMPA_CHECK(data_psi[ipsi++] < 0.0, 1, "negative values in the angular-flux solution");
         }
      }
   }
   
   /* Restore the array for the angular flux: */
   PETSC_CALL(VecRestoreArray(psi, &data_psi));
   
   return 0;
   
}

/* Write the solution to a plain-text file in .vtk format: */
int SNSolver::writeVTK(const std::string& filename) const {
   
   /* Get the number of cells: */
   int num_cells = mesh->getNumCells();
   
   /* Get the number of directions: */
   int num_directions = quadrature.getNumDirections();
   
   /* Get the arrays for the scalar and angular fluxes: */
   PetscScalar *data_phi, *data_psi;
   PETSC_CALL(VecGetArray(phi, &data_phi));
   PETSC_CALL(VecGetArray(psi, &data_psi));
   
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
         file << data_phi[index(i, g, num_groups)] << std::endl;
      file << std::endl;
   }
   
   /* Write the angular flux: */
   for (int g = 0; g < num_groups; g++) {
      for (int m = 0; m < num_directions; m++) {
         file << "SCALARS flux_" << (g+1) << "_" << (m+1) << " double 1" << std::endl;
         file << "LOOKUP_TABLE default" << std::endl;
         for (int i = 0; i < num_cells; i++)
            file << data_psi[index(i, g, m, num_groups, num_directions)] << std::endl;
         file << std::endl;
      }
   }
   
   /* Restore the arrays for the scalar and angular fluxes: */
   PETSC_CALL(VecRestoreArray(phi, &data_phi));
   PETSC_CALL(VecRestoreArray(psi, &data_psi));
   
   return 0;
   
}

/* Write the solution to a binary file in PETSc format: */
int SNSolver::writePETSc(const std::string& filename) const {
   
   /* Write the solution to a binary file: */
   PetscViewer viewer;
   PETSC_CALL(PetscViewerBinaryOpen(MPI_COMM_WORLD, "flux.ptc", FILE_MODE_WRITE, &viewer));
   PETSC_CALL(VecView(psi, viewer));
   PETSC_CALL(PetscViewerDestroy(&viewer));
   
   return 0;
   
}

/* Destroy the solution vectors: */
int SNSolver::destroyVectors() {
   
   /* Destroy the solution vectors: */
   PETSC_CALL(VecDestroy(&phi));
   PETSC_CALL(VecDestroy(&psi));
   
   return 0;
   
}
