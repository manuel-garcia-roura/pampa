#include "Solver.hxx"

/* Initialize: */
int Solver::initialize(int argc, char* argv[], const Model& model) {
   
   /* Get the transport method: */
   const TransportMethod& method = model.getTransportMethod();
   
   /* Initialize SLEPc: */
   static char help[] = "Solver for the generalized eigensystem R*x = (1/keff)*F*x.\n";
   PETSC_CALL(SlepcInitialize(&argc, &argv, (char*)0, help));
   
   /* Create the coefficient matrices: */
   PETSC_CALL(MatCreate(MPI_COMM_WORLD, &R));
   PETSC_CALL(MatSetFromOptions(R));
   PETSC_CALL(MatCreate(MPI_COMM_WORLD, &F));
   PETSC_CALL(MatSetFromOptions(F));
   
   /* Build the coefficient matrices depending on the transport method: */
   switch (method.type) {
      case TM::DIFFUSION : {
         PAMPA_CALL(buildDiffusionMatrices(model), "unable to build the R and F matrices");
         break;
      }
      case TM::SN : {
         PAMPA_CALL(buildSNMatrices(model), "unable to build the R and F matrices");
         break;
      }
   }
   
   /* Create the EPS context: */
   PETSC_CALL(EPSCreate(MPI_COMM_WORLD, &eps));
   PETSC_CALL(EPSSetOperators(eps, R, F));
   PETSC_CALL(EPSSetFromOptions(eps));
   
   /* Get the initial condition, if given: */
   char filename[PETSC_MAX_PATH_LEN];
   PetscBool flag;
   PetscViewer viewer;
   PETSC_CALL(PetscOptionsGetString(NULL, NULL, "-petsc_initial_condition", filename, 
      sizeof(filename), &flag));
   if (flag) {
      PETSC_CALL(PetscViewerBinaryOpen(MPI_COMM_WORLD, filename, FILE_MODE_READ, &viewer));
      switch (method.type) {
         case TM::DIFFUSION : {
            PETSC_CALL(VecLoad(phi_mpi, viewer));
            PETSC_CALL(EPSSetInitialSpace(eps, 1, &phi_mpi));
            break;
         }
         case TM::SN : {
            PETSC_CALL(VecLoad(psi_mpi, viewer));
            PETSC_CALL(EPSSetInitialSpace(eps, 1, &psi_mpi));
            break;
         }
      }
      PETSC_CALL(PetscViewerDestroy(&viewer));
   }
   
   return 0;
   
}

/* Solve the eigensystem to get the neutron flux and the multiplication factor: */
int Solver::solve() {
   
   /* Solve the eigensystem: */
   double t1 = MPI_Wtime();
   PETSC_CALL(EPSSolve(eps));
   double t2 = MPI_Wtime();
   
   /* Get the solver information: */
   ST st;
   KSP ksp;
   EPSType eps_type;
   PetscInt num_eigenvalues, max_eps_iterations, num_eps_iterations, num_ksp_iterations;
   PetscReal eps_tol;
   PETSC_CALL(EPSGetST(eps, &st));
   PETSC_CALL(STGetKSP(st, &ksp));
   PETSC_CALL(EPSGetType(eps, &eps_type));
   PETSC_CALL(EPSGetTolerances(eps, &eps_tol, &max_eps_iterations));
   PETSC_CALL(EPSGetIterationNumber(eps, &num_eps_iterations));
   PETSC_CALL(KSPGetTotalIterations(ksp, &num_ksp_iterations));
   
   /* Print out the solver information: */
   if (mpi::rank == 0) {
      std::cout << "Elapsed time: " << t2-t1 << std::endl;
      std::cout << "Solution method: " << eps_type << std::endl;
      std::cout << "EPS convergence tolerance: " << eps_tol << std::endl;
      std::cout << "Maximum number of EPS iterations: " << max_eps_iterations << std::endl;
      std::cout << "Number of EPS iterations: " << num_eps_iterations << std::endl;
      std::cout << "Number of KSP iterations: " << num_ksp_iterations << std::endl;
   }
   
   return 0;
   
}

/* Output the solution: */
int Solver::output(const std::string& filename, const Model& model) {
   
   /* Get the model data: */
   const TransportMethod& method = model.getTransportMethod();
   const Mesh* mesh = model.getMesh();
   const AngularQuadratureSet& quadrature = model.getAngularQuadratureSet();
   
   /* Get the number of cells: */
   int num_cells = mesh->getNumCells();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Get the number of directions: */
   int num_directions = quadrature.getNumDirections();
   
   /* Get the solution: */
   double lambda;
   switch (method.type) {
      case TM::DIFFUSION : {
         PETSC_CALL(EPSGetEigenpair(eps, 0, &lambda, NULL, phi_mpi, NULL));
         break;
      }
      case TM::SN : {
         PETSC_CALL(EPSGetEigenpair(eps, 0, &lambda, NULL, psi_mpi, NULL));
         break;
      }
   }
   keff = 1.0 / lambda;
   
   /* Gather the solution from all ranks to the master rank: */
   PAMPA_CALL(gatherSolution(model), "unable to gather the flux");
   
   /* Calculate the scalar flux: */
   if (method.type == TM::SN)
      PAMPA_CALL(calculateScalarFlux(model), "unable to calculate the scalar flux");
   
   /* Normalize the flux: */
   PAMPA_CALL(normalizeFlux(model), "unable to normalize the flux");
   
   /* Get the raw data for the scalar and angular fluxes: */
   PetscScalar *data_phi, *data_psi;
   PETSC_CALL(VecGetArray(phi_seq, &data_phi));
   if (method.type == TM::SN) {
      PETSC_CALL(VecGetArray(psi_seq, &data_psi));
   }
   
   /* Check the MPI rank: */
   if (mpi::rank == 0) {
      
      /* Print out the multiplication factor: */
      std::cout << "keff = " << keff << std::endl;
      
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
      if (method.type == TM::SN) {
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
      
   }
   
   /* Restore the array for the scalar and angular fluxes: */
   PETSC_CALL(VecRestoreArray(phi_seq, &data_phi));
   if (method.type == TM::SN) {
      PETSC_CALL(VecRestoreArray(psi_seq, &data_psi));
   }
   
   /* Output the solution in PETSc format: */
   PetscViewer viewer;
   PETSC_CALL(PetscViewerBinaryOpen(MPI_COMM_WORLD, "flux.ptc", FILE_MODE_WRITE, &viewer));
   switch (method.type) {
      case TM::DIFFUSION : {
         PETSC_CALL(VecView(phi_mpi, viewer));
         break;
      }
      case TM::SN : {
         PETSC_CALL(VecView(psi_mpi, viewer));
         break;
      }
   }
   PETSC_CALL(PetscViewerDestroy(&viewer));
   
   return 0;
   
}

/* Finalize: */
int Solver::finalize(const Model& model) {
   
   /* Get the transport method: */
   const TransportMethod& method = model.getTransportMethod();
   
   /* Destroy the EPS context: */
   PETSC_CALL(EPSDestroy(&eps));
   
   /* Destroy the solution vector: */
   PETSC_CALL(VecDestroy(&phi_mpi));
   if (method.type == TM::SN) {
      PETSC_CALL(VecDestroy(&psi_mpi));
   }
   
   /* Destroy the coefficient matrices: */
   PETSC_CALL(MatDestroy(&R));
   PETSC_CALL(MatDestroy(&F));
   
   /* Finalize SLEPc: */
   PETSC_CALL(SlepcFinalize());
   
   return 0;
   
}

/* Build the coefficient matrices for the diffusion method: */
int Solver::buildDiffusionMatrices(const Model& model) {
   
   /* Get the model data: */
   const TransportMethod& method = model.getTransportMethod();
   const Mesh* mesh = model.getMesh();
   const std::vector<Material>& materials = model.getMaterials();
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   const Cells& cells = mesh->getCells();
   const Faces& faces = mesh->getFaces();
   const std::vector<BoundaryCondition>& bcs = mesh->getBoundaryConditions();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Set up the coefficient matrices: */
   int n = num_cells * num_groups;
   PETSC_CALL(MatSetSizes(R, PETSC_DECIDE, PETSC_DECIDE, n, n));
   PETSC_CALL(MatSeqAIJSetPreallocation(R, 6+num_groups, NULL));
   PETSC_CALL(MatMPIAIJSetPreallocation(R, 6+num_groups, NULL, 6+num_groups, NULL));
   PETSC_CALL(MatSetUp(R));
   PETSC_CALL(MatSetSizes(F, PETSC_DECIDE, PETSC_DECIDE, n, n));
   PETSC_CALL(MatSeqAIJSetPreallocation(F, num_groups, NULL));
   PETSC_CALL(MatMPIAIJSetPreallocation(F, num_groups, NULL, num_groups, NULL));
   PETSC_CALL(MatSetUp(F));
   
   /* Get the local ownership range: */
   int r_l1, r_l2, f_l1, f_l2;
   MatGetOwnershipRange(R, &r_l1, &r_l2);
   MatGetOwnershipRange(F, &f_l1, &f_l2);
   PAMPA_CHECK((f_l1 != r_l1) || (f_l2 != r_l2), 1, "wrong local ownership range");
   
   /* Calculate the coefficients for each cell i: */
   for (int i = 0; i < num_cells; i++) {
      
      /* Get the material for cell i: */
      int mat = cells.materials[i];
      
      /* Calculate the coefficients for each group g: */
      for (int g = 0; g < num_groups; g++) {
         
         /* Get the matrix index for cell i and group g: */
         int l = i*num_groups + g;
         
         /* Check if the matrix index is local: */
         if ((l >= r_l2) || (l < r_l1))
            continue;
         
         /* Set the total-reaction term: */
         double r_l_l = materials[mat].sigma_total[g] * cells.volumes[i];
         
         /* Set the group-to-group coupling terms: */
         for (int g2 = 0; g2 < num_groups; g2++) {
            
            /* Get the matrix index for cell i and group g2: */
            int l2 = i*num_groups + g2;
            
            /* Set the (g2 -> g) scattering term: */
            double r_l_l2 = -materials[mat].sigma_scattering[g2][g] * cells.volumes[i];
            if (l2 == l)
               r_l_l += r_l_l2;
            else {
               PETSC_CALL(MatSetValues(R, 1, &l, 1, &l2, &r_l_l2, INSERT_VALUES));
            }
            
            /* Set the (g2 -> g) fission term: */
            double f_l_l2 = materials[mat].chi[g] * materials[mat].nu_sigma_fission[g2] * 
                               cells.volumes[i];
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
                     r_l_l += w * materials[mat].diffusion_coefficient[g] * faces.areas[i][f];
                     
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
               int mat2 = cells.materials[i2];
               
               /* Set the terms for cells with the same materials: */
               if (mat2 == mat) {
                  
                  /* Get the geometrical data: */
                  const std::vector<double>& p_i = cells.centroids[i];
                  const std::vector<double>& p_i2 = cells.centroids[i2];
                  const std::vector<double>& n_i_f = faces.normals[i][f];
                  
                  /* Get the surface leakage factor: */
                  double w = math::surface_leakage_factor(p_i, p_i2, n_i_f);
                  
                  /* Get the leakage term for cell i2: */
                  r_l_l2 = -w * materials[mat].diffusion_coefficient[g] * faces.areas[i][f];
                  
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
                  w_i_i2 *= materials[mat].diffusion_coefficient[g] * faces.areas[i][f];
                  
                  /* Get the surface leakage factor and the weight for cell i2: */
                  double w_i2_i = math::surface_leakage_factor(p_i2, p_f, n_i_f);
                  w_i2_i *= -materials[mat2].diffusion_coefficient[g] * faces.areas[i][f];
                  
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
   PETSC_CALL(MatAssemblyEnd(R, MAT_FINAL_ASSEMBLY));
   PETSC_CALL(MatAssemblyBegin(F, MAT_FINAL_ASSEMBLY));
   PETSC_CALL(MatAssemblyEnd(F, MAT_FINAL_ASSEMBLY));
   
   /* Create the solution vector: */
   PETSC_CALL(MatCreateVecs(R, NULL, &phi_mpi));
   
   return 0;
   
}

/* Build the coefficient matrices for the SN method: */
int Solver::buildSNMatrices(const Model& model) {
   
   /* Get the model data: */
   const TransportMethod& method = model.getTransportMethod();
   const Mesh* mesh = model.getMesh();
   const AngularQuadratureSet& quadrature = model.getAngularQuadratureSet();
   const std::vector<Material>& materials = model.getMaterials();
   
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
   
   /* Set up the coefficient matrices: */
   int n = num_cells * num_groups * num_directions;
   int max_num_columns_r = 6 + num_groups*num_directions;
   int max_num_columns_f = num_groups * num_directions;
   PETSC_CALL(MatSetSizes(R, PETSC_DECIDE, PETSC_DECIDE, n, n));
   PETSC_CALL(MatSeqAIJSetPreallocation(R, max_num_columns_r, NULL));
   PETSC_CALL(MatMPIAIJSetPreallocation(R, max_num_columns_r, NULL, max_num_columns_r, NULL));
   PETSC_CALL(MatSetUp(R));
   PETSC_CALL(MatSetSizes(F, PETSC_DECIDE, PETSC_DECIDE, n, n));
   PETSC_CALL(MatSeqAIJSetPreallocation(F, max_num_columns_f, NULL));
   PETSC_CALL(MatMPIAIJSetPreallocation(F, max_num_columns_f, NULL, max_num_columns_f, NULL));
   PETSC_CALL(MatSetUp(F));
   
   /* Get the local ownership range: */
   int r_l1, r_l2, f_l1, f_l2;
   MatGetOwnershipRange(R, &r_l1, &r_l2);
   MatGetOwnershipRange(F, &f_l1, &f_l2);
   PAMPA_CHECK((f_l1 != r_l1) || (f_l2 != r_l2), 1, "wrong local ownership range");
   
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
            if ((l >= r_l2) || (l < r_l1))
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
   PETSC_CALL(MatAssemblyEnd(R, MAT_FINAL_ASSEMBLY));
   PETSC_CALL(MatAssemblyBegin(F, MAT_FINAL_ASSEMBLY));
   PETSC_CALL(MatAssemblyEnd(F, MAT_FINAL_ASSEMBLY));
   
   /* Create the solution vector: */
   PETSC_CALL(MatCreateVecs(R, NULL, &psi_mpi));
   
   /* Create the vector for the scalar flux: */
   n = num_cells * num_groups;
   PETSC_CALL(VecCreate(MPI_COMM_WORLD, &phi_mpi));
   PETSC_CALL(VecSetSizes(phi_mpi, PETSC_DECIDE, n));
   PETSC_CALL(VecSetFromOptions(phi_mpi));
   
   return 0;
   
}

/* Gather the solution from all ranks to the master rank: */
int Solver::gatherSolution(const Model& model) {
   
   /* Get the model data: */
   const TransportMethod& method = model.getTransportMethod();
   const Mesh* mesh = model.getMesh();
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Gather the scalar flux: */
   VecScatter phi_context;
   PETSC_CALL(VecCreateSeq(MPI_COMM_SELF, num_cells*num_groups, &phi_seq));
   PETSC_CALL(VecZeroEntries(phi_seq));
   PETSC_CALL(VecScatterCreateToAll(phi_mpi, &phi_context, &phi_seq));
   PETSC_CALL(VecScatterBegin(phi_context, phi_mpi, phi_seq, INSERT_VALUES, SCATTER_FORWARD));
   PETSC_CALL(VecScatterEnd(phi_context, phi_mpi, phi_seq, INSERT_VALUES, SCATTER_FORWARD));
   PETSC_CALL(VecScatterDestroy(&phi_context));
   
   /* Gather the angular flux: */
   if (method.type == TM::SN) {
      VecScatter psi_context;
      PETSC_CALL(VecCreateSeq(MPI_COMM_SELF, num_cells*num_groups, &psi_seq));
      PETSC_CALL(VecZeroEntries(psi_seq));
      PETSC_CALL(VecScatterCreateToAll(psi_mpi, &psi_context, &psi_seq));
      PETSC_CALL(VecScatterBegin(psi_context, psi_mpi, psi_seq, INSERT_VALUES, SCATTER_FORWARD));
      PETSC_CALL(VecScatterEnd(psi_context, psi_mpi, psi_seq, INSERT_VALUES, SCATTER_FORWARD));
      PETSC_CALL(VecScatterDestroy(&psi_context));
   }
   
   return 0;
   
}

/* Calculate the scalar flux: */
int Solver::calculateScalarFlux(const Model& model) {
   
   /* Get the model data: */
   const TransportMethod& method = model.getTransportMethod();
   const Mesh* mesh = model.getMesh();
   const AngularQuadratureSet& quadrature = model.getAngularQuadratureSet();
   
   /* Get the number of cells: */
   int num_cells = mesh->getNumCells();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Get the angular quadrature data: */
   int num_directions = quadrature.getNumDirections();
   const std::vector<double>& weights = quadrature.getWeights();
   
   /* Get the raw data for the scalar and angular fluxes: */
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

/* Normalize the flux: */
int Solver::normalizeFlux(const Model& model) {
   
   /* Get the model data: */
   const TransportMethod& method = model.getTransportMethod();
   const Mesh* mesh = model.getMesh();
   const AngularQuadratureSet& quadrature = model.getAngularQuadratureSet();
   const std::vector<Material>& materials = model.getMaterials();
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   const Cells& cells = mesh->getCells();
   
   /* Get the number of energy groups: */
   int num_groups = method.num_groups;
   
   /* Get the angular quadrature data: */
   int num_directions = quadrature.getNumDirections();
   const std::vector<double>& weights = quadrature.getWeights();
   
   /* Get the raw data_phi for the scalar and angular fluxes: */
   PetscScalar *data_phi, *data_psi;
   PETSC_CALL(VecGetArray(phi_seq, &data_phi));
   if (method.type == TM::SN) {
      PETSC_CALL(VecGetArray(psi_seq, &data_psi));
   }
   
   /* Normalize the scalar flux (TODO: normalize correctly with the power!): */
   double vol = 0.0;
   for (int i = 0; i < num_cells; i++)
      if (materials[cells.materials[i]].nu_sigma_fission[1] > 0.0)
         vol += cells.volumes[i];
   double sum = 0.0;
   for (int g = 0; g < num_groups; g++)
      for (int i = 0; i < num_cells; i++)
         sum += data_phi[i*num_groups+g] * materials[cells.materials[i]].nu_sigma_fission[g] * 
                   cells.volumes[i];
   double f = vol / sum;
   for (int g = 0; g < num_groups; g++)
      for (int i = 0; i < num_cells; i++)
         data_phi[i*num_groups+g] *= f;
   if (method.type == TM::SN) {
      f *= 4.0 * M_PI;
      for (int g = 0; g < num_groups; g++)
         for (int i = 0; i < num_cells; i++)
            for (int m = 0; m < num_directions; m++)
               data_psi[i*num_directions*num_groups+g*num_directions+m] *= f;
   }
   
   /* Restore the arrays for the scalar and angular fluxes: */
   PETSC_CALL(VecRestoreArray(phi_seq, &data_phi));
   if (method.type == TM::SN) {
      PETSC_CALL(VecRestoreArray(psi_seq, &data_psi));
   }
   
   return 0;
   
}
