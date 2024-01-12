#include "Solver.hxx"

/* Initialize: */
int Solver::initialize(int argc, char* argv[], const Model &model) {
   
   /* Initialize SLEPc: */
   static char help[] = "Solver for the generalized eigensystem R*phi = (1/keff)*F*phi.\n";
   PETSC_CALL(SlepcInitialize(&argc, &argv, (char*)0, help), "unable to initialize PETSc/SLEPc");
   
   /* Create the coefficient matrices: */
   PETSC_CALL(MatCreate(PETSC_COMM_WORLD, &R), "unable to create the R matrix");
   PETSC_CALL(MatSetFromOptions(R), "unable to set the R options");
   PETSC_CALL(MatCreate(PETSC_COMM_WORLD, &F), "unable to create the F matrix");
   PETSC_CALL(MatSetFromOptions(F), "unable to set the F options");
   
   /* Build the coefficient matrices: */
   PAMPA_CALL(buildMatrices(model), "unable to build the R and F matrices");
   
   /* Create the solution vector: */
   PETSC_CALL(MatCreateVecs(R, NULL, &phi), "unable to create the solution vector");
   
   /* Create the EPS context: */
   PETSC_CALL(EPSCreate(PETSC_COMM_WORLD, &eps), "unable to create the EPS");
   PETSC_CALL(EPSSetOperators(eps, R, F), "unable to set the EPS operators");
   PETSC_CALL(EPSSetFromOptions(eps), "unable to set the EPS options");
   
   /* Get the initial condition, if given: */
   char filename[PETSC_MAX_PATH_LEN];
   PetscBool flag;
   PetscViewer viewer;
   PETSC_CALL(PetscOptionsGetString(NULL, NULL, "-petsc_initial_condition", filename, 
      sizeof(filename), &flag), "unable to get the PETSc initial condition");
   if (flag) {
      PETSC_CALL(PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer), 
         "unable to open the binary file");
      PETSC_CALL(VecDuplicate(phi, &phi0), "unable to create the initial-condition vector");
      PETSC_CALL(VecLoad(phi0, viewer), "unable to read the initial condition");
      PETSC_CALL(PetscViewerDestroy(&viewer), "unable to close the binary file");
      PETSC_CALL(EPSSetInitialSpace(eps, 1, &phi0), "unable to set the initial condition");
   }
   
   return 0;
   
}

/* Solve the eigensystem to get the neutron flux and the multiplication factor: */
int Solver::solve() {
   
   /* Solve the eigensystem: */
   PETSC_CALL(EPSSolve(eps), "unable to solve the eigensystem");
   
   /* Get the solver information: */
   ST st;
   KSP ksp;
   EPSType eps_type;
   PetscInt num_eigenvalues, max_eps_iterations, num_eps_iterations, num_ksp_iterations;
   PetscReal eps_tol;
   PETSC_CALL(EPSGetST(eps, &st), "unable to get the ST context");
   PETSC_CALL(STGetKSP(st, &ksp), "unable to get the KSP context");
   PETSC_CALL(EPSGetType(eps, &eps_type), "unable to get the EPS type");
   PETSC_CALL(EPSGetTolerances(eps, &eps_tol, &max_eps_iterations), 
      "unable to get the EPS convergence tolerance and maximum number of iterations");
   PETSC_CALL(EPSGetIterationNumber(eps, &num_eps_iterations), 
      "unable to get the EPS number of iterations");
   PETSC_CALL(KSPGetTotalIterations(ksp, &num_ksp_iterations), 
      "unable to get the KSP number of iterations");
   
   /* Print out the solver information: */
   std::cout << "Solution method: " << eps_type << std::endl;
   std::cout << "EPS convergence tolerance: " << eps_tol << std::endl;
   std::cout << "Maximum number of EPS iterations: " << max_eps_iterations << std::endl;
   std::cout << "Number of EPS iterations: " << num_eps_iterations << std::endl;
   std::cout << "Number of KSP iterations: " << num_ksp_iterations << std::endl;
   
   return 0;
   
}

/* Output the solution: */
int Solver::output(const std::string &filename, const Model &model) {
   
   /* Get the model data: */
   int num_groups = model.getNumEnergyGroups();
   const Mesh *mesh = model.getMesh();
   const std::vector<Material>& materials = model.getMaterials();
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   const Cells& cells = mesh->getCells();
   
   /* Print out the multiplication factor: */
   double lambda;
   PETSC_CALL(EPSGetEigenpair(eps, 0, &lambda, NULL, phi, NULL), "");
   keff = 1.0 / lambda;
   std::cout << "keff = " << keff << std::endl;
   
   /* Get the raw data to the solution vector: */
   PetscScalar *data;
   PETSC_CALL(VecGetArray(phi, &data), "unable to get the solution array");
   
   /* Normalize the solution to one (TODO: normalize correctly with the power!): */
   double vol = 0.0;
   for (int i = 0; i < num_cells; i++)
      if (materials[cells.materials[i]].nu_sigma_fission[1] > 0.0)
         vol += cells.volumes[i];
   double sum = 0.0;
   for (int g = 0; g < num_groups; g++)
      for (int i = 0; i < num_cells; i++)
         sum += data[i*num_groups+g] * materials[cells.materials[i]].nu_sigma_fission[g] * 
                   cells.volumes[i];
   double f = vol / sum;
   for (int g = 0; g < num_groups; g++)
      for (int i = 0; i < num_cells; i++)
         data[i*num_groups+g] *= f;
   
   /* Write the mesh: */
   PAMPA_CALL(mesh->writeVTK(filename), "unable to write the mesh");
   
   /* Open the output file: */
   std::ofstream file(filename, std::ios_base::app);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Write the solution: */
   for (int g = 0; g < num_groups; g++) {
      file << "SCALARS flux" << (g+1) << " double 1" << std::endl;
      file << "LOOKUP_TABLE default" << std::endl;
      for (int i = 0; i < num_cells; i++)
         file << data[i*num_groups+g] << std::endl;
      file << std::endl;
   }
   
   /* Restore the solution vector: */
   PETSC_CALL(VecRestoreArray(phi, &data), "unable to restore the solution array");
   
   /* Output the solution in PETSc format: */
   PetscViewer viewer;
   PETSC_CALL(PetscViewerBinaryOpen(PETSC_COMM_WORLD, "flux.ptc", FILE_MODE_WRITE, &viewer), 
      "unable to open the binary file");
   PETSC_CALL(VecView(phi, viewer), "unable to write the solution vector");
   PETSC_CALL(PetscViewerDestroy(&viewer), "unable to close the binary file");
   
   return 0;
   
}

/* Finalize: */
int Solver::finalize() {
   
   /* Destroy the EPS context: */
   PETSC_CALL(EPSDestroy(&eps), "unable to destroy the EPS");
   
   /* Destroy the solution vector: */
   PETSC_CALL(VecDestroy(&phi), "unable to destroy the solution vector");
   
   /* Destroy the coefficient matrices: */
   PETSC_CALL(MatDestroy(&R), "unable to destroy the R matrix");
   PETSC_CALL(MatDestroy(&F), "unable to destroy the F matrix");
   
   /* Finalize SLEPc: */
   PETSC_CALL(SlepcFinalize(), "unable to finalize PETSc/SLEPc");
   
   return 0;
   
}

/* Build the coefficient matrices: */
int Solver::buildMatrices(const Model &model) {
   
   /* Get the model data: */
   int num_groups = model.getNumEnergyGroups();
   const Mesh *mesh = model.getMesh();
   const std::vector<Material>& materials = model.getMaterials();
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   const Cells& cells = mesh->getCells();
   const Faces& faces = mesh->getFaces();
   const std::vector<BoundaryCondition> &bcs = mesh->getBoundaryConditions();
   
   /* Set up the coefficient matrices: */
   int n = num_cells * num_groups;
   PETSC_CALL(MatSetSizes(R, PETSC_DECIDE, PETSC_DECIDE, n, n), "unable to set the R size");
   PETSC_CALL(MatSeqAIJSetPreallocation(R, 6+num_groups, NULL), "unable to preallocate R");
   PETSC_CALL(MatSetUp(R), "unable to set up R");
   PETSC_CALL(MatSetSizes(F, PETSC_DECIDE, PETSC_DECIDE, n, n), "unable to set the F size");
   PETSC_CALL(MatSeqAIJSetPreallocation(F, num_groups, NULL), "unable to preallocate F");
   PETSC_CALL(MatSetUp(F), "unable to set up F");
   
   /* Calculate the coefficients for each cell i: */
   for (int i = 0; i < num_cells; i++) {
      
      /* Get the material for cell i: */
      int mat = cells.materials[i];
      
      /* Calculate the coefficients for each group g: */
      for (int g = 0; g < num_groups; g++) {
         
         /* Get the matrix index for cell i and group g: */
         int l = i*num_groups + g;
         
         /* Set the total-reaction term: */
         double r_l_l = materials[mat].sigma_removal[g] * cells.volumes[i];
         
         /* Set the group-to-group coupling terms: */
         for (int g2 = 0; g2 < num_groups; g2++) {
            
            /* Get the matrix index for cell i and group g2: */
            int l2 = i*num_groups + g2;
            
            /* Set the (g2 -> g) scattering term: */
            if (g2 != g) {
               double r_l_l2 = -materials[mat].sigma_scattering[g2][g] * cells.volumes[i];
               PETSC_CALL(MatSetValues(R, 1, &l, 1, &l2, &r_l_l2, INSERT_VALUES), 
                  "unable to set R(l, l2)");
            }
            
            /* Set the (g2 -> g) fission term: */
            double f_l_l2 = materials[mat].chi[g] * materials[mat].nu_sigma_fission[g2] * 
                               cells.volumes[i];
            PETSC_CALL(MatSetValues(F, 1, &l, 1, &l2, &f_l_l2, INSERT_VALUES), 
               "unable to set F(l, l2)");
            
         }
         
         /* Set the cell-to-cell coupling terms: */
         double r_l_l2;
         for (int f = 0; f < faces.neighbours[i].size(); f++) {
            
            /* Get the index for cell i2 (actual cell or boundary condition): */
            /* Note: boundary conditions have negative, 1-based indexes: */
            int i2 = faces.neighbours[i][f];
            
            /* Set boundary conditions: */
            if (i2 < 0) {
               
               /* Check the boundary-condition type: */
               switch(bcs[-i2].type) {
                  
                  /* Set vacuum (zero-flux) boundary conditions: */
                  case bc::VACUUM : {
                     
                     /* Get the geometrical data: */
                     const std::vector<double> &p_i = cells.centroids[i];
                     const std::vector<double> &p_f = faces.centroids[i][f];
                     const std::vector<double> &n_i_f = faces.normals[i][f];
                     
                     /* Get the surface leakage factor: */
                     double w = math::surface_leakage_factor(p_i, p_f, n_i_f);
                     
                     /* Set the leakage term for cell i: */
                     r_l_l += w * materials[mat].diffusion_coefficient[g] * faces.areas[i][f];
                     
                  }
                  
                  /* Set reflective (zero-current) boundary conditions: */
                  case bc::REFLECTIVE : {
                     
                     /* Nothing to be done: */
                     continue;
                     
                  }
                  
                  /* Set Robin boundary conditions: */
                  case bc::ROBIN : {
                     
                     /* Set the leakage term for cell i: */
                     r_l_l -= bcs[-i2].a * faces.areas[i][f];
                     
                  }
                  
               }
               
            }
            
            /* Set cell-to-cell coupling terms depending on the neighbour material: */
            else {
               
               /* Get the matrix index for cell i2 and group g: */
               int l2 = i2*num_groups + g;
               
               /* Get the material for cell i2: */
               int mat2 = cells.materials[i2];
               
               /* Set the terms for cells with the same materials: */
               if (mat2 == mat) {
                  
                  /* Get the geometrical data: */
                  const std::vector<double> &p_i = cells.centroids[i];
                  const std::vector<double> &p_i2 = cells.centroids[i2];
                  const std::vector<double> &n_i_f = faces.normals[i][f];
                  
                  /* Get the surface leakage factor: */
                  double w = math::surface_leakage_factor(p_i, p_i2, n_i_f);
                  
                  /* Get the leakage term for cell i2: */
                  r_l_l2 = -w * materials[mat].diffusion_coefficient[g] * faces.areas[i][f];
                  
               }
               
               /* Set the terms for cells with different materials: */
               else {
                  
                  /* Get the geometrical data: */
                  const std::vector<double> &p_i = cells.centroids[i];
                  const std::vector<double> &p_i2 = cells.centroids[i2];
                  const std::vector<double> &p_f = faces.centroids[i][f];
                  const std::vector<double> &n_i_f = faces.normals[i][f];
                  
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
               PETSC_CALL(MatSetValues(R, 1, &l, 1, &l2, &r_l_l2, INSERT_VALUES), 
                  "unable to set R(l, l2)");
               
            }
            
         }
         
         /* Set the diagonal coefficient: */
         PETSC_CALL(MatSetValues(R, 1, &l, 1, &l, &r_l_l, INSERT_VALUES), "unable to set R(l, l)");
         
      }
   }
   
   /* Assembly the coefficient matrices: */
   PETSC_CALL(MatAssemblyBegin(R, MAT_FINAL_ASSEMBLY), "unable to assembly the R matrix");
   PETSC_CALL(MatAssemblyEnd(R, MAT_FINAL_ASSEMBLY), "unable to assembly the R matrix");
   PETSC_CALL(MatAssemblyBegin(F, MAT_FINAL_ASSEMBLY), "unable to assembly the F matrix");
   PETSC_CALL(MatAssemblyEnd(F, MAT_FINAL_ASSEMBLY), "unable to assembly the F matrix");
   
   return 0;
   
}
