#include "DiffusionSolver.hxx"

/* Check the material data: */
int DiffusionSolver::checkMaterials() {
   
   /* Check the materials: */
   for (int i = 0; i < materials.size(); i++) {
      PAMPA_CHECK(materials(i).num_groups != num_groups, 1, "wrong number of energy groups");
      PAMPA_CHECK(materials(i).sigma_total.empty(), 2, "missing total cross sections");
      PAMPA_CHECK(materials(i).nu_sigma_fission.empty(), 3, "missing nu-fission cross sections");
      PAMPA_CHECK(materials(i).sigma_scattering.empty(), 4, "missing scattering cross sections");
      PAMPA_CHECK(materials(i).diffusion_coefficient.empty(), 5, "missing diffusion coefficients");
      PAMPA_CHECK(materials(i).chi.empty(), 6, "missing fission spectrum");
   }
   
   return 0;
   
}

/* Build the coefficient matrices and the solution vector: */
int DiffusionSolver::build() {
   
   /* Build the coefficient matrices: */
   PAMPA_CALL(buildMatrices(), "unable to build the coefficient matrices");
   
   /* Create the scalar-flux vector: */
   PAMPA_CALL(petsc::create_vector(phi, R), "unable to create the angular-flux vector");
   
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
   
   /* Create, preallocate and set up the coefficient matrices: */
   int size_local = num_cells * num_groups;
   int size_global = num_cells_global * num_groups;
   int num_r_nonzero_max = num_groups + num_faces_max;
   int num_f_nonzero = num_groups;
   PAMPA_CALL(petsc::create_matrix(R, size_local, size_global, num_r_nonzero_max), 
      "unable to create the R coefficient matrix");
   PAMPA_CALL(petsc::create_matrix(F, size_local, size_global, num_f_nonzero), 
      "unable to create the F coefficient matrix");
   
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
         PetscInt l = index(cells.indices(i), g, num_groups);
         int r_i = 1, f_i = 0;
         
         /* Set the total-reaction term: */
         r_l2[0] = l;
         r_l_l2[0] = mat.sigma_total(g) * cells.volumes(i);
         
         /* Set the group-to-group coupling terms: */
         for (int g2 = 0; g2 < num_groups; g2++) {
            
            /* Get the matrix index for cell i and group g2: */
            PetscInt l2 = index(cells.indices(i), g2, num_groups);
            
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
                  
                  /* Other boundary conditions (not implemented): */
                  default : {
                     
                     /* Not implemented: */
                     PAMPA_CHECK(true, 1, "boundary condition not implemented");
                     
                     break;
                     
                  }
                  
               }
               
            }
            
            /* Set the cell-to-cell coupling terms depending on the neighbour material: */
            else {
               
               /* Get the matrix index for cell i2 and group g: */
               PetscInt l2 = index(cells.indices(i2), g, num_groups);
               
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

/* Get the solution after solving the eigensystem: */
int DiffusionSolver::getSolution() {
   
   /* Get the scalar flux from the EPS context: */
   double lambda;
   PETSC_CALL(EPSGetEigenpair(eps, 0, &lambda, NULL, phi, NULL));
   keff = 1.0 / lambda;
   
   /* Normalize the scalar flux: */
   PAMPA_CALL(normalizeScalarFlux(), "unable to normalize the scalar flux");
   
   return 0;
   
}

/* Write the solution to a plain-text file in .vtk format: */
int DiffusionSolver::writeVTK(const std::string& filename) const {
   
   /* Get the number of cells: */
   int num_cells = mesh->getNumCells();
   
   /* Get the array for the scalar flux: */
   PetscScalar* data_phi;
   PETSC_CALL(VecGetArray(phi, &data_phi));
   
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
   
   /* Restore the array for the scalar flux: */
   PETSC_CALL(VecRestoreArray(phi, &data_phi));
   
   return 0;
   
}

/* Write the solution to a binary file in PETSc format: */
int DiffusionSolver::writePETSc(const std::string& filename) const {
   
   /* Write the solution to a binary file: */
   PAMPA_CALL(petsc::write("flux.ptc", phi), "unable to write the solution");
   
   return 0;
   
}

/* Destroy the solution vectors: */
int DiffusionSolver::destroyVectors() {
   
   /* Destroy the scalar-flux vector: */
   PETSC_CALL(VecDestroy(&phi));
   
   return 0;
   
}
