#include "DiffusionSolver.hxx"

/* Read the solver from a plain-text input file: */
int DiffusionSolver::read(std::ifstream& file, Array1D<Solver*>& solvers) {
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty() || line[0] == "}") break;
      
      /* Get the next keyword: */
      if (line[0] == "energy-groups") {
         
         /* Get the number of energy groups: */
         PAMPA_CALL(utils::read(num_energy_groups, 1, INT_MAX, line[1]), 
            "wrong number of energy groups");
         
      }
      else if (line[0] == "power") {
         
         /* Get the total power: */
         PAMPA_CALL(utils::read(power, 0.0, DBL_MAX, line[1]), "wrong power level");
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 1, "unrecognized keyword '" + line[0] + "'");
         
      }
      
   }
   
   return 0;
   
}

/* Build the coefficient matrices and the RHS vector: */
int DiffusionSolver::buildMatrices(int n, double dt) {
   
   /* Get the boundary conditions: */
   if (bcs.empty()) bcs = mesh->getBoundaryConditions();
   
   /* Copy the scalar flux from the previous time step: */
   if (n > n0+1) {
      if (phi0 == 0) {
         PETSC_CALL(VecDuplicate(phi, &phi0));
      }
      PETSC_CALL(VecCopy(phi, phi0));
      n0++;
   }
   
   /* Initialize the matrix rows for R and F: */
   PetscInt r_l2[num_energy_groups+num_faces_max];
   PetscScalar r_l_l2[num_energy_groups+num_faces_max];
   PetscInt f_l2[num_energy_groups];
   PetscScalar f_l_l2[num_energy_groups];
   
   /* Get the arrays with the raw data: */
   PetscScalar *b_data, *S_data, *phi0_data;
   if (n > 0) {
      PETSC_CALL(VecGetArray(b, &b_data));
      PETSC_CALL(VecGetArray(S, &S_data));
      PETSC_CALL(VecGetArray(phi0, &phi0_data));
   }
   
   /* Calculate the coefficients for each cell i: */
   for (int i = 0; i < num_cells; i++) {
      
      /* Get the material for cell i: */
      const Material& mat = materials(cells.materials(i));
      
      /* Calculate the coefficients for each group g: */
      for (int g = 0; g < num_energy_groups; g++) {
         
         /* Get the matrix index for cell i and group g: */
         PetscInt r_i = 1, f_i = 0, l = index(cells.indices(i), g);
         
         /* Set the total-reaction term: */
         r_l2[0] = l;
         r_l_l2[0] = mat.sigma_total(g) * cells.volumes(i);
         
         /* Set the time-derivative term: */
         if (n > 0) {
            
            /* Set the delayed neutron source: */
            b_data[index(i, g)] = mat.chi_delayed(g) * S_data[i];
            
            /* Get the time-derivative term: */
            double d = cells.volumes(i) / (mat.velocity(g)*dt);
            
            /* Set the source term for cell i and group g in the RHS vector: */
            b_data[index(i, g)] += d * phi0_data[index(i, g)];
            
            /* Set the diagonal term for cell i and group g: */
            r_l_l2[0] += d;
            
         }
         
         /* Set the group-to-group coupling terms: */
         for (int g2 = 0; g2 < num_energy_groups; g2++) {
            
            /* Get the matrix index for cell i and group g2: */
            PetscInt l2 = index(cells.indices(i), g2);
            
            /* Set the (g2 -> g) scattering term: */
            if (l2 == l)
               r_l_l2[0] += -mat.sigma_scattering(g2, g) * cells.volumes(i);
            else
               r_l_l2[r_i] = -mat.sigma_scattering(g2, g) * cells.volumes(i);
            
            /* Set the (g2 -> g) fission term: */
            if (n == 0) {
               f_l2[f_i] = l2;
               f_l_l2[f_i++] = mat.chi_eff(g) * mat.nu_sigma_fission(g2) * cells.volumes(i);
            }
            else {
               if (l2 == l)
                  r_l_l2[0] += -(1.0-mat.beta_total) * mat.chi_prompt(g) * 
                                  mat.nu_sigma_fission(g2) * cells.volumes(i) / keff;
               else
                  r_l_l2[r_i] += -(1.0-mat.beta_total) * mat.chi_prompt(g) * 
                                    mat.nu_sigma_fission(g2) * cells.volumes(i) / keff;
            }
            
            /* Keep the index for the R matrix: */
            if (l2 != l)
               r_l2[r_i++] = l2;
            
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
               PetscInt l2 = index(cells.indices(i2), g);
               
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
         if (n == 0) {
            PETSC_CALL(MatSetValues(F, 1, &l, f_i, f_l2, f_l_l2, INSERT_VALUES));
         }
         
      }
      
   }
   
   /* Restore the arrays with the raw data: */
   if (n > 0) {
      PETSC_CALL(VecRestoreArray(b, &b_data));
      PETSC_CALL(VecRestoreArray(S, &S_data));
      PETSC_CALL(VecRestoreArray(phi0, &phi0_data));
   }
   
   /* Assembly the coefficient matrices: */
   PETSC_CALL(MatAssemblyBegin(R, MAT_FINAL_ASSEMBLY));
   PETSC_CALL(MatAssemblyEnd(R, MAT_FINAL_ASSEMBLY));
   if (n == 0) {
      PETSC_CALL(MatAssemblyBegin(F, MAT_FINAL_ASSEMBLY));
      PETSC_CALL(MatAssemblyEnd(F, MAT_FINAL_ASSEMBLY));
   }
   
   return 0;
   
}

/* Solve the linear system and get the solution: */
int DiffusionSolver::getSolution(int n) {
   
   /* Solve the eigen- (R*x = (1/keff)*F*x) or linear (R*x = b) system: */
   if (n == 0) {
      
      /* Solve the eigensystem: */
      PAMPA_CALL(petsc::solve(eps), "unable to solve the eigensystem");
      
      /* Get the scalar flux and the multiplication factor from the EPS context: */
      PetscScalar lambda;
      PETSC_CALL(EPSGetEigenpair(eps, 0, &lambda, NULL, phi, NULL));
      keff = 1.0 / lambda;
      
      /* Normalize the scalar flux: */
      PAMPA_CALL(normalizeScalarFlux(), "unable to normalize the scalar flux");
      
   }
   else {
      
      /* Solve the linear system: */
      PAMPA_CALL(petsc::solve(ksp, b, phi), "unable to solve the linear system");
      
   }
   
   return 0;
   
}

/* Check the material data: */
int DiffusionSolver::checkMaterials(bool transient) {
   
   /* Check the materials: */
   for (int i = 0; i < materials.size(); i++) {
      PAMPA_CHECK(materials(i).num_energy_groups != num_energy_groups, 1, 
         "wrong number of energy groups");
      PAMPA_CHECK(materials(i).sigma_total.empty(), 2, 
         "missing total cross sections");
      PAMPA_CHECK(materials(i).nu_sigma_fission.empty(), 3, 
         "missing nu-fission cross sections");
      PAMPA_CHECK(materials(i).kappa_sigma_fission.empty(), 4, 
         "missing kappa-fission cross sections");
      PAMPA_CHECK(materials(i).sigma_scattering.empty(), 5, 
         "missing scattering cross sections");
      PAMPA_CHECK(materials(i).diffusion_coefficient.empty(), 6, 
         "missing diffusion coefficients");
      PAMPA_CHECK(materials(i).chi_prompt.empty(), 7, 
         "missing prompt fission spectrum");
      if (transient) {
         PAMPA_CHECK(materials(i).chi_delayed.empty(), 8, 
            "missing delayed fission spectrum");
         PAMPA_CHECK(materials(i).velocity.empty(), 9, 
            "missing neutron velocities");
      }
   }
   
   return 0;
   
}

/* Build the coefficient matrices and the solution vector: */
int DiffusionSolver::build() {
   
   /* Create, preallocate and set up the coefficient matrices: */
   int size_local = num_cells * num_energy_groups;
   int size_global = num_cells_global * num_energy_groups;
   PAMPA_CALL(petsc::create(R, size_local, size_global, num_energy_groups+num_faces_max, matrices), 
      "unable to create the R coefficient matrix");
   PAMPA_CALL(petsc::create(F, size_local, size_global, num_energy_groups, matrices), 
      "unable to create the F coefficient matrix");
   
   /* Create the right-hand-side vector: */
   PAMPA_CALL(petsc::create(b, R, vectors), "unable to create the right-hand-side vector");
   
   /* Create the delayed-neutron-source vector: */
   PAMPA_CALL(petsc::create(S, num_cells, num_cells_global, vectors), 
      "unable to create the delayed-neutron-source vector");
   fields.pushBack(Field{"delayed-source", &S, true, false});
   
   /* Create the scalar-flux vectors: */
   PAMPA_CALL(petsc::create(phi, R, vectors), "unable to create the angular-flux vector");
   fields.pushBack(Field{"scalar-flux", &phi, false, true});
   
   /* Create the thermal-power vector: */
   PAMPA_CALL(petsc::create(q, num_cells, num_cells_global, vectors), 
      "unable to create the thermal-power vector");
   fields.pushBack(Field{"power", &q, false, true});
   
   /* Create the production-rate vector: */
   PAMPA_CALL(petsc::create(P, num_cells, num_cells_global, vectors), 
      "unable to create the production-rate vector");
   fields.pushBack(Field{"production-rate", &P, false, true});
   
   return 0;
   
}

/* Write the solution to a plain-text file in .vtk format: */
int DiffusionSolver::writeVTK(const std::string& filename) const {
   
   /* Write the scalar flux in .vtk format: */
   PAMPA_CALL(vtk::write(filename, "flux", phi, num_cells, num_energy_groups), 
      "unable to write the scalar flux");
   
   /* Write the thermal power in .vtk format: */
   PAMPA_CALL(vtk::write(filename, "power", q, num_cells), "unable to write the thermal power");
   
   return 0;
   
}

/* Write the solution to a binary file in PETSc format: */
int DiffusionSolver::writePETSc() const {
   
   /* Write the scalar flux in PETSc format: */
   PAMPA_CALL(petsc::write("flux.ptc", phi), "unable to write the scalar flux");
   
   return 0;
   
}
