#include "DiffusionSolver.hxx"

/* Read the solver from a plain-text input file: */
int DiffusionSolver::read(std::ifstream& file, Array1D<Solver*>& solvers) {
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = input::get_next_line(file);
      if (line.empty() || line[0] == "}") break;
      
      /* Get the next keyword: */
      unsigned int l = 0;
      if (line[l] == "energy-groups") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the number of energy groups: */
         PAMPA_CHECK(input::read(num_energy_groups, 1, INT_MAX, line[++l]), 
            "wrong number of energy groups");
         
      }
      else if (line[l] == "bc") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() < 3, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the mesh boundaries: */
         const Array1D<std::string>& boundaries = mesh->getBoundaries();
         
         /* Initialize the boundary-condition array, if not done yet: */
         if (bcs.empty()) bcs.resize(1+boundaries.size());
         
         /* Get the boundary name and index: */
         std::string name = line[++l];
         int ibc = boundaries.find(name);
         PAMPA_CHECK(ibc < 0, "wrong boundary name");
         
         /* Get the boundary condition (1-based indexed): */
         PAMPA_CHECK(input::read(bcs(ibc+1), line, ++l, file), "wrong boundary condition");
         
      }
      else if (line[l] == "power") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the total power: */
         PAMPA_CHECK(input::read(power, 0.0, DBL_MAX, line[++l]), "wrong power level");
         
      }
      else if (line[l] == "convergence") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 5, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the field name: */
         std::string name = line[++l];
         
         /* Get the convergence error: */
         if (name == "power") {
            PAMPA_CHECK(input::read(dq, name, line, l), "wrong convergence error");
         }
         else if (name == "production-rate") {
            PAMPA_CHECK(input::read(dP, name, line, l), "wrong convergence error");
         }
         else {
            PAMPA_CHECK(true, "wrong field");
         }
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   return 0;
   
}

/* Build the coefficient matrices and the RHS vector: */
int DiffusionSolver::buildMatrices(int n, double dt, double t) {
   
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
   PetscScalar *T_data, *b_data, *S_data, *phi0_data;
   PETSC_CALL(VecGetArray(T, &T_data));
   if (n > 0) {
      PETSC_CALL(VecGetArray(b, &b_data));
      PETSC_CALL(VecGetArray(S, &S_data));
      PETSC_CALL(VecGetArray(phi0, &phi0_data));
   }
   
   /* Calculate the coefficients for each cell i: */
   for (int i = 0; i < num_cells; i++) {
      
      /* Get the material for cell i: */
      const Material* mat = materials(cells.materials(i));
      
      /* Calculate the coefficients for each group g: */
      for (int g = 0; g < num_energy_groups; g++) {
         
         /* Get the matrix index for cell i and group g: */
         PetscInt r_i = 1, f_i = 0, l = index(cells.global_indices(i), g);
         
         /* Set the total-reaction term: */
         r_l2[0] = l;
         r_l_l2[0] = mat->sigmaTotal(g, T_data[i]) * cells.volumes(i);
         
         /* Set the time-derivative term: */
         if (n > 0) {
            
            /* Set the delayed neutron source: */
            b_data[index(i, g)] = mat->chiDelayed(g, T_data[i]) * S_data[i];
            
            /* Get the time-derivative term: */
            double d = cells.volumes(i) / (mat->neutronVelocity(g, T_data[i])*dt);
            
            /* Set the source term for cell i and group g in the RHS vector: */
            b_data[index(i, g)] += d * phi0_data[index(i, g)];
            
            /* Set the diagonal term for cell i and group g: */
            r_l_l2[0] += d;
            
         }
         
         /* Set the group-to-group coupling terms: */
         for (int g2 = 0; g2 < num_energy_groups; g2++) {
            
            /* Get the matrix index for cell i and group g2: */
            PetscInt l2 = index(cells.global_indices(i), g2);
            
            /* Set the (g2 -> g) scattering term: */
            if (l2 == l)
               r_l_l2[0] += -mat->sigmaScattering(g2, g, T_data[i]) * cells.volumes(i);
            else
               r_l_l2[r_i] = -mat->sigmaScattering(g2, g, T_data[i]) * cells.volumes(i);
            
            /* Set the (g2 -> g) fission term: */
            if (n == 0) {
               f_l2[f_i] = l2;
               f_l_l2[f_i++] = mat->chiEffective(g, T_data[i]) * 
                                  mat->sigmaNuFission(g2, T_data[i]) * cells.volumes(i);
            }
            else {
               if (l2 == l)
                  r_l_l2[0] += -(1.0-mat->beta()) * mat->chiPrompt(g, T_data[i]) * 
                                  mat->sigmaNuFission(g2, T_data[i]) * cells.volumes(i) / keff;
               else
                  r_l_l2[r_i] += -(1.0-mat->beta()) * mat->chiPrompt(g, T_data[i]) * 
                                    mat->sigmaNuFission(g2, T_data[i]) * cells.volumes(i) / keff;
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
            int i2 = faces.neighbors(i, f);
            
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
                     r_l_l2[0] += w * mat->diffusionCoefficient(g, T_data[i]) * faces.areas(i, f);
                     
                     break;
                     
                  }
                  
                  /* Set reflective (zero-current) boundary conditions (nothing to be done): */
                  case BC::REFLECTIVE : {
                     
                     break;
                     
                  }
                  
                  /* Set Robin boundary conditions: */
                  case BC::ROBIN : {
                     
                     /* Set the leakage term for cell i: */
                     r_l_l2[0] -= bcs(-i2).f(0)(t) * faces.areas(i, f);
                     
                     break;
                     
                  }
                  
                  /* Other boundary conditions (not implemented): */
                  default : {
                     
                     /* Not implemented: */
                     PAMPA_CHECK(true, "boundary condition not implemented");
                     
                     break;
                     
                  }
                  
               }
               
            }
            
            /* Set the cell-to-cell coupling terms depending on the neighbor material: */
            else {
               
               /* Get the matrix index for cell i2 and group g: */
               PetscInt l2 = index(cells.global_indices(i2), g);
               
               /* Get the material for cell i2: */
               const Material* mat2 = materials(cells.materials(i2));
               
               /* Set the terms for cells with the same materials: */
               if (mat2 == mat) {
                  
                  /* Get the surface leakage factor: */
                  double w = math::surface_leakage_factor(cells.centroids(i), cells.centroids(i2), 
                                faces.normals(i, f));
                  
                  /* Get the leakage term for cell i2: */
                  r = -w * mat->diffusionCoefficient(g, T_data[i]) * faces.areas(i, f);
                  
               }
               
               /* Set the terms for cells with different materials: */
               else {
                  
                  /* Get the surface leakage factor and the weight for cell i: */
                  double w_i_i2 = math::surface_leakage_factor(cells.centroids(i), 
                                     faces.centroids(i, f), faces.normals(i, f));
                  w_i_i2 *= mat->diffusionCoefficient(g, T_data[i]) * faces.areas(i, f);
                  
                  /* Get the surface leakage factor and the weight for cell i2: */
                  double w_i2_i = math::surface_leakage_factor(cells.centroids(i2), 
                                     faces.centroids(i, f), faces.normals(i, f));
                  w_i2_i *= -mat2->diffusionCoefficient(g, T_data[i]) * faces.areas(i, f);
                  
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
   PETSC_CALL(VecRestoreArray(T, &S_data));
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
      PETSC_CALL(EPSSetInitialSpace(eps, 1, &phi));
      PAMPA_CHECK(petsc::solve(eps), "unable to solve the eigensystem");
      
      /* Get the scalar flux and the multiplication factor from the EPS context: */
      PetscScalar lambda;
      PETSC_CALL(EPSGetEigenpair(eps, 0, &lambda, nullptr, phi, nullptr));
      keff = 1.0 / lambda;
      
      /* Normalize the scalar flux: */
      PAMPA_CHECK(normalizeScalarFlux(), "unable to normalize the scalar flux");
      
   }
   else {
      
      /* Solve the linear system: */
      PAMPA_CHECK(petsc::solve(ksp, b, phi), "unable to solve the linear system");
      
   }
   
   return 0;
   
}

/* Check the material data: */
int DiffusionSolver::checkMaterials(bool transient) {
   
   /* Check the materials: */
   for (int i = 0; i < materials.size(); i++) {
      const Material* mat = materials(i);
      PAMPA_CHECK(!(mat->hasNuclearData()), "missing nuclear data");
      PAMPA_CHECK(mat->checkNuclearData(num_energy_groups, true, transient), "wrong nuclear data");
   }
   
   return 0;
   
}

/* Build the coefficient matrices and the solution vector: */
int DiffusionSolver::build() {
   
   /* Create, preallocate and set up the coefficient matrices: */
   int size_local = num_cells * num_energy_groups;
   int size_global = num_cells_global * num_energy_groups;
   int size_cell = num_energy_groups;
   PAMPA_CHECK(petsc::create(R, size_local, size_global, size_cell+num_faces_max, matrices), 
      "unable to create the R coefficient matrix");
   PAMPA_CHECK(petsc::create(F, size_local, size_global, size_cell, matrices), 
      "unable to create the F coefficient matrix");
   
   /* Create the right-hand-side vector: */
   PAMPA_CHECK(petsc::create(b, R, vectors), "unable to create the right-hand-side vector");
   
   /* Create the temperature vector: */
   PAMPA_CHECK(petsc::create(T, num_cells, num_cells_global, vectors), 
      "unable to create the temperature vector");
   fields.pushBack(Field{"temperature", &T, true, false, nullptr});
   
   /* Create the delayed-neutron-source vector: */
   PAMPA_CHECK(petsc::create(S, num_cells, num_cells_global, vectors), 
      "unable to create the delayed-neutron-source vector");
   fields.pushBack(Field{"delayed-source", &S, true, false, nullptr});
   
   /* Create the scalar-flux vector: */
   PAMPA_CHECK(petsc::create(phi, R, vectors), "unable to create the angular-flux vector");
   fields.pushBack(Field{"scalar-flux", &phi, false, false, nullptr});
   
   /* Create the thermal-power vector: */
   PAMPA_CHECK(petsc::create(q, num_cells, num_cells_global, vectors), 
      "unable to create the thermal-power vector");
   fields.pushBack(Field{"power", &q, false, true, &dq});
   
   /* Create the production-rate vector: */
   PAMPA_CHECK(petsc::create(P, num_cells, num_cells_global, vectors), 
      "unable to create the production-rate vector");
   fields.pushBack(Field{"production-rate", &P, false, true, &dP});
   
   return 0;
   
}

/* Write the solution to a plain-text file in .vtk format: */
int DiffusionSolver::writeVTK(const std::string& path, int n) const {
   
   /* Write the scalar flux in .vtk format: */
   PAMPA_CHECK(vtk::write(path + "/output", n, "flux", phi, num_cells, num_energy_groups), 
      "unable to write the scalar flux");
   
   /* Write the thermal power in .vtk format: */
   PAMPA_CHECK(vtk::write(path + "/output", n, "power", q, num_cells), 
      "unable to write the thermal power");
   
   return 0;
   
}

/* Write the solution to a binary file in PETSc format: */
int DiffusionSolver::writePETSc(int n) const {
   
   /* Write the scalar flux in PETSc format: */
   PAMPA_CHECK(petsc::write("scalar_flux", n, phi), "unable to write the scalar flux");
   
   return 0;
   
}
