#include "SNSolver.hxx"

#ifdef WITH_EIGEN
#include <Eigen/Dense>
#endif

/* Read the solver from a plain-text input file: */
int SNSolver::read(std::ifstream& file, Array1D<Solver*>& solvers) {
   
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
      else if (line[l] == "order") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the order (N) of the SN method: */
         PAMPA_CHECK(input::read(order, 1, INT_MAX, line[++l]), "wrong SN order");
         
      }
      else if (line[l] == "mixed-face-interpolation") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the weight between upwind and linear schemes for face interpolation: */
         PAMPA_CHECK(input::read(face_interpolation_delta, 0.0, 1.0, line[++l]), 
            "wrong weight between upwind and linear interpolation");
         
      }
      else if (line[l] == "least-squares-boundary-interpolation") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the switch to use the least-squares gradient for boundary interpolation: */
         PAMPA_CHECK(input::read(boundary_interpolation_ls, line[++l]), 
            "wrong switch for least-squares boundary interpolation");
         
         /* Check if the least-squares gradient is supported: */
         #ifndef WITH_EIGEN
         PAMPA_CHECK(boundary_interpolation_ls, 
            "unsupported least-squares gradient interpolation, Eigen library needed");
         #endif
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   return 0;
   
}

/* Get the mapping and the number of faces for boundary cells: */
int SNSolver::getBoundaryCells(Array1D<int>& num_faces_bc) {
   
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
   
   /* Initialize the coefficients: */
   if (bc) {
      Array1D<int> num_faces_bc;
      PAMPA_CHECK(getBoundaryCells(num_faces_bc), "unable to get the boundary cells");
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
            coefs(ic, f, 0) = (r_i2_f+face_interpolation_delta*r_i_f) / r_i_i2;
            coefs(ic, f, 1) = (1.0-face_interpolation_delta)*r_i_f / r_i_i2;
            
            /* Get the flux weights for the face flux for incoming directions: */
            coefs(ic, f, 2) = (1.0-face_interpolation_delta)*r_i2_f / r_i_i2;
            coefs(ic, f, 3) = (r_i_f+face_interpolation_delta*r_i2_f) / r_i_i2;
            
         }
         
      }
      
   }
   
   return 0;
   
}

/* Build the coefficients for the least-squares gradient-discretization scheme: */
#ifdef WITH_EIGEN
int SNSolver::buildLSGradientScheme(Vector3D<double>& coefs, bool bc) {
   
   /* Get the number of dimensions: */
   int num_dims = mesh->getNumDimensions();
   
   /* Initialize the coefficients: */
   if (bc) {
      Array1D<int> num_faces_bc;
      PAMPA_CHECK(getBoundaryCells(num_faces_bc), "unable to get the boundary cells");
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
#endif

/* Calculate the scalar flux: */
int SNSolver::calculateScalarFlux() {
   
   /* Get the angular quadrature weights: */
   const Array1D<double>& weights = quadrature.getWeights();
   
   /* Get the arrays with the raw data: */
   PetscScalar *phi_data, *psi_data;
   PETSC_CALL(VecGetArray(phi, &phi_data));
   PETSC_CALL(VecGetArray(psi, &psi_data));
   
   /* Integrate the angular flux over all directions: */
   for (int iphi = 0, ipsi = 0, i = 0; i < num_cells; i++) {
      for (int g = 0; g < num_energy_groups; g++) {
         phi_data[iphi] = 0.0;
         for (int m = 0; m < num_directions; m++)
            phi_data[iphi] += weights(m) * psi_data[ipsi++];
         phi_data[iphi] *= 4.0 * M_PI;
         iphi++;
      }
   }
   
   /* Restore the arrays with the raw data: */
   PETSC_CALL(VecRestoreArray(phi, &phi_data));
   PETSC_CALL(VecRestoreArray(psi, &psi_data));
   
   return 0;
   
}

/* Normalize the angular flux: */
int SNSolver::normalizeAngularFlux() {
   
   /* Get the angular quadrature weights: */
   const Array1D<double>& weights = quadrature.getWeights();
   
   /* Get the arrays with the raw data: */
   PetscScalar *T_data, *psi_data;
   PETSC_CALL(VecGetArray(T, &T_data));
   PETSC_CALL(VecGetArray(psi, &psi_data));
   
   /* Get the current power: */
   double p0 = 0.0;
   for (int ipsi = 0, i = 0; i < num_cells; i++) {
      const Material* mat = materials(cells.materials(i));
      for (int g = 0; g < num_energy_groups; g++)
         for (int m = 0; m < num_directions; m++)
            p0 += weights(m) * psi_data[ipsi++] * mat->sigmaKappaFission(g, T_data[i]) * 
                     cells.volumes(i);
   }
   MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, &p0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
   
   /* Normalize the angular flux and check for negative fluxes: */
   double f = power / p0;
   for (int ipsi = 0, i = 0; i < num_cells; i++) {
      for (int g = 0; g < num_energy_groups; g++) {
         for (int m = 0; m < num_directions; m++) {
            psi_data[ipsi] *= f;
            PAMPA_CHECK(psi_data[ipsi] < 0.0, "negative values in the angular-flux solution");
            ipsi++;
         }
      }
   }
   
   /* Restore the arrays with the raw data: */
   PETSC_CALL(VecRestoreArray(T, &T_data));
   PETSC_CALL(VecRestoreArray(psi, &psi_data));
   
   return 0;
   
}

/* Build the coefficient matrices and the RHS vector: */
int SNSolver::buildMatrices(int n, double dt, double t) {
   
   /* Get the boundary conditions: */
   if (bcs.empty()) bcs = mesh->getBoundaryConditions();
   
   /* Copy the angular flux from the previous time step: */
   if (n > n0+1) {
      if (psi0 == 0) {
         PETSC_CALL(VecDuplicate(psi, &psi0));
      }
      PETSC_CALL(VecCopy(psi, psi0));
      n0++;
   }
   
   /* Get the angular quadrature data: */
   const Array2D<double>& directions = quadrature.getDirections();
   const Array1D<double>& weights = quadrature.getWeights();
   const Array2D<double>& axes = quadrature.getAxes();
   const Array2D<int>& reflected_directions = quadrature.getReflectedDirections();
   
   /* Initialize the matrix rows for R and F: */
   PetscInt r_l2[num_energy_groups*num_directions];
   PetscScalar r_l_l2[num_energy_groups*num_directions];
   PetscInt f_l2[num_energy_groups*num_directions];
   PetscScalar f_l_l2[num_energy_groups*num_directions];
   
   /* Get the arrays with the raw data: */
   PetscScalar *T_data, *b_data, *S_data, *psi0_data;
   PETSC_CALL(VecGetArray(T, &T_data));
   if (n > 0) {
      PETSC_CALL(VecGetArray(b, &b_data));
      PETSC_CALL(VecGetArray(S, &S_data));
      PETSC_CALL(VecGetArray(psi0, &psi0_data));
   }
   
   /* Calculate the coefficients for each cell i: */
   for (int i = 0; i < num_cells; i++) {
      
      /* Get the material for cell i: */
      const Material* mat = materials(cells.materials(i));
      
      /* Calculate the coefficients for each group g: */
      for (int g = 0; g < num_energy_groups; g++) {
         
         /* Calculate the coefficients for each direction m: */
         for (int m = 0; m < num_directions; m++) {
            
            /* Get the matrix index for cell i, group g and direction m: */
            PetscInt l = index(cells.global_indices(i), g, m);
            int r_i = 1, f_i = 0;
            
            /* Set the total-reaction term: */
            r_l2[0] = l;
            r_l_l2[0] = mat->sigmaTotal(g, T_data[i]) * cells.volumes(i);
            
            /* Set the time-derivative term: */
            if (n > 0) {
               
               /* Set the delayed neutron source: */
               b_data[index(i, g, m)] = mat->chiDelayed(g, T_data[i]) * weights(m) * S_data[i];
               
               /* Get the time-derivative term: */
               double d = cells.volumes(i) / (mat->neutronVelocity(g, T_data[i])*dt);
               
               /* Set the source term for cell i, group g and direction m in the RHS vector: */
               b_data[index(i, g, m)] += d * psi0_data[index(i, g, m)];
               
               /* Set the diagonal term for cell i, group g and direction m: */
               r_l_l2[0] += d;
               
            }
            
            /* Set the group-to-group coupling terms: */
            for (int g2 = 0; g2 < num_energy_groups; g2++) {
               
               /* Set the direction-to-direction coupling terms: */
               for (int m2 = 0; m2 < num_directions; m2++) {
                  
                  /* Get the matrix index for cell i, group g2 and direction m2: */
                  PetscInt l2 = index(cells.global_indices(i), g2, m2);
                  
                  /* Set the (g2 -> g, m2 -> m) isotropic scattering term: */
                  if (l2 == l)
                     r_l_l2[0] += -mat->sigmaScattering(g2, g, T_data[i]) * weights(m2) * 
                                     cells.volumes(i);
                  else
                     r_l_l2[r_i] = -mat->sigmaScattering(g2, g, T_data[i]) * weights(m2) * 
                                      cells.volumes(i);
                  
                  /* Set the (g2 -> g, m2 -> m) fission term: */
                  if (n == 0) {
                     f_l2[f_i] = l2;
                     f_l_l2[f_i++] = mat->chiEffective(g, T_data[i]) * 
                                        mat->sigmaNuFission(g2, T_data[i]) * weights(m2) * 
                                        cells.volumes(i);
                  }
                  else {
                     if (l2 == l)
                        r_l_l2[0] += -(1.0-mat->beta()) * mat->chiPrompt(g, T_data[i]) * 
                                        mat->sigmaNuFission(g2, T_data[i]) * weights(m2) * 
                                        cells.volumes(i) / keff;
                     else
                        r_l_l2[r_i] += -(1.0-mat->beta()) * mat->chiPrompt(g, T_data[i]) * 
                                          mat->sigmaNuFission(g2, T_data[i]) * weights(m2) * 
                                          cells.volumes(i) / keff;
                  }
                  
                  /* Keep the index for the R matrix: */
                  if (l2 != l)
                     r_l2[r_i++] = l2;
                  
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
                           if (boundary_interpolation_ls) {
							  
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
                                    PetscInt l3 = index(cells.global_indices(i3), g, m);
                                    
                                    /* Get the boundary-cell index: */
                                    int ibc = ic_to_ibc(i);
                                    
                                    /* Set the LS contribution for cell i: */
                                    double w_i3 = math::dot_product(dp, grad_coefs_bc(ibc, f2), 3);
                                    r_l_l2[0] += -w_i3 * w * faces.areas(i, f);
                                    
                                    /* Set the LS contribution for cell i3: */
                                    PetscScalar r = w_i3 * w * faces.areas(i, f);
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
                           PAMPA_CHECK(m2 == -1, "reflected direction not found");
                           
                           /* Get the matrix index for cell i, group g and direction m2: */
                           PetscInt l2 = index(cells.global_indices(i), g, m2);
                           
                           /* Set the leakage term for cell i: */
                           PetscScalar r = w * faces.areas(i, f);
                           PETSC_CALL(MatSetValues(R, 1, &l, 1, &l2, &r, ADD_VALUES));
                           
                        }
                        
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
               
               /* Set the cell-to-cell coupling terms: */
               else {
                  
                  /* Get the matrix index for cell i2, group g and direction m: */
                  PetscInt l2 = index(cells.global_indices(i2), g, m);
                  
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
                  PetscScalar r = w_i2 * w * faces.areas(i, f);
                  PETSC_CALL(MatSetValues(R, 1, &l, 1, &l2, &r, ADD_VALUES));
                  
               }
               
            }
            
            /* Set the matrix rows for R and F: */
            PETSC_CALL(MatSetValues(R, 1, &l, r_i, r_l2, r_l_l2, ADD_VALUES));
            if (n == 0) {
               PETSC_CALL(MatSetValues(F, 1, &l, f_i, f_l2, f_l_l2, INSERT_VALUES));
            }
            
         }
         
      }
      
   }
   
   /* Restore the arrays with the raw data: */
   PETSC_CALL(VecRestoreArray(T, &S_data));
   if (n > 0) {
      PETSC_CALL(VecRestoreArray(b, &b_data));
      PETSC_CALL(VecRestoreArray(S, &S_data));
      PETSC_CALL(VecRestoreArray(psi0, &psi0_data));
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
int SNSolver::getSolution(int n) {
   
   /* Solve the eigen- (R*x = (1/keff)*F*x) or linear (R*x = b) system: */
   if (n == 0) {
      
      /* Solve the eigensystem: */
      PETSC_CALL(EPSSetInitialSpace(eps, 1, &psi));
      PAMPA_CHECK(petsc::solve(eps), "unable to solve the eigensystem");
      
      /* Get the angular flux and the multiplication factor from the EPS context: */
      PetscScalar lambda;
      PETSC_CALL(EPSGetEigenpair(eps, 0, &lambda, nullptr, psi, nullptr));
      keff = 1.0 / lambda;
      
      /* Calculate the scalar flux: */
      PAMPA_CHECK(calculateScalarFlux(), "unable to calculate the scalar flux");
      
      /* Normalize the scalar flux: */
      PAMPA_CHECK(normalizeScalarFlux(), "unable to normalize the scalar flux");
      
      /* Normalize the angular flux: */
      PAMPA_CHECK(normalizeAngularFlux(), "unable to normalize the angular flux");
      
   }
   else {
      
      /* Solve the linear system: */
      PAMPA_CHECK(petsc::solve(ksp, b, phi), "unable to solve the linear system");
      
      /* Calculate the scalar flux: */
      PAMPA_CHECK(calculateScalarFlux(), "unable to calculate the scalar flux");
      
   }
   
   return 0;
   
}

/* Check the material data: */
int SNSolver::checkMaterials(bool transient) {
   
   /* Check the materials: */
   for (int i = 0; i < materials.size(); i++) {
      const Material* mat = materials(i);
      PAMPA_CHECK(!(mat->hasNuclearData()), "missing nuclear data");
      PAMPA_CHECK(mat->checkNuclearData(num_energy_groups, false, transient), "wrong nuclear data");
   }
   
   return 0;
   
}

/* Build the coefficient matrices and the solution vectors: */
int SNSolver::build() {
   
   /* Build the angular quadrature set: */
   quadrature = AngularQuadratureSet(order);
   PAMPA_CHECK(quadrature.build(), "unable to build the angular quadrature set");
   
   /* Get the number of directions: */
   num_directions = quadrature.getNumDirections();
   
   /* Build the cell-to-cell coupling coefficients for the gradient-discretization scheme: */
   PAMPA_CHECK(buildGaussGradientScheme(grad_coefs, false), 
      "unable to build the Gauss gradient discretization");
   if (boundary_interpolation_ls) {
      #ifdef WITH_EIGEN
      PAMPA_CHECK(buildLSGradientScheme(grad_coefs_bc, true), 
         "unable to build the least-squares gradient discretization for boundary cells");
      #endif
   }
   
   /* Create, preallocate and set up the coefficient matrices: */
   int size_local = num_cells * num_energy_groups * num_directions;
   int size_global = num_cells_global * num_energy_groups * num_directions;
   int size_cell = num_energy_groups * num_directions;
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
   size_local = num_cells * num_energy_groups;
   size_global = num_cells_global * num_energy_groups;
   PAMPA_CHECK(petsc::create(phi, size_local, size_global, vectors), 
      "unable to create the scalar-flux vector");
   fields.pushBack(Field{"scalar-flux", &phi, false, false, nullptr});
   
   /* Create the angular-flux vector: */
   PAMPA_CHECK(petsc::create(psi, R, vectors), "unable to create the angular-flux vector");
   fields.pushBack(Field{"angular-flux", &psi, false, false, nullptr});
   
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
int SNSolver::writeVTK(const std::string& path, int n) const {
   
   /* Write the scalar flux in .vtk format: */
   PAMPA_CHECK(vtk::write(path + "/output", n, "flux", phi, num_cells, num_energy_groups), 
      "unable to write the scalar flux");
   
   /* Write the angular flux in .vtk format: */
   PAMPA_CHECK(vtk::write(path + "/output", n, "flux", psi, num_cells, num_energy_groups, 
      num_directions), "unable to write the angular flux");
   
   /* Write the thermal power in .vtk format: */
   PAMPA_CHECK(vtk::write(path + "/output", n, "power", q, num_cells), 
      "unable to write the thermal power");
   
   return 0;
   
}

/* Write the solution to a binary file in PETSc format: */
int SNSolver::writePETSc(int n) const {
   
   /* Write the angular flux in PETSc format: */
   PAMPA_CHECK(petsc::write("angular_flux", n, psi), "unable to write the angular flux");
   
   return 0;
   
}
