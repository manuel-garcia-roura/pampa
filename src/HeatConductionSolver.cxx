#include "HeatConductionSolver.hxx"

/* Solve the linear system to get the solution: */
int HeatConductionSolver::solve(int n, double dt) {
   
   /* Build the coefficient matrix and the RHS vector: */
   PAMPA_CALL(buildMatrix(n, dt), "unable to build the coefficient matrix and the RHS vector");
   
   /* Create the KSP context: */
   if (ksp == 0) {
      PAMPA_CALL(petsc::create(ksp, A), "unable to create the KSP context");
   }
   
   /* Solve the linear system: */
   PAMPA_CALL(petsc::solve(ksp, b, T), "unable to solve the linear system");
   
   /* Get a random volumetric heat source for the next time step: */
   int ip = std::min(n+1, power.size()-1);
   PAMPA_CALL(petsc::random(q), "unable to initialize the volumetric heat source");
   PAMPA_CALL(petsc::normalize(q, power(ip)), "unable to normalize the volumetric heat source");
   
   return 0;
   
}

/* Build the coefficient matrix and the RHS vector: */
int HeatConductionSolver::buildMatrix(int n, double dt) {
   
   /* Get the boundary conditions: */
   if (bcs.empty()) bcs = mesh->getBoundaryConditions();
   
   /* Copy the temperature from the previous time step: */
   if (n > n0+1) {
      PETSC_CALL(VecCopy(T, T0));
      n0++;
   }
   
   /* Initialize the matrix rows for A: */
   PetscInt a_i2[1+num_faces_max];
   PetscScalar a_i_i2[1+num_faces_max];
   
   /* Get the arrays with the raw data: */
   PetscScalar *T0_data, *b_data, *q_data;;
   PETSC_CALL(VecGetArray(T0, &T0_data));
   PETSC_CALL(VecGetArray(b, &b_data));
   PETSC_CALL(VecGetArray(q, &q_data));
   
   /* Calculate the coefficients for each cell i: */
   for (int i = 0; i < num_cells; i++) {
      
      /* Get the material for cell i: */
      const Material& mat = materials(cells.materials(i));
      
      /* Set the volumetric heat source: */
      b_data[i] = q_data[i];
      
      /* Set the time-derivative term: */
      a_i2[0] = cells.indices(i);
      a_i_i2[0] = 0.0;
      if (n > 0) {
         
         /* Get the time-derivative term: */
         double a = mat.rho * mat.cp * cells.volumes(i) / dt;
         
         /* Set the source term for cell i in the RHS vector: */
         b_data[i] += a * T0_data[i];
         
         /* Set the diagonal term for cell i: */
         a_i_i2[0] += a;
         
      }
      
      /* Set the cell-to-cell coupling terms: */
      PetscInt a_i = 1;
      double a;
      for (int f = 0; f < faces.num_faces(i); f++) {
         
         /* Get the index for cell i2 (actual cell or boundary condition): */
         /* Note: boundary conditions have negative, 1-based indexes: */
         int i2 = faces.neighbours(i, f);
         
         /* Set the boundary conditions: */
         if (i2 < 0) {
            
            /* Check the boundary-condition type: */
            switch (bcs(-i2).type) {
               
               /* Set Dirichlet boundary conditions: */
               case BC::DIRICHLET : {
                  
                  /* Get the surface leakage factor: */
                  double w = math::surface_leakage_factor(cells.centroids(i), 
                                faces.centroids(i, f), faces.normals(i, f));
                  
                  /* Set the leakage term for cell i: */
                  a = w * mat.k * faces.areas(i, f);
                  a_i_i2[0] += a;
                  
                  /* Set the leakage term for cell i in the RHS vector: */
                  b_data[i] += a * bcs(-i2).x;
                  
                  break;
                  
               }
               
               /* Set reflective boundary conditions (nothing to be done): */
               case BC::REFLECTIVE : {
                  
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
            
            /* Get the material for cell i2: */
            const Material& mat2 = materials(cells.materials(i2));
            
            /* Set the terms for cells with the same materials: */
            if (&mat2 == &mat) {
               
               /* Get the surface leakage factor: */
               double w = math::surface_leakage_factor(cells.centroids(i), cells.centroids(i2), 
                             faces.normals(i, f));
               
               /* Get the leakage term for cell i2: */
               a = -w * mat.k * faces.areas(i, f);
               
            }
            
            /* Set the terms for cells with different materials: */
            else {
               
               /* Get the surface leakage factor and the weight for cell i: */
               double w_i_i2 = math::surface_leakage_factor(cells.centroids(i), 
                                  faces.centroids(i, f), faces.normals(i, f));
               w_i_i2 *= mat.k * faces.areas(i, f);
               
               /* Get the surface leakage factor and the weight for cell i2: */
               double w_i2_i = math::surface_leakage_factor(cells.centroids(i2), 
                                  faces.centroids(i, f), faces.normals(i, f));
               w_i2_i *= -mat2.k * faces.areas(i, f);
               
               /* Get the leakage term for cell i2: */
               a = -(w_i_i2*w_i2_i) / (w_i_i2+w_i2_i);
               
            }
            
            /* Set the leakage term for cell i: */
            a_i_i2[0] -= a;
            
            /* Set the leakage term for cell i2: */
            a_i2[a_i] = cells.indices(i2);
            a_i_i2[a_i++] = a;
            
         }
         
      }
      
      /* Set the matrix rows for A: */
      PETSC_CALL(MatSetValues(A, 1, &(cells.indices(i)), a_i, a_i2, a_i_i2, INSERT_VALUES));
      
   }
   
   /* Restore the arrays with the raw data: */
   PETSC_CALL(VecRestoreArray(T0, &T0_data));
   PETSC_CALL(VecRestoreArray(b, &b_data));
   PETSC_CALL(VecRestoreArray(q, &q_data));
   
   /* Assembly the coefficient matrix: */
   PETSC_CALL(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
   PETSC_CALL(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
   
   return 0;
   
}

/* Check the material data: */
int HeatConductionSolver::checkMaterials() const {
   
   /* Check the materials: */
   for (int i = 0; i < materials.size(); i++) {
      PAMPA_CHECK(materials(i).rho < 0.0, 1, "missing density");
      PAMPA_CHECK(materials(i).cp < 0.0, 2, "missing specific heat capacity");
      PAMPA_CHECK(materials(i).k < 0.0, 3, "missing thermal conductivity");
   }
   
   return 0;
   
}

/* Build the coefficient matrix, and the solution and RHS vectors: */
int HeatConductionSolver::build() {
   
   /* Create, preallocate and set up the coefficient matrix: */
   PAMPA_CALL(petsc::create(A, num_cells, num_cells_global, 1+num_faces_max, matrices), 
      "unable to create the coefficient matrix");
   
   /* Create the right-hand-side vector: */
   PAMPA_CALL(petsc::create(b, A, vectors), "unable to create the right-hand-side vector");
   
   /* Create the heat-source vector: */
   PAMPA_CALL(petsc::create(q, A, vectors), "unable to create the heat-source vector");
   fields.pushBack(Field{"power", &q, true, false});
   
   /* Create the temperature vectors: */
   PAMPA_CALL(petsc::create(T, A, vectors), "unable to create the temperature vector");
   PAMPA_CALL(petsc::create(T0, A, vectors), "unable to create the temperature vector");
   fields.pushBack(Field{"temperature", &T, false, true});
   
   /* Get a random volumetric heat source for the steady state: */
   PAMPA_CALL(petsc::random(q), "unable to initialize the volumetric heat source");
   PAMPA_CALL(petsc::normalize(q, power(0)), "unable to normalize the volumetric heat source");
   
   return 0;
   
}

/* Print the solution summary to standard output: */
int HeatConductionSolver::printLog(int n) const {
   
   /* Print out the minimum and maximum temperatures: */
   PetscScalar T_min, T_max;
   PETSC_CALL(VecMin(T, NULL, &T_min));
   PETSC_CALL(VecMax(T, NULL, &T_max));
   mpi::print("T_min", T_min);
   mpi::print("T_max", T_max);
   
   return 0;
   
}

/* Write the solution to a plain-text file in .vtk format: */
int HeatConductionSolver::writeVTK(const std::string& filename) const {
   
   /* Write the temperature in .vtk format: */
   PAMPA_CALL(vtk::write(filename, "temperature", T, num_cells), "unable to write the temperature");
   
   return 0;
   
}

/* Write the solution to a binary file in PETSc format: */
int HeatConductionSolver::writePETSc() const {
   
   /* Write the temperature in PETSc format: */
   PAMPA_CALL(petsc::write("temperature.ptc", T), "unable to write the temperature");
   
   return 0;
   
}
