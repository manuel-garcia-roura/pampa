#include "HeatConductionSolver.hxx"

/* Read the solver from a plain-text input file: */
int HeatConductionSolver::read(std::ifstream& file, Array1D<Solver*>& solvers) {
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty() || line[0] == "}") break;
      
      /* Get the next keyword: */
      unsigned int l = 0;
      if (line[l] == "bc") {
         
         /* Initialize the boundary-condition array if not done yet: */
         if (bcs.empty()) bcs.resize(1);
         
         /* Get the boundary condition (1-based indexed): */
         int i;
         BoundaryCondition bc;
         PAMPA_CALL(utils::read(i, bcs.size(), bcs.size(), line[++l]), 
            "wrong boundary condition index");
         PAMPA_CALL(utils::read(bc, line, ++l, file), "wrong boundary condition");
         bcs.pushBack(bc);
         
      }
      else if (line[l] == "fixed") {
         
         /* Get the material and the value: */
         int mat;
         PAMPA_CALL(utils::read(mat, 1, materials.size(), line[++l]), "wrong material index");
         Function temp;
         PAMPA_CALL(utils::read(temp, line, ++l, file), "wrong fixed value");
         
         /* Set the fixed temperature: */
         fixed_temperatures(mat-1) = temp;
         
      }
      else if (line[l] == "power") {
         
         /* Get the total power: */
         PAMPA_CALL(utils::read(power, line, ++l, file), "power data");
         
      }
      else if (line[l] == "convergence") {
         
         /* Get the convergence tolerance and p-norm for nonlinear problems: */
         PAMPA_CALL(utils::read(tol, 0.0, DBL_MAX, line[++l]), "wrong convergence tolerance");
         PAMPA_CALL(utils::read(p, 0.0, DBL_MAX, line[++l]), "wrong convergence p-norm");
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 1, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   return 0;
   
}

/* Solve the linear system to get the solution: */
int HeatConductionSolver::solve(int n, double dt, double t) {
   
   /* Get a random volumetric heat source: */
   if (!(power.empty())) {
      PAMPA_CALL(petsc::random(q), "unable to initialize the volumetric heat source");
      PAMPA_CALL(petsc::normalize(q, power(t)), "unable to normalize the volumetric heat source");
   }
   
   /* Solve the linear system until convegence: */
   bool converged = false;
   while (!converged) {
      
      /* Build the coefficient matrix and the RHS vector: */
      PAMPA_CALL(buildMatrix(n, dt, t), 
         "unable to build the coefficient matrix and the RHS vector");
      
      /* Create the KSP context: */
      if (ksp == 0) {
         PAMPA_CALL(petsc::create(ksp, A), "unable to create the KSP context");
      }
      
      /* Solve the linear system: */
      PAMPA_CALL(petsc::solve(ksp, b, T), "unable to solve the linear system");
      
      /* Evaluate the convergence: */
      converged = true;
      if (nonlinear) {
         if (Tprev == 0) {
            PETSC_CALL(VecDuplicate(T, &Tprev));
            converged = false;
         }
         else {
            double eps;
            PAMPA_CALL(petsc::difference(T, Tprev, p, eps, false), 
               "unable to calculate the convergence error");
            converged = eps < tol;
         }
         PETSC_CALL(VecCopy(T, Tprev));
      }
      
   }
   
   return 0;
   
}

/* Build the coefficient matrix and the RHS vector: */
int HeatConductionSolver::buildMatrix(int n, double dt, double t) {
   
   /* Get the boundary conditions: */
   if (bcs.empty()) bcs = mesh->getBoundaryConditions();
   
   /* Copy the temperature from the previous time step: */
   if (n > n0+1) {
      if (T0 == 0) {
         PETSC_CALL(VecDuplicate(T, &T0));
      }
      PETSC_CALL(VecCopy(T, T0));
      n0++;
   }
   
   /* Initialize the matrix rows for A: */
   PetscInt a_i2[1+num_faces_max];
   PetscScalar a_i_i2[1+num_faces_max];
   
   /* Get the arrays with the raw data: */
   PetscScalar *b_data, *q_data, *T_data, *T0_data;
   PETSC_CALL(VecGetArray(b, &b_data));
   PETSC_CALL(VecGetArray(q, &q_data));
   PETSC_CALL(VecGetArray(T, &T_data));
   if (n > 0) {
      PETSC_CALL(VecGetArray(T0, &T0_data));
   }
   
   /* Calculate the coefficients for each cell i: */
   for (int i = 0; i < num_cells; i++) {
      
      /* Get the material for cell i: */
      const Material& mat = materials(cells.materials(i));
      
      /* Set a fixed temperature: */
      if (!(fixed_temperatures(cells.materials(i)).empty())) {
         b_data[i] = fixed_temperatures(cells.materials(i))(t);
         double a = 1.0;
         PETSC_CALL(MatSetValues(A, 1, &(cells.indices(i)), 1, &(cells.indices(i)), &a, 
            INSERT_VALUES));
         continue;
      }
      
      /* Set the volumetric heat source: */
      b_data[i] = q_data[i];
      
      /* Set the time-derivative term: */
      a_i2[0] = cells.indices(i);
      a_i_i2[0] = 0.0;
      if (n > 0) {
         
         /* Get the time-derivative term: */
         double a = mat.rho(T_data[i]) * mat.cp(T_data[i]) * cells.volumes(i) / dt;
         
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
                  a = w * mat.k(T_data[i]) * faces.areas(i, f);
                  a_i_i2[0] += a;
                  
                  /* Set the leakage term for cell i in the RHS vector: */
                  b_data[i] += a * bcs(-i2).x(t);
                  
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
               a = -w * mat.k(T_data[i]) * faces.areas(i, f);
               
            }
            
            /* Set the terms for cells with different materials: */
            else {
               
               /* Get the surface leakage factor and the weight for cell i: */
               double w_i_i2 = math::surface_leakage_factor(cells.centroids(i), 
                                  faces.centroids(i, f), faces.normals(i, f));
               w_i_i2 *= mat.k(T_data[i]) * faces.areas(i, f);
               
               /* Get the surface leakage factor and the weight for cell i2: */
               double w_i2_i = math::surface_leakage_factor(cells.centroids(i2), 
                                  faces.centroids(i, f), faces.normals(i, f));
               w_i2_i *= -mat2.k(T_data[i]) * faces.areas(i, f);
               
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
   PETSC_CALL(VecRestoreArray(b, &b_data));
   PETSC_CALL(VecRestoreArray(q, &q_data));
   PETSC_CALL(VecRestoreArray(T, &T_data));
   if (n > 0) {
      PETSC_CALL(VecRestoreArray(T0, &T0_data));
   }
   
   /* Assembly the coefficient matrix: */
   PETSC_CALL(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
   PETSC_CALL(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
   
   return 0;
   
}

/* Check the material data: */
int HeatConductionSolver::checkMaterials(bool transient) {
   
   /* Check the materials: */
   for (int i = 0; i < materials.size(); i++) {
      PAMPA_CHECK(materials(i).k0 < 0.0, 1, "missing thermal conductivity");
      if (transient) {
         PAMPA_CHECK(materials(i).rho0 < 0.0, 2, "missing density");
         PAMPA_CHECK(materials(i).cp0 < 0.0, 3, "missing specific heat capacity");
      }
      if (materials(i).thermal_properties != TH::CONSTANT)
         nonlinear = true;
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
   fields.pushBack(Field{"temperature", &T, false, true});
   
   /* Get a random volumetric heat source: */
   if (power.empty()) {
      PAMPA_CALL(petsc::random(q), "unable to initialize the volumetric heat source");
      PAMPA_CALL(petsc::normalize(q, 1.0), "unable to normalize the volumetric heat source");
   }
   
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
