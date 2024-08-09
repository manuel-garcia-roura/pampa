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
         
         /* Initialize the boundary-condition array, if not done yet: */
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
   
   /* Print progress: */
   if (!microscale)
      mpi::print("Run '" + name + "' solver...", true);
   
   /* Get a random volumetric heat source: */
   if (!(power.empty())) {
      if (!microscale && false) {
         PAMPA_CALL(petsc::random(q), "unable to initialize the volumetric heat source");
      }
      else {
         PetscScalar *q_data;
         PETSC_CALL(VecGetArray(q, &q_data));
         for (int i = 0; i < num_cells; i++) {
            if (multiscale)
               q_data[i] = (solvers(cells.materials(i))).getFuelVolume();
            else {
               if ((materials(cells.materials(i)))->isFuel())
                  q_data[i] = 1.0;
               else
                  q_data[i] = 0.0;
            }
         }
         PETSC_CALL(VecRestoreArray(q, &q_data));
      }
      PAMPA_CALL(petsc::normalize(q, power(t)), "unable to normalize the volumetric heat source");
   }
   
   /* Solve the linear system until convegence: */
   int i = 0, imin = 5, imax = 100;
   bool converged = false;
   while (!converged) {
      
      /* Build the coefficient matrix and the RHS vector: */
      PAMPA_CALL(buildMatrix(n, dt, t), 
         "unable to build the coefficient matrix and the RHS vector");
      
      /* Manage the KSP context: */
      if (ksp != 0) {
         PAMPA_CALL(petsc::destroy(ksp), "unable to destroy the KSP context");
      }
      PAMPA_CALL(petsc::create(ksp, A), "unable to create the KSP context");
      
      /* Solve the linear system: */
      PAMPA_CALL(petsc::solve(ksp, b, T), "unable to solve the linear system");
      
      /* Evaluate the convergence: */
      converged = true;
      if (nonlinear) {
         if (Tprev == 0) {
            PETSC_CALL(VecDuplicate(T, &Tprev));
            converged = false;
            if (!microscale)
               mpi::print("Temperature convergence initialized.", true);
         }
         else {
            double eps;
            PAMPA_CALL(petsc::difference(T, Tprev, p, eps, false), 
               "unable to calculate the convergence error");
            converged = eps < tol;
            if (!microscale) {
               mpi::print("Temperature convergence: ", true);
               mpi::print("   - error: " + std::to_string(eps), true);
               mpi::print("   - tolerance: " + std::to_string(tol), true);
               mpi::print("   - converged: " + std::to_string(converged), true);
            }
         }
         PETSC_CALL(VecCopy(T, Tprev));
      }
      
      i++;
      if (i == imax) break;
      if (i < imin) converged = false;
      
   }
   
   /* Print progress: */
   if (!microscale)
      mpi::print("Done.", true);
   
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
   PetscScalar *b_data, *q_data, *qfull_data, *T_data, *Tfull_data, *Tprev_data, *T0_data;
   PETSC_CALL(VecGetArray(b, &b_data));
   PETSC_CALL(VecGetArray(q, &q_data));
   PETSC_CALL(VecGetArray(T, &T_data));
   if (Tprev != 0) {
      PETSC_CALL(VecGetArray(Tprev, &Tprev_data));
   }
   if (multiscale) {
      PETSC_CALL(VecGetArray(qfull, &qfull_data));
      PETSC_CALL(VecGetArray(Tfull, &Tfull_data));
   }
   if (n > 0) {
      PETSC_CALL(VecGetArray(T0, &T0_data));
   }
   
   /* Get the effective thermal properties from the materials or from the microscale problems: */
   Array1D<double> k(num_cells), rho(num_cells), cp(num_cells);
   if (!multiscale) {
      for (int i = 0; i < num_cells; i++) {
         const Material* mat = materials(cells.materials(i));
         k(i) = mat->k(T_data[i]);
         rho(i) = mat->rho(T_data[i]);
         cp(i) = mat->cp(T_data[i]);
      }
   }
   else {
      double w = 0.5;
      if (Tprev != 0)
         for (int i = 0; i < num_cells; i++)
            T_data[i] = w*T_data[i] + (1.0-w)*Tprev_data[i];
      if (s.empty()) {
         s.resize(num_cells, 0.0);
         smicro.resize(num_cells, 0.0);
         smacro.resize(num_cells, 0.0);
         fb.resize(num_cells, 1.0);
      }
      int ic = 0;
      for (int i = 0; i < num_cells; i++) {
         HeatConductionSolver& solver = solvers(cells.materials(i));
         int num_boundaries = (solver.mesh)->getNumBoundaries();
         if (num_boundaries == 1)
            solver.bcs(1).x = Function(fb(i)*T_data[i]);
         else {
            PAMPA_CHECK(faces.num_faces(i) < num_boundaries, 1, 
               "the number of microscale boundaries is larger than the number of macroscale faces");
            for (int f = 0; f < faces.num_faces(i); f++) {
               if (f < num_boundaries) {
                  int i2 = faces.neighbours(i, f);
                  if (i2 < 0)
                     solver.bcs(f+1).x = Function(bcs(-i2).x(t));
                  else
                     solver.bcs(f+1).x = Function(fb(i)*0.5*(T_data[i]+T_data[i2]));
               }
            }
         }
         // double qi = q_data[i] + smicro(i) - smacro(i);
         double qi = q_data[i];
         solver.power = Function(qi);
         PAMPA_CALL(solver.solve(), "unable to solve the microscale problem");
         PAMPA_CALL(solver.getEffectiveThermalProperties(k(i), rho(i), cp(i), s(i), smicro(i), fb(i), t, ic, qfull_data, Tfull_data), 
            "unable to get the effective thermal properties from the microscale problem");
         
         double dTdx = 0.5;
         for (int f = 0; f < faces.num_faces(i); f++) {
            double dT;
            int i2 = faces.neighbours(i, f);
            if (i2 < 0)
               dT = T_data[i] - bcs(-i2).x(t);
            else
               dT = T_data[i] - 0.5*(T_data[i]+T_data[i2]);
            double dx = math::distance(cells.centroids(i), faces.centroids(i, f), 3);
            dTdx += dT/dx * faces.areas(i, f);
         }
         smacro(i) = k(i) * dTdx;
         
         // s(i) = s(i) + smicro(i) - smacro(i);
         
      }
      w = 0.5;
      if (s0.empty())
         s0.resize(num_cells);
      else
         for (int i = 0; i < num_cells; i++)
            s(i) = w*s(i) + (1.0-w)*s0(i);
      for (int i = 0; i < num_cells; i++)
         s0(i) = s(i);
   }
   
   /* Calculate the coefficients for each cell i: */
   for (int i = 0; i < num_cells; i++) {
      
      /* Set a fixed temperature: */
      if (!multiscale) {
         if (!(fixed_temperatures(cells.materials(i)).empty())) {
            b_data[i] = fixed_temperatures(cells.materials(i))(t);
            double a = 1.0;
            PETSC_CALL(MatSetValues(A, 1, &(cells.indices(i)), 1, &(cells.indices(i)), &a, 
               INSERT_VALUES));
            continue;
         }
      }
      
      /* Get the material for cell i: */
      const Material* mat = materials(cells.materials(i));
      
      /* Set the volumetric heat source: */
      b_data[i] = q_data[i];
      if (multiscale)
         b_data[i] -= s(i);
      
      /* Set the time-derivative term: */
      a_i2[0] = cells.indices(i);
      a_i_i2[0] = 0.0;
      if (n > 0) {
         
         /* Get the time-derivative term: */
         double a = rho(i) * cp(i) * cells.volumes(i) / dt;
         
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
                  a = w * k(i) * faces.areas(i, f);
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
                  PAMPA_CHECK(true, 2, "boundary condition not implemented");
                  
                  break;
                  
               }
               
            }
            
         }
         
         /* Set the cell-to-cell coupling terms depending on the neighbour material: */
         else {
            
            /* Treat fixed-temperature materials as Dirichlet boundary conditions: */
            if (!(fixed_temperatures(cells.materials(i2)).empty())) {
               
               /* Get the surface leakage factor: */
               double w = math::surface_leakage_factor(cells.centroids(i), 
                              faces.centroids(i, f), faces.normals(i, f));
               
               /* Set the leakage term for cell i: */
               a = w * k(i) * faces.areas(i, f);
               a_i_i2[0] += a;
               
               /* Set the leakage term for cell i in the RHS vector: */
               b_data[i] += a * fixed_temperatures(cells.materials(i2))(t);
               
               continue;
               
            }
            
            /* Get the material for cell i2: */
            const Material* mat2 = materials(cells.materials(i2));
            
            /* Set the terms for cells with the same materials: */
            if (mat2 == mat && false) {
               
               /* Get the surface leakage factor: */
               double w = math::surface_leakage_factor(cells.centroids(i), cells.centroids(i2), 
                             faces.normals(i, f));
               
               /* Get the leakage term for cell i2: */
               a = -w * k(i) * faces.areas(i, f);
               
            }
            
            /* Set the terms for cells with different materials: */
            else {
               
               /* Get the surface leakage factor and the weight for cell i: */
               double w_i_i2 = math::surface_leakage_factor(cells.centroids(i), 
                                  faces.centroids(i, f), faces.normals(i, f));
               w_i_i2 *= k(i) * faces.areas(i, f);
               
               /* Get the surface leakage factor and the weight for cell i2: */
               double w_i2_i = math::surface_leakage_factor(cells.centroids(i2), 
                                  faces.centroids(i, f), faces.normals(i, f));
               w_i2_i *= -k(i2) * faces.areas(i, f);
               
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
   
   if (!microscale) {
      double qtot = 0.0;
      double ptot = 0.0;
      for (int i = 0; i < num_cells; i++) {
         if (fixed_temperatures(cells.materials(i)).empty()) {
            for (int f = 0; f < faces.num_faces(i); f++) {
               int i2 = faces.neighbours(i, f);
               if (i2 < 0) {
                  double dT = T_data[i] - bcs(-i2).x(t);
                  double dx = math::distance(cells.centroids(i), faces.centroids(i, f), 3);
                  qtot += k(i) * dT/dx * faces.areas(i, f);
               }
               else if (!multiscale && !(fixed_temperatures(cells.materials(i2)).empty())) {
                  double dT = T_data[i] - fixed_temperatures(cells.materials(i2))(t);
                  double dx = math::distance(cells.centroids(i), faces.centroids(i, f), 3);
                  qtot += k(i) * dT/dx * faces.areas(i, f);
               }
            }
            if (multiscale)
               qtot += s(i);
            ptot += q_data[i];
         }
      }
      std::cout << "qtot = " << qtot << ", ptot = " << ptot << std::endl;
   }
   
   /* Restore the arrays with the raw data: */
   PETSC_CALL(VecRestoreArray(b, &b_data));
   PETSC_CALL(VecRestoreArray(q, &q_data));
   PETSC_CALL(VecRestoreArray(T, &T_data));
   if (Tprev != 0) {
      PETSC_CALL(VecRestoreArray(Tprev, &Tprev_data));
   }
   if (multiscale) {
      PETSC_CALL(VecRestoreArray(qfull, &qfull_data));
      PETSC_CALL(VecRestoreArray(Tfull, &Tfull_data));
   }
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
   nonlinear = multiscale;
   for (int i = 0; i < materials.size(); i++) {
      const Material* mat = materials(i);
      PAMPA_CHECK(!(mat->hasThermalProperties()), 1, "missing thermal properties");
      nonlinear |= !(mat->hasConstantThermalProperties());
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
   
   /* Create the temperature vector: */
   PAMPA_CALL(petsc::create(T, A, vectors), "unable to create the temperature vector");
   fields.pushBack(Field{"temperature", &T, false, true});
   VecSet(T, 1200.0);
   
   if (multiscale) {
      int num_cells_full = meshes(meshes.size()-1)->getNumCells();
      int num_cells_global_full = meshes(meshes.size()-1)->getNumCellsGlobal();
      PAMPA_CALL(petsc::create(qfull, num_cells_full, num_cells_global_full, vectors), 
         "unable to create the multiscale heat-source vector");
      PAMPA_CALL(petsc::create(Tfull, num_cells_full, num_cells_global_full, vectors), 
         "unable to create the multiscale temperature vector");
   }
   
   /* Get a random volumetric heat source: */
   if (power.empty()) {
      PAMPA_CALL(petsc::random(q), "unable to initialize the volumetric heat source");
      PAMPA_CALL(petsc::normalize(q, 1.0), "unable to normalize the volumetric heat source");
   }
   
   /* Create the solvers for the microscale problems: */
   for (int i = 1; i < meshes.size()-1; i++) {
      solvers.pushBack(HeatConductionSolver(Array1D<Mesh*>{1, meshes(i)}, materials));
      int num_boundaries = meshes(i)->getNumBoundaries();
      solvers(i-1).bcs = Array1D<BoundaryCondition>{1+num_boundaries, 
         BoundaryCondition{BC::DIRICHLET, 0.0, 0.0}};
      solvers(i-1).fixed_temperatures = fixed_temperatures;
      solvers(i-1).tol = tol;
      solvers(i-1).p = p;
      solvers(i-1).nonlinear = false;
      solvers(i-1).microscale = true;
      PAMPA_CALL(solvers(i-1).build(), "unable to build the microscale-problem solver");
   }
   
   return 0;
   
}

/* Print the solution summary to standard output: */
int HeatConductionSolver::printLog(int n) const {
   
   /* Print out the minimum and maximum temperatures: */
   PetscScalar T_min, T_max;
   PETSC_CALL(VecMin(T, nullptr, &T_min));
   PETSC_CALL(VecMax(T, nullptr, &T_max));
   mpi::print("T_min", T_min);
   mpi::print("T_max", T_max);
   if (multiscale) {
      PETSC_CALL(VecMin(Tfull, nullptr, &T_min));
      PETSC_CALL(VecMax(Tfull, nullptr, &T_max));
      mpi::print("Tfull_min", T_min);
      mpi::print("Tfull_max", T_max);
   }
   
   return 0;
   
}

/* Write the solution to a plain-text file in .vtk format: */
int HeatConductionSolver::writeVTK(const std::string& filename) const {
   
   /* Write the temperature in .vtk format: */
   PAMPA_CALL(vtk::write(filename, "temperature", T, num_cells), "unable to write the temperature");
   
   /* Write the heat source in .vtk format: */
   PAMPA_CALL(vtk::write(filename, "heat-source", q, num_cells), "unable to write the heat source");
   
   /* Write the full temperature field in .vtk format for multiscale problems: */
   if (multiscale) {
      int num_cells_full = meshes(meshes.size()-1)->getNumCells();
      PAMPA_CALL(meshes(meshes.size()-1)->writeVTK("output_full.vtk"), 
         "unable to write the mesh in .vtk format");
      PAMPA_CALL(vtk::write("output_full.vtk", "temperature", Tfull, num_cells_full), 
         "unable to write the temperature");
      PAMPA_CALL(vtk::write("output_full.vtk", "heat-source", qfull, num_cells_full), 
         "unable to write the heat source");
   }
   
   return 0;
   
}

/* Write the solution to a binary file in PETSc format: */
int HeatConductionSolver::writePETSc(int n) const {
   
   /* Write the temperature in PETSc format: */
   std::string filename = "temperature_" + std::to_string(n) + ".ptc";
   PAMPA_CALL(petsc::write(filename, T), "unable to write the temperature");
   
   return 0;
   
}

/* Get the effective thermal properties: */
int HeatConductionSolver::getEffectiveThermalProperties(double& k, double& rho, 
   double& cp, double& s, double& smicro, double& fb, double t, int& ic, PetscScalar* qfull_data, 
   PetscScalar* Tfull_data) const {
   
   /* Get the array with the raw temperature data: */
   PetscScalar *q_data, *T_data;
   PETSC_CALL(VecGetArray(q, &q_data));
   PETSC_CALL(VecGetArray(T, &T_data));
   
   /* Volume-average the properties for all cells: */
   k = 0.0; rho = 0.0; cp = 0.0; s = 0.0; smicro = 0.0;
   double v = 0.0, a = 0.0, Tv = 0.0, Ta = 0.0;
   for (int i = 0; i < num_cells; i++) {
      
      qfull_data[ic] = q_data[i];
      Tfull_data[ic++] = T_data[i];
      
      if (!(fixed_temperatures(cells.materials(i)).empty()))
         continue;
      
      for (int f = 0; f < faces.num_faces(i); f++) {
         int i2 = faces.neighbours(i, f);
         if (i2 < 0) {
            double ki = materials(cells.materials(i))->k(T_data[i]);
            double dT = T_data[i] - bcs(-i2).x(t);
            double dx = math::distance(cells.centroids(i), faces.centroids(i, f), 3);
            smicro += ki * dT/dx * faces.areas(i, f);
            a += faces.areas(i, f);
            Ta += bcs(-i2).x(t) * faces.areas(i, f);
         }
         else if (!(fixed_temperatures(cells.materials(i2)).empty())) {
            double ki = materials(cells.materials(i))->k(T_data[i]);
            double dT = T_data[i] - fixed_temperatures(cells.materials(i2))(t);
            double dx = math::distance(cells.centroids(i), faces.centroids(i, f), 3);
            s += ki * dT/dx * faces.areas(i, f);
         }
      }
      
      /* Get the material for cell i: */
      const Material* mat = materials(cells.materials(i));
      
      /* Get the temperature-dependent contribution for all properties: */
      k += mat->k(T_data[i]) * cells.volumes(i);
      rho += mat->rho(T_data[i]) * cells.volumes(i);
      cp += mat->cp(T_data[i]) * cells.volumes(i);
      
      /* Get the volume: */
      v += cells.volumes(i);
      Tv += T_data[i] * cells.volumes(i);
      
   }
   
   /* Get the final volume-average for all properties: */
   k /= v;
   rho /= v;
   cp /= v;
   fb = (Ta/a) / (Tv/v);
   
   /* Restore the array with the raw temperature data: */
   PETSC_CALL(VecRestoreArray(q, &q_data));
   PETSC_CALL(VecRestoreArray(T, &T_data));
   
   return 0;
   
}

/* Get the volume for fuel materials: */
double HeatConductionSolver::getFuelVolume() const {
   
   double vol = 0.0;
   for (int i = 0; i < num_cells; i++)
      if ((materials(cells.materials(i)))->isFuel())
         vol += cells.volumes(i);
   
   return vol;
   
}
