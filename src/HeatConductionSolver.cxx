#include "HeatConductionSolver.hxx"

/* Read the solver from a plain-text input file: */
int HeatConductionSolver::read(std::ifstream& file, Array1D<Solver*>& solvers) {
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = input::get_next_line(file);
      if (line.empty() || line[0] == "}") break;
      
      /* Get the next keyword: */
      unsigned int l = 0;
      if (line[l] == "bc") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() < 3, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the mesh boundaries: */
         const Array1D<std::string>& boundaries = mesh->getBoundaries();
         
         /* Initialize the boundary-condition array, if not done yet: */
         if (bcs.empty()) bcs.resize(1+boundaries.size());
         
         /* Get the boundary name and index: */
         std::string name = line[++l];
         int i = boundaries.find(name);
         
         /* Get the boundary condition (1-based indexed): */
         if (i >= 0) {
            PAMPA_CHECK(input::read(bcs(i+1), line, ++l, file), "wrong boundary condition");
            continue;
         }
         
         /* Get the material name if this isn't a boundary: */
         PAMPA_CHECK(utils::find(name, materials, i), "wrong material name");
         bcmat_indices(i) = bcs.size();
         
         /* Get the boundary condition: */
         BoundaryCondition bc;
         PAMPA_CHECK(input::read(bc, line, ++l, file), "wrong boundary condition");
         bcs.pushBack(bc);
         
      }
      else if (line[l] == "power") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() < 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the total power: */
         PAMPA_CHECK(input::read(power, 0.0, DBL_MAX, line, ++l, file), "wrong power data");
         
      }
      else if (line[l] == "convergence") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 5, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the field name: */
         std::string name = line[++l];
         
         /* Get the convergence error: */
         if (name == "temperature") {
            PAMPA_CHECK(input::read(dT, name, line, l), "wrong convergence error");
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

/* Solve the linear system to get the solution: */
int HeatConductionSolver::solve(int n, double dt, double t) {
   
   /* Print progress: */
   output::print("Run '" + name + "' solver...", true);
   
   /* Calculate the volumetric heat source from the nodal power: */
   if (mesh_nodal) {
      PAMPA_CHECK(calculateHeatSource(), "unable to calculate the volumetric heat source");
   }
   
   /* Normalize the volumetric heat source: */
   if (!(power.empty())) {
      PAMPA_CHECK(petsc::normalize(q, power(t)), "unable to normalize the volumetric heat source");
   }
   
   /* Solve the linear system until convegence: */
   bool converged = false;
   while (!converged) {
      
      /* Build the coefficient matrix and the RHS vector: */
      PAMPA_CHECK(buildMatrix(n, dt, t), 
         "unable to build the coefficient matrix and the RHS vector");
      
      /* Manage the KSP context: */
      if (ksp != 0) {
         PAMPA_CHECK(petsc::destroy(ksp), "unable to destroy the KSP context");
      }
      PAMPA_CHECK(petsc::create(ksp, A), "unable to create the KSP context");
      
      /* Solve the linear system: */
      PAMPA_CHECK(petsc::solve(ksp, b, T), "unable to solve the linear system");
      
      /* Evaluate the convergence: */
      converged = true;
      if (nonlinear) {
         PAMPA_CHECK(dT.check(T, converged), "unable to check the convergence");
      }
      
   }
   
   /* Calculate the nodal temperatures: */
   if (mesh_nodal) {
      PAMPA_CHECK(calculateNodalTemperatures(), "unable to calculate the nodal temperatures");
   }
   
   /* Calculate the total heat flows in and out of the system: */
   PAMPA_CHECK(calculateHeatFlows(t), "unable to calculate the total heat flows");
   
   /* Print progress: */
   output::print("Done.", true);
   
   return 0;
   
}

/* Initialize the volumetric heat source: */
int HeatConductionSolver::initializeHeatSource() {
   
   /* Initialize the heat sources to zero: */
   PAMPA_CHECK(petsc::set(q, 0.0), "unable to initialize the volumetric heat source");
   if (mesh_nodal) {
      PAMPA_CHECK(petsc::set(qnodal, 0.0), "unable to initialize the nodal volumetric heat source");
   }
   
   /* Get the arrays with the raw data: */
   PetscScalar *q_data, *qnodal_data;
   PETSC_CALL(VecGetArray(q, &q_data));
   if (mesh_nodal) {
      PETSC_CALL(VecGetArray(qnodal, &qnodal_data));
   }
   
   /* Set a uniform volumetric heat source for fuel materials: */
   for (int i = 0; i < num_cells; i++) {
      if (materials(cells.materials(i))->isFuel()) {
         q_data[i] = cells.volumes(i);
         if (mesh_nodal) {
            int in = cells.nodal_indices(i);
            if (in >= 0)
               qnodal_data[in] += cells.volumes(i);
         }
      }
   }
   
   /* Gather the nodal power: */
   if (mesh_nodal) {
      int num_cells_nodal = mesh_nodal->getNumCells();
      MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, qnodal_data, num_cells_nodal, MPI_DOUBLE, MPI_SUM, 
         MPI_COMM_WORLD));
   }
   
   /* Restore the arrays with the raw data: */
   PETSC_CALL(VecRestoreArray(q, &q_data));
   if (mesh_nodal) {
      PETSC_CALL(VecRestoreArray(qnodal, &qnodal_data));
   }
   
   /* Normalize the heat sources to one: */
   PAMPA_CHECK(petsc::normalize(q, 1.0), "unable to normalize the volumetric heat source");
   if (mesh_nodal) {
      PAMPA_CHECK(petsc::normalize(qnodal, 1.0), 
         "unable to normalize the nodal volumetric heat source");
   }
   
   return 0;
   
}

/* Calculate the volumetric heat source from the nodal power: */
int HeatConductionSolver::calculateHeatSource() {
   
   /* Initialize the heat source to zero: */
   PAMPA_CHECK(petsc::set(q, 0.0), "unable to initialize the volumetric heat source");
   
   /* Get the number of nodal cells: */
   int num_cells_nodal = mesh_nodal->getNumCells();
   
   /* Get the arrays with the raw data: */
   PetscScalar *q_data, *qnodal_data;
   PETSC_CALL(VecGetArray(q, &q_data));
   PETSC_CALL(VecGetArray(qnodal, &qnodal_data));
   
   /* Calculate the volumetric heat source for fuel materials: */
   Array1D<double> vol(num_cells_nodal, 0.0);
   for (int i = 0; i < num_cells; i++) {
      if (materials(cells.materials(i))->isFuel()) {
         int in = cells.nodal_indices(i);
         if (in >= 0) {
            q_data[i] = qnodal_data[in] * cells.volumes(i);
            vol(in) += cells.volumes(i);
         }
      }
   }
   
   /* Gather the nodal volumes: */
   MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, &(vol(0)), num_cells_nodal, MPI_DOUBLE, MPI_SUM, 
      MPI_COMM_WORLD));
   
   /* Normalize the heat source with the nodal volumes: */
   for (int i = 0; i < num_cells; i++) {
      if (materials(cells.materials(i))->isFuel()) {
         int in = cells.nodal_indices(i);
         if (in >= 0)
            q_data[i] /= vol(in);
      }
   }
   
   /* Restore the arrays with the raw data: */
   PETSC_CALL(VecRestoreArray(q, &q_data));
   PETSC_CALL(VecRestoreArray(qnodal, &qnodal_data));
   
   return 0;
   
}

/* Calculate the nodal temperatures: */
int HeatConductionSolver::calculateNodalTemperatures() {
   
   /* Get the number of materials and nodal cells: */
   int num_materials = materials.size();
   int num_cells_nodal = mesh_nodal->getNumCells();
   
   /* Get the arrays with the raw data: */
   PetscScalar* T_data;
   Array1D<PetscScalar*> Tnodal_data(num_materials);
   PETSC_CALL(VecGetArray(T, &T_data));
   for (int im = 0; im < num_materials; im++) {
      if (bcmat_indices(im) < 0 && !(materials(im)->isBC())) {
         PETSC_CALL(VecGetArray(Tnodal(im), &(Tnodal_data(im))));
         for (int in = 0; in < num_cells_nodal; in++)
            Tnodal_data(im)[in] = 0.0;
      }
   }
   
   /* Calculate the contributions to the nodal temperatures for each cell: */
   Array2D<double> vol(num_materials, num_cells_nodal, 0.0);
   for (int i = 0; i < num_cells; i++) {
      int im = cells.materials(i);
      if (bcmat_indices(im) < 0 && !(materials(im)->isBC())) {
         int in = cells.nodal_indices(i);
         if (in >= 0) {
            Tnodal_data(im)[in] += T_data[i] * cells.volumes(i);
            vol(im, in) += cells.volumes(i);
         }
      }
   }
   
   /* Gather the nodal volumes: */
   MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, &(vol(0, 0)), num_materials*num_cells_nodal, MPI_DOUBLE, 
      MPI_SUM, MPI_COMM_WORLD));
   
   /* Gather the nodal temperatures: */
   for (int im = 0; im < num_materials; im++) {
      if (bcmat_indices(im) < 0 && !(materials(im)->isBC())) {
         MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, Tnodal_data(im), num_cells_nodal, MPI_DOUBLE, 
            MPI_SUM, MPI_COMM_WORLD));
      }
   }
   
   /* Normalize the temperatures with the nodal volumes: */
   for (int im = 0; im < num_materials; im++)
      if (bcmat_indices(im) < 0 && !(materials(im)->isBC()))
         for (int in = 0; in < num_cells_nodal; in++)
            if (vol(im, in) > 0.0)
               Tnodal_data(im)[in] /= vol(im, in);
   
   /* Restore the arrays with the raw data: */
   PETSC_CALL(VecRestoreArray(T, &T_data));
   for (int im = 0; im < num_materials; im++) {
      if (bcmat_indices(im) < 0 && !(materials(im)->isBC())) {
         PETSC_CALL(VecRestoreArray(Tnodal(im), &(Tnodal_data(im))));
      }
   }
   
   return 0;
   
}

/* Calculate the total heat flows in and out of the system: */
int HeatConductionSolver::calculateHeatFlows(double t) {
   
   /* Get the total heat source: */
   PETSC_CALL(VecSum(q, &qin));
   
   /* Get the arrays with the raw data: */
   PetscScalar *T_data;
   PETSC_CALL(VecGetArray(T, &T_data));
   
   /* Calculate the heat flow out of the boundaries: */
   qout = 0.0;
   for (int i = 0; i < num_cells; i++) {
      if (bcmat_indices(cells.materials(i)) < 0) {
         for (int f = 0; f < faces.num_faces(i); f++) {
            
            /* Get the material for cell i: */
            const Material* mat = materials(cells.materials(i));
            
            /* Get the index for cell i2 (actual cell or boundary condition): */
            int i2 = faces.neighbours(i, f);
            if (i2 >= 0)
               if (bcmat_indices(cells.materials(i2)) >= 0)
                  i2 = -bcmat_indices(cells.materials(i2));
            
            /* Get the heat flow for Dirichlet and convection boundary conditions: */
            if (i2 < 0) {
               switch (bcs(-i2).type) {
                  case BC::DIRICHLET : {
                     double w = math::surface_leakage_factor(cells.centroids(i), 
                                   faces.centroids(i, f), faces.normals(i, f));
                     qout += w * mat->k(T_data[i]) * (T_data[i]-bcs(-i2).f(0)(t)) * 
                                faces.areas(i, f);
                     break;
                  }
                  case BC::CONVECTION : {
                     qout += bcs(-i2).f(0)(t) * (T_data[i]-bcs(-i2).f(1)(t)) * faces.areas(i, f);
                     break;
                  }
                  default : {
                     break;
                  }
               }
            }
            
         }
      }
   }
   
   /* Gather the heat flow out of the boundaries: */
   MPI_CALL(MPI_Allreduce(MPI_IN_PLACE, &qout, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
   
   /* Restore the arrays with the raw data: */
   PETSC_CALL(VecRestoreArray(T, &T_data));
   
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
      
      /* Check for boundary-condition materials: */
      int ibc = bcmat_indices(cells.materials(i));
      if (ibc >= 0) {
         if (bcs(ibc).type == BC::DIRICHLET)
            b_data[i] = bcs(ibc).f(0)(t);
         else
            b_data[i] = 0.0;
         double a = 1.0;
         PETSC_CALL(MatSetValues(A, 1, &(cells.global_indices(i)), 1, &(cells.global_indices(i)), 
            &a, INSERT_VALUES));
         continue;
      }
      
      /* Get the material for cell i: */
      const Material* mat = materials(cells.materials(i));
      
      /* Set the volumetric heat source: */
      b_data[i] = q_data[i];
      
      /* Set the time-derivative term: */
      a_i2[0] = cells.global_indices(i);
      a_i_i2[0] = 0.0;
      if (n > 0) {
         
         /* Get the time-derivative term: */
         double a = mat->rho(T_data[i]) * mat->cp(T_data[i]) * cells.volumes(i) / dt;
         
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
         if (i2 >= 0)
            if (bcmat_indices(cells.materials(i2)) >= 0)
               i2 = -bcmat_indices(cells.materials(i2));
         
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
                  a = w * mat->k(T_data[i]) * faces.areas(i, f);
                  a_i_i2[0] += a;
                  
                  /* Set the leakage term for cell i in the RHS vector: */
                  b_data[i] += a * bcs(-i2).f(0)(t);
                  
                  break;
                  
               }
               
               /* Set reflective boundary conditions (nothing to be done): */
               case BC::REFLECTIVE : {
                  
                  break;
                  
               }
               
               /* Set adiabatic boundary conditions (nothing to be done): */
               case BC::ADIABATIC : {
                  
                  break;
                  
               }
               
               /* Set convection boundary conditions: */
               case BC::CONVECTION : {
                  
                  /* Set the leakage term for cell i: */
                  a = bcs(-i2).f(0)(t) * faces.areas(i, f);
                  a_i_i2[0] += a;
                  
                  /* Set the leakage term for cell i in the RHS vector: */
                  b_data[i] += a * bcs(-i2).f(1)(t);
                  
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
         
         /* Set the cell-to-cell coupling terms depending on the neighbour material: */
         else {
            
            /* Get the material for cell i2: */
            const Material* mat2 = materials(cells.materials(i2));
            
            /* Set the terms for cells with the same materials: */
            if (mat2 == mat && false) {
               
               /* Get the surface leakage factor: */
               double w = math::surface_leakage_factor(cells.centroids(i), cells.centroids(i2), 
                             faces.normals(i, f));
               
               /* Get the leakage term for cell i2: */
               a = -w * mat->k(T_data[i]) * faces.areas(i, f);
               
            }
            
            /* Set the terms for cells with different materials: */
            else {
               
               /* Get the surface leakage factor and the weight for cell i: */
               double w_i_i2 = math::surface_leakage_factor(cells.centroids(i), 
                                  faces.centroids(i, f), faces.normals(i, f));
               w_i_i2 *= mat->k(T_data[i]) * faces.areas(i, f);
               
               /* Get the surface leakage factor and the weight for cell i2: */
               double w_i2_i = math::surface_leakage_factor(cells.centroids(i2), 
                                  faces.centroids(i, f), faces.normals(i, f));
               w_i2_i *= -mat2->k(T_data[i]) * faces.areas(i, f);
               
               /* Get the leakage term for cell i2: */
               a = -(w_i_i2*w_i2_i) / (w_i_i2+w_i2_i);
               
            }
            
            /* Set the leakage term for cell i: */
            a_i_i2[0] -= a;
            
            /* Set the leakage term for cell i2: */
            a_i2[a_i] = cells.global_indices(i2);
            a_i_i2[a_i++] = a;
            
         }
         
      }
      
      /* Set the matrix rows for A: */
      PETSC_CALL(MatSetValues(A, 1, &(cells.global_indices(i)), a_i, a_i2, a_i_i2, INSERT_VALUES));
      
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
      if (bcmat_indices(i) < 0 && !(materials(i)->isBC())) {
         const Material* mat = materials(i);
         PAMPA_CHECK(!(mat->hasThermalProperties()), "missing thermal properties");
         nonlinear = !(mat->hasConstantThermalProperties());
      }
   }
   
   return 0;
   
}

/* Build the coefficient matrix, and the solution and RHS vectors: */
int HeatConductionSolver::build() {
   
   /* Create, preallocate and set up the coefficient matrix: */
   PAMPA_CHECK(petsc::create(A, num_cells, num_cells_global, 1+num_faces_max, matrices), 
      "unable to create the coefficient matrix");
   
   /* Create the right-hand-side vector: */
   PAMPA_CHECK(petsc::create(b, A, vectors), "unable to create the right-hand-side vector");
   
   /* Create the heat-source vector: */
   PAMPA_CHECK(petsc::create(q, A, vectors), "unable to create the heat-source vector");
   fields.pushBack(Field{"power", &q, true, false, nullptr});
   
   /* Create the temperature vector: */
   PAMPA_CHECK(petsc::create(T, A, vectors), "unable to create the temperature vector");
   fields.pushBack(Field{"temperature", &T, false, true, &dT});
   
   /* Create the nodal vectors: */
   if (mesh_nodal) {
      
      /* Get the number of materials and nodal cells: */
      int num_materials = materials.size();
      int num_cells_nodal = mesh_nodal->getNumCells();
      
      /* Create the nodal heat-source vector: */
      PAMPA_CHECK(petsc::create(qnodal, num_cells_nodal, num_cells_nodal, vectors, true), 
         "unable to create the nodal heat-source vector");
      fields.pushBack(Field{"nodal_power", &qnodal, false, false, nullptr});
      
      /* Create the nodal temperature vector for each non-boundary-condition material: */
      Tnodal.resize(num_materials, 0);
      for (int i = 0; i < num_materials; i++) {
         if (bcmat_indices(i) < 0 && !(materials(i)->isBC())) {
            PAMPA_CHECK(petsc::create(Tnodal(i), num_cells_nodal, num_cells_nodal, vectors, true), 
               "unable to create the nodal temperature vector");
            fields.pushBack(Field{materials(i)->name + "_temperature", &Tnodal(i), false, false, 
               nullptr});
         }
      }
      
   }
   
   /* Initialize the volumetric heat source: */
   PAMPA_CHECK(initializeHeatSource(), "unable to initialize the volumetric heat source");
   
   return 0;
   
}

/* Print the solution summary to standard output: */
int HeatConductionSolver::printLog(int n) const {
   
   /* Print out the minimum and maximum temperatures: */
   PetscScalar T_min, T_max;
   PETSC_CALL(VecMin(T, nullptr, &T_min));
   PETSC_CALL(VecMax(T, nullptr, &T_max));
   output::print("T_min", T_min);
   output::print("T_max", T_max);
   
   /* Print out the minimum and maximum nodal temperatures for each material: */
   if (mesh_nodal) {
      for (int i = 0; i < materials.size(); i++) {
         if (bcmat_indices(i) < 0 && !(materials(i)->isBC())) {
            PetscScalar T_min, T_max;
            PETSC_CALL(VecMin(Tnodal(i), nullptr, &T_min));
            PETSC_CALL(VecMax(Tnodal(i), nullptr, &T_max));
            output::print("T_min (" + materials(i)->name + ")", T_min);
            output::print("T_max (" + materials(i)->name + ")", T_max);
         }
      }
   }
   
   /* Print out the total heat flows in and out of the system: */
   output::print("qin", qin);
   output::print("qout", qout);
   
   return 0;
   
}

/* Write the solution to a plain-text file in .vtk format: */
int HeatConductionSolver::writeVTK(const std::string& filename) const {
   
   /* Write the temperature in .vtk format: */
   PAMPA_CHECK(vtk::write(filename, "temperature", T, num_cells), 
      "unable to write the temperature");
   
   /* Write the volumetric heat source in .vtk format: */
   PAMPA_CHECK(vtk::write(filename, "power", q, num_cells), "unable to write the heat source");
   
   /* Write the nodal temperatures in .vtk format: */
   if (mesh_nodal && (mpi::rank == 0)) {
      
      /* Get the number of nodal cells: */
      int num_cells_nodal = mesh_nodal->getNumCells();
      
      /* Write the nodal mesh in .vtk format: */
      PAMPA_CHECK(mesh_nodal->writeVTK("nodal_output.vtk"), 
         "unable to write the nodal mesh in .vtk format");
      
      /* Write the nodal temperature for each non-boundary-condition material in .vtk format: */
      for (int i = 0; i < materials.size(); i++) {
         if (bcmat_indices(i) < 0 && !(materials(i)->isBC())) {
            PAMPA_CHECK(vtk::write("nodal_output.vtk", materials(i)->name + "_temperature", 
               Tnodal(i), num_cells_nodal), "unable to write the nodal temperature");
         }
      }
      
      /* Write the nodal volumetric heat source in .vtk format: */
      PAMPA_CHECK(vtk::write("nodal_output.vtk", "power", qnodal, num_cells_nodal), 
         "unable to write the nodal heat source");
      
   }
   
   return 0;
   
}

/* Write the solution to a binary file in PETSc format: */
int HeatConductionSolver::writePETSc(int n) const {
   
   /* Write the temperature in PETSc format: */
   std::string filename = "temperature_" + std::to_string(n) + ".ptc";
   PAMPA_CHECK(petsc::write(filename, T), "unable to write the temperature");
   
   return 0;
   
}
