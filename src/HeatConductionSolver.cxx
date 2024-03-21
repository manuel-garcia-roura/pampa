#include "HeatConductionSolver.hxx"

/* Initialize: */
int HeatConductionSolver::initialize(int argc, char* argv[]) {
   
   /* Check the material data: */
   PAMPA_CALL(checkMaterials(), "wrong material data");
   
   /* Initialize PETSc: */
   static char help[] = "Solver for the linear system A*x = b.\n";
   PETSC_CALL(PetscInitialize(&argc, &argv, (char*)0, help));
   
   /* Build the coefficient matrix and the solution and RHS vectors: */
   PAMPA_CALL(build(), "unable to build the solver");
   
   /* Create the KSP context: */
   PETSC_CALL(KSPCreate(MPI_COMM_WORLD, &ksp));
   PETSC_CALL(KSPSetOperators(ksp, A, A););
   PETSC_CALL(KSPSetFromOptions(ksp));
   
   return 0;
   
}

/* Solve the linear system to get the solution: */
int HeatConductionSolver::solve() {
   
   /* Solve the linear system: */
   double t1 = MPI_Wtime();
   PETSC_CALL(KSPSolve(ksp, q, T));
   double t2 = MPI_Wtime();
   
   return 0;
   
}

/* Output the solution: */
int HeatConductionSolver::output(const std::string& filename) {
   
   /* Write to a rank directory in parallel runs: */
   std::string path;
   if (mpi::size > 1) {
      std::string dir = std::to_string(mpi::rank);
      path = dir + "/" + filename;
   }
   else
      path = filename;
   
   /* Write the solution to a plain-text file in .vtk format: */
   PAMPA_CALL(writeVTK(path), "unable to output the solution in .vtk format");
   
   /* Write the solution to a binary file in PETSc format: */
   PetscBool write, flag;
   PETSC_CALL(PetscOptionsGetBool(NULL, NULL, "-petsc_write_solution", &write, &flag));
   if (flag && write) {
      PAMPA_CALL(petsc::write("temperature.ptc", T), "unable to write the solution");
   }
   
   return 0;
   
}

/* Finalize: */
int HeatConductionSolver::finalize() {
   
   /* Destroy the KSP context: */
   PETSC_CALL(KSPDestroy(&ksp));
   
   /* Destroy the heat-source vector: */
   PETSC_CALL(VecDestroy(&q));
   
   /* Destroy the temperature vector: */
   PETSC_CALL(VecDestroy(&T));
   
   /* Destroy the coefficient matrix: */
   PETSC_CALL(MatDestroy(&A));
   
   /* Finalize PETSCc: */
   PETSC_CALL(PetscFinalize());
   
   return 0;
   
}

/* Check the material data: */
int HeatConductionSolver::checkMaterials() {
   
   /* Check the materials: */
   for (int i = 0; i < materials.size(); i++) {
      PAMPA_CHECK(materials(i).rho < 0.0, 1, "missing density");
      PAMPA_CHECK(materials(i).cp < 0.0, 1, "missing specific heat capacity");
      PAMPA_CHECK(materials(i).k < 0.0, 1, "missing thermal conductivity");
   }
   
   return 0;
   
}

/* Build the coefficient matrix and the solution and RHS vectors: */
int HeatConductionSolver::build() {
   
   /* Build the coefficient matrix and the RHS vector: */
   PAMPA_CALL(buildMatrix(), "unable to build the coefficient matrix");
   
   /* Create the temperature vector: */
   PAMPA_CALL(petsc::create_vector(T, A), "unable to create the temperature vector");
   
   return 0;
   
}

/* Build the coefficient matrix and the RHS vector: */
int HeatConductionSolver::buildMatrix() {
   
   /* Get the mesh data: */
   int num_cells = mesh->getNumCells();
   int num_cells_global = mesh->getNumCellsGlobal();
   int num_faces_max = mesh->getNumFacesMax();
   const Cells& cells = mesh->getCells();
   const Faces& faces = mesh->getFaces();
   const Array1D<BoundaryCondition>& bcs = mesh->getBoundaryConditions();
   
   /* Create, preallocate and set up the coefficient matrix: */
   PAMPA_CALL(petsc::create_matrix(A, num_cells, num_cells_global, 1+num_faces_max), 
      "unable to create the coefficient matrix");
   
   /* Create the heat-source vector: */
   PAMPA_CALL(petsc::create_vector(q, A, true), "unable to create the heat-source vector");
   
   /* Initialize the matrix rows for A: */
   PetscInt a_i2[1+num_faces_max];
   PetscScalar a_i_i2[1+num_faces_max];
   
   /* Calculate the coefficients for each cell i: */
   for (int i = 0; i < num_cells; i++) {
      
      /* Get the material for cell i: */
      const Material& mat = materials(cells.materials(i));
      
      /* Set the cell-to-cell coupling terms: */
      a_i2[0] = cells.indices(i);
      a_i_i2[0] = 0.0;
      int a_i = 1;
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
                  PETSC_CALL(VecSetValue(q, cells.indices(i), a*bcs(-i2).x, ADD_VALUES));
                  
                  break;
                  
               }
               
               /* Set reflective boundary conditions (nothing to be done): */
               case BC::REFLECTIVE : {
                  
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
   
   /* Assembly the coefficient matrix: */
   PETSC_CALL(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
   PETSC_CALL(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
   
   /* Assembly the heat-source vector: */
   PETSC_CALL(VecAssemblyBegin(q));
   PETSC_CALL(VecAssemblyEnd(q));
   
   return 0;
   
}

/* Write the solution to a plain-text file in .vtk format: */
int HeatConductionSolver::writeVTK(const std::string& filename) const {
   
   /* Get the number of cells: */
   int num_cells = mesh->getNumCells();
   
   /* Get the array for the temperature: */
   PetscScalar* data_T;
   PETSC_CALL(VecGetArray(T, &data_T));
   
   /* Write the mesh: */
   PAMPA_CALL(mesh->writeVTK(filename), "unable to write the mesh");
   
   /* Open the output file: */
   std::ofstream file(filename, std::ios_base::app);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Write the temperature: */
   file << "SCALARS temperature double 1" << std::endl;
   file << "LOOKUP_TABLE default" << std::endl;
   for (int i = 0; i < num_cells; i++)
      file << data_T[i] << std::endl;
   file << std::endl;
   
   /* Restore the array for the temperature: */
   PETSC_CALL(VecRestoreArray(T, &data_T));
   
   return 0;
   
}
