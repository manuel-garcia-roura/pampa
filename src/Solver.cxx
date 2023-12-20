#include "Solver.hxx"

/* The Solver constructor: */
Solver::Solver() {};

/* The Solver destructor: */
Solver::~Solver() {};

/* Initialize: */
int Solver::initialize(int argc, char* argv[], const Model &model) {
   
   /* Initialize SLEPc: */
   static char help[] = "Solver for the generalized eigensystem R*phi = (1/keff)*F*phi.\n";
   PETSC_CALL(SlepcInitialize(&argc, &argv, (char*)0, help), "unable to initialize SLEPc");
   
   /* Create the coefficient matrices: */
   PETSC_CALL(MatCreate(PETSC_COMM_WORLD, &R), "unable to create R");
   PETSC_CALL(MatSetFromOptions(R), "unable to set R options");
   PETSC_CALL(MatCreate(PETSC_COMM_WORLD, &F), "unable to create F");
   PETSC_CALL(MatSetFromOptions(F), "unable to set F options");
   
   /* Build the coefficient matrices: */
   PAMPA_CALL(build_matrices(model), "unable to build R and F matrices");
   
   /* Create the solution vector: */
   PETSC_CALL(MatCreateVecs(R, NULL, &phi), "unable to create phi");
   
   /* Create the EPS context: */
   PETSC_CALL(EPSCreate(PETSC_COMM_WORLD, &eps), "unable to create EPS");
   
   /* Set the operators for a generalized eigensystem: */
   PETSC_CALL(EPSSetOperators(eps, R, F), "unable to set EPS operators");
   
   /* Set the solver parameters: */
   PETSC_CALL(EPSSetFromOptions(eps), "unable to set EPS options");
   
   return 0;
   
};

/* Solve the eigensystem to get the neutron flux and the multiplication factor: */
int Solver::solve() {
   
   /* Solve the eigensystem: */
   PETSC_CALL(EPSSolve(eps), "unable to solve the eigensystem");
   
   /* Get solver information: */
   ST st;
   KSP ksp;
   EPSType eps_type;
   PetscInt num_eigenvalues, max_eps_iterations, num_eps_iterations, num_ksp_iterations;
   PetscReal eps_tol;
   PETSC_CALL(EPSGetST(eps, &st), "unable to get ST");
   PETSC_CALL(STGetKSP(st, &ksp), "unable to get KSP");
   PETSC_CALL(EPSGetType(eps, &eps_type), "unable to get EPS type");
   PETSC_CALL(EPSGetDimensions(eps, &num_eigenvalues, NULL, NULL), "unable to get EPS dimensions");
   PETSC_CALL(EPSGetTolerances(eps, &eps_tol, &max_eps_iterations), "unable to get EPS tolerances");
   PETSC_CALL(EPSGetIterationNumber(eps, &num_eps_iterations), "unable to get EPS iterations");
   PETSC_CALL(KSPGetTotalIterations(ksp, &num_ksp_iterations), "unable to get KSP iterations");
   
   /* Print out solver information: */
   std::cout << "Solution method: " << eps_type << std::endl;
   std::cout << "Number of requested eigenvalues: " << num_eigenvalues << std::endl;
   std::cout << "EPS convergence tolerance: " << eps_tol << std::endl;
   std::cout << "Maximum number of EPS iterations: " << max_eps_iterations << std::endl;
   std::cout << "Number of EPS iterations: " << num_eps_iterations << std::endl;
   std::cout << "Number of KSP iterations: " << num_ksp_iterations << std::endl;
   
   /* Print out the multiplication factor: */
   double lambda;
   PETSC_CALL(EPSGetEigenpair(eps, 0, &lambda, NULL, phi, NULL), "");
   keff = 1.0 / lambda;
   std::cout << "keff = " << keff << std::endl;
   
   return 0;
   
};

/* Output the solution: */
int Solver::output(const std::string &filename, const Model &model) {
   
   /* Open the output file: */
   std::ofstream file(filename, std::ios_base::app);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Get the model mesh: */
   const Mesh *mesh = model.getMesh();
   
   /* Get the number of cells and energy groups: */
   int num_cells = mesh->getNumCells();
   int num_groups = model.getNumEnergyGroups();
   
   /* Get the raw data to the solution vector: */
   PetscScalar *data;
   PETSC_CALL(VecGetArray(phi, &data), "unable to get the solution array");
   
   /* Write the solution: */
   for (int g = 0; g < num_groups; g++) {
      // file << "CELL_DATA " << num_cells << std::endl;
      file << "SCALARS flux" << (g+1) << " double 1" << std::endl;
      file << "LOOKUP_TABLE default" << std::endl;
      for (int i = 0; i < num_cells; i++)
         file << data[i*num_groups+g] << std::endl;
      file << std::endl;
   }
   
   /* Restore the solution vector: */
   PETSC_CALL(VecRestoreArray(phi, &data), "unable to restore the solution array");
   
   return 0;
   
};

/* Finalize: */
int Solver::finalize() {
   
   /* Destroy the EPS context: */
   PETSC_CALL(EPSDestroy(&eps), "unable to destroy EPS");
   
   /* Destroy the solution vector: */
   PETSC_CALL(VecDestroy(&phi), "unable to destroy phi");
   
   /* Destroy the coefficient matrices: */
   PETSC_CALL(MatDestroy(&R), "unable to destroy R");
   PETSC_CALL(MatDestroy(&F), "unable to destroy F");
   
   /* Finalize SLEPc: */
   PETSC_CALL(SlepcFinalize(), "unable to finalize SLEPc");
   
   return 0;
   
};

/* Build the coefficient matrices: */
int Solver::build_matrices(const Model &model) {
   
   /* Get the model mesh: */
   const Mesh *mesh = model.getMesh();
   
   /* Get the number of cells and energy groups: */
   int num_cells = mesh->getNumCells();
   int num_groups = model.getNumEnergyGroups();
   
   /* Set up the coefficient matrices: */
   int n = num_cells * num_groups;
   PETSC_CALL(MatSetSizes(R, PETSC_DECIDE, PETSC_DECIDE, n, n), "unable to set R size");
   PETSC_CALL(MatSetUp(R), "unable to set up R");
   PETSC_CALL(MatSetSizes(F, PETSC_DECIDE, PETSC_DECIDE, n, n), "unable to set F size");
   PETSC_CALL(MatSetUp(F), "unable to set up F");
   
   /* Calculate the coefficients for each cell i and group g: */
   for (int i = 0; i < num_cells; i++) {
      for (int g = 0; g < num_groups; g++) {
         
         /* Get the matrix index for cell i and group g: */
         int l = i*num_groups + g;
         
         /* Set the total-reaction term: */
         double r_l_l = mesh->getVolumeIntegralSigmaRemoval(i, g);
         
         /* Set the group-to-group coupling terms: */
         for (int g2 = 0; g2 < num_groups; g2++) {
            
            /* Get the matrix index for cell i and group g2: */
            int l2 = i*num_groups + g2;
            
            /* Set the (g2 -> g) scattering term: */
            if (g2 != g) {
               double r_l_l2 = -(mesh->getVolumeIntegralSigmaScattering(i, g2, g));
               PETSC_CALL(MatSetValues(R, 1, &l, 1, &l2, &r_l_l2, INSERT_VALUES), 
                  "unable to set R(l, l2)");
            }
            
            /* Set the (g2 -> g) fission term: */
            double f_l_l2 = mesh->getChi(i, g) * mesh->getVolumeIntegralSigmaNuFission(i, g2);
            PETSC_CALL(MatSetValues(F, 1, &l, 1, &l2, &f_l_l2, INSERT_VALUES), 
               "unable to set F(l, l2)");
            
         }
         
         /* Set the cell-to-cell coupling terms: */
         std::vector<int> neighbours = mesh->getFaceNeighbours(i);
         double r_l_l2;
         for (int f = 0; f < neighbours.size(); f++) {
            
            /* Get the index for cell i2 (internal or boundary condition): */
            int i2 = neighbours[f];
            
            /* Set vacuum (zero-flux) boundary conditions: */
            if (i2 == bc::vacuum) {
               
               /* Get the geometrical data: */
               const std::vector<double> &p_i = mesh->getCellCentroid(i);
               const std::vector<double> &p_f = mesh->getFaceCentroid(i, f);
               const std::vector<double> &n_i_f = mesh->getFaceNormal(i, f);
               
               /* Get the surface leakage factor: */
               double w = math::surface_leakage_factor(p_i, p_f, n_i_f);
               
               /* Set the leakage term for cell i: */
               r_l_l += w * mesh->getSurfaceIntegralDiffusionCoefficient(i, f, g);
               
            }
            
            /* Set reflective (zero-current) boundary conditions: */
            else if (i2 == bc::reflective) {
               
               /* Nothing to be done: */
               continue;
               
            }
            
            /* Set Robin boundary conditions (TODO: implement!): */
            else if (i2 == bc::robin) {
               
               /* Not implemented: */
               PAMPA_CHECK(true, 1, "Robin boundary conditions not implemented yet");
               
            }
            
            /* Set cell-to-cell coupling terms depending on the neighbour material: */
            else {
               
               /* Get the matrix index for cell i2 and group g: */
               int l2 = i2*num_groups + g;
               
               /* Get the materials for cell i and i2: */
               int mat = mesh->getMaterial(i);
               int mat2 = mesh->getMaterial(i2);
               
               /* Set the terms for cells with the same materials: */
               if (mat2 == mat) {
                  
                  /* Get the geometrical data: */
                  const std::vector<double> &p_i = mesh->getCellCentroid(i);
                  const std::vector<double> &p_i2 = mesh->getCellCentroid(i2);
                  const std::vector<double> &n_i_f = mesh->getFaceNormal(i, f);
                  
                  /* Get the surface leakage factor: */
                  double w = math::surface_leakage_factor(p_i, p_i2, n_i_f);
                  
                  /* Get the leakage term for cell i2: */
                  r_l_l2 = -w * mesh->getSurfaceIntegralDiffusionCoefficient(i, f, g);
                  
               }
               
               /* Set the terms for cells with different materials: */
               else {
                  
                  /* Get the geometrical data: */
                  const std::vector<double> &p_i = mesh->getCellCentroid(i);
                  const std::vector<double> &p_i2 = mesh->getCellCentroid(i2);
                  const std::vector<double> &p_f = mesh->getFaceCentroid(i, f);
                  const std::vector<double> &n_i_f = mesh->getFaceNormal(i, f);
                  
                  /* Get the surface leakage factor and the weight for cell i: */
                  double w_i_i2 = math::surface_leakage_factor(p_i, p_f, n_i_f);
                  w_i_i2 *= mesh->getSurfaceIntegralDiffusionCoefficient(i, f, g);
                  
                  /* Get the surface leakage factor and the weight for cell i2: */
                  double w_i2_i = math::surface_leakage_factor(p_i2, p_f, n_i_f);
                  w_i2_i *= -mesh->getSurfaceIntegralDiffusionCoefficient(i2, i, f, g);
                  
                  /* Get the leakage term for cell i2: */
                  r_l_l2 = -(w_i_i2*w_i2_i) / (w_i_i2+w_i2_i);
                  
               }
               
               /* Set the leakage term for cell i2: */
               PETSC_CALL(MatSetValues(R, 1, &l, 1, &l2, &r_l_l2, INSERT_VALUES), 
                  "unable to set R(l, l2)");
               
               /* Set the leakage term for cell i: */
               r_l_l -= r_l_l2;
               
            }
            
         }
         
         /* Set the diagonal coefficient: */
         PETSC_CALL(MatSetValues(R, 1, &l, 1, &l, &r_l_l, INSERT_VALUES), "unable to set R(l, l)");
         
      }
   }
   
   /* Assembly the coefficient matrices: */
   PETSC_CALL(MatAssemblyBegin(R, MAT_FINAL_ASSEMBLY), "unable to assembly R");
   PETSC_CALL(MatAssemblyEnd(R, MAT_FINAL_ASSEMBLY), "unable to assembly R");
   PETSC_CALL(MatAssemblyBegin(F, MAT_FINAL_ASSEMBLY), "unable to assembly F");
   PETSC_CALL(MatAssemblyEnd(F, MAT_FINAL_ASSEMBLY), "unable to assembly F");
   
   return 0;
   
};
