#include "CouplingSolver.hxx"

/* Read the solver from a plain-text input file: */
int CouplingSolver::read(std::ifstream& file, Array1D<Solver*>& solvers) {
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty() || line[0] == "}") break;
      
      /* Get the next keyword: */
      unsigned int l = 0;
      if (line[l] == "coupled-solvers") {
         
         /* Get the number of coupled solvers: */
         int num_coupled_solvers;
         PAMPA_CALL(utils::read(num_coupled_solvers, 1, INT_MAX, line[++l]), 
            "wrong number of coupled solvers");
         
         /* Get the coupled solvers: */
         coupled_solvers.resize(num_coupled_solvers, nullptr);
         for (int i = 0; i < num_coupled_solvers; i++) {
            PAMPA_CALL(utils::find(line[++l], solvers, &(coupled_solvers(i))), 
               "unable to find coupled solver");
         }
         
      }
      else if (line[l] == "implicit") {
         
         /* Get the switch to use implicit coupling: */
         PAMPA_CALL(utils::read(implicit, line[++l]), "wrong switch for implicit coupling");
         
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

/* Initialize: */
int CouplingSolver::initialize(bool transient) {
   
   /* Initialize all the coupled solvers: */
   for (int i = 0; i < coupled_solvers.size(); i++) {
      PAMPA_CALL(coupled_solvers(i)->initialize(transient), "unable to initialize the solver");
   }
   
   /* Get the feedback fields: */
   for (int i = 0; i < coupled_solvers.size(); i++) {
      Array1D<Field>& coupled_fields = coupled_solvers(i)->getFields();
      for (int f = 0; f < coupled_fields.size(); f++)
         fields.pushBack(coupled_fields(f));
   }
   
   return 0;
   
}

/* Get the solution: */
int CouplingSolver::solve(int n, double dt, double t) {
   
   /* Print progress: */
   mpi::print("Run '" + name + "' solver...", true);
   
   /* Iterate the solution until convegence: */
   bool converged = false;
   while (!converged) {
      
      /* Reset the convergence flag: */
      converged = true;
      
      /* Get the solution from all the coupled solvers: */
      for (int i = 0; i < coupled_solvers.size(); i++) {
         
         /* Get the solution: */
         PAMPA_CALL(coupled_solvers(i)->solve(n, dt, t), "unable to get the solution");
         
         /* Exchange the output fields calculated by this solver: */
         Array1D<Field>& output_fields = coupled_solvers(i)->getFields();
         for (int f = 0; f < output_fields.size(); f++) {
            if (output_fields(f).output) {
               
               /* Set the field in the coupled solvers that take it as input: */
               for (int i2 = 0; i2 < coupled_solvers.size(); i2++) {
                  if (i2 != i) {
                     Array1D<Field>& input_fields = coupled_solvers(i2)->getFields();
                     for (int f2 = 0; f2 < input_fields.size(); f2++) {
                        if (input_fields(f2).input) {
                           if (input_fields(f2).name == output_fields(f).name) {
                              PETSC_CALL(VecCopy(*(output_fields(f).vec), *(input_fields(f2).vec)));
                              mpi::print("Field exchange: ", true);
                              mpi::print("   - field name: " + output_fields(f).name, true);
                              mpi::print("   - output solver: " + coupled_solvers(i)->name, true);
                              mpi::print("   - input solver: " + coupled_solvers(i2)->name, true);
                           }
                        }
                     }
                  }
               }
               
               /* Evaluate the convergence: */
               if (implicit) {
                  if (output_fields(f).vec0 == nullptr) {
                     output_fields(f).vec0 = new Vec;
                     PETSC_CALL(VecDuplicate(*(output_fields(f).vec), output_fields(f).vec0));
                     converged = false;
                     mpi::print("Field convergence initialized: ", true);
                     mpi::print("   - field name: " + output_fields(f).name, true);
                  }
                  else {
                     double eps;
                     PAMPA_CALL(petsc::difference(*(output_fields(f).vec), 
                        *(output_fields(f).vec0), p, eps, true), 
                        "unable to calculate the convergence error");
                     converged &= eps < tol;
                     mpi::print("Field convergence: ", true);
                     mpi::print("   - field name: " + output_fields(f).name, true);
                     mpi::print("   - error: " + std::to_string(eps), true);
                     mpi::print("   - tolerance: " + std::to_string(tol), true);
                  }
                  PETSC_CALL(VecCopy(*(output_fields(f).vec), *(output_fields(f).vec0)));
               }
               
            }
         }
         
      }
      
      /* Print progress: */
      mpi::print("Coupled solution convergence: ", true);
      mpi::print("   - converged: " + std::to_string(converged), true);
      
   }
   
   /* Print progress: */
   mpi::print("Done.", true);
   
   return 0;
   
}

/* Output the solution: */
int CouplingSolver::output(const std::string& filename, int n, bool write_mesh) const {
   
   /* Write the mesh in .vtk format: */
   if (write_mesh) {
      PAMPA_CALL(mesh->writeVTK(filename), "unable to write the mesh in .vtk format");
   }
   
   /* Output the solution from all the coupled solvers: */
   for (int i = 0; i < coupled_solvers.size(); i++) {
      PAMPA_CALL(coupled_solvers(i)->output(filename, n, false), "unable to output the solution");
   }
   
   return 0;
   
}

/* Finalize: */
int CouplingSolver::finalize() {
   
   /* Destroy the PETSc vectors used to evaluate convergence: */
   for (int i = 0; i < coupled_solvers.size(); i++) {
      Array1D<Field>& output_fields = coupled_solvers(i)->getFields();
      for (int f = 0; f < output_fields.size(); f++) {
         if (output_fields(f).vec0 != nullptr) {
            PETSC_CALL(VecDestroy(output_fields(f).vec0));
            utils::free(&(output_fields(f).vec0));
         }
      }
   }
   
   /* Finalize all the coupled solvers: */
   for (int i = 0; i < coupled_solvers.size(); i++) {
      PAMPA_CALL(coupled_solvers(i)->finalize(), "unable to finalize the solver");
   }
   
   return 0;
   
}
