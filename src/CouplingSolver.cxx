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
         PAMPA_CHECK(utils::read(num_coupled_solvers, 1, INT_MAX, line[++l]), 
            "wrong number of coupled solvers");
         
         /* Get the coupled solvers: */
         coupled_solvers.resize(num_coupled_solvers, nullptr);
         for (int i = 0; i < num_coupled_solvers; i++) {
            PAMPA_CHECK(utils::find(line[++l], solvers, &(coupled_solvers(i))), 
               "unable to find the coupled solver");
         }
         
      }
      else if (line[l] == "implicit") {
         
         /* Get the switch to use implicit coupling: */
         PAMPA_CHECK(utils::read(implicit, line[++l]), "wrong switch for implicit coupling");
         
      }
      else if (line[l] == "convergence") {
         
         /* Get the convergence tolerance and p-norm for nonlinear problems: */
         PAMPA_CHECK(utils::read(tol, 0.0, DBL_MAX, line[++l]), "wrong convergence tolerance");
         PAMPA_CHECK(utils::read(p, 0.0, DBL_MAX, line[++l]), "wrong convergence p-norm");
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   return 0;
   
}

/* Initialize: */
int CouplingSolver::initialize(bool transient) {
   
   /* Initialize all the coupled solvers: */
   for (int i = 0; i < coupled_solvers.size(); i++) {
      PAMPA_CHECK(coupled_solvers(i)->initialize(transient), "unable to initialize the solver");
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
   output::print("Run '" + name + "' solver...", true);
   
   /* Iterate the solution until convegence: */
   bool converged = false;
   while (!converged) {
      
      /* Reset the convergence flag: */
      converged = true;
      
      /* Get the solution from all the coupled solvers: */
      for (int i = 0; i < coupled_solvers.size(); i++) {
         
         /* Get the solution: */
         PAMPA_CHECK(coupled_solvers(i)->solve(n, dt, t), "unable to get the solution");
         
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
                              output::print("Field exchange: ", true);
                              output::print("   - field name: " + output_fields(f).name, true);
                              output::print("   - output solver: " + coupled_solvers(i)->name, true);
                              output::print("   - input solver: " + coupled_solvers(i2)->name, true);
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
                     output::print("Field convergence initialized: ", true);
                     output::print("   - field name: " + output_fields(f).name, true);
                  }
                  else {
                     double eps;
                     PAMPA_CHECK(petsc::difference(*(output_fields(f).vec), 
                        *(output_fields(f).vec0), p, eps, true), 
                        "unable to calculate the convergence error");
                     converged &= eps < tol;
                     output::print("Field convergence: ", true);
                     output::print("   - field name: " + output_fields(f).name, true);
                     output::print("   - error: " + std::to_string(eps), true);
                     output::print("   - tolerance: " + std::to_string(tol), true);
                  }
                  PETSC_CALL(VecCopy(*(output_fields(f).vec), *(output_fields(f).vec0)));
               }
               
            }
         }
         
      }
      
      /* Print progress: */
      output::print("Coupled solution convergence: ", true);
      output::print("   - converged: " + std::to_string(converged), true);
      
   }
   
   /* Print progress: */
   output::print("Done.", true);
   
   return 0;
   
}

/* Output the solution: */
int CouplingSolver::output(const std::string& filename, int n, bool write_mesh) const {
   
   /* Write the mesh in .vtk format: */
   if (write_mesh) {
      PAMPA_CHECK(mesh->writeVTK(filename), "unable to write the mesh in .vtk format");
   }
   
   /* Output the solution from all the coupled solvers: */
   for (int i = 0; i < coupled_solvers.size(); i++) {
      PAMPA_CHECK(coupled_solvers(i)->output(filename, n, false), "unable to output the solution");
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
      PAMPA_CHECK(coupled_solvers(i)->finalize(), "unable to finalize the solver");
   }
   
   return 0;
   
}
