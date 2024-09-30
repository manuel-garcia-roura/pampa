#include "CouplingSolver.hxx"

/* Read the solver from a plain-text input file: */
int CouplingSolver::read(std::ifstream& file, Array1D<Solver*>& solvers) {
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = input::get_next_line(file);
      if (line.empty() || line[0] == "}") break;
      
      /* Get the next keyword: */
      unsigned int l = 0;
      if (line[l] == "coupled-solvers") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() < 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the number of coupled solvers: */
         int num_coupled_solvers;
         PAMPA_CHECK(input::read(num_coupled_solvers, 1, INT_MAX, line[++l]), 
            "wrong number of coupled solvers");
         
         /* Get the coupled solvers: */
         PAMPA_CHECK(line.size() != unsigned(2+num_coupled_solvers), 
            "wrong number of arguments for keyword '" + line[l] + "'");
         coupled_solvers.resize(num_coupled_solvers, nullptr);
         for (int i = 0; i < num_coupled_solvers; i++) {
            PAMPA_CHECK(utils::find(line[++l], solvers, &(coupled_solvers(i))), 
               "unable to find the coupled solver");
         }
         
      }
      else if (line[l] == "implicit") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the switch to use implicit coupling: */
         PAMPA_CHECK(input::read(implicit, line[++l]), "wrong switch for implicit coupling");
         
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
   
   /* Print info: */
   output::print("Run " + name + " solver...", true);
   output::indent(true);
   
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
                              output::print("Feedback performed for " + output_fields(f).name + 
                                 " vector (" + coupled_solvers(i)->name + " -> " + 
                                 coupled_solvers(i2)->name + ").", true);
                           }
                        }
                     }
                  }
               }
               
               /* Evaluate the convergence: */
               if (implicit) {
                  bool conv;
                  PAMPA_CHECK((output_fields(f).delta)->check(*(output_fields(f).vec), conv), 
                     "unable to check the convergence");
                  converged &= conv;
               }
               
            }
         }
         
      }
      
      /* Print info: */
      if (!converged)
         output::print("Coupled solution not converged.", true);
      else
         output::print("Coupled solution converged.", true);
      
   }
   
   /* Print info: */
   output::outdent(true);
   output::print("Done.", true);
   
   return 0;
   
}

/* Output the solution: */
int CouplingSolver::output(const std::string& path, int n, bool write_mesh) const {
   
   /* Write the mesh in .vtk format: */
   if (write_mesh) {
      PAMPA_CHECK(mesh->writeVTK(path + "/output", n), "unable to write the mesh in .vtk format");
   }
   
   /* Output the solution from all the coupled solvers: */
   for (int i = 0; i < coupled_solvers.size(); i++) {
      PAMPA_CHECK(coupled_solvers(i)->output(path, n, false), "unable to output the solution");
   }
   
   return 0;
   
}

/* Finalize: */
int CouplingSolver::finalize() {
   
   /* Finalize all the coupled solvers: */
   for (int i = 0; i < coupled_solvers.size(); i++) {
      PAMPA_CHECK(coupled_solvers(i)->finalize(), "unable to finalize the solver");
   }
   
   return 0;
   
}
