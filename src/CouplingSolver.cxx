#include "CouplingSolver.hxx"

/* Initialize: */
int CouplingSolver::initialize() {
   
   /* Initialize all the solvers: */
   for (int i = 0; i < solvers.size(); i++) {
      PAMPA_CALL(solvers(i)->initialize(), "unable to initialize the solver");
   }
   
   /* Get the feedback fields: */
   for (int i = 0; i < solvers.size(); i++) {
      const Array1D<Field>& coupled_fields = solvers(i)->getFields();
      for (int f = 0; f < coupled_fields.size(); f++)
         fields.pushBack(coupled_fields(f));
   }
   
   return 0;
   
}

/* Get the solution: */
int CouplingSolver::solve(int n, double dt) {
   
   /* Get the solution from all the solvers: */
   for (int i = 0; i < solvers.size(); i++) {
      
      /* Get the solution: */
      PAMPA_CALL(solvers(i)->solve(n, dt), "unable to get the solution from the solver");
      
      /* Exchange the output fields calculated by this solver: */
      const Array1D<Field>& output_fields = solvers(i)->getFields();
      for (int f = 0; f < output_fields.size(); f++) {
         if (output_fields(f).output) {
            
            /* Set the field in the solvers that take it as input: */
            for (int i2 = 0; i2 < solvers.size(); i2++) {
               if (i2 != i) {
                  const Array1D<Field>& input_fields = solvers(i2)->getFields();
                  for (int f2 = 0; f2 < input_fields.size(); f2++) {
                     if (input_fields(f2).input) {
                        if (input_fields(f2).name == output_fields(f).name) {
                           PETSC_CALL(VecCopy(*(output_fields(f).vector), 
                              *(input_fields(f2).vector)));
                        }
                     }
                  }
               }
            }
            
         }
      }
      
   }
   
   return 0;
   
}

/* Output the solution: */
int CouplingSolver::output(const std::string& filename, int n, bool write_mesh) const {
   
   /* Write the mesh in .vtk format: */
   if (write_mesh) {
      PAMPA_CALL(mesh->writeVTK(filename), "unable to write the mesh in .vtk format");
   }
   
   /* Output the solution from all the solvers: */
   for (int i = 0; i < solvers.size(); i++) {
      PAMPA_CALL(solvers(i)->output(filename, n, false), 
         "unable to output the solution from the solver");
   }
   
   return 0;
   
}

/* Finalize: */
int CouplingSolver::finalize() {
   
   /* Finalize all the solvers: */
   for (int i = 0; i < solvers.size(); i++) {
      PAMPA_CALL(solvers(i)->finalize(), "unable to finalize the solver");
   }
   
   return 0;
   
}
