#include "CouplingSolver.hxx"

/* Initialize: */
int CouplingSolver::initialize(bool transient) {
   
   /* Initialize all the solvers: */
   for (int i = 0; i < solvers.size(); i++) {
      PAMPA_CALL(solvers(i)->initialize(transient), "unable to initialize the solver");
   }
   
   /* Get the feedback fields: */
   for (int i = 0; i < solvers.size(); i++) {
      Array1D<Field>& coupled_fields = solvers(i)->getFields();
      for (int f = 0; f < coupled_fields.size(); f++)
         fields.pushBack(coupled_fields(f));
   }
   
   return 0;
   
}

/* Get the solution: */
int CouplingSolver::solve(int n, double dt) {
   
   /* Iterate the solution until convegence: */
   bool converged = false;
   while (!converged) {
      
      /* Reset the convergence flag: */
      converged = true;
      
      /* Get the solution from all the solvers: */
      for (int i = 0; i < solvers.size(); i++) {
         
         /* Get the solution: */
         PAMPA_CALL(solvers(i)->solve(n, dt), "unable to get the solution from the solver");
         
         /* Exchange the output fields calculated by this solver: */
         Array1D<Field>& output_fields = solvers(i)->getFields();
         for (int f = 0; f < output_fields.size(); f++) {
            if (output_fields(f).output) {
               
               /* Set the field in the solvers that take it as input: */
               for (int i2 = 0; i2 < solvers.size(); i2++) {
                  if (i2 != i) {
                     Array1D<Field>& input_fields = solvers(i2)->getFields();
                     for (int f2 = 0; f2 < input_fields.size(); f2++) {
                        if (input_fields(f2).input) {
                           if (input_fields(f2).name == output_fields(f).name) {
                              PETSC_CALL(VecCopy(*(output_fields(f).vec), *(input_fields(f2).vec)));
                           }
                        }
                     }
                  }
               }
               
               /* Evaluate the convergence: */
               if (implicit) {
                  if (output_fields(f).vec0 == NULL) {
                     output_fields(f).vec0 = new Vec;
                     PETSC_CALL(VecDuplicate(*(output_fields(f).vec), output_fields(f).vec0));
                  }
                  else {
                     double eps;
                     PAMPA_CALL(petsc::difference(*(output_fields(f).vec), 
                        *(output_fields(f).vec0), p, eps), 
                        "unable to calculate the convergence error");
                     converged &= eps < tol;
                  }
                  PETSC_CALL(VecCopy(*(output_fields(f).vec), *(output_fields(f).vec0)));
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
      PAMPA_CALL(solvers(i)->output(filename, n, false), "unable to output the solution");
   }
   
   return 0;
   
}

/* Finalize: */
int CouplingSolver::finalize() {
   
   /* Destroy the PETSc vectors used to evaluate convergence: */
   for (int i = 0; i < solvers.size(); i++) {
      Array1D<Field>& output_fields = solvers(i)->getFields();
      for (int f = 0; f < output_fields.size(); f++) {
         if (output_fields(f).vec0 != NULL) {
            PETSC_CALL(VecDestroy(output_fields(f).vec0));
            delete output_fields(f).vec0;
         }
      }
   }
   
   /* Finalize all the solvers: */
   for (int i = 0; i < solvers.size(); i++) {
      PAMPA_CALL(solvers(i)->finalize(), "unable to finalize the solver");
   }
   
   return 0;
   
}
