#pragma once

#include "petsc.hxx"
#include "output.hxx"
#include "utils.hxx"

/* The ConvergenceError class: */
class ConvergenceError {
   
   private:
      
      /* Vector name: */
      std::string name;
      
      /* Convergence norm type: */
      NormType p;
      
      /* Switch to use the relative instead of the absolute error: */
      bool relative;
      
      /* Convergence tolerance: */
      double tol;
      
      /* Previous iteration and vector difference used to evaluate the convergence: */
      Vec vec_prev = 0, vec_diff = 0;
   
   public:
      
      /* The ConvergenceError constructor: */
      ConvergenceError(const std::string& name, NormType p, bool relative, double tol) : 
         name(name), p(p), relative(relative), tol(tol) {}
      
      /* The ConvergenceError destructor: */
      ~ConvergenceError() {}
      
      /* Finalize: */
      int WARN_UNUSED finalize() {
         
         /* Destroy the PETSc vectors: */
         PAMPA_CHECK(petsc::destroy(vec_prev), "unable to destroy the PETSc vector");
         PAMPA_CHECK(petsc::destroy(vec_diff), "unable to destroy the PETSc vector");
         
         return 0;
         
      }
      
      /* Check if the vector has converged: */
      int WARN_UNUSED check(const Vec& vec, bool& converged) {
         
         /* Get the convergence error or initialize the convergence vectors: */
         if (vec_prev == 0) {
            
            /* Create the vectors used to evaluate the convergence: */
            PETSC_CALL(VecDuplicate(vec, &vec_prev));
            PETSC_CALL(VecDuplicate(vec, &vec_diff));
            
            /* Don't converge on the first iteration: */
            converged = false;
            
            /* Print info: */
            output::print("Convergence initialized for " + name + " vector.", true);
            
         }
         else {
            
            /* Get the vector difference: */
            PETSC_CALL(VecWAXPY(vec_diff, -1.0, vec_prev, vec));
            
            /* Get the norm of the difference: */
            double eps;
            PETSC_CALL(VecNorm(vec_diff, p, &eps));
            
            /* Normalize the difference with the L2-norm of the current iteration: */
            if (relative) {
               PetscScalar norm;
               PETSC_CALL(VecNorm(vec, NORM_2, &norm));
               eps /= norm;
            }
            
            /* Check if the convergence error is lower than the tolerance: */
            converged = eps < tol;
            
            /* Print info: */
            std::stringstream message;
            message << std::scientific;
            message << std::setprecision(3);
            message << "Convergence evaluation for " + name + " vector: ";
            if (eps > tol)
               message << eps << " > " << tol << " (not converged).";
            else
               message << eps << " < " << tol << " (converged).";
            output::print(message.str(), true);
            
         }
         
         /* Keep the vector for the next convergence evaluation: */
         PETSC_CALL(VecCopy(vec, vec_prev));
         
         return 0;
         
      }
   
};
