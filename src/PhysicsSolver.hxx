#pragma once

#include "Solver.hxx"
#include "Material.hxx"

/* The PhysicsSolver class: */
class PhysicsSolver : public Solver {
   
   protected:
      
      /* Materials: */
      const Array1D<Material*>& materials;
      
      /* Pointers to the PETSc matrices: */
      Array1D<Mat*> matrices;
      
      /* Pointers to the PETSc vectors: */
      Array1D<Vec*> vectors;
      
      /* Krylov Subspace Solver (KSP) context: */
      KSP ksp = 0;
      
      /* Eigenvalue Problem Solver (EPS) context: */
      EPS eps = 0;
      
      /* Check the material data: */
      virtual int WARN_UNUSED checkMaterials(bool transient = false) {PAMPA_CHECK_VIRTUAL}
      
      /* Build the matrices and vectors: */
      virtual int WARN_UNUSED build() {PAMPA_CHECK_VIRTUAL}
      
      /* Print the solution summary to standard output: */
      virtual int WARN_UNUSED printLog(int n = 0) const {PAMPA_CHECK_VIRTUAL}
      
      /* Write the solution to a plain-text file in .vtk format: */
      virtual int WARN_UNUSED writeVTK(const std::string& filename) const {PAMPA_CHECK_VIRTUAL}
      
      /* Write the solution to a binary file in PETSc format: */
      virtual int WARN_UNUSED writePETSc(int n = 0) const {PAMPA_CHECK_VIRTUAL}
   
   public:
      
      /* The PhysicsSolver constructor: */
      PhysicsSolver(const std::string& name, const Mesh* mesh, 
         const Array1D<Material*>& materials) : Solver(name, mesh), materials(materials) {}
      
      /* The PhysicsSolver destructor: */
      virtual ~PhysicsSolver() {}
      
      /* Initialize: */
      int WARN_UNUSED initialize(bool transient = false);
      
      /* Output the solution: */
      int WARN_UNUSED output(const std::string& filename, int n = 0, bool write_mesh = true) const;
      
      /* Finalize: */
      int WARN_UNUSED finalize();
   
};
