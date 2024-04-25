#pragma once

#include "Solver.hxx"
#include "Material.hxx"

/* The PhysicsSolver class: */
class PhysicsSolver : public Solver {
   
   protected:
      
      /* Mesh dimensions: */
      const int num_cells = -1, num_cells_global = -1, num_faces_max = -1;
      
      /* Mesh cells: */
      const Cells& cells;
      
      /* Mesh faces: */
      const Faces& faces;
      
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
      virtual int WARN_UNUSED writePETSc() const {PAMPA_CHECK_VIRTUAL}
   
   public:
      
      /* The PhysicsSolver constructor: */
      PhysicsSolver(const std::string& name, const Mesh* mesh, 
         const Array1D<Material*>& materials) : Solver(name, mesh), num_cells(mesh->getNumCells()), 
         num_cells_global(mesh->getNumCellsGlobal()), num_faces_max(mesh->getNumFacesMax()), 
         cells(mesh->getCells()), faces(mesh->getFaces()), materials(materials) {}
      
      /* The PhysicsSolver destructor: */
      virtual ~PhysicsSolver() {}
      
      /* Initialize: */
      int WARN_UNUSED initialize(bool transient = false);
      
      /* Output the solution: */
      int WARN_UNUSED output(const std::string& filename, int n = 0, bool write_mesh = true) const;
      
      /* Finalize: */
      int WARN_UNUSED finalize();
   
};
