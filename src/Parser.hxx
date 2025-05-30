#pragma once

#include "Mesh.hxx"
#include "CartesianMesh.hxx"
#include "UnstructuredExtrudedMesh.hxx"
#include "PartitionedMesh.hxx"
#include "Material.hxx"
#include "Solver.hxx"
#include "DiffusionSolver.hxx"
#include "SNSolver.hxx"
#include "PrecursorSolver.hxx"
#include "HeatConductionSolver.hxx"
#include "CouplingSolver.hxx"
#include "input.hxx"
#include "utils.hxx"

/* The Parser class: */
class Parser {
   
   public:
      
      /* The Parser constructor: */
      Parser() {}
      
      /* The Parser destructor: */
      ~Parser() {}
      
      /* Read a plain-text input file: */
      int WARN_UNUSED read(const std::string& filename, Mesh** mesh, Mesh** mesh_nodal, 
         Array1D<Material*>& materials, Array1D<Solver*>& solvers, Array1D<double>& dt);
   
};
