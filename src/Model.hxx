#pragma once

#include <vector>
#include <string>

#include "Mesh.hxx"
#include "Material.hxx"
#include "Solver.hxx"

/* The Model struct: */
struct Model {
   
   /* Model name: */
   std::string model = "default";
   
   /* Number of energy groups: */
   int num_groups = -1;
   
   /* Model mesh: */
   Mesh *mesh = NULL;
   
   /* Model materials: */
   std::vector<Material> materials;
   
   /* Solver: */
   Solver solver;
   
   /* The Model destructor: */
   ~Model() {delete mesh;};
   
};
