#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "Mesh.hxx"
#include "math.hxx"
#include "utils.hxx"

/* The PartitionedMesh class: */
class PartitionedMesh : public Mesh {
   
   public:
      
      /* The PartitionedMesh constructor: */
      PartitionedMesh() {}
      
      /* The PartitionedMesh destructor: */
      ~PartitionedMesh() {}
      
      /* Read the mesh from a plain-text input file: */
      int read(const std::string& filename);
      
      /* Build the mesh (nothing to be done): */
      int build() {return 0;}
   
};
