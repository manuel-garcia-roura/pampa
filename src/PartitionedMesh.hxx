#pragma once

#include "Mesh.hxx"

/* The PartitionedMesh class: */
class PartitionedMesh : public Mesh {
   
   public:
      
      /* The PartitionedMesh constructor: */
      PartitionedMesh() {partitioned = true;}
      
      /* The PartitionedMesh destructor: */
      ~PartitionedMesh() {}
      
      /* Read the mesh from a plain-text input file: */
      int WARN_UNUSED read(const std::string& filename);
      
      /* Build the mesh (nothing to be done): */
      int WARN_UNUSED build() {return 0;}
   
};
