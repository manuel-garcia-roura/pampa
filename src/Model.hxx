#pragma once

#include <vector>
#include <string>

#include "Mesh.hxx"
#include "Material.hxx"

/* The Model class: */
class Model {
   
   private:
      
      /* Model name: */
      std::string name;
      
      /* Number of energy groups: */
      int num_groups;
      
      /* Model mesh: */
      Mesh *mesh;
      
      /* Model materials: */
      std::vector<Material> materials;
   
   public:
      
      /* The Model constructor: */
      Model();
      
      /* The Model destructor: */
      ~Model();
      
      /* Build the model: */
      int build();
      
      /* Output the model: */
      int output(const std::string &filename);
      
      /* Set the number of energy groups: */
      void setNumEnergyGroups(int num_groups);
      
      /* Get the number of energy groups: */
      int getNumEnergyGroups() const;
      
      /* Set the mesh: */
      void setMesh(Mesh *mesh);
      
      /* Get the mesh: */
      const Mesh* getMesh() const;
      
      /* Add a material: */
      void addMaterial(const Material &material);
   
};
