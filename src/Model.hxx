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
      
      /* Build a material: */
      void buildMaterial(Material &material);
   
   public:
      
      /* The Model constructor: */
      Model();
      
      /* The Model destructor: */
      ~Model();
      
      /* Set the number of energy groups: */
      void setNumEnergyGroups(int num_groups);
      
      /* Set the mesh: */
      void setMesh(Mesh *mesh);
      
      /* Add a material: */
      void addMaterial(const Material &material);
      
      /* Get the number of energy groups: */
      int getNumEnergyGroups() const;
      
      /* Get the mesh: */
      const Mesh* getMesh() const;
      
      /* Build the model: */
      int build();
      
      /* Output the model: */
      int output(const std::string &filename);
   
};
