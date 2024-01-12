#pragma once

#include <vector>
#include <string>

#include "Mesh.hxx"
#include "Material.hxx"

/* The Model class: */
class Model {
   
   private:
      
      /* Number of energy groups: */
      int num_groups = -1;
      
      /* Model mesh: */
      Mesh *mesh = NULL;
      
      /* Model materials: */
      std::vector<Material> materials;
      
      /* Build a material: */
      void buildMaterial(Material &material);
   
   public:
      
      /* The Model constructor: */
      Model() {}
      
      /* The Model destructor: */
      ~Model() {delete mesh;}
      
      /* Set the number of energy groups: */
      void setNumEnergyGroups(int num_groups) {this->num_groups = num_groups;}
      
      /* Set the model mesh: */
      void setMesh(Mesh *mesh) {this->mesh = mesh;}
      
      /* Add a model material: */
      void addMaterial(const Material &material) {materials.push_back(material);}
      
      /* Get the number of energy groups: */
      int getNumEnergyGroups() const {return num_groups;}
      
      /* Get the model mesh: */
      const Mesh* getMesh() const {return mesh;}
      
      /* Get the model materials: */
      const std::vector<Material>& getMaterials() const {return materials;}
      
      /* Build the model: */
      int build();
      
   
};
