#pragma once

#include <vector>
#include <string>

#include "Mesh.hxx"
#include "AngularQuadratureSet.hxx"
#include "Material.hxx"

/* The Model class: */
class Model {
   
   private:
      
      /* Transport method: */
      TransportMethod method;
      
      /* Model mesh: */
      Mesh* mesh = NULL;
      
      /* Angular quadrature set: */
      AngularQuadratureSet quadrature;
      
      /* Model materials: */
      std::vector<Material> materials;
   
   public:
      
      /* The Model constructor: */
      Model() {}
      
      /* The Model destructor: */
      ~Model() {delete mesh;}
      
      /* Set the transport method: */
      void setTransportMethod(const TransportMethod& method) {this->method = method;}
      
      /* Set the model mesh: */
      void setMesh(Mesh* mesh) {this->mesh = mesh;}
      
      /* Add a model material: */
      void addMaterial(const Material& material) {materials.push_back(material);}
      
      /* Get the transport method: */
      const TransportMethod& getTransportMethod() const {return method;}
      
      /* Get the model mesh: */
      const Mesh* getMesh() const {return mesh;}
      
      /* Get the angular quadrature set: */
      const AngularQuadratureSet& getAngularQuadratureSet() const {return quadrature;}
      
      /* Get the model materials: */
      const std::vector<Material>& getMaterials() const {return materials;}
      
      /* Build the model: */
      int build();
      
   
};
