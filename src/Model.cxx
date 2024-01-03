#include "Model.hxx"

/* The Model constructor: */
Model::Model() : name("default"), num_groups(-1), mesh(NULL) {};

/* The Model destructor: */
Model::~Model() {delete mesh;};

/* Set the number of energy groups: */
void Model::setNumEnergyGroups(int num_groups) {
   
   this->num_groups = num_groups;
   
};

/* Set the mesh: */
void Model::setMesh(Mesh *mesh) {
   
   this->mesh = mesh;
   
};

/* Add a material: */
void Model::addMaterial(const Material &material) {
   
   materials.push_back(material);
   
};

/* Get the number of energy groups: */
int Model::getNumEnergyGroups() const {
   
   return num_groups;
   
};

/* Get the mesh: */
const Mesh* Model::getMesh() const {
   
   return mesh;
   
};

/* Build the model: */
int Model::build() {
   
   /* Build the mesh: */
   PAMPA_CALL(mesh->build(), "unable to build the mesh");
   
   /* Build the materials: */
   for (int i = 0; i < materials.size(); i++)
      buildMaterial(materials[i]);
   
   /* Set the model materials to the mesh: */
   mesh->setModelMaterials(&materials);
   
   return 0;
   
};

/* Output the model: */
int Model::output(const std::string &filename) {
   
   /* Write the mesh: */
   PAMPA_CALL(mesh->write(filename), "unable to write the mesh");
   
   return 0;
   
};

/* Build a material: */
void Model::buildMaterial(Material &material) {
   
   /* Calculate the removal cross sections: */
   material.sigma_removal.resize(num_groups);
   for (int g = 0; g < num_groups; g++)
      material.sigma_removal[g] = material.sigma_total[g] - material.sigma_scattering[g][g];
   
};
