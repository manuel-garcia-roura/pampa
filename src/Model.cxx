#include "Model.hxx"

/* The Model constructor: */
Model::Model() : name("default"), num_groups(-1), mesh(NULL) {};

/* The Model destructor: */
Model::~Model() {delete mesh;};

/* Build the model: */
int Model::build() {
   
   /* Build the mesh: */
   PAMPA_CALL(mesh->build(), "unable to build the mesh");
   
   /* Set the model materials: */
   mesh->setModelMaterials(&materials);
   
   return 0;
   
};

/* Output the model: */
int Model::output(const std::string &filename) {
   
   /* Write the mesh: */
   PAMPA_CALL(mesh->write(filename), "unable to write the mesh");
   
   return 0;
   
};

/* Set the number of energy groups: */
void Model::setNumEnergyGroups(int num_groups) {
   
   this->num_groups = num_groups;
   
};

/* Get the number of energy groups: */
int Model::getNumEnergyGroups() const {
   
   return num_groups;
   
};

/* Set the mesh: */
void Model::setMesh(Mesh *mesh) {
   
   this->mesh = mesh;
   
};

/* Get the mesh: */
const Mesh* Model::getMesh() const {
   
   return mesh;
   
};

/* Add a material: */
void Model::addMaterial(const Material &material) {
   
   materials.push_back(material);
   
};
