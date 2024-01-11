#include "Model.hxx"

/* Build the model: */
int Model::build() {
   
   /* Build the mesh: */
   PAMPA_CALL(mesh->build(), "unable to build the mesh");
   
   /* Build the materials: */
   for (int i = 0; i < materials.size(); i++)
      buildMaterial(materials[i]);
   
   return 0;
   
}

/* Output the model: */
int Model::output(const std::string &filename) const {
   
   /* Write the mesh: */
   PAMPA_CALL(mesh->write(filename), "unable to write the mesh");
   
   return 0;
   
}

/* Build a material: */
void Model::buildMaterial(Material &material) {
   
   /* Calculate the removal cross sections: */
   material.sigma_removal.resize(num_groups);
   for (int g = 0; g < num_groups; g++)
      material.sigma_removal[g] = material.sigma_total[g] - material.sigma_scattering[g][g];
   
}
