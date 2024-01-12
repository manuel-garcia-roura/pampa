#include "Model.hxx"

/* Build the model: */
int Model::build() {
   
   /* Build the mesh: */
   PAMPA_CALL(mesh->build(), "unable to build the mesh");
   
   /* Build the materials: */
   for (int i = 0; i < materials.size(); i++) {
      PAMPA_CALL(buildMaterial(materials[i]), "unable to build the material");
   }
   
   return 0;
   
}

/* Build a material: */
int Model::buildMaterial(Material &material) {
   
   /* Check the transport method: */
   PAMPA_CHECK(material.method.type != method.type, 1, "wrong transport method");
   PAMPA_CHECK(material.method.num_groups != method.num_groups, 2, "wrong number of energy groups");
   
   /* Calculate the removal cross sections: */
   material.sigma_removal.resize(method.num_groups);
   for (int g = 0; g < method.num_groups; g++)
      material.sigma_removal[g] = material.sigma_total[g] - material.sigma_scattering[g][g];
   
   return 0;
   
}
