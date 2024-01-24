#include "Model.hxx"

/* Build the model: */
int Model::build() {
   
   /* Build the mesh: */
   PAMPA_CALL(mesh->build(), "unable to build the mesh");
   
   /* Build the angular quadrature set: */
   if (method.type == TM::SN) {
      quadrature.setOrder(method.order);
      PAMPA_CALL(quadrature.build(), "unable to build the angular quadrature set");
   }
   
   /* Check the materials: */
   for (int i = 0; i < materials.size(); i++) {
      PAMPA_CHECK(materials[i].method.type != method.type, 1, "wrong transport method");
      PAMPA_CHECK(materials[i].method.num_groups != method.num_groups, 2, 
         "wrong number of energy groups");
   }
   
   return 0;
   
}
