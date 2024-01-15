#include "AngularQuadratureSet.hxx"

/* Build the angular quadrature set: */
int AngularQuadratureSet::build() {
   
   /* Get the number of discrete directions: */
   num_directions = order * (order+2);
   directions.resize(num_directions);
   weights.resize(num_directions);
   
   /* Get the discrete directions and weights in the first octant: */
   switch (order) {
      case 2 : {
         
         /* Get the discrete direction: */
         directions[0] = {1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0)};
         
         /* Get the weight: */
         weights[0] = 1.0;
         
         break;
         
      }
      case 4 : {
         
         /* Get the direction cosines: */
         std::vector<double> mu = std::vector<double>{0.3500212, 0.8688903};
         
         /* Get the discrete directions: */
         directions[0] = {mu[0], mu[0], mu[1]};
         directions[1] = {mu[0], mu[1], mu[0]};
         directions[2] = {mu[1], mu[0], mu[0]};
         
         /* Get the weights: */
         weights[0] = 1.0 / 3.0;
         weights[1] = 1.0 / 3.0;
         weights[2] = 1.0 / 3.0;
         
         break;
         
      }
      case 6 : {
         
         /* Get the direction cosines and weights: */
         std::vector<double> mu = std::vector<double>{0.2666355, 0.6815076, 0.9261808};
         std::vector<double> w = std::vector<double>{0.1761263, 0.1572071};
         
         /* Get the discrete directions: */
         directions[0] = {mu[0], mu[0], mu[2]};
         directions[1] = {mu[0], mu[2], mu[0]};
         directions[2] = {mu[2], mu[0], mu[0]};
         directions[3] = {mu[0], mu[1], mu[1]};
         directions[4] = {mu[1], mu[0], mu[1]};
         directions[5] = {mu[1], mu[1], mu[0]};
         
         /* Get the weights: */
         weights[0] = w[0];
         weights[1] = w[0];
         weights[2] = w[0];
         weights[3] = w[1];
         weights[4] = w[1];
         weights[5] = w[1];
         
         break;
         
      }
      case 8 : {
         
         /* Get the direction cosines and weights: */
         std::vector<double> mu = std::vector<double>{0.2182179, 0.5773503, 0.7867958, 0.9511897};
         std::vector<double> w = std::vector<double>{0.1209877, 0.0907407, 0.0925926};
         
         /* Get the discrete directions: */
         directions[0] = {mu[0], mu[0], mu[3]};
         directions[1] = {mu[0], mu[3], mu[0]};
         directions[2] = {mu[3], mu[0], mu[0]};
         directions[3] = {mu[0], mu[1], mu[2]};
         directions[4] = {mu[0], mu[2], mu[1]};
         directions[5] = {mu[1], mu[0], mu[2]};
         directions[6] = {mu[1], mu[2], mu[0]};
         directions[7] = {mu[2], mu[0], mu[1]};
         directions[8] = {mu[2], mu[1], mu[0]};
         directions[9] = {mu[1], mu[1], mu[1]};
         
         /* Get the weights: */
         weights[0] = w[0];
         weights[1] = w[0];
         weights[2] = w[0];
         weights[3] = w[1];
         weights[4] = w[1];
         weights[5] = w[1];
         weights[6] = w[1];
         weights[7] = w[1];
         weights[8] = w[1];
         weights[9] = w[2];
         
         break;
         
      }
      default : {
         PAMPA_CHECK(true, 1, "SN order not implemented");
      }
   }
   
   /* Get the discrete directions and weights in the other octants: */
   for (int i = 1; i < 8; i++) {
      for (int m = 0; m < num_directions/8; m++) {
         int l = i*num_directions/8 + m;
         directions[l].resize(num_directions);
         directions[l][0] = (i & 1) ? -directions[m][0] : directions[m][0];
         directions[l][1] = (i & 2) ? -directions[m][1] : directions[m][1];
         directions[l][2] = (i & 4) ? -directions[m][2] : directions[m][2];
      }
   }
   
   return 0;
   
}
