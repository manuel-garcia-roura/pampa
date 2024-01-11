#pragma once

#include <vector>
#include <string>

/* The Material struct: */
struct Material {
   
   /* Cross sections: */
   std::vector<double> sigma_total;
   std::vector<double> nu_sigma_fission;
   std::vector<std::vector<double>> sigma_scattering;
   std::vector<double> sigma_removal;
   
   /* Diffusion coefficient: */
   std::vector<double> diffusion_coefficient;
   
   /* Fission spectrum: */
   std::vector<double> chi;
   
};
