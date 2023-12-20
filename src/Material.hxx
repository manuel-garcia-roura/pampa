#pragma once

#include <vector>
#include <string>

/* The Material struct: */
struct Material {
   
   /* Material name: */
   std::string name = "default";
   
   /* Cross sections: */
   std::vector<double> sigma_total;
   std::vector<double> sigma_removal;
   std::vector<double> sigma_absorption;
   std::vector<double> nu_sigma_fission;
   std::vector<std::vector<double>> sigma_scattering;
   
   /* Diffusion coefficient: */
   std::vector<double> diffusion_coefficient;
   
   /* Fission spectrum: */
   std::vector<double> chi;
   
};
