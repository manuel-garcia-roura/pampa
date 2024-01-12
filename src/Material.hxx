#pragma once

#include <vector>
#include <string>

#include "utils.hxx"

/* The Material struct: */
struct Material {
   
   /* Transport method: */
   TransportMethod method;
   
   /* Cross sections: */
   std::vector<double> sigma_total;
   std::vector<double> nu_sigma_fission;
   std::vector<std::vector<double>> sigma_scattering;
   std::vector<std::vector<double>> sigma_scattering_s1;
   std::vector<double> sigma_removal;
   
   /* Diffusion coefficient: */
   std::vector<double> diffusion_coefficient;
   
   /* Fission spectrum: */
   std::vector<double> chi;
   
};
