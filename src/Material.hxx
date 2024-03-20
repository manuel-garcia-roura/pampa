#pragma once

#include "utils.hxx"

/* The Material struct: */
struct Material {
   
   /* Number of energy groups: */
   int num_groups = -1;
   
   /* Cross sections: */
   Array1D<double> sigma_total;
   Array1D<double> nu_sigma_fission;
   Array2D<double> sigma_scattering;
   
   /* Diffusion coefficients: */
   Array1D<double> diffusion_coefficient;
   
   /* Fission spectrum: */
   Array1D<double> chi;
   
   /* Thermal properties: */
   double rho = 0.0, cp = 0.0, k = 0.0;
   
};
