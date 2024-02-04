#pragma once

#include "utils.hxx"

/* The Material struct: */
struct Material {
   
   /* Transport method: */
   TransportMethod method;
   
   /* Cross sections: */
   Array1D<double> sigma_total;
   Array1D<double> nu_sigma_fission;
   Array2D<double> sigma_scattering;
   Array2D<double> sigma_scattering_s1;
   
   /* Diffusion coefficient: */
   Array1D<double> diffusion_coefficient;
   
   /* Fission spectrum: */
   Array1D<double> chi;
   
};
