#pragma once

#include "Pampa.hxx"

/* Calculation driver: */
Pampa pampa;

/* Initialize: */
int WARN_UNUSED initialize(int argc, char* argv[], Array1D<double>& dt);

/* Solve: */
int WARN_UNUSED solve(int n = 0, double dt = 0.0, double t = 0.0);

/* Finalize: */
int WARN_UNUSED finalize();
