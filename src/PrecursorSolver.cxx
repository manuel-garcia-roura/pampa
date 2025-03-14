#include "PrecursorSolver.hxx"

/* Read the solver from a plain-text input file: */
int PrecursorSolver::read(std::ifstream& file, Array1D<Solver*>& solvers) {
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = input::get_next_line(file);
      if (line.empty() || line[0] == "}") break;
      
      /* Get the next keyword: */
      unsigned int l = 0;
      if (line[l] == "precursor-groups") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 2, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the number of delayed-neutron precursor groups: */
         PAMPA_CHECK(input::read(num_precursor_groups, 1, INT_MAX, line[++l]), 
            "wrong number of delayed-neutron precursor groups");
         
      }
      else if (line[l] == "convergence") {
         
         /* Check the number of arguments: */
         PAMPA_CHECK(line.size() != 5, "wrong number of arguments for keyword '" + line[l] + "'");
         
         /* Get the field name: */
         std::string name = line[++l];
         
         /* Get the convergence error: */
         if (name == "precursors") {
            PAMPA_CHECK(input::read(dC, name, line, l), "wrong convergence error");
         }
         else if (name == "delayed-source") {
            PAMPA_CHECK(input::read(dS, name, line, l), "wrong convergence error");
         }
         else {
            PAMPA_CHECK(true, "wrong field");
         }
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, "unrecognized keyword '" + line[l] + "'");
         
      }
      
   }
   
   return 0;
   
}

/* Solve the linear system to get the solution: */
int PrecursorSolver::solve(int n, double dt, double t) {
   
   /* Print info: */
   output::print("Run " + name + " solver...", true);
   output::indent(true);
   
   /* Copy the precursor population from the previous time step: */
   if (n > n0+1) {
      if (C0 == 0) {
         PETSC_CALL(VecDuplicate(C, &C0));
      }
      PETSC_CALL(VecCopy(C, C0));
      n0++;
   }
   
   /* Get the arrays with the raw data: */
   PetscScalar *P_data, *C_data, *S_data, *C0_data;
   PETSC_CALL(VecGetArray(P, &P_data));
   PETSC_CALL(VecGetArray(C, &C_data));
   PETSC_CALL(VecGetArray(S, &S_data));
   if (n > 0) {
      PETSC_CALL(VecGetArray(C0, &C0_data));
   }
   
   /* Get the steady-state or transient precursor population: */
   if (n == 0) {
      for (int i = 0; i < num_cells; i++) {
         const Material* mat = materials(cells.materials(i));
         if (mat->hasPrecursorData())
            for (int g = 0; g < num_precursor_groups; g++)
               C_data[index(i, g)] = (mat->beta(g)/mat->lambda(g)) * P_data[i];
      }
   }
   else {
      for (int i = 0; i < num_cells; i++) {
         const Material* mat = materials(cells.materials(i));
         if (mat->hasPrecursorData())
            for (int g = 0; g < num_precursor_groups; g++)
               C_data[index(i, g)] = (C0_data[index(i, g)]+dt*mat->beta(g)*P_data[i]) / 
                                        (1.0+dt*mat->lambda(g));
      }
   }
   
   /* Get the delayed neutron source: */
   for (int i = 0; i < num_cells; i++) {
      const Material* mat = materials(cells.materials(i));
      S_data[i] = 0.0;
      if (mat->hasPrecursorData())
         for (int g = 0; g < num_precursor_groups; g++)
            S_data[i] += mat->lambda(g) * C_data[index(i, g)];
   }
   
   /* Restore the arrays with the raw data: */
   PETSC_CALL(VecRestoreArray(P, &P_data));
   PETSC_CALL(VecRestoreArray(C, &C_data));
   PETSC_CALL(VecRestoreArray(S, &S_data));
   if (n > 0) {
      PETSC_CALL(VecRestoreArray(C0, &C0_data));
   }
   
   /* Print info: */
   output::outdent(true);
   output::print("Done.", true);
   
   return 0;
   
}

/* Check the material data: */
int PrecursorSolver::checkMaterials(bool transient) {
   
   /* Check the materials: */
   for (int i = 0; i < materials.size(); i++) {
      const Material* mat = materials(i);
      if (mat->hasPrecursorData()) {
         PAMPA_CHECK(mat->checkPrecursorData(num_precursor_groups), "wrong precursor data");
      }
   }
   
   return 0;
   
}

/* Build the solution and source vectors: */
int PrecursorSolver::build() {
   
   /* Create the production-rate vector: */
   PAMPA_CHECK(petsc::create(P, num_cells, num_cells_global, vectors), 
      "unable to create the production-rate vector");
   fields.pushBack(Field{"production-rate", &P, true, false, nullptr});
   
   /* Create the precursor-population vector: */
   int size_local = num_cells * num_precursor_groups;
   int size_global = num_cells_global * num_precursor_groups;
   PAMPA_CHECK(petsc::create(C, size_local, size_global, vectors), 
      "unable to create the precursor-population vector");
   fields.pushBack(Field{"precursors", &C, false, true, &dC});
   
   /* Create the delayed-neutron-source vector: */
   PAMPA_CHECK(petsc::create(S, num_cells, num_cells_global, vectors), 
      "unable to create the delayed-neutron-source vector");
   fields.pushBack(Field{"delayed-source", &S, false, true, &dS});
   
   /* Initialize the production rate: */
   PAMPA_CHECK(petsc::set(P, 1.0), "unable to initialize the production rate");
   PAMPA_CHECK(petsc::normalize(P, 1.0), "unable to normalize the production rate");
   
   return 0;
   
}

/* Print the solution summary to standard output: */
int PrecursorSolver::printLog(int n) const {
   
   /* Print out the minimum and maximum precursor populations: */
   PetscScalar C_min, C_max;
   PETSC_CALL(VecMin(C, nullptr, &C_min));
   PETSC_CALL(VecMax(C, nullptr, &C_max));
   output::print("Precursor population", C_min, C_max, true, 3);
   
   return 0;
   
}

/* Write the solution to a plain-text file in .vtk format: */
int PrecursorSolver::writeVTK(const std::string& path, int n) const {
   
   /* Write the precursor population in .vtk format: */
   PAMPA_CHECK(vtk::write(path + "/output", n, "precursor", C, num_cells, num_precursor_groups), 
      "unable to write the precursor population");
   
   return 0;
   
}

/* Write the solution to a binary file in PETSc format: */
int PrecursorSolver::writePETSc(int n) const {
   
   /* Write the precursor population in PETSc format: */
   PAMPA_CHECK(petsc::write("precursors", n, C), "unable to write the precursor population");
   
   return 0;
   
}
