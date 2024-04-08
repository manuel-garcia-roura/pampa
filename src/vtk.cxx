#include "vtk.hxx"

/* Write a solution vector to a .vtk file: */
int vtk::write(const std::string& filename, const std::string& name, const Vec& v, int num_cells, 
   int num_groups, int num_directions) {
   
   /* Open the output file: */
   std::ofstream file(filename, std::ios_base::app);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Get the array with the raw data: */
   PetscScalar* data;
   PETSC_CALL(VecGetArray(v, &data));
   
   /* Write the solution vector: */
   for (int g = 0; g < num_groups; g++) {
      for (int m = 0; m < num_directions; m++) {
         
         /* Write the section header: */
         file << "SCALARS " << name;
         if (num_groups > 1) file << "_" << (g+1);
         if (num_directions > 1) file << "_" << (m+1);
         file << " double 1" << std::endl;
         file << "LOOKUP_TABLE default" << std::endl;
         
         /* Write the vector data: */
         for (int iv = g*num_directions+m, i = 0; i < num_cells; i++) {
            file << data[iv] << std::endl;
            iv += num_directions*num_groups;
         }
         file << std::endl;
         
      }
   }
   
   /* Restore the array with the raw data: */
   PETSC_CALL(VecRestoreArray(v, &data));
   
   return 0;
   
}
