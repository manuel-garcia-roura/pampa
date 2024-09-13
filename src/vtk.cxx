#include "vtk.hxx"

/* Write a mesh to a .vtk file: */
int vtk::write(const std::string& filename, const Array2D<double>& points, int num_points, 
   const Vector2D<int>& cells, int num_cells, const Array1D<int>& materials, 
   const Array1D<int>& nodal_indices) {
   
   /* Open the output file: */
   std::ofstream file(filename, std::ios_base::out);
   PAMPA_CHECK(!file.is_open(), "unable to open " + filename);
   
   /* Set the precision: */
   file << std::scientific;
   file << std::setprecision(VTK_PRECISION);
   
   /* Write the file header: */
   file << "# vtk DataFile Version 3.0" << std::endl;
   file << "FVM mesh" << std::endl;
   file << "ASCII" << std::endl;
   file << "DATASET UNSTRUCTURED_GRID" << std::endl;
   file << std::endl;
   
   /* Write the point coordinates: */
   file << "POINTS " << num_points << " double" << std::endl;
   for (int i = 0; i < num_points; i++) {
      file << points(i, 0) << " ";
      file << points(i, 1) << " ";
      file << points(i, 2) << std::endl;
   }
   file << std::endl;
   
   /* Get the total number of cell points: */
   int num_cell_points = num_cells;
   for (int i = 0; i < num_cells; i++)
      num_cell_points += cells.size(i);
   
   /* Write the cell points: */
   file << "CELLS " << num_cells << " " << num_cell_points << std::endl;
   for (int i = 0; i < num_cells; i++) {
      num_cell_points = cells.size(i);
      file << num_cell_points;
      for (int j = 0; j < num_cell_points; j++)
         file << " " << cells(i, j);
      file << std::endl;
   }
   file << std::endl;
   
   /* Write the cell types: */
   file << "CELL_TYPES " << num_cells << std::endl;
   for (int i = 0; i < num_cells; i++) {
      num_cell_points = cells.size(i);
      if (num_cell_points == 2)
         file << "3" << std::endl;
      else if (num_cell_points == 3)
         file << "5" << std::endl;
      else if (num_cell_points == 4)
         file << "9" << std::endl;
      else if (num_cell_points == 6)
         file << "13" << std::endl;
      else if (num_cell_points == 8)
         file << "12" << std::endl;
      else if (num_cell_points == 12)
         file << "16" << std::endl;
      else {
         PAMPA_CHECK(true, "wrong cell type");
      }
   }
   file << std::endl;
   
   /* Write the header for the cell data: */
   file << "CELL_DATA " << num_cells << std::endl << std::endl;
   
   /* Write the cell materials: */
   file << "SCALARS materials double 1" << std::endl;
   file << "LOOKUP_TABLE default" << std::endl;
   for (int i = 0; i < num_cells; i++)
      file << materials(i)+1 << std::endl;
   file << std::endl;
   
   /* Write the cell indices in the nodal mesh: */
   if (!(nodal_indices.empty())) {
      file << "SCALARS nodal_indices double 1" << std::endl;
      file << "LOOKUP_TABLE default" << std::endl;
      for (int i = 0; i < num_cells; i++)
         file << nodal_indices(i) << std::endl;
      file << std::endl;
   }
   
   return 0;
   
}

/* Write a solution vector to a .vtk file: */
int vtk::write(const std::string& filename, const std::string& name, const Vec& v, int num_cells, 
   int num_groups, int num_directions) {
   
   /* Open the output file: */
   std::ofstream file(filename, std::ios_base::app);
   PAMPA_CHECK(!file.is_open(), "unable to open " + filename);
   
   /* Set the precision: */
   file << std::scientific;
   file << std::setprecision(VTK_PRECISION);
   
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
         int iv0 = g*num_directions+m, div = num_directions*num_groups;
         for (int iv = iv0, i = 0; i < num_cells; i++) {
            file << data[iv] << std::endl;
            iv += div;
         }
         file << std::endl;
         
      }
   }
   
   /* Restore the array with the raw data: */
   PETSC_CALL(VecRestoreArray(v, &data));
   
   return 0;
   
}
