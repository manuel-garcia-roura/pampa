#pragma once

#include "petsc.hxx"
#include "utils.hxx"

/* The input namespace: */
namespace input {
   
   /* Trim and remove tabs and double spaces from a string: */
   void clean(std::string& s);
   
   /* Get the next line from a file stream: */
   std::vector<std::string> get_next_line(std::ifstream& file);
   
   /* Read an array with n elements of type double from a file stream: */
   int WARN_UNUSED read(Array1D<double>& v, unsigned int n, double x1, double x2, 
      std::ifstream& file);
   
   /* Read an array with n elements of type int from a file stream: */
   int WARN_UNUSED read(Array1D<int>& v, unsigned int n, int x1, int x2, std::ifstream& file);
   
   /* Read an array with (n, m) elements of type int from a file stream: */
   int WARN_UNUSED read(Array2D<double>& v, unsigned int n, unsigned int m, double x1, double x2, 
      std::ifstream& file);
   
   /* Read a vector with n rows of type double and total size nt from a file stream: */
   int WARN_UNUSED read(Vector2D<double>& v, unsigned int n, unsigned int nt, double x1, double x2, 
      std::ifstream& file);
   
   /* Read a vector with n rows of type int and total size nt from a file stream: */
   int WARN_UNUSED read(Vector2D<int>& v, unsigned int n, unsigned int nt, int x1, int x2, 
      std::ifstream& file);
   
   /* Read a vector with (n1, n2, n3) elements of type double from a file stream: */
   int WARN_UNUSED read(Vector3D<double>& v, unsigned int n1, const Array1D<int>& n2, 
      unsigned int n3, double x1, double x2, std::ifstream& file);
   
   /* Read a bool value from a string: */
   int WARN_UNUSED read(bool& q, const std::string& s);
   
   /* Read an int value from a string: */
   int WARN_UNUSED read(int& x, int x1, int x2, const std::string& s);
   
   /* Read a double value from a string: */
   int WARN_UNUSED read(double& x, double x1, double x2, const std::string& s);
   
   /* Read a PETSc norm type from a string: */
   int WARN_UNUSED read(NormType& p, const std::string& s);
   
   /* Read a function from a file stream: */
   int WARN_UNUSED read(Function& f, double x1, double x2, const std::vector<std::string>& line, 
      unsigned int& i, std::ifstream& file);
   
   /* Read multiple functions from a file stream: */
   int WARN_UNUSED read(Array1D<Function>& f, int m, double x1, double x2, 
      const std::vector<std::string>& line, unsigned int& i, std::ifstream& file);
   
   /* Read a boundary condition from a line: */
   int WARN_UNUSED read(BoundaryCondition& bc, const std::vector<std::string>& line, 
      unsigned int& i, std::ifstream& file);
   
}
