#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

/* The utils namespace: */
namespace utils {
   
   /* Trim and remove double spaces from a string: */
   void clean(std::string &s);
   
   /* Get the next line from a file stream: */
   std::vector<std::string> get_next_line(std::ifstream &file);
   
   /* Read a vector with n elements of type double from a file stream: */
   bool read(std::vector<double> &v, int n, std::ifstream &file);
   
   /* Read a vector with n elements of type std::vector<double> of size m from a file stream: */
   bool read(std::vector<std::vector<double>> &v, int n, int m, std::ifstream &file);
   
   /* Read a vector with n elements of type std::vector<int> of unknown size from a file stream: */
   bool read(std::vector<std::vector<int>> &v, int n, std::ifstream &file);
   
}
