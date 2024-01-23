#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "Model.hxx"
#include "Mesh.hxx"
#include "CartesianMesh.hxx"
#include "UnstructuredExtrudedMesh.hxx"
#include "Material.hxx"
#include "utils.hxx"

/* The Parser class: */
class Parser {
   
   public:
      
      /* The Parser constructor: */
      Parser() {}
      
      /* The Parser destructor: */
      ~Parser() {}
      
      /* Read a plain-text input file: */
      int read(const std::string& filename, Model& model);
   
};
