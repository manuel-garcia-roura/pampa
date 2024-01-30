#pragma once

#include <vector>

/* The Array1D class: */
template <class type>
class Array1D {
   
   private:
      
      /* Array dimensions: */
      int n1;
      
      /* Array data: */
      std::vector<type> data;
   
   public:
      
      /* The Array1D default constructor: */
      Array1D() {}
      
      /* The Array1D constructor: */
      Array1D(int n1) : n1(n1) {data.resize(n1);}
      
      /* The Array1D destructor: */
      ~Array1D() {}
      
      /* Subscript operator (write): */
      type& operator() (int i1) {return data[i1];}
      
      /* Subscript operator (read): */
      type operator() (int i1) const {return data[i1];}
   
};
