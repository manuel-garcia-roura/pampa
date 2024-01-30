#pragma once

#include <vector>

/* The Array2D class: */
template <class type>
class Array2D {
   
   private:
      
      /* Array dimensions: */
      int n1, n2;
      
      /* Array data: */
      std::vector<type> data;
   
   public:
      
      /* The Array2D default constructor: */
      Array2D() {}
      
      /* The Array2D constructor: */
      Array2D(int n1, int n2) : n1(n1), n2(n2) {data.resize(n1*n2);}
      
      /* The Array2D destructor: */
      ~Array2D() {}
      
      /* Subscript operator (write): */
      type& operator() (int i1, int i2) {return data[i1*n2+i2];}
      
      /* Subscript operator (read): */
      type operator() (int i1, int i2) const {return data[i1*n2+i2];}
   
};
