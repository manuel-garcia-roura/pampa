#pragma once

#include <vector>

/* The Array3D class: */
template <class type>
class Array3D {
   
   private:
      
      /* Array dimensions: */
      int n1, n2, n3;
      
      /* Array data: */
      std::vector<type> data;
   
   public:
      
      /* The Array3D default constructor: */
      Array3D() {}
      
      /* The Array3D constructor: */
      Array3D(int n1, int n2, int n3) : n1(n1), n2(n2), n3(n3) {data.resize(n1*n2*n3);}
      
      /* The Array3D destructor: */
      ~Array3D() {}
      
      /* Subscript operator (write): */
      type& operator() (int i1, int i2, int i3) {return data[i1*n2*n3+i2*n3+i3];}
      
      /* Subscript operator (read): */
      type operator() (int i1, int i2, int i3) const {return data[i1*n2*n3+i2*n3+i3];}
   
};
