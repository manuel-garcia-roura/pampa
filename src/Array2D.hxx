#pragma once

#include <vector>

/* The Array2D class: */
template <class T>
class Array2D {
   
   private:
      
      /* Array dimensions: */
      int n1, n2;
      
      /* Array data: */
      std::vector<T> v;
   
   public:
      
      /* The Array2D default constructor: */
      Array2D() {}
      
      /* The Array2D constructor: */
      Array2D(int n1, int n2, const T& x0 = T()) : n1(n1), n2(n2) {v.resize(n1*n2, x0);}
      
      /* The Array2D destructor: */
      ~Array2D() {}
      
      /* Subscript operator (write): */
      T& operator() (int i1, int i2) {return v[i1*n2+i2];}
      
      /* Subscript operator (read): */
      T operator() (int i1, int i2) const {return v[i1*n2+i2];}
      
      /* Subscript operator (read): */
      const T* operator() (int i1) const {return &(v[i1*n2]);}
   
};
