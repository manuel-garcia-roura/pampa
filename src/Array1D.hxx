#pragma once

#include <vector>

/* The Array1D class: */
template <class T>
class Array1D {
   
   private:
      
      /* Array dimensions: */
      int n1;
      
      /* Array data: */
      std::vector<T> v;
   
   public:
      
      /* The Array1D default constructor: */
      Array1D() {}
      
      /* The Array1D constructor: */
      Array1D(int n1) : n1(n1) {v.resize(n1);}
      
      /* The Array1D constructor with a default value: */
      Array1D(int n1, T x0) : n1(n1) {v.resize(n1, x0);}
      
      /* The Array1D destructor: */
      ~Array1D() {}
      
      /* Subscript operator (write): */
      T& operator() (int i1) {return v[i1];}
      
      /* Subscript operator (read): */
      T operator() (int i1) const {return v[i1];}
   
};
