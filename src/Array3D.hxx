#pragma once

#include <vector>

/* The Array3D class: */
template <class T>
class Array3D {
   
   private:
      
      /* Array dimensions: */
      int n1 = 0, n2 = 0, n3 = 0;
      
      /* Array data: */
      std::vector<T> v;
   
   public:
      
      /* The Array3D default constructor: */
      Array3D() {}
      
      /* The Array3D constructor: */
      Array3D(int n1, int n2, int n3, const T& x0 = T()) {resize(n1, n2, n3, x0);}
      
      /* The Array3D destructor: */
      ~Array3D() {}
      
      /* Subscript operator (write): */
      T& operator() (int i1, int i2, int i3) {return v[i1*n2*n3+i2*n3+i3];}
      
      /* Subscript operator (read): */
      T operator() (int i1, int i2, int i3) const {return v[i1*n2*n3+i2*n3+i3];}
      
      /* Subscript operator (write): */
      T* operator() (int i1, int i2) {return &(v[i1*n2*n3+i2*n3]);}
      
      /* Subscript operator (read): */
      const T* operator() (int i1, int i2) const {return &(v[i1*n2*n3+i2*n3]);}
      
      /* Resize the array: */
      void resize(int n1, int n2, int n3, const T& x0 = T()) 
         {this->n1 = n1; this->n2 = n2; this->n3 = n3; v.resize(n1*n2*n3, x0);}
   
};
