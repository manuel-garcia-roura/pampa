#pragma once

#include <vector>

#include "Array1D.hxx"

/* The Vector3D class: */
template <class T>
class Vector3D {
   
   private:
      
      /* Vector dimensions: */
      int n1 = 0, n3 = 0, n = 0;
      
      /* First index for each row: */
      std::vector<int> i0{0};
      
      /* Vector data: */
      std::vector<T> v;
      
      /* Build the indexing: */
      void build(const Array1D<int>& n2) {
         
         /* Get the first index for each row: */
         i0.resize(n1+1);
         for (int i = 0; i < n1; i++)
            i0[i+1] = i0[i] + n2(i)*n3;
         
         /* Get the total size: */
         n = i0[n1];
         
      }
   
   public:
      
      /* The Vector3D default constructor: */
      Vector3D() {}
      
      /* The Vector3D constructor: */
      Vector3D(int n1, const Array1D<int>& n2, int n3, const T& x0 = T()) {resize(n1, n2, n3, x0);}
      
      /* The Vector3D destructor: */
      ~Vector3D() {}
      
      /* Subscript operator (write): */
      T& operator() (int i1, int i2, int i3) {return v[i0[i1]+i2*n3+i3];}
      
      /* Subscript operator (read): */
      const T& operator() (int i1, int i2, int i3) const {return v[i0[i1]+i2*n3+i3];}
      
      /* Subscript operator (write): */
      T* operator() (int i1, int i2) {return &(v[i0[i1]+i2*n3]);}
      
      /* Subscript operator (read): */
      const T* operator() (int i1, int i2) const {return &(v[i0[i1]+i2*n3]);}
      
      /* Reserve the vector memory: */
      void reserve(int n) {v.reserve(n);}
      
      /* Resize the vector: */
      void resize(int n1, const Array1D<int>& n2, int n3, const T& x0 = T()) 
         {this->n1 = n1; this->n3 = n3; build(n2); v.resize(n, x0);}
      
      /* Free the vector: */
      void free() {resize(0, Array1D<int>(), 0);}
      
      /* Get the vector size: */
      int size() const {return n;}
      
      /* Get the size of row i1: */
      int size(int i1) const {return (i0[i1+1]-i0[i1]) / n3;}
   
};
