#pragma once

#include <vector>

#include "Array1D.hxx"

/* The Vector2D class: */
template <class T>
class Vector2D {
   
   private:
      
      /* Vector dimensions: */
      int n1 = 0, n = 0;
      
      /* First index for each row: */
      std::vector<int> i0{0};
      
      /* Vector data: */
      std::vector<T> v;
      
      /* Build the indexing: */
      void build(const Array1D<int>& n2) {
         
         /* Get the first index for each row: */
         i0.resize(n1+1);
         for (int i = 0; i < n1; i++)
            i0[i+1] = i0[i] + n2(i);
         
         /* Get the total size: */
         n = i0[n1];
         
      }
   
   public:
      
      /* The Vector2D default constructor: */
      Vector2D() {}
      
      /* The Vector2D constructor: */
      Vector2D(int n1, const Array1D<int>& n2, const T& x0 = T()) {resize(n1, n2, x0);}
      
      /* The Vector2D destructor: */
      ~Vector2D() {}
      
      /* Subscript operator (write): */
      T& operator() (int i1, int i2) {return v[i0[i1]+i2];}
      
      /* Subscript operator (read): */
      const T& operator() (int i1, int i2) const {return v[i0[i1]+i2];}
      
      /* Subscript operator (write): */
      T* operator() (int i1) {return &(v[i0[i1]]);}
      
      /* Subscript operator (read): */
      const T* operator() (int i1) const {return &(v[i0[i1]]);}
      
      /* Resize the array: */
      void resize(int n1, const Array1D<int>& n2, const T& x0 = T()) 
         {this->n1 = n1; build(n2); v.resize(n, x0);}
      
      /* Get the size of row i1: */
      int size(int i1) const {return i0[i1+1]-i0[i1];}
      
      /* Push back a row: */
      void pushBack(const Array1D<T>& u) {
         
         /* Get the new dimensions: */
         n1++;
         n += u.size();
         
         /* Get the first index for the next row: */
         i0.push_back(n);
         
         /* Copy the array elements: */
         v.reserve(n);
         for (int i = 0; i < u.size(); i++)
            v.push_back(u(i));
         
      }
   
};
