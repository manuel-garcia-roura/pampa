#pragma once

#include <vector>

/* The Vector2D class: */
template <class T>
class Vector2D {
   
   private:
      
      /* Array dimensions: */
      int n1, n;
      
      /* First index for each row: */
      std::vector<int> i0;
      
      /* Array data: */
      std::vector<T> v;
      
      /* Build the indexing: */
      void build(const std::vector<int>& n2) {
         
         /* Get the first index for each row: */
         i0.resize(n1+1);
         for (int i = 0; i < n1; i++)
            i0[i+1] = i0[i] + n2[i];
         
         /* Get the total size: */
         n = i0[n1];
         
      }
   
   public:
      
      /* The Vector2D default constructor: */
      Vector2D() {}
      
      /* The Vector2D constructor: */
      Vector2D(int n1, const std::vector<int>& n2, const T& x0 = T()) : n1(n1) 
         {build(n2); v.resize(n, x0);}
      
      /* The Vector2D destructor: */
      ~Vector2D() {}
      
      /* Subscript operator (write): */
      T& operator() (int i1, int i2) {return v[i0[i1]+i2];}
      
      /* Subscript operator (read): */
      T operator() (int i1, int i2) const {return v[i0[i1]+i2];}
      
      /* Subscript operator (read): */
      const T* operator() (int i1) const {return &(v[i0[i1]]);}
      
      /* Get the size of row i1: */
      int size(int i1) const {return i0[i1+1]-i0[i1];}
   
};
