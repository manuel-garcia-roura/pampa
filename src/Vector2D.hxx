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
      void build(const std::vector<int>& n2);
   
   public:
      
      /* The Vector2D default constructor: */
      Vector2D() {}
      
      /* The Vector2D constructor: */
      Vector2D(int n1, const std::vector<int>& n2, T x0 = T()) : n1(n1) 
         {build(n2); v.resize(n, x0);}
      
      /* The Vector2D destructor: */
      ~Vector2D() {}
      
      /* Subscript operator (write): */
      T& operator() (int i1, int i2) {return v[i0[i1]+i2];}
      
      /* Subscript operator (read): */
      T operator() (int i1, int i2) const {return v[i0[i1]+i2];}
      
      /* Subscript operator (read): */
      const T* operator() (int i1) const {return &(v[i0[i1]]);}
   
};
