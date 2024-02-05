#pragma once

#include <vector>

/* The Array1D class: */
template <class T>
class Array1D {
   
   private:
      
      /* Array dimensions: */
      int n1 = 0;
      
      /* Array data: */
      std::vector<T> v;
   
   public:
      
      /* The Array1D default constructor: */
      Array1D() {}
      
      /* The Array1D constructor: */
      Array1D(int n1, const T& x0 = T()) {resize(n1, x0);}
      
      /* The Array1D destructor: */
      ~Array1D() {}
      
      /* Subscript operator (write): */
      T& operator() (int i1) {return v[i1];}
      
      /* Subscript operator (read): */
      T operator() (int i1) const {return v[i1];}
      
      /* Resize the array: */
      void resize(int n1, const T& x0 = T()) {this->n1 = n1; v.resize(n1, x0);}
      
      /* Get the array size: */
      int size() const {return n1;}
      
      /* Check if the array is empty: */
      bool empty() const {return n1 == 0;}
      
      /* Remove all elements by value: */
      void remove(const T& x) {v.erase(std::remove(v.begin(), v.end(), x), v.end()); n1 = v.size();}
   
};
