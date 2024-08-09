#pragma once

#include <vector>
#include <algorithm>

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
      const T& operator() (int i1) const {return v[i1];}
      
      /* Get the last element (write): */
      T& back() {return v.back();}
      
      /* Get the last element (read): */
      const T& back() const {return v.back();}
      
      /* Resize the array: */
      void resize(int n1, const T& x0 = T()) {this->n1 = n1; v.resize(n1, x0);}
      
      /* Free the array: */
      void free() {resize(0);}
      
      /* Get the array size: */
      int size() const {return n1;}
      
      /* Check if the array is empty: */
      bool empty() const {return n1 == 0;}
      
      /* Remove all elements by value: */
      void remove(const T& x) {v.erase(std::remove(v.begin(), v.end(), x), v.end()); n1 = v.size();}
      
      /* Push back an element: */
      void pushBack(const T& x) {v.push_back(x); n1++;}
      
      /* Check if the array contains a value: */
      bool find(const T& x) const {return std::find(v.begin(), v.end(), x) != v.end();}
      
      /* Get the index of the closest value larger than a given value: */
      int lowerBound(const T& x) const {return std::lower_bound(v.begin(), v.end(), x) - v.begin();}
   
};
