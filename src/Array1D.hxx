#pragma once

#include "utils.hxx"

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
      
      /* Set the array values to a single value: */
      void fill(const T& x) {std::fill(v.begin(), v.end(), x);}
      
      /* Get the array size: */
      int size() const {return n1;}
      
      /* Check if the array is empty: */
      bool empty() const {return n1 == 0;}
      
      /* Remove all elements by value: */
      void remove(const T& x) {v.erase(std::remove(v.begin(), v.end(), x), v.end()); n1 = v.size();}
      
      /* Push back an element: */
      void pushBack(const T& x) {v.push_back(x); n1++;}
      
      /* Check if the array contains a value: */
      int find(const T& x) const 
         {auto it = std::find(v.begin(), v.end(), x); return it != v.end() ? it - v.begin() : -1;}
      
      /* Get the index of the closest value larger than a given value: */
      int lowerBound(const T& x) const {return std::lower_bound(v.begin(), v.end(), x) - v.begin();}
      
      /* Get the minimum value: */
      const T& minValue() const {return *(std::min_element(v.begin(), v.end()));}
      
      /* Get the maximum value: */
      const T& maxValue() const {return *(std::max_element(v.begin(), v.end()));}
      
      /* Get the index of the minimum value: */
      int minIndex() const {return std::distance(v.begin(), std::min_element(v.begin(), v.end()));}
      
      /* Get the index of the maximum value: */
      int maxIndex() const {return std::distance(v.begin(), std::max_element(v.begin(), v.end()));}
   
};
