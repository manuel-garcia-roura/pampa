#pragma once

#include "utils.hxx"

/* The Function struct: */
struct Function {
   
   private:
      
      /* Time values: */
      Array1D<double> t0;
      
      /* Function values (y = f(t)): */
      Array1D<double> y0;
   
   public:
      
      /* The Function default constructor: */
      Function() {}
      
      /* The Function constructor for a constant function: */
      Function(double y0) : y0(Array1D<double>{1, y0}) {}
      
      /* The Function constructor for a linear function: */
      Function(const Array1D<double>& t0, const Array1D<double>& y0) : t0(t0), y0(y0) {}
      
      /* The Function destructor: */
      ~Function() {}
      
      /* Check if the function is empty: */
      bool empty() const {return y0.empty();}
      
      /* Call operator: */
      double operator() (double t) const {
         
         /* Use the first value: */
         if (t0.empty() || t < t0(0))
            return y0(0);
         
         /* Use the last value: */
         if (t > t0.back())
            return y0.back();
         
         /* Use a pair of values: */
         int i2 = t0.lowerBound(t);
         int i1 = i2 - 1;
         double f = (t-t0(i1)) / (t0(i2)-t0(i1));
         
         /* Interpolate: */
         return (1.0-f)*y0(i1) + f*y0(i2);
         
      }
   
};
