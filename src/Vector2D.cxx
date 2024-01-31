#include "Vector2D.hxx"

/* Build the indexing: */
template <class T>
void Vector2D<T>::build(const std::vector<int>& n2) {
   
   /* Get the first index for each row: */
   i0.resize(n1);
   for (int i = 0; i < n1-1; i++)
      i0[i+1] = i0[i] + n2[i];
   
   /* Get the total size: */
   n = i0[n1-1] + n2[n1-1];
   
}
