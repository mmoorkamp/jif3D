#ifndef NUMUTIL_H_
#define NUMUTIL_H_

//! Determine the sign of a numeric type by comparing it with zero
template <class NumericType>
int sign (NumericType n)
{
   if (n >= 0) return 1;
   return -1;
}

#endif /*NUMUTIL_H_*/
