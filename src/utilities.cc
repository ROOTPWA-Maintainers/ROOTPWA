//
// package that contains various small utility functions
//


#include "utilities.h"


// create virtual instances of maxPrecisionValue__
template class maxPrecisionValue__<short int>;
template class maxPrecisionValue__<unsigned short int>;
template class maxPrecisionValue__<int>;
template class maxPrecisionValue__<unsigned int>;
template class maxPrecisionValue__<long int>;
template class maxPrecisionValue__<unsigned long int>;

template class maxPrecisionValue__<float>;
template class maxPrecisionValue__<double>;
template class maxPrecisionValue__<long double>;
