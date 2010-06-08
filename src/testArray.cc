//
// tests matrix implementation using vector of vector
// run using:
// g++ -Wall testArray.cxx -o testArray && ./testArray
//


#include <iostream>
#include <complex>

#include "utilities.h"


using namespace std;


int main(int argc, char** argv)
{
  double a[2][3][4];
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      for (unsigned int k = 0; k < 4; ++k)
	a[i][j][k] = i * 100 + j * 10 + k;

  if (1) {
  
    double*** x               = 0;
    const unsigned int dim[3] = {2, 3, 4};
    allocate3DArray(x, dim);
    for (unsigned int i = 0; i < 2; ++i)
      for (unsigned int j = 0; j < 3; ++j)
	for (unsigned int k = 0; k < 4; ++k)
	  x[i][j][k] = i * 100 + j * 10 + k;

    for (unsigned int i = 0; i < 2; ++i)
      for (unsigned int j = 0; j < 3; ++j)
	for (unsigned int k = 0; k < 4; ++k) {
	  if (a[i][j][k] != a[i][j][k])
	    cout << "ERROR!" << endl;
	  cout << "x[" << i << "][" << j << "][" << k << "] = " << x[i][j][k] << " vs. "
	       << "a[" << i << "][" << j << "][" << k << "] = " << a[i][j][k] << endl;
	}

    delete3DArray(x, dim);
    
  }


  if (0) {

    double*            x      = 0;
    const unsigned int dim[3] = {2, 3, 4};
    const unsigned int nmbDim = 3;
    allocatePseudoNdimArray(x, dim, nmbDim);
    for (unsigned int i = 0; i < 2; ++i)
      for (unsigned int j = 0; j < 3; ++j)
	for (unsigned int k = 0; k < 4; ++k) {
	  const unsigned int indices[3] = {i, j, k};
	  const unsigned int offset     = indicesToOffset(indices, dim, nmbDim);
	  x[offset] = i * 100 + j * 10 + k;
	}

    double* y = (double*)a;
    for (unsigned int i = 0; i < 2; ++i)
      for (unsigned int j = 0; j < 3; ++j)
	for (unsigned int k = 0; k < 4; ++k) {
	  const unsigned int indices[3] = {i, j, k};
	  const unsigned int offset     = indicesToOffset(indices, dim, nmbDim);
	  if (x[offset] != a[i][j][k])
	    cout << "ERROR!" << endl;
	  cout << "x[" << offset << "] = " << x[offset] << " vs. "
	       << "a[" << i << "][" << j << "][" << k << "] = " << a[i][j][k] << endl;
	  if (x[offset] != y[offset])
	    cout << "ERROR!" << endl;
	  cout << "x[" << offset << "] = " << x[offset] << " vs. "
	       << "y[" << offset << "] = " << y[offset] << endl;
	  unsigned int testIndices[3];
	  offsetToIndices(offset, dim, nmbDim, testIndices);
	  cout << " indices {";
	  for (unsigned int d = 0; d < nmbDim; ++d) {
	    if (testIndices[d] != indices[d])
	      cout << "ERROR!  ";
	    cout << testIndices[d] << "  ";
	  }
	  cout << "} vs. {" << i << "  " << j << "  " << k << "}" << endl;
	}
    
    delete[] x;
  }
  
  return 0;
}
