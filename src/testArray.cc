//
// tests matrix implementation using vector of vector
// run using:
// g++ -Wall testArray.cxx -o testArray && ./testArray
//


#include <iostream>
#include <complex>

#include "utilities.h"


using namespace std;


template<typename T>
inline
void
offsetToIndices(const T  offset,   // one-dimensional array index
		const T* dim,      // extents of n-dimensional array
		const T  nmbDim,   // number of dimensions
		T*       indices)  // indices to map onto
{
  T index = offset;
  for (T i = nmbDim - 1; i >= 1; --i) {
    indices[i] = index % dim[i];
    index      = index / dim[i];
  }
  indices[0] = index;
}


void allocate(double***&         m,
	      const unsigned int d1 = 2,
	      const unsigned int d2 = 3,
	      const unsigned int d3 = 4)
{
  if (m) {
    for (int i = d1 - 1; i >= 0; --i)
      for (int j = d2 - 1; j >= 0; --j)
	if (m[i][j])
	  delete[] m[i][j];
    for (int i = d1 - 1; i >= 0; --i)
      if (m[i])
	delete[] m[i];
    delete[] m;
  }
    
  m = new double** [d1];
  for (unsigned int i = 0; i < d1; ++i) {
    m[i] = new double*[d2];
    for (unsigned int j = 0; j < d2; ++j) {
      m[i][j] = new double[d3];
      for (unsigned int k = 0; k < d3; ++k)
	m[i][j][k] = i * 100 + j * 10 + k;
    }
  }
}


int main(int argc, char** argv)
{
  double a[2][3][4];
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      for (unsigned int k = 0; k < 4; ++k)
	a[i][j][k] = i * 100 + j * 10 + k;

  if (0) {
  
    for (unsigned int i = 0; i < 2; ++i) {
      cout << i << endl;
      for (unsigned int j = 0; j < 3; ++j) {
	for (unsigned int k = 0; k < 4; ++k)
	  cout << a[i][j][k] << "   ";
	cout << endl;
      }
    }
    cout << endl;

    double* b = a[1][2];
    for (unsigned int k = 0; k < 4; ++k)
      cout << b[k] << "   ";
    cout << endl;

    // does not work
    //   double** c = &(*a[0]);
    //   for (unsigned int j = 0; j < 3; ++j) {
    //     for (unsigned int k = 0; k < 4; ++k)
    //       cout << c[j][k] << "   ";
    //     cout << endl;
    //   }

    double*** d = NULL;
    allocate(d);
    cout << endl << "!!!HERE" << endl;
    allocate(d);
    cout << "!!!HERE2" << endl;
    for (unsigned int i = 0; i < 2; ++i) {
      cout << i << endl;
      for (unsigned int j = 0; j < 3; ++j) {
	for (unsigned int k = 0; k < 4; ++k)
	  cout << d[i][j][k] << "   ";
	cout << endl;
      }
    }
  }

  if (1) {

    double*            x      = 0;
    const unsigned int dim[3] = {2, 3, 4};
    const unsigned int nmbDim = 3;
    allocatePseudoNdimArray<double, unsigned int>(x, dim, nmbDim);
    for (unsigned int i = 0; i < 2; ++i)
      for (unsigned int j = 0; j < 3; ++j)
	for (unsigned int k = 0; k < 4; ++k) {
	  const unsigned int indices[3] = {i, j, k};
	  const unsigned int offset     = indicesToOffset<unsigned int>(indices, dim, nmbDim);
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
	       << " a[" << i << "][" << j << "][" << k << "] = " << a[i][j][k] << endl;
	  if (x[offset] != y[offset])
	    cout << "ERROR!" << endl;
	  cout << "x[" << offset << "] = " << x[offset] << " vs. "
	       << " y[" << offset << "] = " << y[offset] << endl;
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
    
  }
  
  return 0;
}
