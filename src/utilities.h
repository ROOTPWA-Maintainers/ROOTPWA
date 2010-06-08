//
// package that contains various small utility functions
//


#ifndef utilities_h
#define utilities_h


#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <vector>
#include <complex>


// cint has problems parsing glob.h
#ifndef __CINT__
#include <glob.h>
#else
struct glob_t;
int glob(const char *,
	 int,
	 int(*)(const char*, int),
	 glob_t*);
void globfree(glob_t *);
#endif

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TComplex.h"
#include "TMatrixT.h"


//////////////////////////////////////////////////////////////////////////////////
// macros and functions for printout and formatting

// cuts out block "className::methodName" from __PRETTY_FUNCTION__ output
inline
std::string
getClassMethod__(std::string prettyFunction)
{
  size_t pos = prettyFunction.find("(");
  if (pos == std::string::npos)
    return prettyFunction;           // something is not right
  prettyFunction.erase(pos);         // cut away signature
  pos = prettyFunction.rfind(" ");
  if (pos == std::string::npos)
    return prettyFunction;           // something is not right
  prettyFunction.erase(0, pos + 1);  // cut away return type
  return prettyFunction;
}

// macros for printinf errors, warnings, and infos
#define printErr  std::cerr << "!!! " << __PRETTY_FUNCTION__ << " [" << __FILE__ << ":" << __LINE__ << "]: error: "   << std::flush
#define printWarn std::cerr << "??? " << __PRETTY_FUNCTION__ << " [" << __FILE__ << ":" << __LINE__ << "]: warning: " << std::flush
#define printInfo std::cout << ">>> " << getClassMethod__(__PRETTY_FUNCTION__) << "(): info: "  << std::flush


template<typename T> class maxPrecisionValue__;

// output stream manipulator that prints a value with its maximum precision
template<typename T>
inline
maxPrecisionValue__<T>
maxPrecision(const T& value)
{ return maxPrecisionValue__<T>(value); }

// output stream manipulator that prints a value with its maximum precision
// in addition manipulator reserves space so that values will align
template<typename T>
inline
maxPrecisionValue__<T>
maxPrecisionAlign(const T& value)
{ return maxPrecisionValue__<T>(value, maxPrecisionValue__<T>::ALIGN); }

// output stream manipulator that prints a value with maximum precision for double
template<typename T>
inline
maxPrecisionValue__<T>
maxPrecisionDouble(const T& value)
{ return maxPrecisionValue__<T>(value, maxPrecisionValue__<T>::DOUBLE); }

// general helper class that encapsulates a value of type T
template<typename T>
class maxPrecisionValue__ {
public:
  enum modeEnum { PLAIN,
		  ALIGN,
		  DOUBLE};  // forces precision for double
  maxPrecisionValue__(const T&       value,
		      const modeEnum mode = PLAIN)
    : _value(value),
      _mode (mode)
  { }
  std::ostream& print(std::ostream& out) const
  {
    const int nmbDigits = (_mode != DOUBLE) ? std::numeric_limits<T>::digits10 + 1
                                            : std::numeric_limits<double>::digits10 + 1;
    std::ostringstream s;
    s.precision(nmbDigits);
    s.setf(std::ios_base::scientific, std::ios_base::floatfield);
    s << _value;
    switch (_mode) {
    case ALIGN:
      return out << std::setw(nmbDigits + 7) << s.str();  // make space for sign, dot, and exponent
    case PLAIN: case DOUBLE: default:
      return out << s.str();
    }
  }
private:
  const T& _value;
  modeEnum _mode;
};

template<typename T>
inline
std::ostream& operator << (std::ostream&                 out,
			   const maxPrecisionValue__<T>& value)
{ return value.print(out); }


// simple stream operators for some STL classes
template<typename T1, typename T2>
inline
std::ostream&
operator << (std::ostream&            out,
             const std::pair<T1, T2>& pair)
{
  out << "(" << pair.first << ", " << pair.second << ")";
  return out;
}


template<typename T>
inline
std::ostream&
operator << (std::ostream&         out,
	     const std::vector<T>& vec)
{
  out << "{";
  for (unsigned int i = 0; i < (vec.size() - 1); ++i)
    out << vec[i] << ", ";
  out << vec[vec.size() - 1] << "}";
  return out;
}


// simple stream operators for some common ROOT classes
inline
std::ostream&
operator << (std::ostream&   out,
             const TVector3& vec)
{
  out << "(" << vec.X() << ", " << vec.Y() << ", " << vec.Z() << ")";
  return out;
}


inline
std::ostream&
operator << (std::ostream&         out,
             const TLorentzVector& vec)
{
  out << "(" << vec.X() << ", " << vec.Y() << ", " << vec.Z() << "; " << vec.T() << ")";
  return out;
}


inline
std::ostream&
operator << (std::ostream&   out,
             const TComplex& c)
{
  out << "(" << c.Re() << ", " << c.Im() << ")";
  return out;
}


template<typename T>
std::ostream&
operator << (std::ostream&      out,
             const TMatrixT<T>& A)
{
  for (int row = 0; row < A.GetNrows(); ++row) {
    out << "row " << row << " = (";
    for (int col = 0; col < A.GetNcols(); ++col) {
      out << A[row][col];
      if (col < A.GetNcols() - 1)
        out << ", ";
    }
    if (row < A.GetNrows() - 1)
      out << "), " << std::endl;
    else
      out << ")";
  }
  return out;
}


// indents output by offset
inline
void
indent(std::ostream&      out,
       const unsigned int offset)
{
  for (unsigned int i = 0; i < offset; ++i)
    out << " ";
}


// extracts sign of a value
template<typename T>
std::string
sign(const T& val)
{
  if (val < 0)
    return "-";
  if (val == 0)
    return "0";
  return "+";
}


// converts sign character into int
inline
int
sign(const char s)
{
  if (s == '+')
    return +1;
  if (s == '-')
    return -1;
  return 0;
}


// extracts sign from value
template<typename T>
T signum(const T& val)
{
  if (val < 0)
    return -1;
  if (val > 0)
    return +1;
  return 0;
}


// indicates progess by printing percentage complete
inline
std::ostream&
progressIndicator(const long    currentPos,
                  const long    nmbTotal,
                  const int     nmbSteps = 10,
                  std::ostream& out      = std::cout)
{
  const double step = nmbTotal / (double)nmbSteps;
  if ((nmbTotal >= 0) and ((int)(currentPos / step) - (int)((currentPos - 1) / step) != 0))
    out << "    " << std::setw(3) << (int)(currentPos / step) * nmbSteps << " %" << std::endl;
  return out;
}


inline
std::vector<std::string>
tokenizeString(const std::string& in,
	       const std::string& delimiter)
{
  // skip delimiters at beginning.
  std::vector<std::string> tokens;
  std::string::size_type begin = in.find_first_not_of(delimiter, 0);
  std::string::size_type end   = in.find_first_of    (delimiter, begin);
  while ((begin != std::string::npos) or (end != std::string::npos)) {
    tokens.push_back(in.substr(begin, end - begin));
    // skip delimiters
    begin = in.find_first_not_of(delimiter, end);
    end   = in.find_first_of    (delimiter, begin);
  }
  return tokens;
}


//////////////////////////////////////////////////////////////////////////////////
// functions for dynamic allocation of n-dimensional arrays

#ifndef protect__  // macro that converts arguments containing commas into single string; useful for more complex types T
#define protect__(x...) x
#endif
#define vector2(T) std::vector<std::vector<T> >
#define vector3(T) std::vector<std::vector<std::vector<T> > >
#define vector4(T) std::vector<std::vector<std::vector<std::vector<T> > > >


// arrays based on pointer^n variables that can be passed to
// functions without knowing the array size at compile time
// elements can be accessed in the canonical way using []^n operators
// this implementation is, due to the levels of indirection,
// potentially less peformant than pseudo-arrays, which on the other
// hand make element access more cumbersome
// 
// array dimensions are passed as unsigned int arrays, where the order
// reflects the order of the indices:
//
// array[x_0][x_1] ... [x_(n - 1)]
//        |    |            |
//     dim[0]  |            |
//          dim[1] ...   dim[n - 1]
//
// array memory has to be freed in reverser order of its allocation

template<typename D, typename T>
void
delete2DArray(D**&    array,   // two-dimensional array to delete
              const T dim[2])  // extents of two-dimensional array
{
  for (T i = 0; i < dim[0]; ++i)
    if (array[i])
      delete[] array[i];
  delete[] array;
  array = NULL;
}

template<typename D, typename T>
void
allocate2DArray(D**&     array,           // two-dimensional array to create
                const T  dim[2],          // extents of two-dimensional array
                const D* defaultVal = 0)  // optional default value
{
  if (array)
    delete2DArray<D, T>(array, dim);
  array = new D*[dim[0]];
  for (T i = 0; i < dim[0]; ++i) {
    array[i] = new D[dim[1]];
    if (defaultVal)
      for (T j = 0; j < dim[1]; ++j)
	array[i][j] = *defaultVal;
  }
}


template<typename D, typename T>
void
delete3DArray(D***&   array,   // three-dimensional array to delete
              const T dim[3])  // extents of three-dimensional array
{
  for (T i = 0; i < dim[0]; ++i)
    if (array[i])
      delete2DArray<D, T>(array[i], &dim[1]);      
  delete[] array;
}

template<typename D, typename T>
void
allocate3DArray(D***&    array,           // three-dimensional array to create
                const T  dim[3],          // extents of three-dimensional array
                const D* defaultVal = 0)  // optional default value
{
  if (array)
    delete3DArray<D, T>(array, dim);
  array = new D**[dim[0]];
  for (T i = 0; i < dim[0]; ++i)
    allocate2DArray<D, T>(array[i], &dim[1], defaultVal);
}


// pseudo-n-dimensional arrays map the dimensions onto a one-dimensional array
//
// here index is (like in C++ multi-dimensional arrays) row-major, that
// is values for last index i_(n - 1) are consecutive
//
// array[x_0][x_1] ... [x_(n - 1)]
//        |    |            |
//     dim[0]  |            |
//          dim[1] ...   dim[n - 1]
//
// no range checks whatsoever are performed

template<typename D, typename T>
inline
void
allocatePseudoNdimArray(D*&      array,           // pointer to one-dimensional array
			const T* dim,             // extents of n-dimensional array
			const T  nmbDim,          // number of dimensions
			const D* defaultVal = 0)  // optional default value
{
  T nmbElements = dim[0];
  for (T i = 1; i < nmbDim; ++i)
    nmbElements *= dim[i];
  array = new D[nmbElements];
  if (defaultVal)
    for (T i = 0; i < nmbElements; ++i)
      array[i] = *defaultVal;
}


template<typename T>
inline
T
indicesToOffset(const T* indices,  // indices to map to one-dimensional array index
		const T* dim,      // extents of n-dimensional array
		const T  nmbDim)   // number of dimensions
{
  T offset = indices[0];
  for (T i = 1; i < nmbDim; ++i)
    offset = offset * dim[i] + indices[i];
  return offset;
}

template<typename T>
inline
T
indicesToOffset(const std::vector<T>& indices,  // indices to map to one-dimensional array
		const std::vector<T>& dim)      // extents of n-dimensional array
{
  T offset = indices[0];
  for (T i = 1; i < dim.size(); ++i)
    offset = offset * dim[i] + indices[i];
  return offset;
}


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

template<typename T>
inline
void
offsetToIndices(const T               offset,   // one-dimensional array index
		const std::vector<T>& dim,      // extents of n-dimensional array
		std::vector<T>&       indices)  // indices to map onto
{
  T index = offset;
  for (T i = dim.size() - 1; i >= 1; --i) {
    indices[i] = index % dim[i];
    index      = index / dim[i];
  }
  indices[0] = index;
}



//////////////////////////////////////////////////////////////////////////////////
// conversion functions

// converts any class that supports << into a string
template<typename T>
inline
std::string
toString(const T& fromValue)
{
  std::ostringstream to;
  to << fromValue;
  return to.str();

  // much cleaner with BOOST
  // #include <boost/lexical_cast.hpp>
  // double d = 453.23;
  // string str = boost::lexical_cast<string>(d);
}


//////////////////////////////////////////////////////////////////////////////////
// math helper functions and constants

// mathematical constants
const double pi     = 2 * asin((double)1);
const double piHalf = pi / 2;
const double twoPi  = 2 * pi;
const double fourPi = 4 * pi;

const std::complex<double> imag(0,1);

// computes n!
inline
unsigned int
factorial(const unsigned int n)
{
 unsigned int fac = 1;
 for (unsigned int i = 1; i <= n; ++i)
   fac *= i;
 return fac;
}


//////////////////////////////////////////////////////////////////////////////////
// physics helper functions

// computes breakup momentum of 2-body decay
inline
double
breakupMomentum(const double M,   // mass of mother particle
		const double m1,  // mass of daughter particle 1
		const double m2)  // mass of daughter particle 2
{
  if (M < m1 + m2)
    return 0;
  return sqrt((M - m1 - m2) * (M + m1 + m2) * (M - m1 + m2) * (M + m1 - m2)) / (2 * M);
}


// kinematic border in Dalitz plot; PDG 2008 eq. 38.22a, b
// for decay M -> m0 m1 m2
inline
double
dalitzKinematicBorder(const double  mass_2,      // 2-body mass squared on x-axis
		      const double  M,           // 3-body mass
		      const double* m,           // array with the 3 daughter masses
		      const bool    min = true)  // switches between curves for minimum and maximum mass squared on y-axis
{
  if (mass_2 < 0)
    return 0;
  const double  mass   = sqrt(mass_2);
  const double  M_2    = M * M;                                    // 3-body mass squared
  const double  m_2[3] = {m[0] * m[0], m[1] * m[1], m[2] * m[2]};  // daughter masses squared

  // calculate energies of particles 1 and 2 in m01 RF
  const double E1 = (mass_2 - m_2[0] + m_2[1]) / (2 * mass);
  const double E2 = (M_2    - mass_2 - m_2[2]) / (2 * mass);
  const double E1_2  = E1 * E1;
  const double E2_2  = E2 * E2;
  if ((E1_2 < m_2[1]) or (E2_2 < m_2[2]))
    return 0;

  // calculate m12^2
  const double p1     = sqrt(E1_2 - m_2[1]);
  const double p2     = sqrt(E2_2 - m_2[2]);
  const double Esum_2 = (E1 + E2) * (E1 + E2);
  if (min)
    return Esum_2 - (p1 + p2) * (p1 + p2);
  else
    return Esum_2 - (p1 - p2) * (p1 - p2);
}


//////////////////////////////////////////////////////////////////////////////////
// file system helper functions

// expands glob pattern into list of file names
inline
std::vector<std::string>
globFileList(const std::string& globPattern)
{
  std::vector<std::string> fileList;
  glob_t globBuffer;
  glob(globPattern.c_str(), GLOB_NOSORT, NULL, &globBuffer);
  for (unsigned int i = 0; i < globBuffer.gl_pathc; ++i)
    fileList.push_back(globBuffer.gl_pathv[i]);
  globfree(&globBuffer);
  return fileList;
}


#endif  // utilities_h
