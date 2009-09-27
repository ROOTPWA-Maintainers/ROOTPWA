//
// package that contains various small utility functions
//


#ifndef utilities_h
#define utilities_h


#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>

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


template <typename T> class maxPrecisionValue__;

// output stream manipulator that prints a value with its maximum precision
template <typename T>
inline
maxPrecisionValue__<T>
maxPrecision(T value)
{ return maxPrecisionValue__<T>(value); }

// output stream manipulator that prints a value with its maximum precision
// in addition manipulator reserves space so that values will align
template <typename T>
inline
maxPrecisionValue__<T>
maxPrecisionAlign(T value)
{ return maxPrecisionValue__<T>(value, maxPrecisionValue__<T>::ALIGN); }

// general helper class that encapsulates a value of type T
template <typename T>
class maxPrecisionValue__ {
public:
  enum modeEnum { PLAIN,
		  ALIGN };
  maxPrecisionValue__(const T        value,
		      const modeEnum mode = PLAIN)
    : value_(value),
      mode_ (mode)
  { }
  std::ostream& print(std::ostream&  out) const
  {
    const int nmbDigits = std::numeric_limits<double>::digits10 + 1;
    std::ostringstream s;
    s.precision(nmbDigits);
    s.setf(std::ios_base::scientific, std::ios_base::floatfield);
    s << value_;
    switch (mode_) {
    case ALIGN:
      return out << std::setw(nmbDigits + 7) << s.str();  // make space for sign, dot, and exponent
    case PLAIN: default:
      return out << s.str();
    }
  }
private:
  T        value_;
  modeEnum mode_;
};

template <typename T>
inline
std::ostream& operator << (std::ostream&                 out,
			   const maxPrecisionValue__<T>& value)
{ return value.print(out); }


// simple stream operators for some common ROOT classes
inline
std::ostream&
operator << (std::ostream&   o,
             const TVector3& vec)
{
  o << "(" << vec.X() << ", " << vec.Y() << ", " << vec.Z() << ")";
  return o;
}

inline
std::ostream&
operator << (std::ostream&         o,
             const TLorentzVector& vec)
{
  o << "(" << vec.X() << ", " << vec.Y() << ", " << vec.Z() << "; " << vec.T() << ")";
  return o;
}

inline
std::ostream&
operator << (std::ostream&   o,
             const TComplex& c)
{
  o << "(" << c.Re() << ", " << c.Im() << ")";
  return o;
}

template <typename T>
std::ostream&
operator << (std::ostream&      o,
             const TMatrixT<T>& A)
{
  o << "(";
  for (int row = 0; row < A.GetNrows(); ++row) {
    o << "(";
    for (int col = 0; col < A.GetNcols(); ++col) {
      o << A[row][col];
      if (col < A.GetNcols() - 1)
        o << ", ";
    }
    if (row < A.GetNrows() - 1)
      o << "), ";
    else
      o << ")";
  }
  o << ")";
  return o;
}


#endif  // utilities_h
