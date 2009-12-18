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

#include "TCMatrix.h"


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


// simple stream operators for some SLT classes
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


template <typename T>
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


inline
std::ostream&
operator << (std::ostream&   out,
             const TCMatrix& A)
{
  for (int row = 0; row < A.nrows(); ++row) {
    out << "row " << row << " = (";
    for (int col = 0; col < A.ncols(); ++col) {
      out << A(row, col);
      if (col < A.ncols() - 1)
        out << ", ";
    }
    if (row < A.nrows() - 1)
      out << "), " << std::endl;
    else
      out << ")";
  }
  return out;
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
  if ((E1_2 < m_2[1]) || (E2_2 < m_2[2]))
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


#endif  // utilities_h
