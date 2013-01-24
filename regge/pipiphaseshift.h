//-----------------------------------------------------------
//
// Description:
//      Pion-Pion phase shift parameterizations and elasticities
//
//      Literature:
//      Pelaez, Yndurain              Phys.Rev.D71 (2005) 074016
//      Kaminski, Pelaez, Yndurain    Phys.Rev.D74 (2006) 014001
//      Kaminski, Pelaez, Yndurain    Phys.Rev.D77 (2008) 054015
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef PIPIPHASESHIFT_HH
#define PIPIPHASESHIFT_HH

// Base Class Headers ----------------

// Collaborating Class Headers -------
#include <complex>
#include <boost/numeric/ublas/vector.hpp>
 #include <boost/numeric/ublas/vector_proxy.hpp>
 #include <boost/numeric/ublas/matrix.hpp>
 #include <boost/numeric/ublas/triangular.hpp>
 #include <boost/numeric/ublas/lu.hpp>
 #include <boost/numeric/ublas/io.hpp>

 namespace ublas = boost::numeric::ublas;
 typedef ublas::matrix<double> rmatrix;

// Collaborating Class Declarations --

static const double mpic=0.13957018;
static const double mpic2=mpic*mpic;
static const double mKc=0.493677;
static const double mKc2=mKc*mKc;

// two-particle threshold function (breakup momentum for identical particles)
std::complex<double> kThr(std::complex<double> s, double m2){
  return 0.5*sqrt(s-4.0*m2);
}


/// interface class for Phaseshift functions
class absPhaseshift {
public:
 // Constructors/Destructors ---------
 absPhaseshift(double s0=1):_s0(s0){;}
  virtual ~absPhaseshift(){;}

 // Operators
 // Accessors -----------------------
 virtual double cotDelta(std::complex<double> s) const = 0;
 virtual double eta(std::complex<double> s) const      = 0;
 // Modifiers -----------------------
 void setS0(double s0){_s0=s0;}

 // Operations ----------------------

 /// Conformal variable
 std::complex<double> w(std::complex<double> s) const;

 protected:
 double _s0; /// scale variable for conformal transform

};

///////////////////////////////////////////////////////////
///// PHASE SHIFT PARAMETERIZATIONS ala Pelaez et al //////
///////////////////////////////////////////////////////////


/// S0 wave (I=0 l=0) ///////////////////////////////////
class s0wave : public absPhaseshift {
 public:
  s0wave(double z02=0.01947983515,
	 double B0=4.3, double B1=-26.7, double B2=-14.1);
  virtual ~s0wave(){;}

  virtual double cotDelta(std::complex<double> s) const;
  virtual double eta(std::complex<double> s) const ;
  
  void setParameters(double z0, double B0, double B1, double B2)
  {_z02=z0*z0;_B0=B0;_B1=B1;_B2=B2;}

 private:
  double _z02;  /// adler zero is at 1/2 z02;
  double _B0,_B1,_B2; /// expansion parameters

  /// K-Matrix parameters
  rmatrix _gamma; // 
  double _alpha[2];
  double _beta[2];
  double _mu,_invmu;
  double _M12,_M22;
  
 };








#endif
