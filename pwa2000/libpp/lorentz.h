#ifndef __LORENTZ_H_
#define __LORENTZ_H_
	

#include <cmath>
#include <Vec.h>
#include <matrix.h>


class rotation : public matrix<double> {

public:

  rotation()
    : matrix<double>(3, 3)
  { }
  rotation(double, double, double);
  ~rotation() { }
		
  rotation set(double, double, double);
  rotation set(const threeVec&);
		
  friend threeVec operator *= (threeVec&, const rotation&);
  //		friend particle operator *= (particle&, const rotation&);

};


class lorentzTransform : public matrix<double> {

  friend class fourVec;
	
private:

  double _beta;
  double _gamma;

public:

  lorentzTransform()
    : matrix<double>(4, 4),
      _beta(0.),
      _gamma(1.)
  { }
  lorentzTransform(const matrix<double>& m)
    : matrix<double>(m),
      _beta(0.),
      _gamma(1.)
  { }	
  lorentzTransform(double, double, double);
  lorentzTransform(const threeVec& beta);
  lorentzTransform(const fourVec&);
  lorentzTransform(const rotation&);
  ~lorentzTransform() { }

  lorentzTransform set(double, double, double); // rotation
  lorentzTransform set(const threeVec& beta); // boost
  lorentzTransform set(const fourVec&); //boost to rest frame
  lorentzTransform set(const rotation&);

  double beta() { return _beta;}

  //		lorentzTransform operator *= (const lorentzTransform&);
  friend fourVec operator *= (fourVec&, const lorentzTransform&);

};


#endif
