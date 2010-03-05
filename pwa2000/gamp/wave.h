#ifndef WAVE_H
#define WAVE_H


#include <complex>
#include <string>

#include "lorentz.h"
#include "particle.h"


#define getsign(x) (((x)>0)?"+":"-")


class wave : public particle {

public:

  wave();
  wave(const wave& wv);
  virtual ~wave();

  wave& operator =  (const wave&             wv);
  wave& operator *= (const lorentzTransform& L);

  wave setM    (const int    m);
  wave setSlope(const double b);
  wave setT    (const double t);
  wave channel (const char*  ch);

  fourVec     getBeam  () const { return _beam;    } 
  fourVec     getTarget() const { return _target;  } 
  int         M        () const { return _m;       }
  std::string channel  () const { return _channel; }

  wave fill(const event& e,
	    const int    debug = 0);

  wave& setupFrames(const int debug = 0);

  std::complex<double> decayAmp(const int debug = 0);

  std::string sprint     (const std::string& space = " ") const;
  void        print      () const;
  void        printFrames() const;

private:

  int         _m;
  int         _epsilon;
  fourVec     _beam;
  fourVec     _target;
  std::string _channel;
  double      _b;
  double      _t;

};


#endif
