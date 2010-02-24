#ifndef __PARTICLE_H_
#define __PARTICLE_H_


#include <list>
#include <string>
#include <complex>

#include "Vec.h"
#include "lorentz.h"
#include "particleData.h"

	
class massDep;
class event;
class decay;


class particle : public particleData {
	
public:
			
  particle();
  particle(const particle& p);
  particle(const particleData& data, const int charge);
  virtual ~particle();

  particle& operator =  (const particle&         p);
  int       operator == (const particle&         p) { return Mass() == p.Mass(); }
  int       operator != (const particle&         p) { return Mass() != p.Mass(); }
  int       operator <  (const particle&         p) { return Mass() <  p.Mass(); }
  particle& operator *= (const lorentzTransform& L);
  friend particle operator * (const lorentzTransform& L,
			      const particle&         p);

  particle& setCharge (const int       charge);
  particle& set4P     (const fourVec&  p4);
  particle& set3P     (const threeVec& p3);
  particle& Index     (const int       i);
  particle& setDecay  (const decay&    d);
  particle& setMassDep(massDep*        md);

  decay*   Decay()  const { return _decay;           }
  int      Stable() const { return (_decay) ? 0 : 1; }
  fourVec  get4P()  const { return _p;               }
  threeVec get3P()  const { return _p.V();           }
  int      Index()  const { return _index;           }
  int      Charge() const { return _charge;          }
  fourVec* get4P(particle* p,
		 const int debug = 0);

  std::list<int>& helicities ();
  particle&       addHelicity(const int lam);

  int is(const std::string& name) const;

  fourVec setupFrames(const int debug = 0);

  double q () const;
  double q0() const;

  std::complex<double> breitWigner() const;

  std::complex<double> decayAmp(const int lambda,
				const int debug = 0);

  std::string sprint     (const std::string& space = " ") const;
  void        print      () const;
  void        printFrames() const;

  void debug(const int d = 1) { _particle_debug = d; }

private:
			
  static int     _particle_debug;
  int            _lambda;
  int            _charge;
  fourVec        _p;
  decay*         _decay;
  int            _index;
  std::list<int> _helicities;
  int            _inRestFrame;
  massDep*       _massDep;

};

	
class decay {
	
public:

  decay();
  decay(const decay&);
  virtual ~decay();
		
  decay& addChild  (const particle& p);
  decay& setL      (const int       l);
  decay& setS      (const int       s);
  decay& calculateS();

  int L() const { return _l; }
  int S() const { return _s; }
  fourVec* get4P(particle* part,
		 const int debug = 0);

  decay& operator =  (const decay&            d);
  decay& operator *= (const lorentzTransform& L);

  fourVec fill(const event& e,
	       const int    debug = 0);

  decay& setupFrames(const lorentzTransform& T,
		     const int               debug = 0);

  std::complex<double> expt_amp(const double b,
				const double t,
				const int    debug = 0);

  std::complex<double> amp(const int j,
			   const int m,
			   const int debug = 0);
		
  void print      () const;
  void printFrames() const;

  void debug(const int d = 1) { _decay_debug = d; }

  std::list<particle> _children;

private:

  void _init(const std::list<particle>& children,
	     const int                  l,
	     const int                  s,
	     const double               mass);

  static int          _decay_debug;
  std::list<particle> _childrenInFrames;
  int                 _l;
  int                 _s;
  double              _mass;
	
};


#define _PARTICLE_H
#endif
