#ifndef EVENT_H
#define EVENT_H


#include <iostream>
#include <string>
#include <list>

#include "Vec.h"
#include "lorentz.h"
#include "particle.h"
	

class event {

public:
			
  event();
  event(const event&);
  virtual ~event();

  event& operator = (const event& e);
  friend event operator * (const lorentzTransform&,
			   const event&);

  event& addfinal  (const particle&);
  event& addinitial(const particle& p);
  event& erase();

  particle beam  () const { return *_beam;   }
  particle target() const { return *_target; }
  event&   beam  (const particle&);
  event&   target(const particle&);

  int OK(const double epsilon = 1e-6) const;	

  fourVec getPartPFinal(const std::string& name,
			const int          charge,
			const int          index,
			const int          debug = 0) const;
  fourVec getPartPInitial(const std::string& name,
			  const int          charge,
			  const int          index) const;

  int                 f_charge   () const;
  double              f_mass     () const;
  std::list<particle> f_mesons   () const;
  std::list<particle> f_baryons  () const;
  std::list<particle> f_particles() const { return _final; }
  particle            f_particle (const std::string& name,
				  const int          charge,
				  const int          index) const;
  void set_f_mesons(const std::list<particle>& l);

  int                 i_charge   () const;
  std::list<particle> i_mesons   () const;
  std::list<particle> i_baryons  () const;
  std::list<particle> i_particles() const { return _initial; }
  particle            i_particle (const std::string& name,
				  const int          charge,
				  const int          index) const;

  threeVec mesonPlane () const;
  threeVec baryonPlane() const;

  void print() const;
  friend std::ostream& operator << (std::ostream& os, const event& e);
  std::ostream& write1(std::ostream& os) const;
  std::ostream& write2(std::ostream& os) const;

  friend std::istream& operator >> (std::istream& is, event& e);
  std::istream& read1(std::istream& is);
  std::istream& read2(std::istream& is);

  event& setIOVersion(const int ver);

protected:

  std::list<particle> _initial;
  std::list<particle> _final;
  particle*           _beam;
  particle*           _target;
  int                 _ioversion;

};

	
#endif  // EVENT_H
