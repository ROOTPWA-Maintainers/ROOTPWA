#line 321 "event.nw"
#ifndef EVENT_H
#define EVENT_H

#include <iostream>
#include <string>
#include <list>
#include <pputil.h>
#include <Vec.h>
#include <lorentz.h>
#include <particle.h>
	
	extern particleDataTable PDGtable;
	
	class particle;
	class event {
		protected:
			std::list<particle> _initial;
			std::list<particle> _final;
			particle* _beam;
			particle* _target;
			int _ioversion;
		public:
			
#line 38 "event.nw"
	
#line 52 "event.nw"
	event();

#line 58 "event.nw"
	event(const event&);

#line 64 "event.nw"
	~event();

#line 70 "event.nw"
	event& operator=(const event& e);

#line 39 "event.nw"
	
#line 82 "event.nw"
	event& addfinal(const particle&);	

#line 88 "event.nw"
	event& addinitial(const particle& p);	

#line 94 "event.nw"
	event& erase();

#line 100 "event.nw"
	event& beam(const particle&);	

#line 106 "event.nw"
	event& target(const particle&);	


#line 40 "event.nw"
	
#line 121 "event.nw"
	int OK(double epsilon) const;	

#line 127 "event.nw"
	particle beam() const;	

#line 133 "event.nw"
	particle target() const;	

#line 141 "event.nw"
	fourVec getPartPFinal(std::string name,int charge,int index,int debug=0) const;

#line 147 "event.nw"
	fourVec getPartPInitial(std::string name,int charge,int index) const;

#line 153 "event.nw"
	int f_charge() const;

	  double f_mass() const;


#line 159 "event.nw"
	std::list<particle> f_mesons() const;

#line 165 "event.nw"
	std::list<particle> f_baryons() const;

#line 171 "event.nw"
	std::list<particle> f_particles() const;

#line 177 "event.nw"
	particle f_particle(const std::string& name, int charge, int index) const;

#line 183 "event.nw"
	int i_charge() const;

#line 189 "event.nw"
	std::list<particle> i_mesons() const;

#line 195 "event.nw"
	void set_f_mesons(const std::list<particle>& l);

#line 201 "event.nw"
	std::list<particle> i_baryons() const;

#line 207 "event.nw"
	std::list<particle> i_particles() const;

#line 213 "event.nw"
	particle i_particle(const std::string& name, int charge, int index) const;

#line 219 "event.nw"
	threeVec mesonPlane() const;

#line 225 "event.nw"
	threeVec baryonPlane() const;

#line 41 "event.nw"
	
#line 236 "event.nw"
	void print() const;

#line 241 "event.nw"
	friend std::istream& operator>>(std::istream& is, event& e);	

#line 246 "event.nw"
	friend std::ostream& operator<<(std::ostream& os, event& e);	

#line 262 "event.nw"
	std::istream& read1(std::istream& is);	

#line 268 "event.nw"
	std::ostream& write1(std::ostream& os);	

#line 284 "event.nw"
	std::istream& read2(std::istream& is);

#line 290 "event.nw"
	std::ostream& write2(std::ostream& os);

#line 296 "event.nw"
	event& setIOVersion(int ver);

#line 42 "event.nw"
	
#line 312 "event.nw"
	friend event operator*(const lorentzTransform&,const event&);
	 
#line 344 "event.nw"
	};
	
#endif

