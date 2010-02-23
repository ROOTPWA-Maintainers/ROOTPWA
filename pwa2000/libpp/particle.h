#ifndef __PARTICLE_H_
#define __PARTICLE_H_

#include <list>
#include <string>
#include <Vec.h>
#include <matrix.h>
#include <lorentz.h>
#include <particleData.h>
#include <pputil.h>
#include <massDep.h>

	
	class decay;
	class massDep;

	class particle:public particleData {
		private:
			
		static int _particle_debug;
		int _lambda;
		int _charge;
		fourVec _p;
		decay *_decay;
		int _index;
		std::list<int> _helicities;
		int _inRestFrame;
		massDep* _massDep;
	
		public:
			
		particle();
		particle(const particleData&,int);
		particle(const particle&);
		~particle();
		particle& operator=(const particle&);
	
		particle& setCharge(int);
		particle& set4P(const fourVec&);
		particle& set3P(const threeVec&);
		particle& Index(int);
		particle& setDecay(const decay&);
		particle& setMassDep(massDep* md) {_massDep = md; return *this;}
		friend particle operator*(const lorentzTransform&, const particle&);
		particle& operator*=(const lorentzTransform&);

		decay* Decay() const;
		int Stable() const;
		int operator==(const particle&);
		int operator!=(const particle&);
		int operator<(const particle&);

		fourVec get4P() const;
		fourVec* get4P(particle* p, int debug = 0) const;
		threeVec get3P() const;
		int Index() const;
		int Charge() const;
		std::list<int>& helicities();
		particle& addHelicity(int lam);
		int is(std::string) const;
		fourVec setupFrames(int debug = 0);
		std::string sprint(std::string space = " ");
		double q() const;
		double q0() const;
		std::complex<double> breitWigner() const;

		std::complex<double> decayAmp(int, int debug = 0);

		void print() const;
		void printFrames() const;

		void debug(int d = 1) {
			_particle_debug = d;
		}

	};

	
class event;
class decay {
	
	private:
		// list<particle> _children;
		static int _decay_debug;
		std::list<particle> _childrenInFrames;
		int _l;
		int _s;
		double _mass;
		void _init(const std::list<particle>&,int,int,double);
	
	
 	public:
		std::list<particle> _children;
		decay();
		decay(const decay&);
		~decay();
		
		decay& addChild(const particle&);
		decay& setL(int);
		decay& setS(int);
		int L() const;
		int S() const;

		decay& calculateS();
		fourVec fill(const event&, int debug=0);
		fourVec* get4P(particle* part, int debug=0);

		decay& setupFrames(lorentzTransform T, int debug=0);
		std::complex<double> amp(int,int,int debug=0) const;
		std::complex<double> expt_amp(double b, double t, int debug=0) const;
		
		decay& operator*=(const lorentzTransform& L);
		decay& operator=(const decay&);
		
		void print() const;
		void printFrames() const;

		void debug(int d = 1) {
			_decay_debug = d;
		}

};


#define _PARTICLE_H
#endif
