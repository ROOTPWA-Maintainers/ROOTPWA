#ifndef MASSDEP_H
#define MASSDEP_H

#include <iostream>
#include <complex>
#include <vector>
#include <particle.h>
#include <matrix.h>

class particle;

class massDep {
	public:
		massDep() { }
		virtual ~massDep() { }
		virtual massDep* create() const = 0;
		virtual massDep* clone() const = 0;

		virtual void print() { std::cout << "massDep"; }
		virtual std::complex<double> val(particle&) = 0;

};


class breitWigner:public massDep {
	public:
		breitWigner() {;}
		~breitWigner() {;}
		breitWigner(const breitWigner&) {;}
		// breitWigner& operator=(const breitWigner) {;}
		breitWigner* create() const {return new breitWigner();}
		breitWigner* clone() const {return new breitWigner(*this);}

		virtual void print() { std::cout << "breitWigner"; }
		std::complex<double> val(particle&);

};

class flat:public massDep {
	public:
		flat() {;}
		~flat() {;}
		flat(const flat&) {;}
		// flat& operator=(const flat) {;}
		flat* create() const {return new flat();}
		flat* clone() const {return new flat(*this);}

		virtual void print() { std::cout << "flat"; }
		std::complex<double> val(particle&);

};


class AMP_M:public massDep {

	private:
		int _Pmax;
		int _Nmax;
		matrix<std::complex<double> > _rho;
		matrix<std::complex<double> > _M;
		matrix<std::complex<double> > _T;
		matrix<std::complex<double> > _f;
		std::vector<matrix<std::complex<double> > > _a;
		std::vector<matrix<std::complex<double> > > _c;
		matrix<double> _sP;

	public:
		int ves_sheet;

		AMP_M();
		~AMP_M() {;}
		AMP_M(const AMP_M&) {;}
		// AMP_M& operator=(const AMP_M) {;}
		virtual massDep* create() const {return new AMP_M();}
		virtual massDep* clone() const {return new AMP_M(*this);}

		virtual void print() { std::cout << "AMP_M"; }
		std::complex<double> val(particle&);


};

class AMP_ves:public AMP_M {
	public:
		AMP_ves():AMP_M() {ves_sheet = 1;}
		~AMP_ves() {;}
		AMP_ves(const AMP_ves&) {;}
		// AMP_ves& operator=(const AMP_ves) {;}
		virtual massDep* create() const {return new AMP_ves();}
		virtual massDep* clone() const {return new AMP_ves(*this);}

		virtual void print() { std::cout << "AMP_ves"; }
		std::complex<double> val(particle&);
};

// class AMP_sigma:public massDep {

// 	private:
// 		int _Pmax;
// 		int _Nmax;
// 		matrix<std::complex<double> > _rho;
// 		matrix<std::complex<double> > _M;
// 		matrix<std::complex<double> > _T;
// 		matrix<std::complex<double> > _f;
// 		std::vector<matrix<std::complex<double> > > _a;
// 		std::vector<matrix<std::complex<double> > > _c;
// 		matrix<double> _sP;

// 	public:
// 		int ves_sheet;

// 		AMP_M();
// 		~AMP_M() {;}
// 		AMP_M(const AMP_M&) {;}
// 		// AMP_M& operator=(const AMP_M) {;}
// 		virtual massDep* create() const {return new AMP_M();}
// 		virtual massDep* clone() const {return new AMP_M(*this);}

// 		virtual void print() { std::cout << "AMP_M"; }
// 		std::complex<double> val(particle&);


// };

#endif


