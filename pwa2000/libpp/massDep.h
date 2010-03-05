#ifndef MASSDEP_H
#define MASSDEP_H


#include <iostream>
#include <complex>
#include <vector>

#include "particle.h"
#include "matrix.h"


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
		virtual ~breitWigner() {;}
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
		virtual ~flat() {;}
		flat(const flat&) {;}
		// flat& operator=(const flat) {;}
		flat* create() const {return new flat();}
		flat* clone() const {return new flat(*this);}

		virtual void print() { std::cout << "flat"; }
		std::complex<double> val(particle&);

};

/** @brief AMP parameterization of pipi s-wave
 *  
 *  We have introduced a small modification by setting the off-diagonal 
 *  elements of the M-matrix to zero.
 */
class AMP_M:public massDep {

	protected:
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
		virtual ~AMP_M() {;}
		AMP_M(const AMP_M&) {;}
		// AMP_M& operator=(const AMP_M) {;}
		virtual massDep* create() const {return new AMP_M();}
		virtual massDep* clone() const {return new AMP_M(*this);}

		virtual void print() { std::cout << "AMP_M"; }
		std::complex<double> val(particle&);


};

/** @brief old VES parameterization
 *
 *  Brute force subtraction of the f0(980)
 */ 
class AMP_ves:public AMP_M {
	public:
		AMP_ves():AMP_M() {ves_sheet = 1;}
		virtual ~AMP_ves() {;}
		AMP_ves(const AMP_ves&) {;}
		// AMP_ves& operator=(const AMP_ves) {;}
		virtual massDep* create() const {return new AMP_ves();}
		virtual massDep* clone() const {return new AMP_ves(*this);}

		virtual void print() { std::cout << "AMP_ves"; }
		std::complex<double> val(particle&);
};


/** @brief Kachaev's version of the AMP parameterization
 * 
 * From the original fortran code:
 * Source: K.L.Au et al, Phys.Rev. D35, P 1633. M solution.
 * 04-Mar-2003 See eps_k1.for for description.
 * Here matrix M=K^{-1} is parametrized with one pole.
 * Misprint in the article (other than in K1--K3 solutions)
 * was corrected.
 *
 * 14-Mar-2003 Nice amplitude for pi-pi S-wave without f0(975).
 * It is smooth and nicely tends to zero after approx 1.5 GeV.
 * f0(975) pole excluded; coupling to KK zeroed; set C411=C422=0.
 * The largest effect from C411, zeroing of C422 looks insignificant.
 */ 

class AMP_kach:public AMP_M {
	public:
		AMP_kach();
		virtual ~AMP_kach() {;}
		AMP_kach(const AMP_kach&) {;}
		// AMP_ves& operator=(const AMP_ves) {;}
		virtual massDep* create() const {return new AMP_kach();}
		virtual massDep* clone() const {return new AMP_kach(*this);}

		virtual void print() { std::cout << "AMP_kach"; }
  //std::complex<double> val(particle&);
};


#endif
