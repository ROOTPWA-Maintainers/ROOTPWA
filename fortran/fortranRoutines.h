//
// collection of definitions for external Fortran routines
//


#ifndef FORTRANROUTINES_HH
#define FORTRANROUTINES_HH


#include <complex>


#ifdef __cplusplus
extern "C" {

	// various Breit-Wigner parametrizations taken from Dima's program
	// Blatt-Weisskopf barrier factor
	void blwa_(float* f,   // square of Blatt-Weisskopf barrier factor
	           float* p,   // q * r dimensionless breakup momentum
	           int*   l);  // relative orbital angular momentum of decay daughters

	// Breit-Wigner with mass-dependent width
	void bwl_(std::complex<float>* bw,     // calculated amplitude
	          int*                 l,      // relative orbital angular momentum of decay daughters
	          float*               w,      // di-meson mass [GeV/c^2]
	          float*               w0,     // nominal resonance mass [GeV/c^2]
	          float*               g0,     // nominal resonance width [GeV/c^2]
	          float*               p,      // di-meson breakup momentum [GeV/c]
	          float*               p0,     // nominal di-meson breakup momentum [GeV/c]
	          float*               rval);  // interaction radius [fm^{-1}]

	// s-wave Breit-Wigner with constant width
	void bws_(std::complex<float>* bw,   // calculated amplitude
	          float*               w,    // di-meson mass [GeV/c^2]
	          float*               w0,   // nominal resonance mass [GeV/c^2]
	          float*               g0);  // nominal resonance width [GeV/c^2]

	// rho(770) Breit-Wigner
	void bwrhoo_(std::complex<float>* bw,   // calculated amplitude
	             float*               w,    // di-meson mass [GeV/c^2]
	             float*               w0,   // nominal resonance mass [GeV/c^2]
	             float*               g0);  // nominal resonance width [GeV/c^2]

	// f_0(980) Breit-Wigner
	void bwf0_(std::complex<float>* bw,   // calculated amplitude
	           float*               w,    // di-meson mass [GeV/c^2]
	           float*               w0,   // nominal resonance mass [GeV/c^2]
	           float*               g0);  // nominal resonance width [GeV/c^2]

	// a_0(980) Flatte
	void bwf0_(std::complex<float>* bw,   // calculated amplitude
	           float*               w,    // di-meson mass [GeV/c^2]
	           float*               x1,   // not used
	           float*               x2,   // not used
	           float*               x3);  // not used


	// various parametrizations of the pi pi s-wave taken from Dima's program
	// author I. A. Kachaev (IHEP, Protvino)
	// based on PRD 35, 1633 (1987)
	//
	// K1-solution
	void epsk1_(std::complex<float>* bw,     // calculated amplitude
	            float*               w,      // di-meson mass [GeV/c^2]
	            int*                 n,      // switch between pi pi (1) and K K (2) final state
	            float*               alfa);  // mixing coefficient; usually 0

	// K1'-solution
	void epsk1p_(std::complex<float>* bw,     // calculated amplitude
	             float*               w,      // di-meson mass [GeV/c^2]
	             int*                 n,      // switch between pi pi (1) and K K (2) final state
	             float*               alfa);  // mixing coefficient; usually 0

	// M-solution
	void epsm_(std::complex<float>* bw,     // calculated amplitude
	           float*               w,      // di-meson mass [GeV/c^2]
	           int*                 n,      // switch between pi pi (1) and K K (2) final state
	           float*               alfa);  // mixing coefficient; usually 0

	// M-solution with f_0(980) removed
	void epsmx_(std::complex<float>* bw,     // calculated amplitude
	            float*               w,      // di-meson mass [GeV/c^2]
	            int*                 n,      // switch between pi pi (1) and K K (2) final state
	            float*               alfa);  // mixing coefficient; usually 0

}
#endif  // __cplusplus


#endif  // FORTRANROUTINES_HH
