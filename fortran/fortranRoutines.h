//
// example collection of definitions for external Fortran routines
//


#ifndef FORTRANROUTINES_HH
#define FORTRANROUTINES_HH


#include <complex>


#ifdef __cplusplus
extern "C" {


	//////////////////////////////////////////////////////////////////////////////
	// functions defined in bw_example.f
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


}
#endif  // __cplusplus


#endif  // FORTRANROUTINES_HH
