//
// collection of definitions for external Fortran routines
//


#ifndef FORTRANROUTINES_HH
#define FORTRANROUTINES_HH


#include <complex>


#ifdef __cplusplus
extern "C" {

  // various parameteriaztions of the pi pi s-wave taken from Dima's program
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
