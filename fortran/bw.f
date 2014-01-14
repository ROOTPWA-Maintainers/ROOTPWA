C*
C* PWA library, Breigth-Wigner functions.
C* Note that generally all functions from this file should return
C* invariant matrix elements for use on Dalitz plot in PWA,
C* while well-known BW functions are "partial waves" or "Argand plot
C* amplitudes" and should be multiplied by M/P to get matrix elements.
C*
      COMPLEX FUNCTION BWL(L,W,W0,G0,P,P0,RVAL)
C
C 31-JUN-1995. New version of BWL function.
C BWL(M) is Lorentz-invariant (not Argand plot) amplitude.
C BWL = M0*G0*BL(q,q0)/(S0-S - i*M0*G(M)),
C G(M)= G0*(Rho/Rho0)*[BL(q,q0)]^2, BL(q,q0) = FL(q*R)/FL(q0*R),
C FL(Z) are Blatt-Weisskopf formfactors (See BLWA(Z,L)),
C R is range constant about 1 fermi = 1./0.1973 Gev/c = 5.068 Gev^-1
C This routine returns GENERAL L-wave Breight-Wigner function
C for L = 0,1,2,3,4 including Blatt-Weisskopf form factors.
C Source: S.U.Chung, PWA in K-matrix formalizm,
C F.v.Hippel and C.Quigg, Phys.Rev.5,624(1972)
C
      IMPLICIT NONE
      INTEGER L			! Spin of resonance.
      REAL W, W0		! Dimeson mass current/in resonance, Gev.
      REAL P, P0		! Decay momentum --"--
      REAL G, G0		! Width --"--
      REAL RVAL,R, R0		! Barrier factor given/used/default, Gev^-1
      PARAMETER (R0 = 5.068)	! Default barrier factor = 1./1 Fermi
      REAL BLWA, BL, BWR, BWI, C, DENOM

      BWL = (0.,0.)
      IF (P.LE.0.) RETURN
      R = RVAL
      IF (R.LE.0.) R = R0
C
C Compute variable width.
C Use constant width if P0.LE.0 (res. under threshold)
C
      IF (P0.GT.0.) GOTO 10
	BL = 1.
	G  = G0
	GOTO 20
   10 CONTINUE
	BL = BLWA(P*R,L)/BLWA(P0*R,L)
	G  = G0*BL*(P/W)*(W0/P0)
   20 CONTINUE
C--
      BWR = W0**2 - W**2	! REAL PART
      BWI = W0*G		! IMAGE PART
      C   = W0*G0*SQRT(BL)	! numerator
      DENOM = C / (BWR**2+BWI**2)
      BWL = CMPLX(BWR*DENOM, BWI*DENOM)	! BWL = C/(BWR-i*BWI)
      RETURN
      END

      REAL FUNCTION BLWA (P, L)
C
C Returns SQUARE of the Blatt-Weisskopf barrier factor.
C P = Q*R, dimensionless decay momentum. R is about 1 Fermi.
C L - orbital momentum.
C
      IMPLICIT NONE
      REAL P, Z
      INTEGER L

      Z = P*P
      GOTO (10,20,30,40), L
	BLWA = 1.			! L = 0 or unknown
	GOTO 99
   10 CONTINUE				! L = 1
	BLWA = Z/(Z+1.)			! 2.*
	GOTO 99
   20 CONTINUE				! L = 2
	BLWA = Z*Z/((Z+3.)*Z+9.)	! 13.*
	GOTO 99
   30 CONTINUE					! L = 3
	BLWA = Z*Z*Z/(((Z+6.)*Z+45.)*Z+225.)	! 277.*
	GOTO 99
   40 CONTINUE					! L = 4, 12746.*
	BLWA = (Z*Z)**2/((((Z+10.)*Z+135.)*Z+1575.)*Z+11025.)
   99 RETURN
      END

      COMPLEX FUNCTION BWS(W,W0,G0)
C
C This routine returns SIMPLE BREIGHT-WIGNER function with CONSTANT width.
C BWS = W0*G0/(W0*W0-W*W - i*W0*G0) . Source: XXXXXX
C
      REAL W, W0, G0		! MASS CURRENT,MASS AND WIDTH OF RESONANCE.

      BWR = W0*W0 - W*W
      BWI = W0*G0
      DENOM = BWI / (BWR**2+BWI**2)
      BWS   = CMPLX(BWR*DENOM, BWI*DENOM)
      END

      COMPLEX FUNCTION BWRHOO(W, W0, G0)
C
C This routine returns RHO --> PI-PI amplitude, just as BWL function
C for L=1, but G=G0*W0/W*(P/P0)*(2*P**2)/(P**2+P0**2), i.e.
C Blatt-Weisskopf factor is replaced by (2*P**2)/(P**2+P0**2)
C Source: D.Bisello et al, Phys.Rev. D39, P 701.
C
      IMPLICIT NONE
      REAL W,  G, S		! Current mass, width (Gev)
      REAL W0, G0		! Resonance mass,width (Gev)
      REAL P,  P0		! Decay momentum for current/rho mass
      REAL WPI, BL, BWR, BWI, C, DENOM
      PARAMETER (WPI = 0.1395675)	! PI+- MASS
C
      S  = W*W
      BWRHOO = (0.,0.)
      IF (S.LE.(4.*WPI**2) ) RETURN
      P  = SQRT(S - 4.*WPI**2)		! CURR MOMENTUM *2
      P0 = SQRT(W0**2 - 4.*WPI**2)	! RHO --"--
      BL = 2.*P**2/(P0**2+P**2)
      G  = G0*(P/W)*(W0/P0)*BL
C--
      BWR = W0**2 - S
      BWI = W0*G
      C   = W0*G0*SQRT(BL)
      DENOM = C / (BWR**2+BWI**2)
      BWRHOO = CMPLX(BWR*DENOM, BWI*DENOM)	! BW = C/(BWR-i*BWI)
      RETURN
      END

      COMPLEX FUNCTION BWF0 (W, W0, G0)
C
C RETURNS S-WAVE BW FUNCTION FOR F0 -> PI+PI-
C This should remains the same for all BWL modifications,
C it is used in EPS* Au, Morgan,Pennington "cut" parametrizations.
C Probably this isn't correct S-wave BW form !
C
      IMPLICIT NONE
      REAL W, W0, G0, P, P0, R0, WPI, RHOPI, WM, G, BWR, BWI, C, DENOM
      PARAMETER  (R0 = 5.07, WPI = 0.1395675)	! PI+- MASS
      RHOPI(WM) = 0.5*SQRT(WM**2 - 4.*WPI**2)

      BWF0 = (0.,0.)
      IF (W .LE.(2.*WPI)) GOTO 990
      IF (W0.LE.(2.*WPI)) GOTO 990
      IF (G0.LE.0.) GOTO 990
      P  = RHOPI (W)
      P0 = RHOPI (W0)
C     BW  = BWL (0, W, W0, G0, P, P0, 0.)
C--
      G = G0*P/P0		! Should be phase_spase/ph_sp_0
      BWR = W0**2 - W**2
      BWI = W0*G
      C   = BWI * W/P
      DENOM = C / (BWR**2+BWI**2)
      BWF0 = CMPLX(BWR*DENOM, BWI*DENOM)	! BW = C/(BWR-i*BWI)
  990 RETURN
      END

!
! Check correctness !!! For PWA the function should return
! Argand plot (or "partial wave" or Breit-Wigner) amplitude divided
! by 2-particle phase space k/m. This is just invariant matrix element
! for use on Dalitz plot.
!
!
	COMPLEX FUNCTION BWA0 (W, X1, X2, X3)
C
C This routine returns the corrected ETA-PI Breit-Wigner S-wave
C amplitude. Source: S. Flatte, Phys. Lett. 63B (1976) 224-227.
C 27.03.95. Add function parameters: WR - mass, G0 - full width,
C GG - G(K)/G(ETA). Additional parameters can be given
C via COMMON/ZDATA/.
C
C
C     IMPLICIT NONE
      Parameter (WET = 0.5488)		! ETA  MASS
      Parameter (WPI = 0.1395675)	! PI+- MASS
      Parameter (WK  = 0.495)		! K    MASS average
      REAL X1, X2, X3, RDATA
      COMMON /ZDATA/ RDATA(100)
C
C  S. Flatte best fit parameters
C
      Parameter (WR = 0.969)		! Resonace mass
      Parameter (G0 = 0.082)		! "Total width"
      Parameter (GG = 2.0)		! G(K)/G(Eta)
      Parameter (CE = 0.5466394)	! See text below
      Parameter (CK = 0.5466394)	! ---"--- (D)
C--
      Real    W, WW		! Di-meson mass, mass squared (Gev**2)
      Complex	GK		! -I * (K-K channel width)
C--
      IF ( Inited.ne.0 ) go to 1000
	SHD2K	= (2.0*WK)**2		! 2K threshold	
	SHDEP	= (WET+WPI)**2		! ETA-PI threshold	
	EMP	= (WET-WPI)**2		! --"-- pseudo-threshold
	WR2	= WR**2
!	CK = CE*RDATA(13)		!(D)
!	Type *,'BWA0:	CK=',CK		!(D)
	Inited	= 1
c
c	  Computation of the CE,CK parameters according
c	  to the  S. Flatte definitions
c
c	CE = WR*G0/(SQRT((WR2-SHDEP)*(WR2-EMP))/(2.*WR))
c	CK = 0.5*GG*CE
c	TYPE *,CE,CK
c	STOP

 1000 Continue
      BWA0 = (0.,0.)
      WW  = W*W
      IF ( WW .lt. SHDEP ) Return

! [EDT]
      P = SQRT((WW-SHDEP)*(WW-EMP))/(2*W)
      IF ( WW.ge.SHD2K ) GK = CMPLX(0.0,-CK*SQRT(WW - SHD2K))
      IF ( WW.lt.SHD2K ) GK = CMPLX(CK*SQRT(SHD2K - WW),0.0)
      BWA0 = (1.,0.) / ( CMPLX(WR2-WW,-CE*P) + GK )
      Return
      End
