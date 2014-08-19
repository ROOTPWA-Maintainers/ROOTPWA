	SUBROUTINE EPSK1P(BW, W, N, ALFA)	! K1-prime solution
C	SUBROUTINE EPSAMP(N, S, BW)		! Cut at 1.7 Gev.
C
C This routine returns I=0 S-wave amplitude, multyplied
C by 2-body phase space, for reactions
C  N = 1  Pi+Pi- --> Pi+Pi-
C  N = 2  Pi+Pi- --> K+K-
C Source: K.L.Au et al, Phys.Rev. D35, P 1633.
C         K1-PRIME solution is used.
C This function is based on fit, which parametrization
C implicitly impoused coupled-channel unitarity. In this
C parametrization
C  F = K*(1-i*Rho*K)^-1 == (M-i*Rho)^-1 , M==K^-1
C  F(i,j) - matrix of amplitudes for initial & final states:
C    1 - PI+Pi-
C    2 - K anti-K
C  Rho - diagonal matrix with diag.elem. phase space Rho1,Rho2
C  K - real symmetric matrix and M - it's inverse.
C Input  - S == M(PiPi,KK)^2  -- real
C Output - BW = Rho1*F(1,1) for N=1 (PiPi),
C          BW = Rho2*F(1,2) for N=2 (KK)     -- complex
C  MPi, MK - masses of Pi+- and average kaon mass K+-/K0 (!)
C  Rho(i) = 2*k(i)/M , k - final-state three-momentum.
C Notation: UPPER indexes in source are THE FIRST indexes here.
C Parametrization of K-matrix elements:
C  K(i,j) = (S-S0)/(4*MK^2)*{SUM over P}[fc(p,i)*fc(p,j)/(S(p)-S)/(S(p)-S0)] +
C         + {SUM over N=0,N}[C(n,i,j)*(S/(4*MK^2) -1)^N]
C Note: probably brackets are missing in the paper, and
C       polinomial background must be scaled by (S-S0)/(4*MK^2) !
C Note that there are not used here coupling coefficients for the
C "production process" (c-it's index):
C  Sigma(c,i) = coupl(c)*Rho(i)*|Fprod(c,i)|^2
C  Fprod(c,i) = 1/(S-S0)*{SUM over J}[Alf(c,j)*F(j,i)}
C  Alf(c,j) = {SUM over N}[Alfa(n,i)*(S/4/MK^2)^N] ! c omitted
C
c	IMPLICIT NONE
	PARAMETER  ( S0   = -0.0141,
     +		   S1   =  0.9226 )
	PARAMETER  ( FC11 = -0.2334,
     +		   FC12 =  0.5969 )
	PARAMETER  ( C011 =  0.7871,
     +		   C111 = -0.5610,
     +		   C211 =  1.6987,
     +		   C311 = -2.0451,
     +		   C411 =  0.6361 )
	PARAMETER  ( C012 = -3.3270,
     +		   C112 = -0.4788,
     +		   C212 =  1.1362,
     +		   C312 = -1.0623,
     +		   C412 =  0.6290 )
	PARAMETER  ( C022 = -2.7914,
     +		   C122 =  7.5952,
     +		   C222 =  4.5612,
     +		   C322 = -0.9356,
     +		   C422 =  0 )
	PARAMETER  ( ALFA01 = -0.4700,
     +		   ALFA11 =  0.6593,
     +		   ALFA21 =  0.2036 )
	PARAMETER  ( ALFA02 =  3.542,
     +		   ALFA12 = -3.824,
     +		   ALFA22 =  1.284 )
C--
	PARAMETER  ( WPI = 0.1395675 )	! PI+- MASS
	PARAMETER  ( WKC = 0.493646  )	! K+-  MASS
	PARAMETER  ( WK0 = 0.497671  )	! K0   MASS
	PARAMETER  ( WK = 0.5*(WKC+WK0)  )	! K MEAN MASS
C--
	PARAMETER  ( WCUT = 1.7 )	! MASS CUT
C--
	INTEGER N		! CHANNEL NUMBER
	REAL    W, S		! DIMESON MASS, MASS SQUARED (GEV**2)
	COMPLEX BW		! AMPLITUDE
	REAL ALFA, X		! COUPLING COEFF, NORMALIZED S.
	COMPLEX F(2,2), DETF	! F MATRIX & DETERM.
	REAL    K(2,2), M(2,2)	! K,M MATRIXES
	REAL DETK, DETMX, DETMY
	REAL R(2)		! PHASE SPASE
	REAL RKC,RK0, RK, TEMP	! --"-- FOR K+-,K0,K_MEAN
	REAL RR2, RI2		! REAL&IMAGE PART OF PHASE SPACE FOR K+-
	REAL ALF(2)		! COUPLING COEFFICIENTS (???)
	COMPLEX FPRIM(2)	! INITIAL PROC. AMPLITUDE
C
C CALCULATE PHASE SPACE FOR ALL CHANNELS.
C FOR THE CALCULATIONS OF F "MEAN" K IS USED.
C NOTE THAT BELOW THRESHOLD  i*R(2) ---> -ABS(R(2))
C
	S  = W*W
	BW = (0.,0.)
	IF ( N.EQ.2 .AND. W.GT.WCUT)	     RETURN
	IF ( N.EQ.1 .AND. S.LE.(4.*WPI**2) ) RETURN
	IF ( N.EQ.2 .AND. S.LE.(4.*WKC**2) ) RETURN
	R(1) = SQRT(1. - 4.*WPI**2/S)		! PI+-
	RR2  = 0.
	RI2  = 0.
	RKC  = SQRT(ABS(1. - 4.*WKC**2/S))	! K+-
	RK0  = SQRT(ABS(1. - 4.*WK0**2/S))	! K0
	IF (S.GT.(4.*WKC**2)) THEN
	  RI2  = RKC
	ELSE
	  RR2  = - RKC
	  RKC  = 0.				! PHYSICAL PHASE SPACE
	END IF
	IF (S.GT.(4.*WK0**2) ) THEN
	  RI2  = RI2 + RK0
	ELSE
	  RR2  = RR2 - RK0
	  RK0  = 0.
	END IF
	RR2  = 0.5*RR2			! K MEAN
	RI2  = 0.5*RI2
C
C FIRST COMPUTE SUM OF POLES (ONE POLE FOR K1 SOLUTION).
C
	K(1,1) = FC11*FC11/((S1-S)*(S1-S0))
	K(1,2) = FC11*FC12/((S1-S)*(S1-S0))
	K(2,2) = FC12*FC12/((S1-S)*(S1-S0))
C-- NOTE ! IN INITIAL VERSION ADDITION & SCALING
C          WAS IN OTHER ORDER !
C-- ADD POLINOMIAL BACKGROUND.
	X    = S/(4.*WK**2) -1.
	K(1,1) = K(1,1) + ((((C411*X+C311)*X)+C211)*X+C111)*X+C011
	K(1,2) = K(1,2) + ((((C412*X+C312)*X)+C212)*X+C112)*X+C012
	K(2,2) = K(2,2) + ((((C422*X+C322)*X)+C222)*X+C122)*X+C022
C-- SCALE K-MATRIX
	TEMP = (S-S0)/(4.*WK**2)
	K(1,1) = K(1,1)*TEMP
	K(1,2) = K(1,2)*TEMP
	K(2,2) = K(2,2)*TEMP
C
C COMPUTE M = K^-1 : (A  B)^-1		    ( C  -B)
C		     (B  C)    = 1/(AC-B^2)*(-B   A)
C
	DETK = 1.0/(K(1,1)*K(2,2)-K(1,2)*K(1,2))	! 1./DET(K)
	M(1,1) =   K(2,2)*DETK
	M(1,2) = - K(1,2)*DETK
	M(2,2) =   K(1,1)*DETK
C
C ADD REAL PART OF K PHASE SPACE TO THE M-MATRIX BELOW THRESHOLD.
C
	M(2,2) = M(2,2) + RR2
	R(2)   = RI2
C
C COMPUTE F = (M - i*R)^-1
C
	DETMX = M(1,1)*M(2,2) - M(1,2)*M(1,2) - R(1)*R(2)
	DETMY = -M(1,1)*R(2) - M(2,2)*R(1)
	DETF  = CMPLX(DETMX/(DETMX**2+DETMY**2),-DETMY/(DETMX**2+DETMY**2))
	F(1,1) = CMPLX(M(2,2),-R(2)) * DETF
	F(1,2) = CMPLX(-M(1,2),  0.) * DETF
	F(2,2) = CMPLX(M(1,1),-R(1)) * DETF
C
C F MATRIX READY, COMPUTE COUPLING COEFFICIENTS (?)
C
c	X = X + 1.
c	ALF(1) = ((ALFA21*X)+ALFA11)*X+ALFA01
c	ALF(2) = ((ALFA22*X)+ALFA12)*X+ALFA02
C
C  COMPUTE AMPLITUDES FOR INITIAL CHANNEL.
C  MULTIPLY F BY PHASE SPACE.
C
c	FPRIM(1) = 1./(S-S0)*(ALF(1)*F(1,1)+ALF(2)*F(1,2))
c	FPRIM(2) = 1./(S-S0)*(ALF(1)*F(1,2)+ALF(2)*F(2,2))
C--
C	R(2) = RKC			! K+-
	IF (N .EQ. 1) BW = (1.-ALFA)*F(1,1) + ALFA*F(2,1)
	IF (N .EQ. 2) BW = (1.-ALFA)*F(1,2) + ALFA*F(2,2)
	RETURN
	END
