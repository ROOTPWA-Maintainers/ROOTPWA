	SUBROUTINE EPSK1 (BW, W, N, ALFA)	! K1 solution
C	SUBROUTINE EPSAMP(N, S, BW)
C
C *** Modified version. Change EPS_K1P.FOR, EPS_K3.FOR, EPS_M.FOR
C *** accordingly to _use_ them.
C
C This routine returns I=0 S-wave amplitude for reactions
C  N = 1  Pi+Pi- --> Pi+Pi-
C  N = 2  Pi+Pi- --> K+K-
C Source: K.L.Au et al, Phys.Rev. D35, P 1633.
C         K1 solution is used.
C This function is based on fit, which parametrization
C implicitly impoused coupled-channel unitarity. In this
C parametrization
C  F = K*(1-i*Rho*K)^-1 == (M-i*Rho)^-1 , M==K^-1
C  F(i,j) - (symm.)matrix of amplitudes for initial & final states:
C    1 - PI+Pi-,    2 - K anti-K(neutral & charged)
C  Rho - diagonal matrix with diag.elem. phase space Rho1,Rho2
C  K - real symmetric matrix and M - it's inverse.
C Input  - W == M(PiPi,KK) -- real, N -- integer, ALFA -- real
C Output - BW = (1.-ALFA)*F(1,1) + ALFA*F(2,1)  for N=1 (PiPi),
C          BW = ALFA*F(1,2) + (1.-ALFA)*F(2,2)  for N=2 (KK)  -- complex
C  ALFA = 0. means pure PiPi initial state.
C  MPi, MK - masses of Pi+- and average kaon mass K+-/K0 (!)
C  Rho(i) = 2*k(i)/M , k - final-state three-momentum.
C Notation: UPPER indexes in source are THE FIRST indexes here.
C Parametrization of K-matrix elements:
C  K(i,j) = (S-S0)/(4*MK^2)*{SUM over P}[fc(p,i)*fc(p,j)/(S(p)-S)/(S(p)-S0)] +
C         + {SUM over N=0,N}[C(n,i,j)*(S/(4*MK^2) -1)^N]
C Note: probably brackets are missing in the paper, and
C       polinomial background must be scaled by (S-S0)/(4*MK^2) !
C Note that fitted coupling coefficients for the
C "production process" (c-it's index) are not used here:
C  Sigma(c,i) = coupl(c)*Rho(i)*|Fprod(c,i)|^2
C  Fprod(c,i) = 1/(S-S0)*{SUM over J}[Alf(c,j)*F(j,i)}
C  Alf(c,j) = {SUM over N}[Alfa(n,j)*(S/4/MK^2)^N] ! c omitted
C
	PARAMETER  ( S0 = -0.0110,
     +		   S1 =  0.9247 )
	PARAMETER  ( FC11 = -0.2242,
     +		   FC12 =  0.5829  )
	PARAMETER  ( C011 =  0.7347,
     +		   C111 = -0.5266,
     +		   C211 =  2.6151,
     +		   C311 = -1.7747,
     +		   C411 =  0.8031 )
	PARAMETER  ( C012 = -3.2762,
     +		   C112 = -0.6662,
     +		   C212 =  0.8778,
     +		   C312 = -2.1190,
     +		   C412 =  0.2319 )
	PARAMETER  ( C022 = -2.6785,
     +		   C122 =  7.9951,
     +		   C222 =  5.5763,
     +		   C322 = -1.4956,
     +		   C422 =  0 )
	PARAMETER  ( ALFA01 = -0.4012,
     +		   ALFA11 =  0.5468,
     +		   ALFA21 =  0.2440 )
	PARAMETER  ( ALFA02 =  3.273,
     +		   ALFA12 = -3.483,
     +		   ALFA22 =  1.183 )
C--
	PARAMETER  ( WPI = 0.1395675 )	! PI+- MASS
	PARAMETER  ( WKC = 0.493646  )	! K+-  MASS
	PARAMETER  ( WK0 = 0.497671  )	! K0   MASS
	PARAMETER  ( WK = 0.5*(WKC+WK0)  )	! K MEAN MASS
C--
	COMPLEX BW		! AMPLITUDE
	REAL    W, S		! DIMESON MASS, MASS SQUARED (GEV**2)
	INTEGER N		! CHANNEL NUMBER
	REAL ALFA		! COUPLING (OR MIXING) COEFFICIENT
	REAL X			! NORMALIZED S
	COMPLEX F(2,2), DETF	! F MATRIX & DETERM.
	REAL    K(2,2), M(2,2)	! K,M MATRIXES
	REAL DETK, DETMX, DETMY
	REAL R(2)		! PHASE SPASE
	REAL RKC,RK0, RK, TEMP	! --"-- FOR K+-,K0,K_MEAN
	REAL RR2, RI2		! REAL&IMAGE PART OF PHASE SPACE FOR K+-
C	REAL ALF(2)		! COUPLING COEFFICIENTS (???)
C	COMPLEX FPRIM(2)	! INITIAL PROC. AMPLITUDE
C
C CALCULATE PHASE SPACE FOR ALL CHANNELS.
C FOR THE CALCULATIONS OF F "MEAN" K IS USED.
C NOTE THAT BELOW THRESHOLD  i*R(2) ---> -ABS(R(2))
C
C	PARAMETER (W_MAX = 1.7)				! BAD FIT AFTER IT !

	BW = (0.,0.)
	S  = W*W
C	IF ( N.EQ.1 .AND. S.GT.( W_MAX**2) ) RETURN	! HIGH MASS CUT
	IF ( N.EQ.1 .AND. S.LE.(4.*WPI**2) ) RETURN
	IF ( N.EQ.2 .AND. S.LE.(4.*WKC**2) ) RETURN
	IF (S .EQ. S1) RETURN				! DIVIDE BY ZERO.
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
C NOTE ! IN INITIAL VERSION ADDITION & SCALING WAS IN OTHER ORDER !
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
	DETK = K(1,1)*K(2,2)-K(1,2)*K(1,2)	! DET(K)
	IF (DETK .EQ. 0.) DETK = 1.E-6
	DETK = 1.0/ DETK
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
	F(2,1) = F(1,2)
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
