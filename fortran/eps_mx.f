      SUBROUTINE EPSMX (BW, W, N, ALFA)  ! M solution - no KK coupling.
C
C Source: K.L.Au et al, Phys.Rev. D35, P 1633. M solution.
C 04-Mar-2003 See eps_k1.for for description.
C Here matrix M=K^{-1} is parametrized with one pole.
C Misprint in the article (other than in K1--K3 solutions)
C was corrected.
C
C 14-Mar-2003 Nice amplitude for pi-pi S-wave without f0(975).
C It is smooth and nicely tends to zero after approx 1.5 GeV.
C f0(975) pole excluded; coupling to KK zeroed; set C411=C422=0.
C The largest effect from C411, zeroing of C422 looks insignificant.
C
      IMPLICIT NONE
      REAL S0, S1, FC11,FC12, A11,A12,A22, C011,C111,C211,C311,C411,
     &     C012,C112,C212,C312,C412, C022,C122,C222,C322,C422,
     &     WPI,WKC,WK0,WK
      PARAMETER (S0   = -0.0074, S1   =  0.9828)
      PARAMETER (FC11 =  0.1968, FC12 = -0.0154)
      PARAMETER (A11  =  0.1131, A12  =  0.0150, A22  = -0.3216)
      PARAMETER (C011 =  0.0337, C111 = -0.3185, C211 = -0.0942,
     &           C311 = -0.5927, C411 =  0.)
C    &           C311 = -0.5927, C411 =  0.1957)
      PARAMETER (C012 = -0.2826, C112 =  0.0918, C212 =  0.1669,
     &           C312 = -0.2082, C412 = -0.1386)
      PARAMETER (C022 =  0.3010, C122 = -0.5140, C222 =  0.1176,
     &           C322 =  0.5204, C422 =  0.)
C    &           C322 =  0.5204, C422 = -0.3977)
C     REAL ALFA01,ALFA11,ALFA21, ALFA02,ALFA12,ALFA22,
C     PARAMETER (ALFA01 =  0.1393, ALFA11 = -0.02775, ALFA21 =  0.3952)
C     PARAMETER (ALFA02 =  3.241,  ALFA12 = -3.432,   ALFA22 =  1.141 )
C--
      PARAMETER (WPI = 0.1395675, WKC = 0.493646)     ! PI+- MASS, K+- MASS
      PARAMETER (WK0 = 0.497671,  WK = 0.5*(WKC+WK0)) ! K0 MASS, K MEAN MASS
C--
      COMPLEX BW              ! AMPLITUDE
      REAL W, S               ! DIMESON MASS, MASS SQUARED (GEV**2)
      INTEGER N               ! CHANNEL NUMBER
      REAL ALFA               ! COUPLING (OR MIXING) COEFFICIENT
      REAL X                  ! NORMALIZED S
      COMPLEX F(2,2), DETF    ! F MATRIX & DETERM.
      REAL    M(2,2)          ! M MATRIX
      REAL DETMX, DETMY
      REAL R(2)               ! PHASE SPASE
      REAL RKC, RK0           ! --"-- FOR K+-,K0
      REAL RR2, RI2           ! REAL&IMAGE PART OF PHASE SPACE FOR K+-
C     REAL ALF(2)             ! COUPLING COEFFICIENTS (???)
C     COMPLEX FPRIM(2)        ! INITIAL PROC. AMPLITUDE
C
C CALCULATE PHASE SPACE FOR ALL CHANNELS.
C FOR THE CALCULATIONS OF F "MEAN" K IS USED.
C NOTE THAT BELOW THRESHOLD  i*R(2) ---> -ABS(R(2))
C
      BW = (0.,0.)
      S  = W*W
      IF ( N.EQ.1 .AND. S.LE.(4.*WPI**2) ) RETURN
      IF ( N.EQ.2 .AND. S.LE.(4.*WKC**2) ) RETURN
      IF (S .EQ. S1) RETURN
      R(1) = SQRT(1. - 4.*WPI**2/S)           ! PI+-
      RR2  = 0.
      RI2  = 0. 
      RKC  = SQRT(ABS(1. - 4.*WKC**2/S))      ! K+-
      RK0  = SQRT(ABS(1. - 4.*WK0**2/S))      ! K0
      IF (S.GT.(4.*WKC**2)) THEN
        RI2  = RKC
      ELSE
        RR2  = - RKC
        RKC  = 0.                       ! PHYSICAL PHASE SPACE
      END IF
      IF (S.GT.(4.*WK0**2) ) THEN
        RI2  = RI2 + RK0
      ELSE
        RR2  = RR2 - RK0
        RK0  = 0.
      END IF
      RR2  = 0.5*RR2              ! K MEAN
      RI2  = 0.5*RI2
C
C First compute sum of poles (one pole for M solution).
C Misprint in the sign was corrected !!!
C
C     M(1,1) = A11/(S-S0) + FC11*FC11/(S1-S)
C     M(1,2) = A12/(S-S0) + FC11*FC12/(S1-S)
C     M(2,2) = A22/(S-S0) + FC12*FC12/(S1-S)
C     M(1,1) = A11/(S-S0) - FC11*FC11/(S1-S)
C     M(1,2) = A12/(S-S0) - FC11*FC12/(S1-S)
C     M(2,2) = A22/(S-S0) - FC12*FC12/(S1-S)
      M(1,1) = A11/(S-S0)
      M(1,2) = 0.
      M(2,2) = A22/(S-S0)
C
C-- Add polinomial background.
C
      X    = S/(4.*WK**2) -1.
      M(1,1) = M(1,1) + ((((C411*X+C311)*X)+C211)*X+C111)*X+C011
C     M(1,2) = M(1,2) + ((((C412*X+C312)*X)+C212)*X+C112)*X+C012
      M(2,2) = M(2,2) + ((((C422*X+C322)*X)+C222)*X+C122)*X+C022
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
      DETF  = CMPLX(DETMX/(DETMX**2+DETMY**2),
     &  -DETMY/(DETMX**2+DETMY**2))
      F(1,1) = CMPLX(M(2,2),-R(2)) * DETF
      F(1,2) = CMPLX(-M(1,2),  0.) * DETF
      F(2,2) = CMPLX(M(1,1),-R(1)) * DETF
      F(2,1) = F(1,2)
C
C F MATRIX READY, COMPUTE COUPLING COEFFICIENTS (?)
C
C     X = X + 1.
C     ALF(1) = ((ALFA21*X)+ALFA11)*X+ALFA01
C     ALF(2) = ((ALFA22*X)+ALFA12)*X+ALFA02
C
C COMPUTE AMPLITUDES FOR INITIAL CHANNEL.
C MULTIPLY F BY PHASE SPACE.
C
C     FPRIM(1) = 1./(S-S0)*(ALF(1)*F(1,1)+ALF(2)*F(1,2))
C     FPRIM(2) = 1./(S-S0)*(ALF(1)*F(1,2)+ALF(2)*F(2,2))
C--
C     R(2) = RKC                  ! K+-
      IF (N .EQ. 1) BW = (1.-ALFA)*F(1,1) + ALFA*F(2,1)
      IF (N .EQ. 2) BW = (1.-ALFA)*F(1,2) + ALFA*F(2,2)
      RETURN
      END
