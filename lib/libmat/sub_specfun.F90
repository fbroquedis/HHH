!     BEGIN OF SUBROUTINE *** RJBESL ***                                
      SUBROUTINE RJBESL(XINPUT,ALPHA,NB,B,NCALC) 
!                                                                       
!     ##################################################################
!                                                                       
!     *** T. WIEDER, 05.03.1998 ***                                     
!                                                                       
!     SOURCE FOR SUBROUTINE *** RJBESL ***:                             
!                                                                       
!     SUBPROGRAM LIBRARY SPECFUN                                        
!                                                                       
!     AVAILABILITIY: PUBLIC DOMAIN                                      
!                                                                       
!     DEVELOPER: W.J. CODY AND L. STOLTZ, APPLIED MATHEMATICS DIVISION, 
!     ARGONNE NATIONAL LABORATORY, ARGONNE, IL 60439.                   
!                                                                       
!     DISTRIBUTOR: NETLIB                                               
!                                                                       
!     ##################################################################
!                                                                       
!---------------------------------------------------------------------  
! This routine calculates Bessel functions J sub(N+ALPHA) (X)           
!   for non-negative argument X, and non-negative order N+ALPHA.        
!                                                                       
!                                                                       
!  Explanation of variables in the calling sequence.                    
!                                                                       
!   XIN = X *** T. WIEDER, 15.04.98 ***                                 
!   X     - working precision non-negative real argument for which      
!           J's are to be calculated.                                   
!   ALPHA - working precision fractional part of order for which        
!           J's or exponentially scaled J'r (J*exp(X)) are              
!           to be calculated.  0 <= ALPHA < 1.0.                        
!   NB  - integer number of functions to be calculated, NB > 0.         
!           The first function calculated is of order ALPHA, and the    
!           last is of order (NB - 1 + ALPHA).                          
!   B  - working precision output vector of length NB.  If RJBESL       
!           terminates normally (NCALC=NB), the vector B contains the   
!           functions J/ALPHA/(X) through J/NB-1+ALPHA/(X), or the      
!           corresponding exponentially scaled functions.               
!   NCALC - integer output variable indicating possible errors.         
!           Before using the vector B, the user should check that       
!           NCALC=NB, i.e., all orders have been calculated to          
!           the desired accuracy.  See Error Returns below.             
!                                                                       
!                                                                       
!*******************************************************************    
!*******************************************************************    
!                                                                       
!  Explanation of machine-dependent constants                           
!                                                                       
!   it     = Number of bits in the mantissa of a working precision      
!            variable                                                   
!   NSIG   = Decimal significance desired.  Should be set to            
!            INT(LOG10(2)*it+1).  Setting NSIG lower will result        
!            in decreased accuracy while setting NSIG higher will       
!            increase CPU time without increasing accuracy.  The        
!            truncation error is limited to a relative error of         
!            T=.5*10**(-NSIG).                                          
!   ENTEN  = 10.0 ** K, where K is the largest integer such that        
!            ENTEN is machine-representable in working precision        
!   ENSIG  = 10.0 ** NSIG                                               
!   RTNSIG = 10.0 ** (-K) for the smallest integer K such that          
!            K .GE. NSIG/4                                              
!   ENMTEN = Smallest ABS(X) such that X/4 does not underflow           
!   XLARGE = Upper limit on the magnitude of X.  If ABS(X)=N,           
!            then at least N iterations of the backward recursion       
!            will be executed.  The value of 10.0 ** 4 is used on       
!            every machine.                                             
!                                                                       
!                                                                       
!     Approximate values for some important machines are:               
!                                                                       
!                                                                       
!                            it    NSIG    ENTEN       ENSIG            
!                                                                       
!   CRAY-1        (S.P.)     48     15    1.0E+2465   1.0E+15           
!   Cyber 180/855                                                       
!     under NOS   (S.P.)     48     15    1.0E+322    1.0E+15           
!   IEEE (IBM/XT,                                                       
!     SUN, etc.)  (S.P.)     24      8    1.0E+38     1.0E+8            
!   IEEE (IBM/XT,                                                       
!     SUN, etc.)  (D.P.)     53     16    1.0D+308    1.0D+16           
!   IBM 3033      (D.P.)     14      5    1.0D+75     1.0D+5            
!   VAX           (S.P.)     24      8    1.0E+38     1.0E+8            
!   VAX D-Format  (D.P.)     56     17    1.0D+38     1.0D+17           
!   VAX G-Format  (D.P.)     53     16    1.0D+307    1.0D+16           
!                                                                       
!                                                                       
!                           RTNSIG      ENMTEN      XLARGE              
!                                                                       
!   CRAY-1        (S.P.)    1.0E-4    1.84E-2466   1.0E+4               
!   Cyber 180/855                                                       
!     under NOS   (S.P.)    1.0E-4    1.25E-293    1.0E+4               
!   IEEE (IBM/XT,                                                       
!     SUN, etc.)  (S.P.)    1.0E-2    4.70E-38     1.0E+4               
!   IEEE (IBM/XT,                                                       
!     SUN, etc.)  (D.P.)    1.0E-4    8.90D-308    1.0D+4               
!   IBM 3033      (D.P.)    1.0E-2    2.16D-78     1.0D+4               
!   VAX           (S.P.)    1.0E-2    1.17E-38     1.0E+4               
!   VAX D-Format  (D.P.)    1.0E-5    1.17D-38     1.0D+4               
!   VAX G-Format  (D.P.)    1.0E-4    2.22D-308    1.0D+4               
!                                                                       
!*******************************************************************    
!*******************************************************************    
!                                                                       
!  Error returns                                                        
!                                                                       
!    In case of an error,  NCALC .NE. NB, and not all J's are           
!    calculated to the desired accuracy.                                
!                                                                       
!    NCALC .LT. 0:  An argument is out of range. For example,           
!       NBES .LE. 0, ALPHA .LT. 0 or .GT. 1, or X is too large.         
!       In this case, B(1) is set to zero, the remainder of the         
!       B-vector is not calculated, and NCALC is set to                 
!       MIN(NB,0)-1 so that NCALC .NE. NB.                              
!                                                                       
!    NB .GT. NCALC .GT. 0: Not all requested function values could      
!       be calculated accurately.  This usually occurs because NB is    
!       much larger than ABS(X).  In this case, B(N) is calculated      
!       to the desired accuracy for N .LE. NCALC, but precision         
!       is lost for NCALC .LT. N .LE. NB.  If B(N) does not vanish      
!       for N .GT. NCALC (because it is too small to be represented),   
!       and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K        
!       significant figures of B(N) can be trusted.                     
!                                                                       
!                                                                       
!  Intrinsic and other functions required are:                          
!                                                                       
!     ABS, AINT, COS, DBLE, GAMMA (or DGAMMA), INT, MAX, MIN,           
!                                                                       
!     REAL, SIN, SQRT                                                   
!                                                                       
!                                                                       
!  Acknowledgement                                                      
!                                                                       
!   This program is based on a program written by David J. Sookne       
!   (2) that computes values of the Bessel functions J or I of real     
!   argument and integer order.  Modifications include the restriction  
!   of the computation to the J Bessel function of non-negative real    
!   argument, the extension of the computation to arbitrary positive    
!   order, and the elimination of most underflow.                       
!                                                                       
!  References: "A Note on Backward Recurrence Algorithms," Olver,       
!               F. W. J., and Sookne, D. J., Math. Comp. 26, 1972,      
!               pp 941-947.                                             
!                                                                       
!              "Bessel Functions of Real Argument and Integer Order,"   
!               Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp       
!               125-132.                                                
!                                                                       
!  Latest modification: March 19, 1990                                  
!                                                                       
!  Author: W. J. Cody                                                   
!          Applied Mathematics Division                                 
!          Argonne National Laboratory                                  
!          Argonne, IL  60439                                           
!                                                                       
!---------------------------------------------------------------------  
!S    REAL               GAMMA,                                         
!     .. Scalar Arguments ..                                            
      DOUBLE PRECISION ALPHA,XINPUT 
      INTEGER NB,NCALC 
!     ..                                                                
!     .. Array Arguments ..                                             
      DOUBLE PRECISION B(NB) 
!     ..                                                                
!     .. Local Scalars ..                                               
      DOUBLE PRECISION ALP2EM,ALPEM,CAPP,CAPQ,EIGHTH,EM,EN,ENMTEN,ENSIG,&
     &                 ENTEN,FOUR,GNU,HALF,HALFX,ONE,ONE30,P,PI2,PLAST, &
     &                 POLD,PSAVE,PSAVEL,RTNSIG,S,SOMME,T,T1,TEMPA,TEMPB, &
     &                 TEMPC,TEST,THREE,THREE5,TOVER,TWO,TWOFIV,TWOPI1, &
     &                 TWOPI2,VCOS,VSIN,X,XC,XIN,XK,XLARGE,XM,Z,ZERO    
      INTEGER I,J,K,L,M,MAGX,N,NBMX,NEND,NSTART 
!     ..                                                                
!     .. Local Arrays ..                                                
      DOUBLE PRECISION FACT(25) 
!     ..                                                                
!     .. External Functions ..                                          
      DOUBLE PRECISION DGAMMA 
      EXTERNAL DGAMMA 
!     ..                                                                
!     .. Intrinsic Functions ..                                         
      INTRINSIC ABS,AINT,COS,DBLE,INT,MAX,MIN,SIN,SQRT 
!     ..                                                                
!     .. Statement Functions ..                                         
      DOUBLE PRECISION CONV,FUNC 
!     ..                                                                
!     .. Data statements ..                                             
!---------------------------------------------------------------------  
!  Mathematical constants                                               
!                                                                       
!   PI2    - 2 / PI                                                     
!   TWOPI1 - first few significant digits of 2 * PI                     
!   TWOPI2 - (2*PI - TWOPI) to working precision, i.e.,                 
!            TWOPI1 + TWOPI2 = 2 * PI to extra precision.               
!---------------------------------------------------------------------  
!    DATA PI2, TWOPI1, TWOPI2 /0.636619772367581343075535E0,6.28125E0,  
!   1 1.935307179586476925286767E-3/                                    
!    DATA ZERO, EIGHTH, HALF, ONE /0.0E0,0.125E0,0.5E0,1.0E0/           
!    DATA TWO, THREE, FOUR, TWOFIV /2.0E0,3.0E0,4.0E0,25.0E0/           
!    DATA ONE30, THREE5 /130.0E0,35.0E0/                                
!---------------------------------------------------------------------  
!  Machine-dependent parameters                                         
!---------------------------------------------------------------------  
!S    DATA ENTEN, ENSIG, RTNSIG /1.0E38,1.0E8,1.0E-2/                   
!S    DATA ENMTEN, XLARGE /1.2E-37,1.0E4/                               
!---------------------------------------------------------------------  
!     Factorial(N)                                                      
!---------------------------------------------------------------------  
!S    DATA FACT /1.0E0,1.0E0,2.0E0,6.0E0,24.0E0,1.2E2,7.2E2,5.04E3,     
!S   1 4.032E4,3.6288E5,3.6288E6,3.99168E7,4.790016E8,6.2270208E9,      
!S   2 8.71782912E10,1.307674368E12,2.0922789888E13,3.55687428096E14,   
!S   3 6.402373705728E15,1.21645100408832E17,2.43290200817664E18,       
!S   4 5.109094217170944E19,1.12400072777760768E21,                     
!S   5 2.585201673888497664E22,6.2044840173323943936E23/                
      DATA PI2,TWOPI1,TWOPI2/0.636619772367581343075535D0,6.28125D0,    &
     &     1.935307179586476925286767D-3/                               
      DATA ZERO,EIGHTH,HALF,ONE/0.0D0,0.125D0,0.5D0,1.0D0/ 
      DATA TWO,THREE,FOUR,TWOFIV/2.0D0,3.0D0,4.0D0,25.0D0/ 
      DATA ONE30,THREE5/130.0D0,35.0D0/ 
      DATA ENTEN,ENSIG,RTNSIG/1.0D38,1.0D17,1.0D-4/ 
      DATA ENMTEN,XLARGE/1.2D-37,1.0D4/ 
      DATA FACT/1.0D0,1.0D0,2.0D0,6.0D0,24.0D0,1.2D2,7.2D2,5.04D3,      &
     &     4.032D4,3.6288D5,3.6288D6,3.99168D7,4.790016D8,6.2270208D9,  &
     &     8.71782912D10,1.307674368D12,2.0922789888D13,                &
     &     3.55687428096D14,6.402373705728D15,1.21645100408832D17,      &
     &     2.43290200817664D18,5.109094217170944D19,                    &
     &     1.12400072777760768D21,2.585201673888497664D22,              &
     &     6.2044840173323943936D23/                                    
!     ..                                                                
!     .. Statement Function definitions ..                              
!                                                                       
!---------------------------------------------------------------------  
! Statement functions for conversion and the gamma function.            
!---------------------------------------------------------------------  
!S    CONV(I) = REAL(I)                                                 
!S    FUNC(X) = GAMMA(X)                                                
      CONV(I) = DBLE(I) 
      FUNC(X) = DGAMMA(X) 
!     ..                                                                
!                                                                       
!     *** T. WIEDER, 15.04.1998 ***                                     
!     *** AMENDMENT SUGGESTED BY DR. T.R. HOPKINS, 15.04.1998 ***       
      X = XINPUT 

!                                                                       
!---------------------------------------------------------------------  
! Check for out of range arguments.                                     
!---------------------------------------------------------------------  
      MAGX = INT(X) 
      IF ((NB.GT.0) .AND. (X.GE.ZERO) .AND. (X.LE.XLARGE) .AND.         &
     &    (ALPHA.GE.ZERO) .AND. (ALPHA.LT.ONE)) THEN                    
!---------------------------------------------------------------------  
! Initialize result array to zero.                                      
!---------------------------------------------------------------------  
          NCALC = NB 
          DO 10 I = 1,NB 
              B(I) = ZERO 
   10     CONTINUE 
!---------------------------------------------------------------------  
! Branch to use 2-term ascending series for small X and asymptotic      
! form for large X when NB is not too large.                            
!---------------------------------------------------------------------  
          IF (X.LT.RTNSIG) THEN 
!---------------------------------------------------------------------  
! Two-term ascending series for small X.                                
!---------------------------------------------------------------------  
              TEMPA = ONE 
              ALPEM = ONE + ALPHA 
              HALFX = ZERO 
              IF (X.GT.ENMTEN) HALFX = HALF*X 
              IF (ALPHA.NE.ZERO) TEMPA = HALFX**ALPHA/                  &
     &                                   (ALPHA*FUNC(ALPHA))            
              TEMPB = ZERO 
              IF ((X+ONE).GT.ONE) TEMPB = -HALFX*HALFX 
              B(1) = TEMPA + TEMPA*TEMPB/ALPEM 
              IF ((X.NE.ZERO) .AND. (B(1).EQ.ZERO)) NCALC = 0 
              IF (NB.NE.1) THEN 
                  IF (X.LE.ZERO) THEN 
                      DO 20 N = 2,NB 
                          B(N) = ZERO 
   20                 CONTINUE 
                                                                        
                  ELSE 
!---------------------------------------------------------------------  
! Calculate higher order functions.                                     
!---------------------------------------------------------------------  
                      TEMPC = HALFX 
                      TOVER = (ENMTEN+ENMTEN)/X 
                      IF (TEMPB.NE.ZERO) TOVER = ENMTEN/TEMPB 
                      DO 30 N = 2,NB 
                          TEMPA = TEMPA/ALPEM 
                          ALPEM = ALPEM + ONE 
                          TEMPA = TEMPA*TEMPC 
                          IF (TEMPA.LE.TOVER*ALPEM) TEMPA = ZERO 
                          B(N) = TEMPA + TEMPA*TEMPB/ALPEM 
                          IF ((B(N).EQ.ZERO) .AND.                      &
     &                        (NCALC.GT.N)) NCALC = N - 1               
   30                 CONTINUE 
                  END IF 
                                                                        
              END IF 
                                                                        
          ELSE IF ((X.GT.TWOFIV) .AND. (NB.LE.MAGX+1)) THEN 
!---------------------------------------------------------------------  
! Asymptotic series for X .GT. 21.0.                                    
!---------------------------------------------------------------------  
              XC = SQRT(PI2/X) 
              XIN = (EIGHTH/X)**2 
              M = 11 
              IF (X.GE.THREE5) M = 8 
              IF (X.GE.ONE30) M = 4 
              XM = FOUR*CONV(M) 
!---------------------------------------------------------------------  
! Argument reduction for SIN and COS routines.                          
!---------------------------------------------------------------------  
              T = AINT(X/ (TWOPI1+TWOPI2)+HALF) 
              Z = ((X-T*TWOPI1)-T*TWOPI2) - (ALPHA+HALF)/PI2 
              VSIN = SIN(Z) 
              VCOS = COS(Z) 
              GNU = ALPHA + ALPHA 
              DO 50 I = 1,2 
                  S = ((XM-ONE)-GNU)* ((XM-ONE)+GNU)*XIN*HALF 
                  T = (GNU- (XM-THREE))* (GNU+ (XM-THREE)) 
                  CAPP = S*T/FACT(2*M+1) 
                  T1 = (GNU- (XM+ONE))* (GNU+ (XM+ONE)) 
                  CAPQ = S*T1/FACT(2*M+2) 
                  XK = XM 
                  K = M + M 
                  T1 = T 
                  DO 40 J = 2,M 
                      XK = XK - FOUR 
                      S = ((XK-ONE)-GNU)* ((XK-ONE)+GNU) 
                      T = (GNU- (XK-THREE))* (GNU+ (XK-THREE)) 
                      CAPP = (CAPP+ONE/FACT(K-1))*S*T*XIN 
                      CAPQ = (CAPQ+ONE/FACT(K))*S*T1*XIN 
                      K = K - 2 
                      T1 = T 
   40             CONTINUE 
                  CAPP = CAPP + ONE 
                  CAPQ = (CAPQ+ONE)* (GNU*GNU-ONE)* (EIGHTH/X) 
                  B(I) = XC* (CAPP*VCOS-CAPQ*VSIN) 
                  IF (NB.EQ.1) GO TO 180 
                  T = VSIN 
                  VSIN = -VCOS 
                  VCOS = T 
                  GNU = GNU + TWO 
   50         CONTINUE 
!---------------------------------------------------------------------  
! If  NB .GT. 2, compute J(X,ORDER+I)  I = 2, NB-1                      
!---------------------------------------------------------------------  
              IF (NB.GT.2) THEN 
                  GNU = ALPHA + ALPHA + TWO 
                  DO 60 J = 3,NB 
                      B(J) = GNU*B(J-1)/X - B(J-2) 
                      GNU = GNU + TWO 
   60             CONTINUE 
              END IF 
!---------------------------------------------------------------------  
! Use recurrence to generate results.  First initialize the             
! calculation of P*S.                                                   
!---------------------------------------------------------------------  
          ELSE 
              NBMX = NB - MAGX 
              N = MAGX + 1 
              EN = CONV(N+N) + (ALPHA+ALPHA) 
              PLAST = ONE 
              P = EN/X 
!---------------------------------------------------------------------  
! Calculate general significance test.                                  
!---------------------------------------------------------------------  
              TEST = ENSIG + ENSIG 
              IF (NBMX.GE.3) THEN 
!---------------------------------------------------------------------  
! Calculate P*S until N = NB-1.  Check for possible overflow.           
!---------------------------------------------------------------------  
                  TOVER = ENTEN/ENSIG 
                  NSTART = MAGX + 2 
                  NEND = NB - 1 
                  EN = CONV(NSTART+NSTART) - TWO + (ALPHA+ALPHA) 
                  DO 90 K = NSTART,NEND 
                      N = K 
                      EN = EN + TWO 
                      POLD = PLAST 
                      PLAST = P 
                      P = EN*PLAST/X - POLD 
                      IF (P.GT.TOVER) THEN 
!---------------------------------------------------------------------  
! To avoid overflow, divide P*S by TOVER.  Calculate P*S until          
! ABS(P) .GT. 1.                                                        
!---------------------------------------------------------------------  
                          TOVER = ENTEN 
                          P = P/TOVER 
                          PLAST = PLAST/TOVER 
                          PSAVE = P 
                          PSAVEL = PLAST 
                          NSTART = N + 1 
   70                     N = N + 1 
                          EN = EN + TWO 
                          POLD = PLAST 
                          PLAST = P 
                          P = EN*PLAST/X - POLD 
                          IF (P.LE.ONE) GO TO 70 
                          TEMPB = EN/X 
!---------------------------------------------------------------------  
! Calculate backward test and find NCALC, the highest N such that       
! the test is passed.                                                   
!---------------------------------------------------------------------  
                          TEST = POLD*PLAST* (HALF-HALF/ (TEMPB*TEMPB)) 
                          TEST = TEST/ENSIG 
                          P = PLAST*TOVER 
                          N = N - 1 
                          EN = EN - TWO 
                          NEND = MIN(NB,N) 
                          DO 80 L = NSTART,NEND 
                              POLD = PSAVEL 
                              PSAVEL = PSAVE 
                              PSAVE = EN*PSAVEL/X - POLD 
                              IF (PSAVE*PSAVEL.GT.TEST) THEN 
                                  NCALC = L - 1 
                                  GO TO 110 
                                                                        
                              END IF 
                                                                        
   80                     CONTINUE 
                          NCALC = NEND 
                          GO TO 110 
                                                                        
                      END IF 
                                                                        
   90             CONTINUE 
                  N = NEND 
                  EN = CONV(N+N) + (ALPHA+ALPHA) 
!---------------------------------------------------------------------  
! Calculate special significance test for NBMX .GT. 2.                  
!---------------------------------------------------------------------  
                  TEST = MAX(TEST,SQRT(PLAST*ENSIG)*SQRT(P+P)) 
              END IF 
!---------------------------------------------------------------------  
! Calculate P*S until significance test passes.                         
!---------------------------------------------------------------------  
  100         N = N + 1 
              EN = EN + TWO 
              POLD = PLAST 
              PLAST = P 
              P = EN*PLAST/X - POLD 
              IF (P.LT.TEST) GO TO 100 
!---------------------------------------------------------------------  
! Initialize the backward recursion and the normalization SOMME.          
!---------------------------------------------------------------------  
  110         N = N + 1 
              EN = EN + TWO 
              TEMPB = ZERO 
              TEMPA = ONE/P 
              M = 2*N - 4* (N/2) 
              SOMME = ZERO 
              EM = CONV(N/2) 
              ALPEM = (EM-ONE) + ALPHA 
              ALP2EM = (EM+EM) + ALPHA 
              IF (M.NE.0) SOMME = TEMPA*ALPEM*ALP2EM/EM 
              NEND = N - NB 
              IF (NEND.GT.0) THEN 
!---------------------------------------------------------------------  
! Recur backward via difference equation, calculating (but not          
! storing) B(N), until N = NB.                                          
!---------------------------------------------------------------------  
                  DO 120 L = 1,NEND 
                      N = N - 1 
                      EN = EN - TWO 
                      TEMPC = TEMPB 
                      TEMPB = TEMPA 
                      TEMPA = (EN*TEMPB)/X - TEMPC 
                      M = 2 - M 
                      IF (M.NE.0) THEN 
                          EM = EM - ONE 
                          ALP2EM = (EM+EM) + ALPHA 
                          IF (N.EQ.1) GO TO 130 
                          ALPEM = (EM-ONE) + ALPHA 
                          IF (ALPEM.EQ.ZERO) ALPEM = ONE 
                          SOMME = (SOMME+TEMPA*ALP2EM)*ALPEM/EM 
                      END IF 
                                                                        
  120             CONTINUE 
              END IF 
!---------------------------------------------------------------------  
! Store B(NB).                                                          
!---------------------------------------------------------------------  
  130         B(N) = TEMPA 
              IF (NEND.GE.0) THEN 
                  IF (NB.LE.1) THEN 
                      ALP2EM = ALPHA 
                      IF ((ALPHA+ONE).EQ.ONE) ALP2EM = ONE 
                      SOMME = SOMME + B(1)*ALP2EM 
                      GO TO 160 
                                                                        
                  ELSE 
!---------------------------------------------------------------------  
! Calculate and store B(NB-1).                                          
!---------------------------------------------------------------------  
                      N = N - 1 
                      EN = EN - TWO 
                      B(N) = (EN*TEMPA)/X - TEMPB 
                      IF (N.EQ.1) GO TO 150 
                      M = 2 - M 
                      IF (M.NE.0) THEN 
                          EM = EM - ONE 
                          ALP2EM = (EM+EM) + ALPHA 
                          ALPEM = (EM-ONE) + ALPHA 
                          IF (ALPEM.EQ.ZERO) ALPEM = ONE 
                          SOMME = (SOMME+B(N)*ALP2EM)*ALPEM/EM 
                      END IF 
                                                                        
                  END IF 
                                                                        
              END IF 
                                                                        
              NEND = N - 2 
              IF (NEND.NE.0) THEN 
!---------------------------------------------------------------------  
! Calculate via difference equation and store B(N), until N = 2.        
!---------------------------------------------------------------------  
                  DO 140 L = 1,NEND 
                      N = N - 1 
                      EN = EN - TWO 
                      B(N) = (EN*B(N+1))/X - B(N+2) 
                      M = 2 - M 
                      IF (M.NE.0) THEN 
                          EM = EM - ONE 
                          ALP2EM = (EM+EM) + ALPHA 
                          ALPEM = (EM-ONE) + ALPHA 
                          IF (ALPEM.EQ.ZERO) ALPEM = ONE 
                          SOMME = (SOMME+B(N)*ALP2EM)*ALPEM/EM 
                      END IF 
                                                                        
  140             CONTINUE 
              END IF 
!---------------------------------------------------------------------  
! Calculate B(1).                                                       
!---------------------------------------------------------------------  
              B(1) = TWO* (ALPHA+ONE)*B(2)/X - B(3) 
  150         EM = EM - ONE 
              ALP2EM = (EM+EM) + ALPHA 
              IF (ALP2EM.EQ.ZERO) ALP2EM = ONE 
              SOMME = SOMME + B(1)*ALP2EM 
              
!---------------------------------------------------------------------  
! Normalize.  Divide all B(N) by SOMME.                                   
!---------------------------------------------------------------------  
  160         IF ((ALPHA+ONE).NE.ONE) SOMME = SOMME*FUNC(ALPHA)*            &
     &            (X*HALF)** (-ALPHA)                                   
              TEMPA = ENMTEN 
              IF (SOMME.GT.ONE) TEMPA = TEMPA*SOMME 
              DO 170 N = 1,NB 
                  IF (ABS(B(N)).LT.TEMPA) B(N) = ZERO 
                  
                  B(N) = B(N)/SOMME 
                 
  170         CONTINUE 
          END IF 
!---------------------------------------------------------------------  
! Error return -- X, NB, or ALPHA is out of range.                      
!---------------------------------------------------------------------  
      ELSE 
          B(1) = ZERO 
          NCALC = MIN(NB,0) - 1 
      END IF 
!---------------------------------------------------------------------  
! Exit                                                                  
!---------------------------------------------------------------------  
  180 RETURN 
! ---------- Last line of RJBESL ----------                             
      END                                           
!     END OF SUBROUTINE *** RJBESL ***                                  
                                                                      
!     BEGIN OF SUBROUTINE *** RYBESL ***                                
      SUBROUTINE RYBESL(X,ALPHA,NB,BY,NCALC) 
!                                                                       
!     ##################################################################
!                                                                       
!     *** T. WIEDER, 05.03.1998 ***                                     
!                                                                       
!     SOURCE FOR *** SUBROUTINE RYBESL ***:                             
!                                                                       
!     SUBPROGRAM LIBRARY SPECFUN                                        
!                                                                       
!     AVAILABILITIY: PUBLIC DOMAIN                                      
!                                                                       
!     DEVELOPER: W.J. CODY AND L. STOLTZ, APPLIED MATHEMATICS DIVISION, 
!     ARGONNE NATIONAL LABORATORY, ARGONNE, IL 60439.                   
!                                                                       
!     DISTRIBUTOR: NETLIB                                               
!                                                                       
!     ##################################################################
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
!  This routine calculates Bessel functions Y SUB(N+ALPHA) (X)          
!  for non-negative argument X, and non-negative order N+ALPHA.         
!                                                                       
!                                                                       
! Explanation of variables in the calling sequence                      
!                                                                       
! X     - Working precision non-negative real argument for which        
!         Y's are to be calculated.                                     
! ALPHA - Working precision fractional part of order for which          
!         Y's are to be calculated.  0 .LE. ALPHA .LT. 1.0.             
! NB    - Integer number of functions to be calculated, NB .GT. 0.      
!         The first function calculated is of order ALPHA, and the      
!         last is of order (NB - 1 + ALPHA).                            
! BY    - Working precision output vector of length NB.  If the         
!         routine terminates normally (NCALC=NB), the vector BY         
!         contains the functions Y(ALPHA,X), ... , Y(NB-1+ALPHA,X),     
!         If (0 .LT. NCALC .LT. NB), BY(I) contains correct function    
!         values for I .LE. NCALC, and contains the ratios              
!         Y(ALPHA+I-1,X)/Y(ALPHA+I-2,X) for the rest of the array.      
! NCALC - Integer output variable indicating possible errors.           
!         Before using the vector BY, the user should check that        
!         NCALC=NB, i.e., all orders have been calculated to            
!         the desired accuracy.  See error returns below.               
!                                                                       
!                                                                       
!*******************************************************************    
!*******************************************************************    
!                                                                       
! Explanation of machine-dependent constants                            
!                                                                       
!   beta   = Radix for the floating-point system                        
!   p      = Number of significant base-beta digits in the              
!            significand of a floating-point number                     
!   minexp = Smallest representable power of beta                       
!   maxexp = Smallest power of beta that overflows                      
!   EPS    = beta ** (-p)                                               
!   DEL    = Machine number below which sin(x)/x = 1; approximately     
!            SQRT(EPS).                                                 
!   XMIN   = Smallest acceptable argument for RBESY; approximately      
!            max(2*beta**minexp,2/XINF), rounded up                     
!   XINF   = Largest positive machine number; approximately             
!            beta**maxexp                                               
!   THRESH = Lower bound for use of the asymptotic form; approximately  
!            AINT(-LOG10(EPS/2.0))+1.0                                  
!   XLARGE = Upper bound on X; approximately 1/DEL, because the sine    
!            and cosine functions have lost about half of their         
!            precision at that point.                                   
!                                                                       
!                                                                       
!     Approximate values for some important machines are:               
!                                                                       
!                        beta    p     minexp      maxexp      EPS      
!                                                                       
!  CRAY-1        (S.P.)    2    48     -8193        8191    3.55E-15    
!  Cyber 180/185                                                        
!    under NOS   (S.P.)    2    48      -975        1070    3.55E-15    
!  IEEE (IBM/XT,                                                        
!    SUN, etc.)  (S.P.)    2    24      -126         128    5.96E-8     
!  IEEE (IBM/XT,                                                        
!    SUN, etc.)  (D.P.)    2    53     -1022        1024    1.11D-16    
!  IBM 3033      (D.P.)   16    14       -65          63    1.39D-17    
!  VAX           (S.P.)    2    24      -128         127    5.96E-8     
!  VAX D-Format  (D.P.)    2    56      -128         127    1.39D-17    
!  VAX G-Format  (D.P.)    2    53     -1024        1023    1.11D-16    
!                                                                       
!                                                                       
!                         DEL      XMIN      XINF     THRESH  XLARGE    
!                                                                       
! CRAY-1        (S.P.)  5.0E-8  3.67E-2466 5.45E+2465  15.0E0  2.0E7    
! Cyber 180/855                                                         
!   under NOS   (S.P.)  5.0E-8  6.28E-294  1.26E+322   15.0E0  2.0E7    
! IEEE (IBM/XT,                                                         
!   SUN, etc.)  (S.P.)  1.0E-4  2.36E-38   3.40E+38     8.0E0  1.0E4    
! IEEE (IBM/XT,                                                         
!   SUN, etc.)  (D.P.)  1.0D-8  4.46D-308  1.79D+308   16.0D0  1.0D8    
! IBM 3033      (D.P.)  1.0D-8  2.77D-76   7.23D+75    17.0D0  1.0D8    
! VAX           (S.P.)  1.0E-4  1.18E-38   1.70E+38     8.0E0  1.0E4    
! VAX D-Format  (D.P.)  1.0D-9  1.18D-38   1.70D+38    17.0D0  1.0D9    
! VAX G-Format  (D.P.)  1.0D-8  2.23D-308  8.98D+307   16.0D0  1.0D8    
!                                                                       
!*******************************************************************    
!*******************************************************************    
!                                                                       
! Error returns                                                         
!                                                                       
!  In case of an error, NCALC .NE. NB, and not all Y's are              
!  calculated to the desired accuracy.                                  
!                                                                       
!  NCALC .LT. -1:  An argument is out of range. For example,            
!       NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE.          
!       XMAX.  In this case, BY(1) = 0.0, the remainder of the          
!       BY-vector is not calculated, and NCALC is set to                
!       MIN0(NB,0)-2  so that NCALC .NE. NB.                            
!  NCALC = -1:  Y(ALPHA,X) .GE. XINF.  The requested function           
!       values are set to 0.0.                                          
!  1 .LT. NCALC .LT. NB: Not all requested function values could        
!       be calculated accurately.  BY(I) contains correct function      
!       values for I .LE. NCALC, and and the remaining NB-NCALC         
!       array elements contain 0.0.                                     
!                                                                       
!                                                                       
! Intrinsic functions required are:                                     
!                                                                       
!     DBLE, EXP, INT, MAX, MIN, REAL, SQRT                              
!                                                                       
!                                                                       
! Acknowledgement                                                       
!                                                                       
!  This program draws heavily on Temme's Algol program for Y(a,x)       
!  and Y(a+1,x) and on Campbell's programs for Y_nu(x).  Temme's        
!  scheme is used for  x < THRESH, and Campbell's scheme is used        
!  in the asymptotic region.  Segments of code from both sources        
!  have been translated into Fortran 77, merged, and heavily modified.  
!  Modifications include parameterization of machine dependencies,      
!  use of a new approximation for ln(gamma(x)), and built-in            
!  protection against over/underflow.                                   
!                                                                       
! References: "Bessel functions J_nu(x) and Y_nu(x) of real             
!              order and real argument," Campbell, J. B.,               
!              Comp. Phy. Comm. 18, 1979, pp. 133-142.                  
!                                                                       
!             "On the numerical evaluation of the ordinary              
!              Bessel function of the second kind," Temme,              
!              N. M., J. Comput. Phys. 21, 1976, pp. 343-350.           
!                                                                       
!  Latest modification: March 19, 1990                                  
!                                                                       
!  Modified by: W. J. Cody                                              
!               Applied Mathematics Division                            
!               Argonne National Laboratory                             
!               Argonne, IL  60439                                      
!                                                                       
!---------------------------------------------------------------------- 
!S    REAL                                                              
!      DIMENSION BY(NB), CH(21)                                         
!     *** T. WIEDER, 07.02.1999 ***                                     
!     .. Scalar Arguments ..                                            
      DOUBLE PRECISION ALPHA,X 
      INTEGER NB,NCALC 
!     ..                                                                
!     .. Array Arguments ..                                             
      DOUBLE PRECISION BY(NB+1) 
!     ..                                                                
!     .. Local Scalars ..                                               
      DOUBLE PRECISION ALFA,AYE,B,C,COSMU,D,D1,D2,DDIV,DEL,DEN,DIV,DMU, &
     &                 E,EIGHT,EN,EN1,ENU,EPS,EVEN,EX,F,FIVPI,G,GAMMA,H,&
     &                 HALF,ODD,ONBPI,ONE,ONE5,P,PA,PA1,PI,PIBY2,PIM5,Q,&
     &                 Q0,QA,QA1,R,S,SINMU,SQ2BPI,TEN9,TERM,THREE,      &
     &                 THRESH,TWO,TWOBYX,X2,XINF,XLARGE,XMIN,XNA,YA,YA1,&
     &                 ZERO                                             
      INTEGER I,K,NA 
!     ..                                                                
!     .. Local Arrays ..                                                
      DOUBLE PRECISION CH(21) 
!     ..                                                                
!     .. Intrinsic Functions ..                                         
      INTRINSIC ABS,AINT,COS,INT,LOG,MIN,SIN,SQRT 
!     ..                                                                
!     .. Data statements ..                                             
!---------------------------------------------------------------------- 
!  Mathematical constants                                               
!    FIVPI = 5*PI                                                       
!    PIM5 = 5*PI - 15                                                   
!    ONBPI = 1/PI                                                       
!    PIBY2 = PI/2                                                       
!    SQ2BPI = SQUARE ROOT OF 2/PI                                       
!---------------------------------------------------------------------- 
!S    DATA ZERO,HALF,ONE,TWO,THREE/0.0E0,0.5E0,1.0E0,2.0E0,3.0E0/       
!S    DATA EIGHT,ONE5,TEN9/8.0E0,15.0E0,1.9E1/                          
!S    DATA FIVPI,PIBY2/1.5707963267948966192E1,1.5707963267948966192E0/ 
!S    DATA PI,SQ2BPI/3.1415926535897932385E0,7.9788456080286535588E-1/  
!S    DATA PIM5,ONBPI/7.0796326794896619231E-1,3.1830988618379067154E-1/
!---------------------------------------------------------------------- 
!  Machine-dependent constants                                          
!---------------------------------------------------------------------- 
!S    DATA DEL,XMIN,XINF,EPS/1.0E-4,2.36E-38,3.40E38,5.96E-8/           
!S    DATA THRESH,XLARGE/8.0E0,1.0E4/                                   
!     *** T. WIEDER, 11.03.1998 ***                                     
!      DATA DEL,XMIN,XINF,EPS/1.0D-8,4.46D-308,1.79D308,1.11D-16/       
!---------------------------------------------------------------------- 
!  Coefficients for Chebyshev polynomial expansion of                   
!         1/gamma(1-x), abs(x) .le. .5                                  
!---------------------------------------------------------------------- 
!S    DATA CH/-0.67735241822398840964E-23,-0.61455180116049879894E-22,  
!S   1         0.29017595056104745456E-20, 0.13639417919073099464E-18,  
!S   2         0.23826220476859635824E-17,-0.90642907957550702534E-17,  
!S   3        -0.14943667065169001769E-14,-0.33919078305362211264E-13,  
!S   4        -0.17023776642512729175E-12, 0.91609750938768647911E-11,  
!S   5         0.24230957900482704055E-09, 0.17451364971382984243E-08,  
!S   6        -0.33126119768180852711E-07,-0.86592079961391259661E-06,  
!S   7        -0.49717367041957398581E-05, 0.76309597585908126618E-04,  
!S   8         0.12719271366545622927E-02, 0.17063050710955562222E-02,  
!S   9        -0.76852840844786673690E-01,-0.28387654227602353814E+00,  
!S   A         0.92187029365045265648E+00/                              
      DATA ZERO,HALF,ONE,TWO,THREE/0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/ 
      DATA EIGHT,ONE5,TEN9/8.0D0,15.0D0,1.9D1/ 
      DATA FIVPI,PIBY2/1.5707963267948966192D1,1.5707963267948966192D0/ 
      DATA PI,SQ2BPI/3.1415926535897932385D0,7.9788456080286535588D-1/ 
      DATA PIM5,ONBPI/7.0796326794896619231D-1,3.1830988618379067154D-1/ 
      DATA DEL,XMIN,XINF,EPS/1.0D-8,2.23D-307,1.79D308,1.11D-16/ 
      DATA THRESH,XLARGE/16.0D0,1.0D8/ 
      DATA CH/-0.67735241822398840964D-23,-0.61455180116049879894D-22,  &
     &     0.29017595056104745456D-20,0.13639417919073099464D-18,       &
     &     0.23826220476859635824D-17,-0.90642907957550702534D-17,      &
     &     -0.14943667065169001769D-14,-0.33919078305362211264D-13,     &
     &     -0.17023776642512729175D-12,0.91609750938768647911D-11,      &
     &     0.24230957900482704055D-09,0.17451364971382984243D-08,       &
     &     -0.33126119768180852711D-07,-0.86592079961391259661D-06,     &
     &     -0.49717367041957398581D-05,0.76309597585908126618D-04,      &
     &     0.12719271366545622927D-02,0.17063050710955562222D-02,       &
     &     -0.76852840844786673690D-01,-0.28387654227602353814D+00,     &
     &     0.92187029365045265648D+00/                                  
!     ..                                                                
!---------------------------------------------------------------------- 
      EX = X 
      ENU = ALPHA 
      IF ((NB.GT.0) .AND. (X.GE.XMIN) .AND. (EX.LT.XLARGE) .AND.        &
     &    (ENU.GE.ZERO) .AND. (ENU.LT.ONE)) THEN                        
          XNA = AINT(ENU+HALF) 
          NA = INT(XNA) 
          IF (NA.EQ.1) ENU = ENU - XNA 
          IF (ENU.EQ.-HALF) THEN 
              P = SQ2BPI/SQRT(EX) 
              YA = P*SIN(EX) 
              YA1 = -P*COS(EX) 
                                                                        
          ELSE IF (EX.LT.THREE) THEN 
!---------------------------------------------------------------------- 
!  Use Temme's scheme for small X                                       
!---------------------------------------------------------------------- 
              B = EX*HALF 
              D = -LOG(B) 
              F = ENU*D 
              E = B** (-ENU) 
              IF (ABS(ENU).LT.DEL) THEN 
                  C = ONBPI 
                                                                        
              ELSE 
                  C = ENU/SIN(ENU*PI) 
              END IF 
!---------------------------------------------------------------------- 
!  Computation of sinh(f)/f                                             
!---------------------------------------------------------------------- 
              IF (ABS(F).LT.ONE) THEN 
                  X2 = F*F 
                  EN = TEN9 
                  S = ONE 
                  DO 10 I = 1,9 
                      S = S*X2/EN/ (EN-ONE) + ONE 
                      EN = EN - TWO 
   10             CONTINUE 
                                                                        
              ELSE 
                  S = (E-ONE/E)*HALF/F 
              END IF 
!---------------------------------------------------------------------- 
!  Computation of 1/gamma(1-a) using Chebyshev polynomials              
!---------------------------------------------------------------------- 
              X2 = ENU*ENU*EIGHT 
              AYE = CH(1) 
              EVEN = ZERO 
              ALFA = CH(2) 
              ODD = ZERO 
              DO 20 I = 3,19,2 
                  EVEN = - (AYE+AYE+EVEN) 
                  AYE = -EVEN*X2 - AYE + CH(I) 
                  ODD = - (ALFA+ALFA+ODD) 
                  ALFA = -ODD*X2 - ALFA + CH(I+1) 
   20         CONTINUE 
              EVEN = (EVEN*HALF+AYE)*X2 - AYE + CH(21) 
              ODD = (ODD+ALFA)*TWO 
              GAMMA = ODD*ENU + EVEN 
!---------------------------------------------------------------------- 
!  End of computation of 1/gamma(1-a)                                   
!---------------------------------------------------------------------- 
              G = E*GAMMA 
              E = (E+ONE/E)*HALF 
              F = TWO*C* (ODD*E+EVEN*S*D) 
              E = ENU*ENU 
              P = G*C 
              Q = ONBPI/G 
              C = ENU*PIBY2 
              IF (ABS(C).LT.DEL) THEN 
                  R = ONE 
                                                                        
              ELSE 
                  R = SIN(C)/C 
              END IF 
                                                                        
              R = PI*C*R*R 
              C = ONE 
              D = -B*B 
              H = ZERO 
              YA = F + R*Q 
              YA1 = P 
              EN = ZERO 
   30         EN = EN + ONE 
              IF (ABS(G/ (ONE+ABS(YA)))+ABS(H/ (ONE+ABS(YA1))).GT.      &
     &            EPS) THEN                                             
                  F = (F*EN+P+Q)/ (EN*EN-E) 
                  C = C*D/EN 
                  P = P/ (EN-ENU) 
                  Q = Q/ (EN+ENU) 
                  G = C* (F+R*Q) 
                  H = C*P - EN*G 
                  YA = YA + G 
                  YA1 = YA1 + H 
                  GO TO 30 
                                                                        
              END IF 
                                                                        
              YA = -YA 
              YA1 = -YA1/B 
                                                                        
          ELSE IF (EX.LT.THRESH) THEN 
!---------------------------------------------------------------------- 
!  Use Temme's scheme for moderate X                                    
!---------------------------------------------------------------------- 
              C = (HALF-ENU)* (HALF+ENU) 
              B = EX + EX 
              E = (EX*ONBPI*COS(ENU*PI)/EPS) 
              E = E*E 
              P = ONE 
              Q = -EX 
              R = ONE + EX*EX 
              S = R 
              EN = TWO 
   40         IF (R*EN*EN.LT.E) THEN 
                  EN1 = EN + ONE 
                  D = (EN-ONE+C/EN)/S 
                  P = (EN+EN-P*D)/EN1 
                  Q = (-B+Q*D)/EN1 
                  S = P*P + Q*Q 
                  R = R*S 
                  EN = EN1 
                  GO TO 40 
                                                                        
              END IF 
                                                                        
              F = P/S 
              P = F 
              G = -Q/S 
              Q = G 
   50         EN = EN - ONE 
              IF (EN.GT.ZERO) THEN 
                  R = EN1* (TWO-P) - TWO 
                  S = B + EN1*Q 
                  D = (EN-ONE+C/EN)/ (R*R+S*S) 
                  P = D*R 
                  Q = D*S 
                  E = F + ONE 
                  F = P*E - G*Q 
                  G = Q*E + P*G 
                  EN1 = EN 
                  GO TO 50 
                                                                        
              END IF 
                                                                        
              F = ONE + F 
              D = F*F + G*G 
              PA = F/D 
              QA = -G/D 
              D = ENU + HALF - P 
              Q = Q + EX 
              PA1 = (PA*Q-QA*D)/EX 
              QA1 = (QA*Q+PA*D)/EX 
              B = EX - PIBY2* (ENU+HALF) 
              C = COS(B) 
              S = SIN(B) 
              D = SQ2BPI/SQRT(EX) 
              YA = D* (PA*S+QA*C) 
              YA1 = D* (QA1*S-PA1*C) 
                                                                        
          ELSE 
!---------------------------------------------------------------------- 
!  Use Campbell's asymptotic scheme.                                    
!---------------------------------------------------------------------- 
              NA = 0 
              D1 = AINT(EX/FIVPI) 
              I = INT(D1) 
              DMU = ((EX-ONE5*D1)-D1*PIM5) - (ALPHA+HALF)*PIBY2 
              IF (I-2* (I/2).EQ.0) THEN 
                  COSMU = COS(DMU) 
                  SINMU = SIN(DMU) 
                                                                        
              ELSE 
                  COSMU = -COS(DMU) 
                  SINMU = -SIN(DMU) 
              END IF 
                                                                        
              DDIV = EIGHT*EX 
              DMU = ALPHA 
              DEN = SQRT(EX) 
              DO 80 K = 1,2 
                  P = COSMU 
                  COSMU = SINMU 
                  SINMU = -P 
                  D1 = (TWO*DMU-ONE)* (TWO*DMU+ONE) 
                  D2 = ZERO 
                  DIV = DDIV 
                  P = ZERO 
                  Q = ZERO 
                  Q0 = D1/DIV 
                  TERM = Q0 
                  DO 60 I = 2,20 
                      D2 = D2 + EIGHT 
                      D1 = D1 - D2 
                      DIV = DIV + DDIV 
                      TERM = -TERM*D1/DIV 
                      P = P + TERM 
                      D2 = D2 + EIGHT 
                      D1 = D1 - D2 
                      DIV = DIV + DDIV 
                      TERM = TERM*D1/DIV 
                      Q = Q + TERM 
                      IF (ABS(TERM).LE.EPS) GO TO 70 
   60             CONTINUE 
   70             P = P + ONE 
                  Q = Q + Q0 
                  IF (K.EQ.1) THEN 
                      YA = SQ2BPI* (P*COSMU-Q*SINMU)/DEN 
                                                                        
                  ELSE 
                      YA1 = SQ2BPI* (P*COSMU-Q*SINMU)/DEN 
                  END IF 
                                                                        
                  DMU = DMU + ONE 
   80         CONTINUE 
          END IF 
                                                                        
          IF (NA.EQ.1) THEN 
              H = TWO* (ENU+ONE)/EX 
              IF (H.GT.ONE) THEN 
                  IF (ABS(YA1).GT.XINF/H) THEN 
                      H = ZERO 
                      YA = ZERO 
                  END IF 
                                                                        
              END IF 
                                                                        
              H = H*YA1 - YA 
              YA = YA1 
              YA1 = H 
          END IF 
!---------------------------------------------------------------------- 
!  Now have first one or two Y's                                        
!---------------------------------------------------------------------- 
          BY(1) = YA 
          BY(2) = YA1 
          IF (YA1.EQ.ZERO) THEN 
              NCALC = 1 
                                                                        
          ELSE 
              AYE = ONE + ALPHA 
              TWOBYX = TWO/EX 
              NCALC = 2 
              DO 90 I = 3,NB 
                  IF (TWOBYX.LT.ONE) THEN 
                      IF (ABS(BY(I-1))*TWOBYX.GE.XINF/AYE) GO TO 100 
                                                                        
                  ELSE 
                      IF (ABS(BY(I-1)).GE.XINF/AYE/TWOBYX) GO TO 100 
                  END IF 
                                                                        
                  BY(I) = TWOBYX*AYE*BY(I-1) - BY(I-2) 
                  AYE = AYE + ONE 
                  NCALC = NCALC + 1 
   90         CONTINUE 
          END IF 
                                                                        
  100     DO 110 I = NCALC + 1,NB 
              BY(I) = ZERO 
  110     CONTINUE 
                                                                        
      ELSE 
          BY(1) = ZERO 
          NCALC = MIN(NB,0) - 1 
      END IF 
                                                                        
      RETURN 
!---------- Last line of RYBESL ----------                              
      END                                           
!     END OF SUBROUTINE *** RYBESL ***                                  
! !                                                                       
!     BEGIN OF FUNCTION *** DLGAMA ***                                  
!S    REAL FUNCTION ALGAMA(X)                                           
      DOUBLE PRECISION FUNCTION DLGAMA(X) 
!                                                                       
!     ##################################################################
!                                                                       
!     *** T. WIEDER, 05.03.1998 ***                                     
!                                                                       
!     SOURCE FOR *** FUNCTION DLGAMA ***:                               
!                                                                       
!     SUBPROGRAM LIBRARY SPECFUN                                        
!                                                                       
!     AVAILABILITIY: PUBLIC DOMAIN                                      
!                                                                       
!     DEVELOPER: W.J. CODY AND L. STOLTZ, APPLED MATHEMATICS DIVISION,  
!     ARGONNE NATIONAL LABORATORY, ARGONNE, IL 60439.                   
!                                                                       
!                                                                       
!     DISTRIBUTOR: NETLIB                                               
!                                                                       
!     ##################################################################
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
! This routine calculates the LOG(GAMMA) function for a positive real   
!   argument X.  Computation is based on an algorithm outlined in       
!   references 1 and 2.  The program uses rational functions that       
!   theoretically approximate LOG(GAMMA) to at least 18 significant     
!   decimal digits.  The approximation for X > 12 is from reference     
!   3, while approximations for X < 12.0 are similar to those in        
!   reference 1, but are unpublished.  The accuracy achieved depends    
!   on the arithmetic system, the compiler, the intrinsic functions,    
!   and proper selection of the machine-dependent constants.            
!                                                                       
!                                                                       
!*********************************************************************  
!*********************************************************************  
!                                                                       
! Explanation of machine-dependent constants                            
!                                                                       
! beta   - radix for the floating-point representation                  
! maxexp - the smallest positive power of beta that overflows           
! XBIG   - largest argument for which LN(GAMMA(X)) is representable     
!          in the machine, i.e., the solution to the equation           
!                  LN(GAMMA(XBIG)) = beta**maxexp                       
! XINF   - largest machine representable floating-point number;         
!          approximately beta**maxexp.                                  
! EPS    - The smallest positive floating-point number such that        
!          1.0+EPS .GT. 1.0                                             
! FRTBIG - Rough estimate of the fourth root of XBIG                    
!                                                                       
!                                                                       
!     Approximate values for some important machines are:               
!                                                                       
!                           beta      maxexp         XBIG               
!                                                                       
! CRAY-1        (S.P.)        2        8191       9.62E+2461            
! Cyber 180/855                                                         
!   under NOS   (S.P.)        2        1070       1.72E+319             
! IEEE (IBM/XT,                                                         
!   SUN, etc.)  (S.P.)        2         128       4.08E+36              
! IEEE (IBM/XT,                                                         
!   SUN, etc.)  (D.P.)        2        1024       2.55D+305             
! IBM 3033      (D.P.)       16          63       4.29D+73              
! VAX D-Format  (D.P.)        2         127       2.05D+36              
! VAX G-Format  (D.P.)        2        1023       1.28D+305             
!                                                                       
!                                                                       
!                           XINF        EPS        FRTBIG               
!                                                                       
! CRAY-1        (S.P.)   5.45E+2465   7.11E-15    3.13E+615             
! Cyber 180/855                                                         
!   under NOS   (S.P.)   1.26E+322    3.55E-15    6.44E+79              
! IEEE (IBM/XT,                                                         
!   SUN, etc.)  (S.P.)   3.40E+38     1.19E-7     1.42E+9               
! IEEE (IBM/XT,                                                         
!   SUN, etc.)  (D.P.)   1.79D+308    2.22D-16    2.25D+76              
! IBM 3033      (D.P.)   7.23D+75     2.22D-16    2.56D+18              
! VAX D-Format  (D.P.)   1.70D+38     1.39D-17    1.20D+9               
! VAX G-Format  (D.P.)   8.98D+307    1.11D-16    1.89D+76              
!                                                                       
!**************************************************************         
!**************************************************************         
!                                                                       
! Error returns                                                         
!                                                                       
!  The program returns the value XINF for X .LE. 0.0 or when            
!     overflow would occur.  The computation is believed to             
!     be free of underflow and overflow.                                
!                                                                       
!                                                                       
! Intrinsic functions required are:                                     
!                                                                       
!      LOG                                                              
!                                                                       
!                                                                       
! References:                                                           
!                                                                       
!  1) W. J. Cody and K. E. Hillstrom, 'Chebyshev Approximations for     
!     the Natural Logarithm of the Gamma Function,' Math. Comp. 21,     
!     1967, pp. 198-203.                                                
!                                                                       
!  2) K. E. Hillstrom, ANL/AMD Program ANLC366S, DGAMMA/DLGAMA, May,    
!     1969.                                                             
!                                                                       
!  3) Hart, Et. Al., Computer Approximations, Wiley and sons, New       
!     York, 1968.                                                       
!                                                                       
!                                                                       
!  Authors: W. J. Cody and L. Stoltz                                    
!           Argonne National Laboratory                                 
!                                                                       
!  Latest modification: June 16, 1988                                   
!                                                                       
!---------------------------------------------------------------------- 
!S    REAL                                                              
!     .. Scalar Arguments ..                                            
      DOUBLE PRECISION X 
!     ..                                                                
!     .. Local Scalars ..                                               
      DOUBLE PRECISION CORR,D1,D2,D4,EPS,FOUR,FRTBIG,HALF,ONE,PNT68,RES,&
     &                 SQRTPI,THRHAL,TWELVE,TWO,XBIG,XDEN,XINF,XM1,XM2, &
     &                 XM4,XNUM,Y,YSQ,ZERO                              
      INTEGER I 
!     ..                                                                
!     .. Local Arrays ..                                                
      DOUBLE PRECISION C(7),P1(8),P2(8),P4(8),Q1(8),Q2(8),Q4(8) 
!     ..                                                                
!     .. Intrinsic Functions ..                                         
      INTRINSIC LOG 
!     ..                                                                
!     .. Data statements ..                                             
!---------------------------------------------------------------------- 
!  Mathematical constants                                               
!---------------------------------------------------------------------- 
!S    DATA ONE,HALF,TWELVE,ZERO/1.0E0,0.5E0,12.0E0,0.0E0/,              
!S   1     FOUR,THRHAL,TWO,PNT68/4.0E0,1.5E0,2.0E0,0.6796875E0/,        
!S   2     SQRTPI/0.9189385332046727417803297E0/                        
!---------------------------------------------------------------------- 
!  Machine dependent parameters                                         
!---------------------------------------------------------------------- 
!S    DATA XBIG,XINF,EPS,FRTBIG/4.08E36,3.401E38,1.19E-7,1.42E9/        
!---------------------------------------------------------------------- 
!  Numerator and denominator coefficients for rational minimax          
!     approximation over (0.5,1.5).                                     
!---------------------------------------------------------------------- 
!S    DATA D1/-5.772156649015328605195174E-1/                           
!S    DATA P1/4.945235359296727046734888E0,2.018112620856775083915565E2,
!S   1        2.290838373831346393026739E3,1.131967205903380828685045E4,
!S   2        2.855724635671635335736389E4,3.848496228443793359990269E4,
!S   3        2.637748787624195437963534E4,7.225813979700288197698961E3/
!S    DATA Q1/6.748212550303777196073036E1,1.113332393857199323513008E3,
!S   1        7.738757056935398733233834E3,2.763987074403340708898585E4,
!S   2        5.499310206226157329794414E4,6.161122180066002127833352E4,
!S   3        3.635127591501940507276287E4,8.785536302431013170870835E3/
!---------------------------------------------------------------------- 
!  Numerator and denominator coefficients for rational minimax          
!     Approximation over (1.5,4.0).                                     
!---------------------------------------------------------------------- 
!S    DATA D2/4.227843350984671393993777E-1/                            
!S    DATA P2/4.974607845568932035012064E0,5.424138599891070494101986E2,
!S   1        1.550693864978364947665077E4,1.847932904445632425417223E5,
!S   2        1.088204769468828767498470E6,3.338152967987029735917223E6,
!S   3        5.106661678927352456275255E6,3.074109054850539556250927E6/
!S    DATA Q2/1.830328399370592604055942E2,7.765049321445005871323047E3,
!S   1        1.331903827966074194402448E5,1.136705821321969608938755E6,
!S   2        5.267964117437946917577538E6,1.346701454311101692290052E7,
!S   3        1.782736530353274213975932E7,9.533095591844353613395747E6/
!---------------------------------------------------------------------- 
!  Numerator and denominator coefficients for rational minimax          
!     Approximation over (4.0,12.0).                                    
!---------------------------------------------------------------------- 
!S    DATA D4/1.791759469228055000094023E0/                             
!S    DATA P4/1.474502166059939948905062E4,2.426813369486704502836312E6,
!S   1        1.214755574045093227939592E8,2.663432449630976949898078E9,
!S   2      2.940378956634553899906876E10,1.702665737765398868392998E11,
!S   3      4.926125793377430887588120E11,5.606251856223951465078242E11/
!S    DATA Q4/2.690530175870899333379843E3,6.393885654300092398984238E5,
!S   2        4.135599930241388052042842E7,1.120872109616147941376570E9,
!S   3      1.488613728678813811542398E10,1.016803586272438228077304E11,
!S   4      3.417476345507377132798597E11,4.463158187419713286462081E11/
!---------------------------------------------------------------------- 
!  Coefficients for minimax approximation over (12, INF).               
!---------------------------------------------------------------------- 
!S    DATA C/-1.910444077728E-03,8.4171387781295E-04,                   
!S   1     -5.952379913043012E-04,7.93650793500350248E-04,              
!S   2     -2.777777777777681622553E-03,8.333333333333333331554247E-02, 
!S   3      5.7083835261E-03/                                           
      DATA ONE,HALF,TWELVE,ZERO/1.0D0,0.5D0,12.0D0,0.0D0/,FOUR,THRHAL,  &
     &     TWO,PNT68/4.0D0,1.5D0,2.0D0,0.6796875D0/,                    &
     &     SQRTPI/0.9189385332046727417803297D0/                        
      DATA XBIG,XINF,EPS,FRTBIG/2.55D305,1.79D308,2.22D-16,2.25D76/ 
      DATA D1/-5.772156649015328605195174D-1/ 
      DATA P1/4.945235359296727046734888D0,2.018112620856775083915565D2,&
     &     2.290838373831346393026739D3,1.131967205903380828685045D4,   &
     &     2.855724635671635335736389D4,3.848496228443793359990269D4,   &
     &     2.637748787624195437963534D4,7.225813979700288197698961D3/   
      DATA Q1/6.748212550303777196073036D1,1.113332393857199323513008D3,&
     &     7.738757056935398733233834D3,2.763987074403340708898585D4,   &
     &     5.499310206226157329794414D4,6.161122180066002127833352D4,   &
     &     3.635127591501940507276287D4,8.785536302431013170870835D3/   
      DATA D2/4.227843350984671393993777D-1/ 
      DATA P2/4.974607845568932035012064D0,5.424138599891070494101986D2,&
     &     1.550693864978364947665077D4,1.847932904445632425417223D5,   &
     &     1.088204769468828767498470D6,3.338152967987029735917223D6,   &
     &     5.106661678927352456275255D6,3.074109054850539556250927D6/   
      DATA Q2/1.830328399370592604055942D2,7.765049321445005871323047D3,&
     &     1.331903827966074194402448D5,1.136705821321969608938755D6,   &
     &     5.267964117437946917577538D6,1.346701454311101692290052D7,   &
     &     1.782736530353274213975932D7,9.533095591844353613395747D6/   
      DATA D4/1.791759469228055000094023D0/ 
      DATA P4/1.474502166059939948905062D4,2.426813369486704502836312D6,&
     &     1.214755574045093227939592D8,2.663432449630976949898078D9,   &
     &     2.940378956634553899906876D10,1.702665737765398868392998D11, &
     &     4.926125793377430887588120D11,5.606251856223951465078242D11/ 
      DATA Q4/2.690530175870899333379843D3,6.393885654300092398984238D5,&
     &     4.135599930241388052042842D7,1.120872109616147941376570D9,   &
     &     1.488613728678813811542398D10,1.016803586272438228077304D11, &
     &     3.417476345507377132798597D11,4.463158187419713286462081D11/ 
      DATA C/-1.910444077728D-03,8.4171387781295D-04,                   &
     &     -5.952379913043012D-04,7.93650793500350248D-04,              &
     &     -2.777777777777681622553D-03,8.333333333333333331554247D-02, &
     &     5.7083835261D-03/                                           
!     ..                                                                
!---------------------------------------------------------------------- 
      Y = X 
      IF ((Y.GT.ZERO) .AND. (Y.LE.XBIG)) THEN 
          IF (Y.LE.EPS) THEN 
              RES = -LOG(Y) 
                                                                        
          ELSE IF (Y.LE.THRHAL) THEN 
!---------------------------------------------------------------------- 
!  EPS .LT. X .LE. 1.5                                                  
!---------------------------------------------------------------------- 
              IF (Y.LT.PNT68) THEN 
                  CORR = -LOG(Y) 
                  XM1 = Y 
                                                                        
              ELSE 
                  CORR = ZERO 
                  XM1 = (Y-HALF) - HALF 
              END IF 
                                                                        
              IF ((Y.LE.HALF) .OR. (Y.GE.PNT68)) THEN 
                  XDEN = ONE 
                  XNUM = ZERO 
                  DO 10 I = 1,8 
                      XNUM = XNUM*XM1 + P1(I) 
                      XDEN = XDEN*XM1 + Q1(I) 
   10             CONTINUE 
                  RES = CORR + (XM1* (D1+XM1* (XNUM/XDEN))) 
                                                                        
              ELSE 
                  XM2 = (Y-HALF) - HALF 
                  XDEN = ONE 
                  XNUM = ZERO 
                  DO 20 I = 1,8 
                      XNUM = XNUM*XM2 + P2(I) 
                      XDEN = XDEN*XM2 + Q2(I) 
   20             CONTINUE 
                  RES = CORR + XM2* (D2+XM2* (XNUM/XDEN)) 
              END IF 
                                                                        
          ELSE IF (Y.LE.FOUR) THEN 
!---------------------------------------------------------------------- 
!  1.5 .LT. X .LE. 4.0                                                  
!---------------------------------------------------------------------- 
              XM2 = Y - TWO 
              XDEN = ONE 
              XNUM = ZERO 
              DO 30 I = 1,8 
                  XNUM = XNUM*XM2 + P2(I) 
                  XDEN = XDEN*XM2 + Q2(I) 
   30         CONTINUE 
              RES = XM2* (D2+XM2* (XNUM/XDEN)) 
                                                                        
          ELSE IF (Y.LE.TWELVE) THEN 
!---------------------------------------------------------------------- 
!  4.0 .LT. X .LE. 12.0                                                 
!---------------------------------------------------------------------- 
              XM4 = Y - FOUR 
              XDEN = -ONE 
              XNUM = ZERO 
              DO 40 I = 1,8 
                  XNUM = XNUM*XM4 + P4(I) 
                  XDEN = XDEN*XM4 + Q4(I) 
   40         CONTINUE 
              RES = D4 + XM4* (XNUM/XDEN) 
                                                                        
          ELSE 
!---------------------------------------------------------------------- 
!  Evaluate for argument .GE. 12.0,                                     
!---------------------------------------------------------------------- 
              RES = ZERO 
              IF (Y.LE.FRTBIG) THEN 
                  RES = C(7) 
                  YSQ = Y*Y 
                  DO 50 I = 1,6 
                      RES = RES/YSQ + C(I) 
   50             CONTINUE 
              END IF 
                                                                        
              RES = RES/Y 
              CORR = LOG(Y) 
              RES = RES + SQRTPI - HALF*CORR 
              RES = RES + Y* (CORR-ONE) 
          END IF 
                                                                        
      ELSE 
!---------------------------------------------------------------------- 
!  Return for bad arguments                                             
!---------------------------------------------------------------------- 
          RES = XINF 
      END IF 
!---------------------------------------------------------------------- 
!  Final adjustments and return                                         
!---------------------------------------------------------------------- 
!S    ALGAMA = RES                                                      
      DLGAMA = RES 
      RETURN 
! ---------- Last line of DLGAMA ----------                             
      END                                           
!     END OF FUNCTION *** DLGAMA ***                                    

