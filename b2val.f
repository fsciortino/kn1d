      REAL FUNCTION B2VAL(XVAL,YVAL,IDX,IDY,TX,TY,NX,NY,
     *  KX,KY,BCOEF,WORK)
C***BEGIN PROLOGUE  B2VAL
C***DATE WRITTEN   25 MAY 1982
C***REVISION DATE  25 MAY 1982
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  E1A
C***KEYWORDS  INTERPOLATION, TWO-DIMENSIONS, GRIDDED DATA, SPLINES,
C             PIECEWISE POLYNOMIALS
C***AUTHOR  BOISVERT, RONALD, NBS
C             SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS
C             WASHINGTON, DC 20234
C***PURPOSE  B2VAL EVALUATES THE PIECEWISE POLYNOMIAL INTERPOLATING
C            FUNCTION CONSTRUCTED BY THE ROUTINE B2INK OR ONE OF ITS
C            PARTIAL DERIVATIVES.
C***DESCRIPTION
C
C   B2VAL evaluates the tensor product piecewise polynomial interpolant
C   constructed by the routine B2INK or one of its derivatives  at  the
C   point (XVAL,YVAL). To evaluate the interpolant itself, set IDX=IDY=
C   0, to evaluate the first partial with respect to x, set  IDX=1,IDY=
C   0, and so on.
C
C   B2VAL returns 0.0E0 if (XVAL,YVAL) is out of range. That is, if
C            XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR.
C            YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY)
C   If the knots  TX  and  TY  were  chosen  by  B2INK,  then  this  is
C   equivalent to
C            XVAL.LT.X(1) .OR. XVAL.GT.X(NX)+EPSX .OR.
C            YVAL.LT.Y(1) .OR. YVAL.GT.Y(NY)+EPSY
C   where EPSX = 0.1*(X(NX)-X(NX-1)) and EPSY = 0.1*(Y(NY)-Y(NY-1)).
C
C   The input quantities TX, TY, NX, NY, KX, KY, and  BCOEF  should  be
C   unchanged since the last call of B2INK.
C
C
C   I N P U T
C   ---------
C
C   XVAL    Real scalar
C           X coordinate of evaluation point.
C
C   YVAL    Real scalar
C           Y coordinate of evaluation point.
C
C   IDX     Integer scalar
C           X derivative of piecewise polynomial to evaluate.
C
C   IDY     Integer scalar
C           Y derivative of piecewise polynomial to evaluate.
C
C   TX      Real 1D array (size NX+KX)
C           Sequence of knots defining the piecewise polynomial in
C           the x direction.  (Same as in last call to B2INK.)
C
C   TY      Real 1D array (size NY+KY)
C           Sequence of knots defining the piecewise polynomial in
C           the y direction.  (Same as in last call to B2INK.)
C
C   NX      Integer scalar
C           The number of interpolation points in x.
C           (Same as in last call to B2INK.)
C
C   NY      Integer scalar
C           The number of interpolation points in y.
C           (Same as in last call to B2INK.)
C
C   KX      Integer scalar
C           Order of polynomial pieces in x.
C           (Same as in last call to B2INK.)
C
C   KY      Integer scalar
C           Order of polynomial pieces in y.
C           (Same as in last call to B2INK.)
C
C   BCOEF   Real 2D array (size NX by NY)
C           The B-spline coefficients computed by B2INK.
C
C   WORK    Real 1D array (size 3*max(KX,KY) + KY)
C           A working storage array.
C
C***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES,
C                 SPRINGER-VERLAG, NEW YORK, 1978.
C***ROUTINES CALLED  INTRV,BVALU
C***END PROLOGUE  B2VAL
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   MODIFICATION
C   ------------
C
C   ADDED CHECK TO SEE IF X OR Y IS OUT OF RANGE, IF SO, RETURN 0.0
C
C   R.F. BOISVERT, NIST
C   22 FEB 00
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  PARAMETERS
C
      INTEGER
     *        IDX, IDY, NX, NY, KX, KY
      REAL
     *     XVAL, YVAL, TX(*), TY(*), BCOEF(NX,NY), WORK(*)
C
C  LOCAL VARIABLES
C
      INTEGER
     *        ILOY, INBVX, INBV, K, LEFTY, MFLAG, KCOL, IW
      REAL
     *     BVALU
C
      DATA ILOY /1/,  INBVX /1/
C     SAVE ILOY    ,  INBVX
C
C
C***FIRST EXECUTABLE STATEMENT
      B2VAL = 0.0E0
C  NEXT STATEMENT - RFB MOD
      IF (XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR.
     +    YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY))      GO TO 100
      CALL INTRV(TY,NY+KY,YVAL,ILOY,LEFTY,MFLAG)
      IF (MFLAG .NE. 0)  GO TO 100
         IW = KY + 1
         KCOL = LEFTY - KY
         DO 50 K=1,KY
            KCOL = KCOL + 1
            WORK(K) = BVALU(TX,BCOEF(1,KCOL),NX,KX,IDX,XVAL,INBVX,
     *                      WORK(IW))
   50    CONTINUE
         INBV = 1
         KCOL = LEFTY - KY + 1
         B2VAL = BVALU(TY(KCOL),WORK,KY,KY,IDY,YVAL,INBV,WORK(IW))
  100 CONTINUE
      RETURN
      END
      FUNCTION BVALU(T,A,N,K,IDERIV,X,INBV,WORK)
C***BEGIN PROLOGUE  BVALU
C***DATE WRITTEN   800901   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  E3,K6
C***KEYWORDS  B-SPLINE,DATA FITTING,INTERPOLATION,SPLINE
C***AUTHOR  AMOS, D. E., (SNLA)
C***PURPOSE  Evaluates the B-representation of a B-spline at X for the
C            function value or any of its derivatives.
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Reference
C         SIAM J. Numerical Analysis, 14, No. 3, June, 1977, pp.441-472.
C
C     Abstract
C         BVALU is the BVALUE function of the reference.
C
C         BVALU evaluates the B-representation (T,A,N,K) of a B-spline
C         at X for the function value on IDERIV = 0 or any of its
C         derivatives on IDERIV = 1,2,...,K-1.  Right limiting values
C         (right derivatives) are returned except at the right end
C         point X=T(N+1) where left limiting values are computed.  The
C         spline is defined on T(K) .LE. X .LE. T(N+1).  BVALU returns
C         a fatal error message when X is outside of this interval.
C
C         To compute left derivatives or left limiting values at a
C         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
C
C         BVALU calls INTRV
C
C     Description of Arguments
C         Input
C          T       - knot vector of length N+K
C          A       - B-spline coefficient vector of length N
C          N       - number of B-spline coefficients
C                    N = sum of knot multiplicities-K
C          K       - order of the B-spline, K .GE. 1
C          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
C                    IDERIV=0 returns the B-spline value
C          X       - argument, T(K) .LE. X .LE. T(N+1)
C          INBV    - an initialization parameter which must be set
C                    to 1 the first time BVALU is called.
C
C         Output
C          INBV    - INBV contains information for efficient process-
C                    ing after the initial call and INBV must not
C                    be changed by the user.  Distinct splines require
C                    distinct INBV parameters.
C          WORK    - work vector of length 3*K.
C          BVALU   - value of the IDERIV-th derivative at X
C
C     Error Conditions
C         An improper input is a fatal error
C***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
C                 SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
C                 JUNE 1977, PP. 441-472.
C***ROUTINES CALLED  INTRV,XERROR
C***END PROLOGUE  BVALU
C
C
      INTEGER I,IDERIV,IDERP1,IHI,IHMKMJ,ILO,IMK,IMKPJ, INBV, IPJ,
     1 IP1, IP1MJ, J, JJ, J1, J2, K, KMIDER, KMJ, KM1, KPK, MFLAG, N
      REAL A, FKMJ, T, WORK, X
C     DIMENSION T(N+K), WORK(3*K)
      DIMENSION T(*), A(N), WORK(*)
C***FIRST EXECUTABLE STATEMENT  BVALU
      BVALU = 0.0E0
      IF(K.LT.1) GO TO 102
      IF(N.LT.K) GO TO 101
      IF(IDERIV.LT.0 .OR. IDERIV.GE.K) GO TO 110
      KMIDER = K - IDERIV
C
C *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
C     (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
      KM1 = K - 1
      CALL INTRV(T, N+1, X, INBV, I, MFLAG)
      IF (X.LT.T(K)) GO TO 120
      IF (MFLAG.EQ.0) GO TO 20
      IF (X.GT.T(I)) GO TO 130
   10 IF (I.EQ.K) GO TO 140
      I = I - 1
      IF (X.EQ.T(I)) GO TO 10
C
C *** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES
C     WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K
C
   20 IMK = I - K
      DO 30 J=1,K
        IMKPJ = IMK + J
        WORK(J) = A(IMKPJ)
   30 CONTINUE
      IF (IDERIV.EQ.0) GO TO 60
      DO 50 J=1,IDERIV
        KMJ = K - J
        FKMJ = FLOAT(KMJ)
        DO 40 JJ=1,KMJ
          IHI = I + JJ
          IHMKMJ = IHI - KMJ
          WORK(JJ) = (WORK(JJ+1)-WORK(JJ))/(T(IHI)-T(IHMKMJ))*FKMJ
   40   CONTINUE
   50 CONTINUE
C
C *** COMPUTE VALUE AT *X* IN (T(I),(T(I+1)) OF IDERIV-TH DERIVATIVE,
C     GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV).
   60 IF (IDERIV.EQ.KM1) GO TO 100
      IP1 = I + 1
      KPK = K + K
      J1 = K + 1
      J2 = KPK + 1
      DO 70 J=1,KMIDER
        IPJ = I + J
        WORK(J1) = T(IPJ) - X
        IP1MJ = IP1 - J
        WORK(J2) = X - T(IP1MJ)
        J1 = J1 + 1
        J2 = J2 + 1
   70 CONTINUE
      IDERP1 = IDERIV + 1
      DO 90 J=IDERP1,KM1
        KMJ = K - J
        ILO = KMJ
        DO 80 JJ=1,KMJ
          WORK(JJ) = (WORK(JJ+1)*WORK(KPK+ILO)+WORK(JJ)
     1              *WORK(K+JJ))/(WORK(KPK+ILO)+WORK(K+JJ))
          ILO = ILO - 1
   80   CONTINUE
   90 CONTINUE
  100 BVALU = WORK(1)
      RETURN
C
C
  101 CONTINUE
      CALL XERROR( ' BVALU,  N DOES NOT SATISFY N.GE.K',34,2,1)
      RETURN
  102 CONTINUE
      CALL XERROR( ' BVALU,  K DOES NOT SATISFY K.GE.1',34,2,1)
      RETURN
  110 CONTINUE
      CALL XERROR( ' BVALU,  IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K',
     1 49, 2, 1)
      RETURN
  120 CONTINUE
      CALL XERROR( ' BVALU,  X IS N0T GREATER THAN OR EQUAL TO T(K)',
     1 47, 2, 1)
      RETURN
  130 CONTINUE
      CALL XERROR( ' BVALU,  X IS NOT LESS THAN OR EQUAL TO T(N+1)',
     1 46, 2, 1)
      RETURN
  140 CONTINUE
      CALL XERROR( ' BVALU,  A LEFT LIMITING VALUE CANN0T BE OBTAINED AT
     1 T(K)',     57, 2, 1)
      RETURN
      END
      SUBROUTINE INTRV(XT,LXT,X,ILO,ILEFT,MFLAG)
C***BEGIN PROLOGUE  INTRV
C***DATE WRITTEN   800901   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  E3,K6
C***KEYWORDS  B-SPLINE,DATA FITTING,INTERPOLATION,SPLINE
C***AUTHOR  AMOS, D. E., (SNLA)
C***PURPOSE  Computes the largest integer ILEFT in 1.LE.ILEFT.LE.LXT
C            such that XT(ILEFT).LE.X where XT(*) is a subdivision
C            of the X interval.
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Reference
C         SIAM J. Numerical Analysis, 14, No. 3, June 1977, pp. 441-472.
C
C     Abstract
C         INTRV is the INTERV routine of the reference.
C
C         INTRV computes the largest integer ILEFT in 1 .LE. ILEFT .LE.
C         LXT such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
C         the X interval.  Precisely,
C
C                      X .LT. XT(1)                1         -1
C         if  XT(I) .LE. X .LT. XT(I+1)  then  ILEFT=I  , MFLAG=0
C           XT(LXT) .LE. X                         LXT        1,
C
C         That is, when multiplicities are present in the break point
C         to the left of X, the largest index is taken for ILEFT.
C
C     Description of Arguments
C         Input
C          XT      - XT is a knot or break point vector of length LXT
C          LXT     - length of the XT vector
C          X       - argument
C          ILO     - an initialization parameter which must be set
C                    to 1 the first time the spline array XT is
C                    processed by INTRV.
C
C         Output
C          ILO     - ILO contains information for efficient process-
C                    ing after the initial call, and ILO must not be
C                    changed by the user.  Distinct splines require
C                    distinct ILO parameters.
C          ILEFT   - largest integer satisfying XT(ILEFT) .LE. X
C          MFLAG   - signals when X lies out of bounds
C
C     Error Conditions
C         None
C***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
C                 SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
C                 JUNE 1977, PP. 441-472.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  INTRV
C
C
      INTEGER IHI, ILEFT, ILO, ISTEP, LXT, MFLAG, MIDDLE
      REAL X, XT
      DIMENSION XT(LXT)
C***FIRST EXECUTABLE STATEMENT  INTRV
      IHI = ILO + 1
      IF (IHI.LT.LXT) GO TO 10
      IF (X.GE.XT(LXT)) GO TO 110
      IF (LXT.LE.1) GO TO 90
      ILO = LXT - 1
      IHI = LXT
C
   10 IF (X.GE.XT(IHI)) GO TO 40
      IF (X.GE.XT(ILO)) GO TO 100
C
C *** NOW X .LT. XT(IHI) . FIND LOWER BOUND
      ISTEP = 1
   20 IHI = ILO
      ILO = IHI - ISTEP
      IF (ILO.LE.1) GO TO 30
      IF (X.GE.XT(ILO)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 20
   30 ILO = 1
      IF (X.LT.XT(1)) GO TO 90
      GO TO 70
C *** NOW X .GE. XT(ILO) . FIND UPPER BOUND
   40 ISTEP = 1
   50 ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI.GE.LXT) GO TO 60
      IF (X.LT.XT(IHI)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 50
   60 IF (X.GE.XT(LXT)) GO TO 110
      IHI = LXT
C
C *** NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL
   70 MIDDLE = (ILO+IHI)/2
      IF (MIDDLE.EQ.ILO) GO TO 100
C     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1
      IF (X.LT.XT(MIDDLE)) GO TO 80
      ILO = MIDDLE
      GO TO 70
   80 IHI = MIDDLE
      GO TO 70
C *** SET OUTPUT AND RETURN
   90 MFLAG = -1
      ILEFT = 1
      RETURN
  100 MFLAG = 0
      ILEFT = ILO
      RETURN
  110 MFLAG = 1
      ILEFT = LXT
      RETURN
      END
