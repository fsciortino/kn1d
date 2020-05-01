 
      SUBROUTINE B2INK(X,NX,Y,NY,FCN,LDF,KX,KY,TX,TY,BCOEF,WORK,IFLAG)
C***BEGIN PROLOGUE  B2INK
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
C***PURPOSE  B2INK DETERMINES A PIECEWISE POLYNOMIAL FUNCTION THAT
C            INTERPOLATES TWO-DIMENSIONAL GRIDDED DATA. USERS SPECIFY
C            THE POLYNOMIAL ORDER (DEGREE+1) OF THE INTERPOLANT AND
C            (OPTIONALLY) THE KNOT SEQUENCE.
C***DESCRIPTION
C
C   B2INK determines the parameters of a function that interpolates the
C   two-dimensional gridded data (X(i),Y(j),FCN(i,j)) for i=1,..,NX and
C   j=1,..,NY. The  interpolating  function  and  its  derivatives  may
C   subsequently be evaluated by the function B2VAL.
C
C   The interpolating  function  is  a  piecewise  polynomial  function
C   represented as a tensor product of one-dimensional  B-splines.  The
C   form of this function is
C
C                          NX   NY
C              S(x,y)  =  SUM  SUM  a   U (x) V (y)
C                         i=1  j=1   ij  i     j
C
C   where the functions U(i)  and  V(j)  are  one-dimensional  B-spline
C   basis functions. The coefficients a(i,j) are chosen so that
C
C         S(X(i),Y(j)) = FCN(i,j)   for i=1,..,NX and j=1,..,NY
C
C   Note that  for  each  fixed  value  of  y  S(x,y)  is  a  piecewise
C   polynomial function of x alone, and for each fixed value of x  S(x,
C   y) is a piecewise polynomial function of y alone. In one  dimension
C   a piecewise polynomial may  be  created  by  partitioning  a  given
C   interval into subintervals and defining a distinct polynomial piece
C   on each one. The points where adjacent subintervals meet are called
C   knots. Each of the functions U(i) and V(j)  above  is  a  piecewise
C   polynomial.
C
C   Users of B2INK choose the order (degree+1) of the polynomial pieces
C   used to define the piecewise polynomial in each  of  the  x  and  y
C   directions (KX and KY).  Users  also  may  define  their  own  knot
C   sequence in x and y separately (TX and TY).  If  IFLAG=0,  however,
C   B2INK will choose sequences of knots that  result  in  a  piecewise
C   polynomial interpolant with KX-2 continuous partial derivatives  in
C   x and KY-2 continuous partial derivatives in y. (KX knots are taken
C   near each endpoint, not-a-knot end conditions  are  used,  and  the
C   remaining knots are placed at data points  if  KX  is  even  or  at
C   midpoints between data points if KX is  odd.  The  y  direction  is
C   treated similarly.)
C
C   After a call to B2INK, all  information  necessary  to  define  the
C   interpolating function are contained in the parameters NX, NY,  KX,
C   KY, TX, TY, and BCOEF. These quantities should not be altered until
C   after the last call of the evaluation routine B2VAL.
C
C
C   I N P U T
C   ---------
C
C   X       Real 1D array (size NX)
C           Array of x abcissae. Must be strictly increasing.
C
C   NX      Integer scalar (.GE. 3)
C           Number of x abcissae.
C
C   Y       Real 1D array (size NY)
C           Array of y abcissae. Must be strictly increasing.
C
C   NY      Integer scalar (.GE. 3)
C           Number of y abcissae.
C
C   FCN     Real 2D array (size LDF by NY)
C           Array of function values to interpolate. FCN(I,J) should
C           contain the function value at the point (X(I),Y(J))
C
C   LDF     Integer scalar (.GE. NX)
C           The actual leading dimension of FCN used in the calling
C           calling program.
C
C   KX      Integer scalar (.GE. 2, .LT. NX)
C           The order of spline pieces in x.
C           (Order = polynomial degree + 1)
C
C   KY      Integer scalar (.GE. 2, .LT. NY)
C           The order of spline pieces in y.
C           (Order = polynomial degree + 1)
C
C
C   I N P U T   O R   O U T P U T
C   -----------------------------
C
C   TX      Real 1D array (size NX+KX)
C           The knots in the x direction for the spline interpolant.
C           If IFLAG=0 these are chosen by B2INK.
C           If IFLAG=1 these are specified by the user.
C                      (Must be non-decreasing.)
C
C   TY      Real 1D array (size NY+KY)
C           The knots in the y direction for the spline interpolant.
C           If IFLAG=0 these are chosen by B2INK.
C           If IFLAG=1 these are specified by the user.
C                      (Must be non-decreasing.)
C
C
C   O U T P U T
C   -----------
C
C   BCOEF   Real 2D array (size NX by NY)
C           Array of coefficients of the B-spline interpolant.
C           This may be the same array as FCN.
C
C
C   M I S C E L L A N E O U S
C   -------------------------
C
C   WORK    Real 1D array (size NX*NY + max( 2*KX*(NX+1),
C                                  2*KY*(NY+1) ))
C           Array of working storage.
C
C   IFLAG   Integer scalar.
C           On input:  0 == knot sequence chosen by B2INK
C                      1 == knot sequence chosen by user.
C           On output: 1 == successful execution
C                      2 == IFLAG out of range
C                      3 == NX out of range
C                      4 == KX out of range
C                      5 == X not strictly increasing
C                      6 == TX not non-decreasing
C                      7 == NY out of range
C                      8 == KY out of range
C                      9 == Y not strictly increasing
C                     10 == TY not non-decreasing
C
C***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES,
C                 SPRINGER-VERLAG, NEW YORK, 1978.
C               CARL DE BOOR, EFFICIENT COMPUTER MANIPULATION OF TENSOR
C                 PRODUCTS, ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C                 VOL. 5 (1979), PP. 173-182.
C***ROUTINES CALLED  BTPCF,BKNOT
C***END PROLOGUE  B2INK
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  PARAMETERS
C
      INTEGER
     *        NX, NY, LDF, KX, KY, IFLAG
      REAL
     *     X(NX), Y(NY), FCN(LDF,NY), TX(*), TY(*), BCOEF(NX,NY),
     *     WORK(*)
C
C  LOCAL VARIABLES
C
      INTEGER
     *        I, IW, NPK
C
C  -----------------------
C  CHECK VALIDITY OF INPUT
C  -----------------------
C
C***FIRST EXECUTABLE STATEMENT
      IF ((IFLAG .LT. 0) .OR. (IFLAG .GT. 1))  GO TO 920
      IF (NX .LT. 3)  GO TO 930
      IF (NY .LT. 3)  GO TO 970
      IF ((KX .LT. 2) .OR. (KX .GE. NX))  GO TO 940
      IF ((KY .LT. 2) .OR. (KY .GE. NY))  GO TO 980
      DO 10 I=2,NX
         IF (X(I) .LE. X(I-1))  GO TO 950
   10 CONTINUE
      DO 20 I=2,NY
         IF (Y(I) .LE. Y(I-1))  GO TO 990
   20 CONTINUE
      IF (IFLAG .EQ. 0)  GO TO 50
         NPK = NX + KX
         DO 30 I=2,NPK
            IF (TX(I) .LT. TX(I-1))  GO TO 960
   30    CONTINUE
         NPK = NY + KY
         DO 40 I=2,NPK
            IF (TY(I) .LT. TY(I-1))  GO TO 1000
   40    CONTINUE
   50 CONTINUE
C
C  ------------
C  CHOOSE KNOTS
C  ------------
C
      IF (IFLAG .NE. 0)  GO TO 100
         CALL BKNOT(X,NX,KX,TX)
         CALL BKNOT(Y,NY,KY,TY)
  100 CONTINUE
C
C  -------------------------------
C  CONSTRUCT B-SPLINE COEFFICIENTS
C  -------------------------------
C
      IFLAG = 1
      IW = NX*NY + 1
      CALL BTPCF(X,NX,FCN,LDF,NY,TX,KX,WORK,WORK(IW))
      CALL BTPCF(Y,NY,WORK,NY,NX,TY,KY,BCOEF,WORK(IW))
      GO TO 9999
C
C  -----
C  EXITS
C  -----
C
  920 CONTINUE
      CALL XERRWV('B2INK -   IFLAG=I1 IS OUT OF RANGE.',
     *            35,2,1,1,IFLAG,I2,0,R1,R2)
      IFLAG = 2
      GO TO 9999
C
  930 CONTINUE
      IFLAG = 3
      CALL XERRWV('B2INK -   NX=I1 IS OUT OF RANGE.',
     *            32,IFLAG,1,1,NX,I2,0,R1,R2)
      GO TO 9999
C
  940 CONTINUE
      IFLAG = 4
      CALL XERRWV('B2INK -   KX=I1 IS OUT OF RANGE.',
     *            32,IFLAG,1,1,KX,I2,0,R1,R2)
      GO TO 9999
C
  950 CONTINUE
      IFLAG = 5
      CALL XERRWV('B2INK -   X ARRAY MUST BE STRICTLY INCREASING.',
     *            46,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
C
  960 CONTINUE
      IFLAG = 6
      CALL XERRWV('B2INK -   TX ARRAY MUST BE NON-DECREASING.',
     *            42,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
C
  970 CONTINUE
      IFLAG = 7
      CALL XERRWV('B2INK -   NY=I1 IS OUT OF RANGE.',
     *            32,IFLAG,1,1,NY,I2,0,R1,R2)
      GO TO 9999
C
  980 CONTINUE
      IFLAG = 8
      CALL XERRWV('B2INK -   KY=I1 IS OUT OF RANGE.',
     *            32,IFLAG,1,1,KY,I2,0,R1,R2)
      GO TO 9999
C
  990 CONTINUE
      IFLAG = 9
      CALL XERRWV('B2INK -   Y ARRAY MUST BE STRICTLY INCREASING.',
     *            46,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
C
 1000 CONTINUE
      IFLAG = 10
      CALL XERRWV('B2INK -   TY ARRAY MUST BE NON-DECREASING.',
     *            42,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
C
 9999 CONTINUE
      RETURN
      END
      SUBROUTINE BKNOT(X,N,K,T)
C***BEGIN PROLOGUE  BKNOT
C***REFER TO  B2INK,B3INK
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***END PROLOGUE  BKNOT
C
C  --------------------------------------------------------------------
C  BKNOT CHOOSES A KNOT SEQUENCE FOR INTERPOLATION OF ORDER K AT THE
C  DATA POINTS X(I), I=1,..,N.  THE N+K KNOTS ARE PLACED IN THE ARRAY
C  T.  K KNOTS ARE PLACED AT EACH ENDPOINT AND NOT-A-KNOT END
C  CONDITIONS ARE USED.  THE REMAINING KNOTS ARE PLACED AT DATA POINTS
C  IF N IS EVEN AND BETWEEN DATA POINTS IF N IS ODD.  THE RIGHTMOST
C  KNOT IS SHIFTED SLIGHTLY TO THE RIGHT TO INSURE PROPER INTERPOLATION
C  AT X(N) (SEE PAGE 350 OF THE REFERENCE).
C  --------------------------------------------------------------------
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  PARAMETERS
C
      INTEGER
     *        N, K
      REAL
     *     X(N), T(*)
C
C  LOCAL VARIABLES
C
      INTEGER
     *        I, J, IPJ, NPJ, IP1
      REAL
     *     RNOT
C
C
C  ----------------------------
C  PUT K KNOTS AT EACH ENDPOINT
C  ----------------------------
C
C     (SHIFT RIGHT ENPOINTS SLIGHTLY -- SEE PG 350 OF REFERENCE)
      RNOT = X(N) + 0.10E0*( X(N)-X(N-1) )
      DO 110 J=1,K
         T(J) = X(1)
         NPJ = N + J
         T(NPJ) = RNOT
  110 CONTINUE
C
C  --------------------------
C  DISTRIBUTE REMAINING KNOTS
C  --------------------------
C
      IF (MOD(K,2) .EQ. 1)  GO TO 150
C
C     CASE OF EVEN K --  KNOTS AT DATA POINTS
C
      I = (K/2) - K
      JSTRT = K+1
      DO 120 J=JSTRT,N
         IPJ = I + J
         T(J) = X(IPJ)
  120 CONTINUE
      GO TO 200
C
C     CASE OF ODD K --  KNOTS BETWEEN DATA POINTS
C
  150 CONTINUE
      I = (K-1)/2 - K
      IP1 = I + 1
      JSTRT = K + 1
      DO 160 J=JSTRT,N
         IPJ = I + J
         T(J) = 0.50E0*( X(IPJ) + X(IPJ+1) )
  160 CONTINUE
  200 CONTINUE
C
      RETURN
      END
      SUBROUTINE BTPCF(X,N,FCN,LDF,NF,T,K,BCOEF,WORK)
C***BEGIN PROLOGUE  BTPCF
C***REFER TO  B2INK,B3INK
C***ROUTINES CALLED  BINTK,BNSLV
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***END PROLOGUE  BTPCF
C
C  -----------------------------------------------------------------
C  BTPCF COMPUTES B-SPLINE INTERPOLATION COEFFICIENTS FOR NF SETS
C  OF DATA STORED IN THE COLUMNS OF THE ARRAY FCN. THE B-SPLINE
C  COEFFICIENTS ARE STORED IN THE ROWS OF BCOEF HOWEVER.
C  EACH INTERPOLATION IS BASED ON THE N ABCISSA STORED IN THE
C  ARRAY X, AND THE N+K KNOTS STORED IN THE ARRAY T. THE ORDER
C  OF EACH INTERPOLATION IS K. THE WORK ARRAY MUST BE OF LENGTH
C  AT LEAST 2*K*(N+1).
C  -----------------------------------------------------------------
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  PARAMETERS
C
      INTEGER
     *        N, LDF, K
      REAL
     *     X(N), FCN(LDF,NF), T(*), BCOEF(NF,N), WORK(*)
C
C  LOCAL VARIABLES
C
      INTEGER
     *        I, J, K1, K2, IQ, IW
C
C  ---------------------------------------------
C  CHECK FOR NULL INPUT AND PARTITION WORK ARRAY
C  ---------------------------------------------
C
C***FIRST EXECUTABLE STATEMENT
      IF (NF .LE. 0)  GO TO 500
      K1 = K - 1
      K2 = K1 + K
      IQ = 1 + N
      IW = IQ + K2*N+1
C
C  -----------------------------
C  COMPUTE B-SPLINE COEFFICIENTS
C  -----------------------------
C
C
C   FIRST DATA SET
C
      CALL BINTK(X,FCN,T,N,K,WORK,WORK(IQ),WORK(IW))
      DO 20 I=1,N
         BCOEF(1,I) = WORK(I)
   20 CONTINUE
C
C  ALL REMAINING DATA SETS BY BACK-SUBSTITUTION
C
      IF (NF .EQ. 1)  GO TO 500
      DO 100 J=2,NF
         DO 50 I=1,N
            WORK(I) = FCN(I,J)
   50    CONTINUE
         CALL BNSLV(WORK(IQ),K2,N,K1,K1,WORK)
         DO 60 I=1,N
            BCOEF(J,I) = WORK(I)
   60    CONTINUE
  100 CONTINUE
C
C  ----
C  EXIT
C  ----
C
  500 CONTINUE
      RETURN
      END

      SUBROUTINE BINTK(X,Y,T,N,K,BCOEF,Q,WORK)
C***BEGIN PROLOGUE  BINTK
C***DATE WRITTEN   800901   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  E1A
C***KEYWORDS  B-SPLINE,DATA FITTING,INTERPOLATION,SPLINE
C***AUTHOR  AMOS, D. E., (SNLA)
C***PURPOSE  Produces the B-spline coefficients, BCOEF, of the
C            B-spline of order K with knots T(I), I=1,...,N+K, which
C            takes on the value Y(I) at X(I), I=1,...,N.
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     References
C
C          A Practical Guide to Splines by C. de Boor, Applied
C          Mathematics Series 27, Springer, 1979.
C
C     Abstract
C
C         BINTK is the SPLINT routine of the reference.
C
C         BINTK produces the B-spline coefficients, BCOEF, of the
C         B-spline of order K with knots T(I), I=1,...,N+K, which
C         takes on the value Y(I) at X(I), I=1,...,N.  The spline or
C         any of its derivatives can be evaluated by calls to BVALU.
C         The I-th equation of the linear system A*BCOEF = B for the
C         coefficients of the interpolant enforces interpolation at
C         X(I)), I=1,...,N.  Hence, B(I) = Y(I), all I, and A is
C         a band matrix with 2K-1 bands if A is invertible. The matrix
C         A is generated row by row and stored, diagonal by diagonal,
C         in the rows of Q, with the main diagonal going into row K.
C         The banded system is then solved by a call to BNFAC (which
C         constructs the triangular factorization for A and stores it
C         again in Q), followed by a call to BNSLV (which then
C         obtains the solution BCOEF by substitution). BNFAC does no
C         pivoting, since the total positivity of the matrix A makes
C         this unnecessary.  The linear system to be solved is
C         (theoretically) invertible if and only if
C                 T(I) .LT. X(I)) .LT. T(I+K),        all I.
C         Equality is permitted on the left for I=1 and on the right
C         for I=N when K knots are used at X(1) or X(N).  Otherwise,
C         violation of this condition is certain to lead to an error.
C
C         BINTK calls BSPVN, BNFAC, BNSLV, XERROR
C
C     Description of Arguments
C         Input
C           X       - vector of length N containing data point abscissa
C                     in strictly increasing order.
C           Y       - corresponding vector of length N containing data
C                     point ordinates.
C           T       - knot vector of length N+K
C                     since T(1),..,T(K) .LE. X(1) and T(N+1),..,T(N+K)
C                     .GE. X(N), this leaves only N-K knots (not nec-
C                     essarily X(I)) values) interior to (X(1),X(N))
C           N       - number of data points, N .GE. K
C           K       - order of the spline, K .GE. 1
C
C         Output
C           BCOEF   - a vector of length N containing the B-spline
C                     coefficients
C           Q       - a work vector of length (2*K-1)*N, containing
C                     the triangular factorization of the coefficient
C                     matrix of the linear system being solved.  The
C                     coefficients for the interpolant of an
C                     additional data set (X(I)),YY(I)), I=1,...,N
C                     with the same abscissa can be obtained by loading
C                     YY into BCOEF and then executing
C                         call BNSLV(Q,2K-1,N,K-1,K-1,BCOEF)
C           WORK    - work vector of length 2*K
C
C     Error Conditions
C         Improper  input is a fatal error
C         Singular system of equations is a fatal error
C***REFERENCES  D.E. AMOS, *COMPUTATION WITH SPLINES AND B-SPLINES*,
C                 SAND78-1968,SANDIA LABORATORIES,MARCH,1979.
C               C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
C                 SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
C                 JUNE 1977, PP. 441-472.
C               C. DE BOOR, *A PRACTICAL GUIDE TO SPLINES*, APPLIED
C                 MATHEMATICS SERIES 27, SPRINGER, 1979.
C***ROUTINES CALLED  BNFAC,BNSLV,BSPVN,XERROR
C***END PROLOGUE  BINTK
C
C
      INTEGER IFLAG, IWORK, K, N, I, ILP1MX, J, JJ, KM1, KPKM2, LEFT,
     1 LENQ, NP1
      REAL BCOEF(N), Y(N), Q(*), T(*), X(N), XI, WORK(*)
C     DIMENSION Q(2*K-1,N), T(N+K)
C***FIRST EXECUTABLE STATEMENT  BINTK
      IF(K.LT.1) GO TO 100
      IF(N.LT.K) GO TO 105
      JJ = N - 1
      IF(JJ.EQ.0) GO TO 6
      DO 5 I=1,JJ
      IF(X(I).GE.X(I+1)) GO TO 110
    5 CONTINUE
    6 CONTINUE
      NP1 = N + 1
      KM1 = K - 1
      KPKM2 = 2*KM1
      LEFT = K
C                ZERO OUT ALL ENTRIES OF Q
      LENQ = N*(K+KM1)
      DO 10 I=1,LENQ
        Q(I) = 0.0E0
   10 CONTINUE
C
C  ***   LOOP OVER I TO CONSTRUCT THE  N  INTERPOLATION EQUATIONS
      DO 50 I=1,N
        XI = X(I)
        ILP1MX = MIN0(I+K,NP1)
C        *** FIND  LEFT  IN THE CLOSED INTERVAL (I,I+K-1) SUCH THAT
C                T(LEFT) .LE. X(I) .LT. T(LEFT+1)
C        MATRIX IS SINGULAR IF THIS IS NOT POSSIBLE
        LEFT = MAX0(LEFT,I)
        IF (XI.LT.T(LEFT)) GO TO 80
   20   IF (XI.LT.T(LEFT+1)) GO TO 30
        LEFT = LEFT + 1
        IF (LEFT.LT.ILP1MX) GO TO 20
        LEFT = LEFT - 1
        IF (XI.GT.T(LEFT+1)) GO TO 80
C        *** THE I-TH EQUATION ENFORCES INTERPOLATION AT XI, HENCE
C        A(I,J) = B(J,K,T)(XI), ALL J. ONLY THE  K  ENTRIES WITH  J =
C        LEFT-K+1,...,LEFT ACTUALLY MIGHT BE NONZERO. THESE  K  NUMBERS
C        ARE RETURNED, IN  BCOEF (USED FOR TEMP.STORAGE HERE), BY THE
C        FOLLOWING
   30   CALL BSPVN(T, K, K, 1, XI, LEFT, BCOEF, WORK, IWORK)
C        WE THEREFORE WANT  BCOEF(J) = B(LEFT-K+J)(XI) TO GO INTO
C        A(I,LEFT-K+J), I.E., INTO  Q(I-(LEFT+J)+2*K,(LEFT+J)-K) SINCE
C        A(I+J,J)  IS TO GO INTO  Q(I+K,J), ALL I,J,  IF WE CONSIDER  Q
C        AS A TWO-DIM. ARRAY , WITH  2*K-1  ROWS (SEE COMMENTS IN
C        BNFAC). IN THE PRESENT PROGRAM, WE TREAT  Q  AS AN EQUIVALENT
C        ONE-DIMENSIONAL ARRAY (BECAUSE OF FORTRAN RESTRICTIONS ON
C        DIMENSION STATEMENTS) . WE THEREFORE WANT  BCOEF(J) TO GO INTO
C        ENTRY
C            I -(LEFT+J) + 2*K + ((LEFT+J) - K-1)*(2*K-1)
C                   =  I-LEFT+1 + (LEFT -K)*(2*K-1) + (2*K-2)*J
C        OF  Q .
        JJ = I - LEFT + 1 + (LEFT-K)*(K+KM1)
        DO 40 J=1,K
          JJ = JJ + KPKM2
          Q(JJ) = BCOEF(J)
   40   CONTINUE
   50 CONTINUE
C
C     ***OBTAIN FACTORIZATION OF  A  , STORED AGAIN IN  Q.
      CALL BNFAC(Q, K+KM1, N, KM1, KM1, IFLAG)
      GO TO (60, 90), IFLAG
C     *** SOLVE  A*BCOEF = Y  BY BACKSUBSTITUTION
   60 DO 70 I=1,N
        BCOEF(I) = Y(I)
   70 CONTINUE
      CALL BNSLV(Q, K+KM1, N, KM1, KM1, BCOEF)
      RETURN
C
C
   80 CONTINUE
      CALL XERROR(  'BINTK,  SOME ABSCISSA WAS NOT IN THE SUPPORT OF THE
     1 CORRESPONDING BASIS FUNCTION AND THE SYSTEM IS SINGULAR.',108,2,1
     2)
      RETURN
   90 CONTINUE
      CALL XERROR( 'BINTK,  THE SYSTEM OF SOLVER DETECTS A SINGULAR SYST
     1EM ALTHOUGH THE THEORETICAL CONDITIONS FOR A SOLUTION WERE SATISFI
     2ED.',121,8,1)
      RETURN
  100 CONTINUE
      CALL XERROR( ' BINTK,  K DOES NOT SATISFY K.GE.1', 34, 2, 1)
      RETURN
  105 CONTINUE
      CALL XERROR( ' BINTK,  N DOES NOT SATISFY N.GE.K', 34, 2, 1)
      RETURN
  110 CONTINUE
      CALL XERROR( ' BINTK,  X(I) DOES NOT SATISFY X(I).LT.X(I+1) FOR SO
     1ME I', 56, 2, 1)
      RETURN
      END
      SUBROUTINE BNFAC(W,NROWW,NROW,NBANDL,NBANDU,IFLAG)
C***BEGIN PROLOGUE  BNFAC
C***REFER TO  BINT4,BINTK
C
C  BNFAC is the BANFAC routine from
C        * A Practical Guide to Splines *  by C. de Boor
C
C  Returns in  W  the lu-factorization (without pivoting) of the banded
C  matrix  A  of order  NROW  with  (NBANDL + 1 + NBANDU) bands or diag-
C  onals in the work array  W .
C
C *****  I N P U T  ******
C  W.....Work array of size  (NROWW,NROW)  containing the interesting
C        part of a banded matrix  A , with the diagonals or bands of  A
C        stored in the rows of  W , while columns of  A  correspond to
C        columns of  W . This is the storage mode used in  LINPACK  and
C        results in efficient innermost loops.
C           Explicitly,  A  has  NBANDL  bands below the diagonal
C                            +     1     (main) diagonal
C                            +   NBANDU  bands above the diagonal
C        and thus, with    MIDDLE = NBANDU + 1,
C          A(I+J,J)  is in  W(I+MIDDLE,J)  for I=-NBANDU,...,NBANDL
C                                              J=1,...,NROW .
C        For example, the interesting entries of A (1,2)-banded matrix
C        of order  9  would appear in the first  1+1+2 = 4  rows of  W
C        as follows.
C                          13 24 35 46 57 68 79
C                       12 23 34 45 56 67 78 89
C                    11 22 33 44 55 66 77 88 99
C                    21 32 43 54 65 76 87 98
C
C        All other entries of  W  not identified in this way with an en-
C        try of  A  are never referenced .
C  NROWW.....Row dimension of the work array  W .
C        must be  .GE.  NBANDL + 1 + NBANDU  .
C  NBANDL.....Number of bands of  A  below the main diagonal
C  NBANDU.....Number of bands of  A  above the main diagonal .
C
C *****  O U T P U T  ******
C  IFLAG.....Integer indicating success( = 1) or failure ( = 2) .
C     If  IFLAG = 1, then
C  W.....contains the LU-factorization of  A  into a unit lower triangu-
C        lar matrix  L  and an upper triangular matrix  U (both banded)
C        and stored in customary fashion over the corresponding entries
C        of  A . This makes it possible to solve any particular linear
C        system  A*X = B  for  X  by A
C              CALL BNSLV ( W, NROWW, NROW, NBANDL, NBANDU, B )
C        with the solution X  contained in  B  on return .
C     If  IFLAG = 2, then
C        one of  NROW-1, NBANDL,NBANDU failed to be nonnegative, or else
C        one of the potential pivots was found to be zero indicating
C        that  A  does not have an LU-factorization. This implies that
C        A  is singular in case it is totally positive .
C
C *****  M E T H O D  ******
C     Gauss elimination  W I T H O U T  pivoting is used. The routine is
C  intended for use with matrices  A  which do not require row inter-
C  changes during factorization, especially for the  T O T A L L Y
C  P O S I T I V E  matrices which occur in spline calculations.
C     The routine should not be used for an arbitrary banded matrix.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  BNFAC
C
      INTEGER IFLAG, NBANDL, NBANDU, NROW, NROWW, I, IPK, J, JMAX, K,
     1 KMAX, MIDDLE, MIDMK, NROWM1
      REAL W(NROWW,NROW), FACTOR, PIVOT
C
C***FIRST EXECUTABLE STATEMENT  BNFAC
      IFLAG = 1
      MIDDLE = NBANDU + 1
C                         W(MIDDLE,.) CONTAINS THE MAIN DIAGONAL OF  A .
      NROWM1 = NROW - 1
      IF (NROWM1) 120, 110, 10
   10 IF (NBANDL.GT.0) GO TO 30
C                A IS UPPER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO .
      DO 20 I=1,NROWM1
        IF (W(MIDDLE,I).EQ.0.0E0) GO TO 120
   20 CONTINUE
      GO TO 110
   30 IF (NBANDU.GT.0) GO TO 60
C              A IS LOWER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO AND
C                 DIVIDE EACH COLUMN BY ITS DIAGONAL .
      DO 50 I=1,NROWM1
        PIVOT = W(MIDDLE,I)
        IF (PIVOT.EQ.0.0E0) GO TO 120
        JMAX = MIN0(NBANDL,NROW-I)
        DO 40 J=1,JMAX
          W(MIDDLE+J,I) = W(MIDDLE+J,I)/PIVOT
   40   CONTINUE
   50 CONTINUE
      RETURN
C
C        A  IS NOT JUST A TRIANGULAR MATRIX. CONSTRUCT LU FACTORIZATION
   60 DO 100 I=1,NROWM1
C                                  W(MIDDLE,I)  IS PIVOT FOR I-TH STEP .
        PIVOT = W(MIDDLE,I)
        IF (PIVOT.EQ.0.0E0) GO TO 120
C                 JMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN COLUMN  I
C                     BELOW THE DIAGONAL .
        JMAX = MIN0(NBANDL,NROW-I)
C              DIVIDE EACH ENTRY IN COLUMN  I  BELOW DIAGONAL BY PIVOT .
        DO 70 J=1,JMAX
          W(MIDDLE+J,I) = W(MIDDLE+J,I)/PIVOT
   70   CONTINUE
C                 KMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN ROW  I  TO
C                     THE RIGHT OF THE DIAGONAL .
        KMAX = MIN0(NBANDU,NROW-I)
C                  SUBTRACT  A(I,I+K)*(I-TH COLUMN) FROM (I+K)-TH COLUMN
C                  (BELOW ROW  I ) .
        DO 90 K=1,KMAX
          IPK = I + K
          MIDMK = MIDDLE - K
          FACTOR = W(MIDMK,IPK)
          DO 80 J=1,JMAX
            W(MIDMK+J,IPK) = W(MIDMK+J,IPK) - W(MIDDLE+J,I)*FACTOR
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
C                                       CHECK THE LAST DIAGONAL ENTRY .
  110 IF (W(MIDDLE,NROW).NE.0.0E0) RETURN
  120 IFLAG = 2
      RETURN
      END
      SUBROUTINE BNSLV(W,NROWW,NROW,NBANDL,NBANDU,B)
C***BEGIN PROLOGUE  BNSLV
C***REFER TO  BINT4,BINTK
C
C  BNSLV is the BANSLV routine from
C        * A Practical Guide to Splines *  by C. de Boor
C
C  Companion routine to  BNFAC . It returns the solution  X  of the
C  linear system  A*X = B  in place of  B , given the LU-factorization
C  for  A  in the work array  W from BNFAC.
C
C *****  I N P U T  ******
C  W, NROWW,NROW,NBANDL,NBANDU.....describe the LU-factorization of a
C        banded matrix  A  of order  NROW  as constructed in  BNFAC .
C        For details, see  BNFAC .
C  B.....Right side of the system to be solved .
C
C *****  O U T P U T  ******
C  B.....Contains the solution  X , of order  NROW .
C
C *****  M E T H O D  ******
C     (With  A = L*U, as stored in  W,) the unit lower triangular system
C  L(U*X) = B  is solved for  Y = U*X, and  Y  stored in  B . Then the
C  upper triangular system  U*X = Y  is solved for  X  . The calcul-
C  ations are so arranged that the innermost loops stay within columns.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  BNSLV
C
      INTEGER NBANDL, NBANDU, NROW, NROWW, I, J, JMAX, MIDDLE, NROWM1
      REAL W(NROWW,NROW), B(NROW)
C***FIRST EXECUTABLE STATEMENT  BNSLV
      MIDDLE = NBANDU + 1
      IF (NROW.EQ.1) GO TO 80
      NROWM1 = NROW - 1
      IF (NBANDL.EQ.0) GO TO 30
C                                 FORWARD PASS
C            FOR I=1,2,...,NROW-1, SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN
C            OF  L )  FROM RIGHT SIDE  (BELOW I-TH ROW) .
      DO 20 I=1,NROWM1
        JMAX = MIN0(NBANDL,NROW-I)
        DO 10 J=1,JMAX
          B(I+J) = B(I+J) - B(I)*W(MIDDLE+J,I)
   10   CONTINUE
   20 CONTINUE
C                                 BACKWARD PASS
C            FOR I=NROW,NROW-1,...,1, DIVIDE RIGHT SIDE(I) BY I-TH DIAG-
C            ONAL ENTRY OF  U, THEN SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN
C            OF  U)  FROM RIGHT SIDE  (ABOVE I-TH ROW).
   30 IF (NBANDU.GT.0) GO TO 50
C                                A  IS LOWER TRIANGULAR .
      DO 40 I=1,NROW
        B(I) = B(I)/W(1,I)
   40 CONTINUE
      RETURN
   50 I = NROW
   60 B(I) = B(I)/W(MIDDLE,I)
      JMAX = MIN0(NBANDU,I-1)
      DO 70 J=1,JMAX
        B(I-J) = B(I-J) - B(I)*W(MIDDLE-J,I)
   70 CONTINUE
      I = I - 1
      IF (I.GT.1) GO TO 60
   80 B(1) = B(1)/W(MIDDLE,1)
      RETURN
      END
      SUBROUTINE BSPVN(T,JHIGH,K,INDEX,X,ILEFT,VNIKX,WORK,IWORK)
C***BEGIN PROLOGUE  BSPVN
C***DATE WRITTEN   800901   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  E3,K6
C***KEYWORDS  B-SPLINE,DATA FITTING,INTERPOLATION,SPLINE
C***AUTHOR  AMOS, D. E., (SNLA)
C***PURPOSE  Calculates the value of all (possibly) nonzero basis
C            functions at X.
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Reference
C         SIAM J. Numerical Analysis, 14, No. 3, June, 1977, pp.441-472.
C
C     Abstract
C         BSPVN is the BSPLVN routine of the reference.
C
C         BSPVN calculates the value of all (possibly) nonzero basis
C         functions at X of order MAX(JHIGH,(J+1)*(INDEX-1)), where
C         T(K) .LE. X .LE. T(N+1) and J=IWORK is set inside the routine
C         on the first call when INDEX=1.  ILEFT is such that T(ILEFT)
C         .LE. X .LT. T(ILEFT+1).  A call to INTRV(T,N+1,X,ILO,ILEFT,
C         MFLAG) produces the proper ILEFT.  BSPVN calculates using the
C         basic algorithm needed in BSPVD.  If only basis functions are
C         desired, setting JHIGH=K and INDEX=1 can be faster than
C         calling BSPVD, but extra coding is required for derivatives
C         (INDEX=2) and BSPVD is set up for this purpose.
C
C         Left limiting values are set up as described in BSPVD.
C
C     Description of Arguments
C         Input
C          T       - knot vector of length N+K, where
C                    N = number of B-spline basis functions
C                    N = sum of knot multiplicities-K
C          JHIGH   - order of B-spline, 1 .LE. JHIGH .LE. K
C          K       - highest possible order
C          INDEX   - INDEX = 1 gives basis functions of order JHIGH
C                          = 2 denotes previous entry with WORK, IWORK
C                              values saved for subsequent calls to
C                              BSPVN.
C          X       - argument of basis functions,
C                    T(K) .LE. X .LE. T(N+1)
C          ILEFT   - largest integer such that
C                    T(ILEFT) .LE. X .LT. T(ILEFT+1)
C
C         Output
C          VNIKX   - vector of length K for spline values.
C          WORK    - a work vector of length 2*K
C          IWORK   - a work parameter.  Both WORK and IWORK contain
C                    information necessary to continue for INDEX = 2.
C                    When INDEX = 1 exclusively, these are scratch
C                    variables and can be used for other purposes.
C
C     Error Conditions
C         Improper input is a fatal error.
C***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
C                 SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
C                 JUNE 1977, PP. 441-472.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  BSPVN
C
C
      INTEGER ILEFT, IMJP1, INDEX, IPJ, IWORK, JHIGH, JP1, JP1ML, K, L
      REAL T, VM, VMPREV, VNIKX, WORK, X
C     DIMENSION T(ILEFT+JHIGH)
      DIMENSION T(*), VNIKX(K), WORK(*)
C     CONTENT OF J, DELTAM, DELTAP IS EXPECTED UNCHANGED BETWEEN CALLS.
C     WORK(I) = DELTAP(I), WORK(K+I) = DELTAM(I), I = 1,K
C***FIRST EXECUTABLE STATEMENT  BSPVN
      IF(K.LT.1) GO TO 90
      IF(JHIGH.GT.K .OR. JHIGH.LT.1) GO TO 100
      IF(INDEX.LT.1 .OR. INDEX.GT.2) GO TO 105
      IF(X.LT.T(ILEFT) .OR. X.GT.T(ILEFT+1)) GO TO 110
      GO TO (10, 20), INDEX
   10 IWORK = 1
      VNIKX(1) = 1.0E0
      IF (IWORK.GE.JHIGH) GO TO 40
C
   20 IPJ = ILEFT + IWORK
      WORK(IWORK) = T(IPJ) - X
      IMJP1 = ILEFT - IWORK + 1
      WORK(K+IWORK) = X - T(IMJP1)
      VMPREV = 0.0E0
      JP1 = IWORK + 1
      DO 30 L=1,IWORK
        JP1ML = JP1 - L
        VM = VNIKX(L)/(WORK(L)+WORK(K+JP1ML))
        VNIKX(L) = VM*WORK(L) + VMPREV
        VMPREV = VM*WORK(K+JP1ML)
   30 CONTINUE
      VNIKX(JP1) = VMPREV
      IWORK = JP1
      IF (IWORK.LT.JHIGH) GO TO 20
C
   40 RETURN
C
C
   90 CONTINUE
      CALL XERROR( ' BSPVN,  K DOES NOT SATISFY K.GE.1', 34, 2, 1)
      RETURN
  100 CONTINUE
      CALL XERROR( ' BSPVN,  JHIGH DOES NOT SATISFY 1.LE.JHIGH.LE.K',
     1 47, 2, 1)
      RETURN
  105 CONTINUE
      CALL XERROR( ' BSPVN,  INDEX IS NOT 1 OR 2',28,2,1)
      RETURN
  110 CONTINUE
      CALL XERROR( ' BSPVN,  X DOES NOT SATISFY T(ILEFT).LE.X.LE.T(ILEFT
     1+1)', 55, 2, 1)
      RETURN
      END
