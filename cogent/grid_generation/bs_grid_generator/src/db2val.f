CC+------------------------------------------------------------------+CC
CC+------------------------------------------------------------------+CC
CC+                                                                  +CC
CC+              Lawrence Livermore National Laboratory              +CC
CC+                                                                  +CC
CC+        Livermore Computing   Mathematical Software Library       +CC
CC+           Mathematical Software Service Library (MSSL)           +CC
CC+                                                                  +CC
CC+------------------------------------------------------------------+CC
CC+                                                                  +CC
CC+  Each MSSL routine has been tested to some extent and            +CC
CC+  meets reasonable documentation and programming standards.       +CC
CC+                                                                  +CC
CC+  These routines are distributed exclusively for use in support   +CC
CC+  of LLNL programs.  Check with the LLNL Code Release Center or   +CC
CC+  the LC Client Services HotLine, (925)422-4531, before moving    +CC
CC+  this source code to a non-LLNL system.                          +CC
CC+                                                                  +CC
CC+  +------------------------------------------------------------+  +CC
CC+  +                        N O T I C E                         +  +CC
CC+  +------------------------------------------------------------+  +CC
CC+  +  This report was prepared as an account of work sponsored  +  +CC
CC+  +  by the United States government.  Neither the United      +  +CC
CC+  +  States government nor any of their employees, nor any of  +  +CC
CC+  +  their contractors, subcontractors, or their employees,    +  +CC
CC+  +  makes any warranty, expressed or implied, or assumes any  +  +CC
CC+  +  legal liability or responsibility for the accuracy,       +  +CC
CC+  +  completeness or usefulness of any information, apparatus, +  +CC
CC+  +  product or process disclosed, or represents that its use  +  +CC
CC+  +  would not infringe privately-owned rights.                +  +CC
CC+  +------------------------------------------------------------+  +CC
CC+                                                                  +CC
CC+  Please report any suspected errors in this routine to the LC    +CC
CC+  Client Services Hotline, (925)422-4531.                         +CC
CC+                                                                  +CC
CC+------------------------------------------------------------------+CC
CC+------------------------------------------------------------------+CC
*DECK DB2VAL
      DOUBLE PRECISION FUNCTION DB2VAL (XVAL, YVAL, IDX, IDY, TX, TY,
     *   NX, NY, KX, KY, FCN, LDF, ICONT, IWORK, WORK, IFLAG)
C***BEGIN PROLOGUE  DB2VAL
C***PURPOSE  DB2VAL evaluates the two-dimensional piecewise polynomial
C            interpolating function constructed by the routine DB2INT.
C            Either function values or partial derivative values may be
C            be requested.
C***LIBRARY   MSSL
C     >> Revised or new BSPLINE routine, provisionally accepted for SLATEC 5.0.
C***CATEGORY  E3, K6
C***TYPE      DOUBLE PRECISION (B2VAL-S, DB2VAL-D)
C***KEYWORDS  TWO DIMENSIONS, GRIDDED DATA, SPLINES,
C             PIECEWISE POLYNOMIALS, EVALUATION, DIFFERENTIATION
C***AUTHOR  Boisvert, R. F., (NIST)
C             Computing and Applied Mathematics Laboratory
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   DB2VAL evaluates the tensor product piecewise polynomial interpolant
C   constructed by the routine DB2INT or one of its derivatives at the
C   point (XVAL,YVAL).  The variables IDX and IDY indicate which
C   function is to be evaluated.  If B(x,y) is the interpolant,
C
C                            IDX+IDY
C                           d
C             DB2VAL  =   ----------- B (XVAL,YVAL)
C                           IDX   IDY
C                         dx    dy
C
C   Thus, to evaluate the interpolant itself, set IDX=IDY=0. To get the
C   first partial derivative with respect to x, set IDX=1 and IDY=0, and
C   so on.
C
C   Since B is a piecewise polynomial of degree KX-1 in x and KY-1 in y,
C   DB2VAL returns zero whenever IDX.GE.KX or IDY.GE.KY.  DB2VAL also
C   returns zero if (XVAL,YVAL) is out of range, that is, if
C
C            XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR.
C            YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY).
C
C   If the knots TX and TY were chosen by DB2INT, then this is
C   equivalent to
C
C              XVAL.LT.X(1) .OR. XVAL.GT.X(NX) .OR.
C              YVAL.LT.Y(1) .OR. YVAL.GT.Y(NY).
C
C   The input quantities TX, TY, NX, NY, KX, KY, and FCN should be
C   unchanged since the call to DB2INT which produced FCN.
C
C   Note that the derivative functions computed by DB2VAL may be
C   discontinuous when x or y correspond to knots.  (When DB2INT selects
C   knots this occurs only when IDX=KX-1 or IDY=KY-1.)  In these cases
C   DB2VAL returns right limiting values (right derivatives), except at
C   the rightmost knot where left limiting values are returned.
C
C
C   I N P U T
C   ---------
C
C   XVAL    Double precision scalar
C           X coordinate of evaluation point.
C
C   YVAL    Double precision scalar
C           Y coordinate of evaluation point.
C
C   IDX     Integer scalar  (IDX .GE. 0)
C           Indicates the x derivative of piecewise polynomial to
C           evaluate: IDX=J for the Jth partial derivative with respect
C           to x. DB2VAL will return 0 if IDX.GE.KX.
C
C   IDY     Integer scalar  (IDY .GE. 0)
C           Indicates the y derivative of piecewise polynomial to
C           evaluate: IDY=J for the Jth partial derivative with respect
C           to y. DB2VAL will return 0 if IDY.GE.KY.
C
C   TX      Double precision 1D array (size NX+KX)
C           Sequence of knots defining the piecewise polynomial in
C           the x direction.  (Unchanged since the call to DB2INT which
C           produced FCN.)
C
C   TY      Double precision 1D array (size NY+KY)
C           Sequence of knots defining the piecewise polynomial in
C           the y direction.  (Unchanged since the call to DB2INT which
C           produced FCN.)
C
C   NX      Integer scalar (NX .GE. KX)
C           The number of interpolation points in x. (Unchanged since
C           the call to DB2INT which produced FCN.)
C
C   NY      Integer scalar (NY .GE. KY)
C           The number of interpolation points in y. (Unchanged since
C           the call to DB2INT which produced FCN.)
C
C   KX      Integer scalar (KX .GE. 2)
C           Order of polynomial pieces in x. (Unchanged since the call
C           to DB2INT which produced FCN.)
C
C   KY      Integer scalar (KY .GE. 2)
C           Order of polynomial pieces in y. (Unchanged since the call
C           to DB2INT which produced FCN.)
C
C   FCN     Double precision 2D array (size LDF by NY)
C           The B-spline coefficients computed by DB2INT.
C
C   LDF     Integer scalar (LDF .GE. NX)
C           The actual leading dimension of FCN used in the calling
C           program.
C
C
C   I N P U T / O U T P U T
C   -----------------------
C
C   ICONT   Integer scalar
C           An input flag.  Set ICONT=0 on the first call. Subsequent
C           calls evaluating the same piecewise polynomial function
C           (i.e., the same TX, TY, NX, NY, KX, KY, and FCN) should
C           have ICONT=1.  (See notes on efficiency below.)  As a
C           convenience, DB2VAL sets ICONT=1 on output.
C
C   IWORK   Integer 1D array (size 6)
C           A working storage array.  NOTE: This array is used for
C           communication between successive calls of DB2VAL.  It must
C           not be modified between successive calls when ICONT=1.
C           Different interpolants require different work arrays.
C
C   WORK    Double precision 1D array (size 3*max(KX,KY) + KY + 1)
C           A working storage array.  NOTE: This array is used for
C           communication between successive calls of DB2VAL.  It must
C           not be modified between successive calls when ICONT=1.
C           Different interpolants require different work arrays.
C
C
C   O U T P U T
C   -----------
C
C   IFLAG   Integer scalar.
C           An output flag.  Possible values are
C           0 : successful execution
C           1 : KX out of range
C           2 : NX or LDF out of range
C           3 : KY out of range
C           4 : NY out of range
C           5 : IDX or IDY out of range
C
C           Note: Error conditions are checked in the order shown here.
C                 Only the first error detected is reported.
C
C
C   N O T E S   O N   E F F I C I E N C Y
C   -------------------------------------
C
C   DB2VAL is designed so that it will be reasonably efficient when it is
C   being used to evaluate the fitted function on a grid.  The most
C   favorable situation occurs when DB2VAL is called repeatedly with the
C   same XVAL, the same IDX, and with slowly increasing YVAL.  This
C   implies that when calling DB2VAL in a double loop over x and y :
C
C     (a) vary x in the OUTER loop, and
C     (b) use separate loops or distinct ICONT,IWORK,WORK for
C         each value of IDX
C
C   A typical example of this usage is as follows.
C
C      ICONT = 0
C      DO 20 I=1,NX
C         DO 10 J=1,NY
C            G(I,J) = DB2VAL(X(I),Y(J),IDX,IDY,TX,TY,NX,NY,KX,KY,
C     +                     FCN,LDF,ICONT,IWORK,WORK,IFLAG)
C            IF (IFLAG .NE. 0) GO TO 9999
C   10    CONTINUE
C   20 CONTINUE
C
C   Note that ICONT=0 initially.  This signals DB2VAL that it is starting
C   a new function.  DB2VAL sets ICONT=1 on output, so ICONT=1 on all
C   subsequent calls, which tells DB2VAL to attempt to reuse information
C   from the previous call to speed the evaluation process.
C
C
C   B A C K G R O U N D
C   -------------------
C
C   DB2VAL evaluates the following function or its partial derivatives
C
C                        NX   NY
C            B(x,y)  =  SUM  SUM  FCN(i,j) U (x) V (y)
C                       i=1  j=1            i     j
C   where
C
C     U (x)  is the ith (one-dimensional) B-spline basis function
C      i     defined by NX, KX, and TX, and
C
C     V (y)  is the jth (one-dimensional) B-spline basis function
C      j     defined by NY, KY, and TY.
C
C   See (de Boor, 1978) for a description of the B-spline basis.
C
C   The summation above can be written as
C
C                                 NY
C   (1)               B(x,y)  =  SUM  b  V (y),
C                                j=1   j  j
C   where
C                          NX
C   (2)             b  =  SUM  FCN(i,j) U (x),   j = 1,..,NY.
C                    j    i=1            i
C
C   Note that each summation is the evaluation of a one-dimensional
C   B-spline, which can be done by calls to DB1VAL.  At most KY basis
C   functions in (1) are nonzero.  The indices of these functions are
C   determined, and only the coefficients b(j) which multiply nonzero
C   basis functions are computed using (2).
C
C
C   A C K N O W L E D G E M E N T
C   -----------------------------
C
C   Thanks to Fred Fritsch of Lawrence Livermore National Laboratory
C   for his critical review of this code which lead to many
C   improvements.
C
C***SEE ALSO  DB2INT
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C               Carl de Boor, Efficient computer manipulation of tensor
C                 products, ACM Transactions on Mathematical Software 5,
C                 (1979), pp. 173-182.
C***ROUTINES CALLED  DB1VAL, DINTRV, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   820525  DATE WRITTEN
C   900910  Modified to meet SLATEC standards.
C   930510  Eliminated multiple calls to XERMSG.  (FNF)
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C   940202  Substantial revision.  (RFB)
C           The revision allows for increased efficiency when evaluating
C           on a fine grid.  This is accomplished by reusing information
C           computed on previous calls.  Changes visible to users are :
C           (a) new arguments (ICONT and IWORK), and
C           (b) increase in size of WORK.
C   951130  Changed library to MSSL. (FNF)
C   000503  Changed library to MSSL in error message calls. (ACH)
C   000504  Removed warning comment about SLATEC externals; it no longer
C           applies, as all needed subordinates are now in MSSL. (ACH)
C***END PROLOGUE  DB2VAL
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          IDX, IDY, NX, NY, KX, KY, LDF, ICONT, IWORK(6),
     +                 IFLAG
      DOUBLE PRECISION XVAL, YVAL, TX(NX+KX), TY(NY+KY), FCN(LDF,NY),
     +                 WORK(*)
C
C  .. LOCAL VARIABLES
C
      CHARACTER*6      LIBRAR, SUBROU
      PARAMETER        (LIBRAR = 'MSSL', SUBROU = 'DB2VAL')
      CHARACTER*50     MESSAG
      LOGICAL          NEEDCO
      INTEGER          IDLAST, IERR, IHAVE, INX, ILOY, IY, J, JMIN,
     +                 JMAX, ILOY1, IPOS, IW, LEFTY, LYLAST, MFLAG
      DOUBLE PRECISION DB1VAL, XLAST
C
      EXTERNAL         DINTRV, DB1VAL, XERMSG
C
C***FIRST EXECUTABLE STATEMENT  DB2VAL
      DB2VAL = 0
C
C  ------------------------------
C  RESET STATE FROM PREVIOUS CALL
C  ------------------------------
C
C  These state variables indicate state of DB2VAL at end of last call.
C  They tell us the range of validity of the coefficients b(j) in (2)
C  computed in the previous call and saved in WORK.  We want to
C  avoid recomputing these if possible.
C
C  IWORK(1) == IHAVE  == 1 if coefficients b(j) in (1) are available
C  IWORK(2) == INX    == knot interval of XVAL (for DB1VAL)
C  IWORK(3) == ILOY   == knot interval of YVAL (for DINTRV)
C  IWORK(4) == LYLAST == value of LEFTY from DINTRV for YVAL
C  IWORK(5) == IDLAST == value of IDX
C  IWORK(6) == JMIN   == smallest index j of nonzero basis b(j) at YVAL
C   WORK(1) == XLAST  == value of x variable
C
      IF ( ICONT .NE. 1 ) THEN
         IHAVE = 0
         INX = 1
         ILOY = 1
         JMIN = 0
         IDLAST = -1
         LYLAST = -1
         XLAST = 1.0D30
      ELSE
         IHAVE  = IWORK(1)
         INX    = IWORK(2)
         ILOY   = IWORK(3)
         LYLAST = IWORK(4)
         IDLAST = IWORK(5)
         JMIN   = IWORK(6)
         XLAST  = WORK(1)
      ENDIF
C
C     ... set LEFTY so that it has a value when state is saved after
C         premature termination
C
      LEFTY = 0
C
C  -------------
C  SPECIAL CASES
C  -------------
C
C     ... CHECK INPUT FOR ERRORS
C
      IF (KX. LT. 1) THEN
         IFLAG = 1
         MESSAG = 'KX IS OUT OF RANGE'
      ELSE IF ((NX .LT. KX) .OR. (NX .GT. LDF)) THEN
         IFLAG = 2
         MESSAG = 'NX OR LDF IS OUT OF RANGE'
      ELSE IF (KY. LT. 1) THEN
         IFLAG = 3
         MESSAG = 'KY IS OUT OF RANGE'
      ELSE IF (NY .LT. KY) THEN
         IFLAG = 4
         MESSAG = 'NY IS OUT OF RANGE'
      ELSE IF ((IDX .LT. 0) .OR. (IDY .LT. 0)) THEN
         IFLAG = 5
         MESSAG = 'IDX OR IDY IS OUT OF RANGE'
      ELSE
         IFLAG = 0
      ENDIF
C
      IF ( IFLAG .NE. 0 ) THEN
         CALL XERMSG(LIBRAR,SUBROU,MESSAG,IFLAG,1)
         GO TO 9999
      ENDIF
C
C     ... QUICK RETURN WHEN SPLINE IDENTICALLY ZERO
C
      IF ( (IDX  .GE. KX       ) .OR.
     +     (IDY  .GE. KY       ) .OR.
     +     (XVAL .LT. TX(1)    ) .OR.
     +     (XVAL .GT. TX(NX+KX)) .OR.
     +     (YVAL .LT. TY(1)    ) .OR.
     +     (YVAL .GT. TY(NY+KY))      ) GO TO 9999
C
C  ---------------------
C  EVALUATE THE B-SPLINE
C  ---------------------
C
C     ... PARTITION WORK ARRAY
C
      IY =  2
      IW = IY + KY
C
C     ... FIND KNOT INTERVAL CONTAINING YVAL
C
      CALL DINTRV(TY,NY+KY,YVAL,ILOY,LEFTY,MFLAG)
      IF (MFLAG .NE. 0) THEN
C
C        ... YVAL .EQ. T(NY+KY),  ADJUST TO GET LEFT LIMITING VALUE
C
   10    CONTINUE
         LEFTY = LEFTY - 1
         IF (YVAL. EQ. TY(LEFTY)) GO TO 10
      ENDIF
C
C     ... CHECK IF WE ALREADY HAVE COEFFICIENTS b(j) IN (2)
C
      NEEDCO = (IHAVE .EQ. 0     ) .OR.
     +         (IDX   .NE. IDLAST) .OR.
     +         (LEFTY .NE. LYLAST) .OR.
     +         (XVAL  .NE. XLAST )
C
C     ... GENERATE COEFFICIENTS IF NECESSARY
C
      IF ( NEEDCO ) THEN
C
C        ... FIND RANGE OF INDICES OF NONZERO BASIS FUNCTIONS IN (1)
C            (JMIN,JMAX) = (SMALLEST,LARGEST)
C
         IF (LEFTY .LT. KY) THEN
            JMIN = 1
            JMAX = KY
         ELSEIF (LEFTY .GT. NY) THEN
            JMIN = NY - KY + 1
            JMAX = NY
         ELSE
            JMIN = LEFTY - KY + 1
            JMAX = LEFTY
         ENDIF
C
C        ... GET COEFFICIENTS OF NONZERO BASIS b(j) IN (2)
C
         IPOS = IY - 1
         DO 30 J=JMIN,JMAX
            IPOS = IPOS + 1
            WORK(IPOS) = DB1VAL(XVAL,IDX,TX,NX,KX,FCN(1,J),INX,WORK(IW),
     +                         IERR)
   30    CONTINUE
         IHAVE = 1
C
      ENDIF
C
C     ... EVALUATE INTERPOLANT USING (1)
C
      ILOY1 = KY - 1
      DB2VAL = DB1VAL(YVAL,IDY,TY(JMIN),KY,KY,WORK(IY),ILOY1,WORK(IW),
     +              IERR)
C
C  -------------------
C  SAVE STATE AND EXIT
C  -------------------
C
 9999 CONTINUE
      IWORK(1) = IHAVE
      IWORK(2) = INX
      IWORK(3) = ILOY
      IWORK(4) = LEFTY
      IWORK(5) = IDX
      IWORK(6) = JMIN
      WORK(1) = XVAL
      ICONT = 1
C
      RETURN
      END
*DECK DB1VAL
      DOUBLE PRECISION FUNCTION DB1VAL (X, IDERIV, T, N, K, A, INBV,
     +   WORK, IFLAG)
C***BEGIN PROLOGUE  DB1VAL
C***PURPOSE  Evaluates the B-representation of a spline at X for the
C            function value or any of its derivatives.
C***LIBRARY   MSSL
C     >> Revised or new BSPLINE routine, provisionally accepted for SLATEC 5.0.
C***CATEGORY  E3, K6
C***TYPE      DOUBLE PRECISION (B1VAL-S, DB1VAL-D)
C***KEYWORDS  B-SPLINE, DIFFERENTIATION, EVALUATION, SPLINES
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   DB1VAL evaluates the B-representation (T,N,K,A) of a spline or one
C   of its derivatives at X.  Variable IDERIV indicates which function
C   is to be evaluated.  If B(x) is the spline,
C
C                              IDERIV
C                             d
C                 DB1VAL  =   -------- B(X)
C                               IDERIV
C                             dx
C
C   Thus, to evaluate the spline itself, set IDERIV=0, to get the first
C   derivative with respect to x, set IDERIV=1, and so on.
C
C   Since B is a piecewise polynomial of degree K-1 DB1VAL returns zero
C   whenever IDERIV.GE.K.  DB1VAL also returns zero if X is out of
C   range, that is, if
C
C                    X.LT.T(1) .OR. X.GT.T(N+K).
C
C   Note that the derivative functions computed by DB1VAL may be
C   discontinuous when X corresponds to a knot.   In these cases DB1VAL
C   returns right limiting values (right derivatives), except at T(N+K)
C   where left limiting values are returned.
C
C
C   I N P U T
C   ---------
C
C   X       Double precision scalar
C           The point at which the B-spline is to be evaluated.
C
C   IDERIV  Integer scalar  (IDERIV .GE. 0)
C           Indicates the derivative of the B-spline to evaluate:
C           IDERIV=0 for the function itself, IDERIV=J for the Jth
C           derivative. DB1VAL will return 0 if IDERIV.GE.K.
C
C   T       Double precision 1D array (size N+K)
C           Sequence of knots defining the B-spline.
C
C   N       Integer scalar
C           The number of B-spline coefficients. (N = sum of knot
C           multiplicities - K.)
C
C   K       Integer scalar
C           Order of the B-spline.
C
C   A       Double precision 1D array (size N)
C           The B-spline coefficients.
C
C
C   I N P U T   A N D   O U T P U T
C   -------------------------------
C
C   INBV    Integer scalar
C           An initialization parameter which must be set to 1 the first
C           time DB1VAL is called. On output, contains information for
C           efficient processing after the initial call and must not be
C           changed by the user.  Distinct splines require distinct INBV
C           parameters.
C
C
C   O T H E R
C   ---------
C
C   WORK    Double precision 1D array (size 3*K)
C           Workspace.
C
C   IFLAG   Integer scalar.
C           On return IFLAG indicates the status of the output.
C
C           Output:  0 : successful execution
C                    1 : K out of range
C                    2 : N out of range
C                    3 : IDERIV out of range
C                    4 : Unable to find left-limiting value at T(N+K)
C
C
C  DB1VAL is a version of the routine BVALUE written by Carl de Boor and
C  given in the reference.
C
C  This routine replaces the former SLATEC routine BVALU.  It differs
C  from BVALU in that
C
C    (1) the order of arguments in the calling sequence has been
C        changed and a new parameter, IFLAG, has been added,
C    (2) X may be any real number,
C    (2) IDERIV may be any non-negative integer,
C    (3) right-limiting values are returned when X=T(N+1) and
C        T(N+K).GT.T(N+1), and
C    (4) the prescription given in BVALU for obtaining left-limiting
C        values at interior knots does not work.
C
C
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C***ROUTINES CALLED  DINTRV, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   900726  DATE WRITTEN
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C   951130  Changed library to MSSL. (FNF)
C   000503  Changed library to MSSL in error message calls. (ACH)
C   000504  Removed warning comment about SLATEC externals; it no longer
C           applies, as all needed subordinates are now in MSSL. (ACH)
C***END PROLOGUE  DB1VAL
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          N, K, IDERIV, INBV, IFLAG
      DOUBLE PRECISION T(N+K), A(N), X, WORK(3*K)
C
C  .. LOCAL VARIABLES
C
      CHARACTER*6      LIBRAR, SUBROU
      PARAMETER        (LIBRAR = 'MSSL', SUBROU = 'DB1VAL')
      CHARACTER*50     MESSAG
      INTEGER          I, IDL, IDR, ILO, IMK, IP1, J, JJ, JMIN, JMAX,
     +                 KM1, KMJ, MFLAG, NPK
      DOUBLE PRECISION FKMJ
C
      EXTERNAL         DINTRV, XERMSG
C
C***FIRST EXECUTABLE STATEMENT  DB1VAL
      DB1VAL = 0
      NPK = N + K
C
C  -------------
C  SPECIAL CASES
C  -------------
C
C     ... CHECK INPUT FOR ERRORS
C
      IF (K. LT. 1)  GO TO 9001
      IF (N .LT. K)  GO TO 9002
      IF (IDERIV .LT. 0)  GO TO 9003
C
C     ... QUICK RETURN WHEN SPLINE EVALUATES TO ZERO
C
      IF (IDERIV .GE. K)  GO TO 9999
      IF ((X .LT. T(1)) .OR. (X .GT. T(NPK)))  GO TO 9999
C
C  -------------------------------
C  FIND KNOT INTERVAL CONTAINING X
C  -------------------------------
C
C     ... FIND LARGEST I IN (1,N+K) SUCH THAT T(I) .LE. X .LT. T(I+1)
C
      CALL DINTRV(T, NPK, X, INBV, I, MFLAG)
      IF (MFLAG .NE. 0) THEN
C
C        ... X .EQ. T(N+K),  USE LEFT LIMITING VALUE
C
   10    CONTINUE
         IF (I .EQ. 1)  GO TO 9004
         I = I - 1
         IF (X. EQ. T(I)) GO TO 10
      ENDIF
C
C  ------------------------------------------------------
C  STORE COEFFICIENTS OF NON-ZERO BASIS FUNCTIONS IN WORK
C  ------------------------------------------------------
C
C     ... IF THERE ARE LESS THAN K BASIS FUNCTIONS, WE ASSUME FICTITIOUS
C         ONES WITH ZERO COEFFICIENTS AND KNOTS AT T(1) OR T(N+K)
C
      IMK = I - K
      JMIN = MAX(1-IMK,1)
      JMAX = MIN(NPK-I,K)
      DO 20 J=1,JMIN-1
         WORK(J) = 0
   20 CONTINUE
      DO 30 J=JMIN,JMAX
         WORK(J) = A(IMK+J)
   30 CONTINUE
      DO 40 J=JMAX+1,K
         WORK(J) = 0
   40 CONTINUE
C
C  ----------------------------
C  COMPUTE AUXILIARY QUANTITIES
C  ----------------------------
C
C     ... DL(J) = X - T(I+1-J) = WORK(IDL+J), J=1,..,K-1
C     ... DR(J) = T(I+J) - X   = WORK(IDR+J), J=1,..,K-1
C
      IDL = K
      IDR = IDL + K
      KM1 = K - 1
      IP1 = I + 1
      JMAX = MIN(I,KM1)
      DO 50 J=1,JMAX
         WORK(IDL+J) = X - T(IP1-J)
   50 CONTINUE
      DO 60 J=JMAX+1,KM1
         WORK(IDL+J) = WORK(IDL+JMAX)
   60 CONTINUE
      JMAX = MIN(NPK-I,KM1)
      DO 70 J=1,JMAX
         WORK(IDR+J) = T(I+J) - X
   70 CONTINUE
      DO 80 J=JMAX+1,KM1
         WORK(IDR+J) = WORK(IDR+JMAX)
   80 CONTINUE
C
C  ----------------------------------------
C  DIFFERENCE THE COEFFICIENTS IDERIV TIMES
C  ----------------------------------------
C
      DO 100 J=1,IDERIV
        KMJ = K - J
        FKMJ = KMJ
        ILO = KMJ
        DO 90 JJ=1,KMJ
           WORK(JJ) = FKMJ*(WORK(JJ+1)-WORK(JJ))
     +                /(WORK(IDL+ILO)+WORK(IDR+JJ))
           ILO = ILO - 1
   90   CONTINUE
  100 CONTINUE
C
C  ---------------------------------
C  COMPUTE IDERIV-TH DERIVATIVE AT X
C  ---------------------------------
C
C     ... B-SPLINE COEFFICIENTS ARE IN WORK(1),..,WORK(K-IDERIV)
C
      DO 120 J=IDERIV+1,KM1
         KMJ = K - J
         ILO = KMJ
         DO 110 JJ=1,KMJ
            WORK(JJ) = (WORK(JJ+1)*WORK(IDL+ILO)+WORK(JJ)*WORK(IDR+JJ))
     +                 /(WORK(IDL+ILO)+WORK(IDR+JJ))
            ILO = ILO - 1
  110    CONTINUE
  120 CONTINUE
      DB1VAL = WORK(1)
      IFLAG = 0
      GO TO 9999
C
C  -----------
C  ERROR EXITS
C  -----------
C
 9001 CONTINUE
      IFLAG = 1
      MESSAG = 'K DOES NOT SATISFY K.GE.1'
      GO TO 9900
C
 9002 CONTINUE
      IFLAG = 2
      MESSAG = 'N DOES NOT SATISFY N.GE.K'
      GO TO 9900
C
 9003 CONTINUE
      IFLAG = 3
      MESSAG = 'IDERIV IS LESS THAN ZERO'
      GO TO 9900
C
 9004 CONTINUE
      IFLAG = 4
      MESSAG = 'NO LEFT LIMITING VALUE AT X=T(N+K)'
      GO TO 9900
C
 9900 CONTINUE
      CALL XERMSG(LIBRAR,SUBROU,MESSAG,IFLAG,1)
C
C  ----
C  EXIT
C  ----
C
 9999 CONTINUE
      RETURN
      END
*DECK DINTRV
      SUBROUTINE DINTRV (XT, LXT, X, ILO, ILEFT, MFLAG)
C***BEGIN PROLOGUE  DINTRV
C***PURPOSE  Compute the largest integer ILEFT in 1 .LE. ILEFT .LE. LXT
C            such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
C            the X interval.
C***LIBRARY   MSSL
C***CATEGORY  E3, K6
C***TYPE      DOUBLE PRECISION (INTRV-S, DINTRV-D)
C***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Abstract    **** a double precision routine ****
C         DINTRV is the INTERV routine of the reference.
C
C         DINTRV computes the largest integer ILEFT in 1 .LE. ILEFT .LE.
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
C
C         Input      XT,X are double precision
C          XT      - XT is a knot or break point vector of length LXT
C          LXT     - length of the XT vector
C          X       - argument
C          ILO     - an initialization parameter which must be set
C                    to 1 the first time the spline array XT is
C                    processed by DINTRV.
C
C         Output
C          ILO     - ILO contains information for efficient process-
C                    ing after the initial call and ILO must not be
C                    changed by the user.  Distinct splines require
C                    distinct ILO parameters.
C          ILEFT   - largest integer satisfying XT(ILEFT) .LE. X
C          MFLAG   - signals when X lies out of bounds
C
C     Error Conditions
C         None
C
C***REFERENCES  Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   000501  Changed library to MSSL (ACH)
C***END PROLOGUE  DINTRV
C
      INTEGER IHI, ILEFT, ILO, ISTEP, LXT, MFLAG, MIDDLE
      DOUBLE PRECISION X, XT
      DIMENSION XT(*)
C***FIRST EXECUTABLE STATEMENT  DINTRV
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
*DECK XERMSG
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
C-----------------------------------------------------------------------
C Subroutines XERMSG, XSETF, XSETUN, and the function routine IXSAV, as
C given here, constitute a simplified version of the SLATEC error 
C handling package.  Written by A. C. Hindmarsh, 18 November 1992.
C
C All arguments are input arguments.
C LIBRAR = Library name (character array).  Prefixed to message.
C SUBROU = Routine name (character array).  Prefixed to message.
C MESSG  = The message (character array).
C NERR   = Integer error number.  Prefixed to message.
C LEVEL  = The error level..
C          0 or 1 means recoverable (control returns to caller).
C          2 means fatal (run is aborted--see note below).
C
C Note..  This routine has been simplified in the following ways..
C 1. A single prefix line is printed with NERR, SUBROU, and LIBRAR.
C 2. The message in MESSG is printed, unaltered, on lines of up to 72
C    characters each using a format of (A).
C 3. If LEVEL = 2, control passes to the statement   STOP
C    to abort the run.  This statement may be machine-dependent.
C
C For a different default logical unit number, IXSAV (or a subsidiary
C routine that it calls) will need to be modified.
C For a different run-abort command, change the statement following
C statement 100 at the end.
C-----------------------------------------------------------------------
C Subroutines called by XERMSG.. None
C Function routines called by XERMSG.. IXSAV
C Intrinsic function used by XERMSG.. LEN
C-----------------------------------------------------------------------
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      INTEGER NERR, LEVEL
      INTEGER I1, I2, IL, IXSAV, LENMSG, LLEN, LUNIT, MESFLG, NLINES
      PARAMETER (LLEN = 72)
C
C Get message print flag and logical unit number. ----------------------
      LUNIT = IXSAV (1, 0, .FALSE.)
      MESFLG = IXSAV (2, 0, .FALSE.)
      IF (MESFLG .EQ. 0) GO TO 100
C Write NERR, SUBROU, and LIBRAR. --------------------------------------
      I1 = LEN(SUBROU)
      I2 = LEN(LIBRAR)
      WRITE (LUNIT, 10) NERR, SUBROU(1:I1), LIBRAR(1:I2)
  10  FORMAT(/,'***Error number ',I6,' from ',A,' in library ',A,'***')
C Write the message. ---------------------------------------------------
      LENMSG = LEN(MESSG)
      NLINES = ( (LENMSG - 1)/LLEN ) + 1
      DO 20 IL = 1,NLINES
        I1 = 1 + (IL - 1)*LLEN
        I2 = MIN(IL*LLEN,LENMSG)
        WRITE (LUNIT,'(A)') MESSG(I1:I2)
  20    CONTINUE
C Abort the run if LEVEL = 2. ------------------------------------------
 100  IF (LEVEL .NE. 2) RETURN
      STOP
C----------------------- End of Subroutine XERMSG ----------------------
      END
*DECK IXSAV
      INTEGER FUNCTION IXSAV (IPAR, IVALUE, ISET)
C***BEGIN PROLOGUE  IXSAV
C***SUBSIDIARY
C***PURPOSE  Save and recall error message control parameters.
C***CATEGORY  R3C
C***TYPE      ALL (IXSAV-A)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  IXSAV saves and recalls one of two error message parameters:
C    LUNIT, the logical unit number to which messages are printed, and
C    MESFLG, the message print flag.
C  This is a modification of the SLATEC library routine J4SAVE.
C
C  Saved local variables..
C   LUNIT  = Logical unit number for messages.  The default is obtained
C            by a call to IUMACH (may be machine-dependent).
C   MESFLG = Print control flag..
C            1 means print all messages (the default).
C            0 means no printing.
C
C  On input..
C    IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
C    IVALUE = The value to be set for the parameter, if ISET = .TRUE.
C    ISET   = Logical flag to indicate whether to read or write.
C             If ISET = .TRUE., the parameter will be given
C             the value IVALUE.  If ISET = .FALSE., the parameter
C             will be unchanged, and IVALUE is a dummy argument.
C
C  On return..
C    IXSAV = The (old) value of the parameter.
C
C***SEE ALSO  XERRWD, XERRWV
C***ROUTINES CALLED  IUMACH
C***REVISION HISTORY  (YYMMDD)
C   921118  DATE WRITTEN
C   930329  Modified prologue to SLATEC format. (FNF)
C   930915  Added IUMACH call to get default output unit.  (ACH)
C   930922  Minor cosmetic changes. (FNF)
C   010425  Type declaration for IUMACH added. (ACH)
C***END PROLOGUE  IXSAV
C
C Subroutines called by IXSAV.. None
C Function routine called by IXSAV.. IUMACH
C-----------------------------------------------------------------------
C**End
      LOGICAL ISET
      INTEGER IPAR, IVALUE
C-----------------------------------------------------------------------
      INTEGER IUMACH, LUNIT, MESFLG
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this routine.
C-----------------------------------------------------------------------
      SAVE LUNIT, MESFLG
      DATA LUNIT/-1/, MESFLG/1/
C
C***FIRST EXECUTABLE STATEMENT  IXSAV
      IF (IPAR .EQ. 1) THEN
        IF (LUNIT .EQ. -1) LUNIT = IUMACH()
        IXSAV = LUNIT
        IF (ISET) LUNIT = IVALUE
        ENDIF
C
      IF (IPAR .EQ. 2) THEN
        IXSAV = MESFLG
        IF (ISET) MESFLG = IVALUE
        ENDIF
C
      RETURN
C----------------------- End of Function IXSAV -------------------------
      END
*DECK IUMACH
      INTEGER FUNCTION IUMACH()
C***BEGIN PROLOGUE  IUMACH
C***PURPOSE  Provide standard output unit number.
C***CATEGORY  R1
C***TYPE      INTEGER (IUMACH-I)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C *Usage:
C        INTEGER  LOUT, IUMACH
C        LOUT = IUMACH()
C
C *Function Return Values:
C     LOUT : the standard logical unit for Fortran output.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   930915  DATE WRITTEN
C   930922  Made user-callable, and other cosmetic changes. (FNF)
C***END PROLOGUE  IUMACH
C
C*Internal Notes:
C  The built-in value of 6 is standard on a wide range of Fortran
C  systems.  This may be machine-dependent.
C**End
C***FIRST EXECUTABLE STATEMENT  IUMACH
      IUMACH = 6
C
      RETURN
C----------------------- End of Function IUMACH ------------------------
      END
*DECK XSETF
      SUBROUTINE XSETF (MFLAG)
C***BEGIN PROLOGUE  XSETF
C***PURPOSE  Reset the error print control flag.
C***CATEGORY  R3A
C***TYPE      ALL (XSETF-A)
C***KEYWORDS  ERROR CONTROL
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C   XSETF sets the error print control flag to MFLAG:
C      MFLAG=1 means print all messages (the default).
C      MFLAG=0 means no printing.
C
C***SEE ALSO  XERRWD, XERRWV
C***REFERENCES  (NONE)
C***ROUTINES CALLED  IXSAV
C***REVISION HISTORY  (YYMMDD)
C   921118  DATE WRITTEN
C   930329  Added SLATEC format prologue. (FNF)
C   930407  Corrected SEE ALSO section. (FNF)
C   930922  Made user-callable, and other cosmetic changes. (FNF)
C***END PROLOGUE  XSETF
C
C Subroutines called by XSETF.. None
C Function routine called by XSETF.. IXSAV
C-----------------------------------------------------------------------
C**End
      INTEGER MFLAG, JUNK, IXSAV
C
C***FIRST EXECUTABLE STATEMENT  XSETF
      IF (MFLAG .EQ. 0 .OR. MFLAG .EQ. 1) JUNK = IXSAV (2,MFLAG,.TRUE.)
      RETURN
C----------------------- End of Subroutine XSETF -----------------------
      END
*DECK XSETUN
      SUBROUTINE XSETUN (LUN)
C***BEGIN PROLOGUE  XSETUN
C***PURPOSE  Reset the logical unit number for error messages.
C***CATEGORY  R3B
C***TYPE      ALL (XSETUN-A)
C***KEYWORDS  ERROR CONTROL
C***DESCRIPTION
C
C   XSETUN sets the logical unit number for error messages to LUN.
C
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***SEE ALSO  XERRWD, XERRWV
C***REFERENCES  (NONE)
C***ROUTINES CALLED  IXSAV
C***REVISION HISTORY  (YYMMDD)
C   921118  DATE WRITTEN
C   930329  Added SLATEC format prologue. (FNF)
C   930407  Corrected SEE ALSO section. (FNF)
C   930922  Made user-callable, and other cosmetic changes. (FNF)
C***END PROLOGUE  XSETUN
C
C Subroutines called by XSETUN.. None
C Function routine called by XSETUN.. IXSAV
C-----------------------------------------------------------------------
C**End
      INTEGER LUN, JUNK, IXSAV
C
C***FIRST EXECUTABLE STATEMENT  XSETUN
      IF (LUN .GT. 0) JUNK = IXSAV (1,LUN,.TRUE.)
      RETURN
C----------------------- End of Subroutine XSETUN ----------------------
      END
