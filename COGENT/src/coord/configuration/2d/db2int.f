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
*DECK DB2INT
      SUBROUTINE DB2INT (X, NX, Y, NY, KX, KY, TX, TY, FCN, LDF, WORK,
     +   IFLAG)
C***BEGIN PROLOGUE  DB2INT
C***PURPOSE  DB2INT determines a piecewise polynomial function that
C            interpolates two-dimensional gridded data. Users specify
C            the polynomial order (degree+1) of the interpolant and
C            (optionally) the knot sequence.
C***LIBRARY   MSSL
C     >> Revised or new BSPLINE routine, provisionally accepted for SLATEC 5.0.
C***CATEGORY  E2A
C***TYPE      DOUBLE PRECISION (B2INT-S, DB2INT-D)
C***KEYWORDS  INTERPOLATION, TWO DIMENSIONS, GRIDDED DATA, SPLINES,
C             PIECEWISE POLYNOMIALS
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   DB2INT determines the parameters of a function that interpolates the
C   two-dimensional gridded data (X(i),Y(j),FCN(i,j)) for i=1,..,NX and
C   j=1,..,NY.  The interpolating function and its derivatives may be
C   subsequently evaluated by the function DB2VAL.
C
C   The interpolating function is a piecewise polynomial (pp) function
C   represented as a tensor product of one-dimensional B-splines.  The
C   form of this function is
C
C                          NX   NY
C              S(x,y)  =  SUM  SUM  a   U (x) V (y)
C                         i=1  j=1   ij  i     j
C
C   where the functions U(i) and V(j) are one-dimensional B-spline
C   basis functions. The coefficients a(i,j) are chosen so that
C
C         S(X(i),Y(j)) = FCN(i,j)   for i=1,..,NX and j=1,..,NY.
C
C   Note that for each fixed value of y S(x,y) is a pp function of x
C   alone, and for each fixed x S(x,y) is a pp function of y alone.  In
C   one dimension a piecewise polynomial may be created by partitioning
C   a given interval into subintervals and defining a distinct poly-
C   nomial on each one. The points where adjacent subintervals meet are
C   called  knots.  Each of the functions U(i) and V(j) above is a
C   piecewise polynomial.
C
C   Users of DB2INT choose the order (degree+1) of the polynomial pieces
C   used to define the interpolant in each of the x and y directions
C   (KX and KY). Users also may define their own knot sequence in x and
C   y separately (TX and TY).  If IFLAG=0, however, DB2INT will choose
C   knots that result in a pp interpolant with KX-2 continuous partial
C   derivatives in x and KY-2 continuous partial derivatives in y.  The
C   interpolating function is identically zero outside the rectangular
C   region defined by the knots. See below for more information on knot
C   selection.
C
C   After a call to DB2INT, all information necessary to define the
C   interpolating function is contained in the parameters NX, NY, KX,
C   KY, TX, TY, and FCN.  These quantities should not be altered until
C   after the last call of the evaluation routine DB2VAL.
C
C
C   I N P U T
C   ---------
C
C   X       Double precision 1D array (size NX)
C           Array of x abscissae. Must be strictly increasing.
C
C   NX      Integer scalar (2 .LE. NX .LE. LDF)
C           Number of x abscissae.
C
C   Y       Double precision 1D array (size NY)
C           Array of y abscissae. Must be strictly increasing.
C
C   NY      Integer scalar (NY .GE. 2)
C           Number of y abscissae.
C
C   KX      Integer scalar (2 .LE. KX .LE. NX)
C           The order (degree + 1) of polynomial pieces in x.
C
C   KY      Integer scalar (2 .LE. KY .LE. NY)
C           The order (degree + 1) of polynomial pieces in y.
C
C
C   I N P U T   O R   O U T P U T
C   -----------------------------
C
C   TX      Double precision 1D array (size NX+KX)
C           The knots in the x direction. If IFLAG=0 these are chosen
C           by DB2INT. If IFLAG=1 they are specified by the user as
C           input. The knots must be non-decreasing and also satisfy
C                  TX(i).LT.X(i).LT.TX(i+KX),  i=1,..,NX.
C           Equality on the left is permitted for i=1 and equality on
C           the right is permitted for i=NX. (See information below on
C           knot placement.)
C
C   TY      Double precision 1D array (size NY+KY)
C           The knots in the y direction. If IFLAG=0 these are chosen
C           by DB2INT. If IFLAG=1 they are specified by the user as
C           input. The knots must be non-decreasing and also satisfy
C                  TY(i).LT.Y(i).LT.TY(i+KY),  i=1,..,NY.
C           Equality on the left is permitted for i=1 and equality on
C           the right is permitted for i=NY. (See information below on
C           knot placement.)
C
C
C   I N P U T   A N D   O U T P U T
C   -------------------------------
C
C   FCN     Double precision 2D array (size LDF by NY)
C           Input : Array of function values to interpolate. FCN(I,J)
C                   should contain the function value at the point
C                   (X(I),Y(J)).
C           Output: Array of coefficients of the B-spline interpolant.
C
C   LDF     Integer scalar (LDF .GE. NX)
C           The actual leading dimension of FCN used in the calling
C           program.
C
C
C   O T H E R
C   ---------
C
C   WORK    Double precision 1D array
C           Array of working storage. Must be dimensioned of length
C           at least NX*NY + 2*max(KX*(NX+1),KY*(NY+1)).
C
C   IFLAG   Integer scalar.
C           Must be set by the user before DB2INT is called.
C           On return IFLAG indicates the status of the output.
C
C           Input :  0 : knot sequence chosen by user
C                    1 : knot sequence chosen by DB2INT
C
C           Output:  0 : successful execution
C                    2 : IFLAG out of range
C                    3 : NX or LDF out of range
C                    4 : KX out of range
C                    5 : X not strictly increasing
C                    6 : TX is an illegal knot sequence
C                    7 : NY out of range
C                    8 : KY out of range
C                    9 : Y not strictly increasing
C                   10 : TY is an illegal knot sequence
C
C           Note: Error conditions are checked in the order shown here.
C                 Only the first error detected is reported.
C
C           After a successful return, DB2INT may be recalled without
C           resetting IFLAG, provided that only the array FCN has
C           changed.
C
C
C   K N O T   S E L E C T I O N
C   ---------------------------
C
C   In this section we describe the relationship between data points
C   and knots.  Users who choose to let DB2INT select the knot sequences
C   TX and TY (by setting IFLAG=1) may skip this discussion.
C
C   To describe the relationship between knots and data points we  first
C   consider the simpler case of one-dimensional interpolation; in
C   particular, we consider interpolating the data X(i), i=1,..,N with
C   a pp function of order K.
C
C   Knots are the points where the individual polynomial pieces join up,
C   and hence where the pp function may suffer a loss of smoothness.  To
C   define a pp function, one needs N+K knots. If the knots are distinct
C   the interpolant will be as smooth as possible (continuous, with K-2
C   continuous derivatives).  If two adjacent knots come together, the
C   smoothness of the function is reduced at that point.  In general, if
C   M knots are assigned to the same point, the pp function will have
C   K-M-1 continuous derivatives there. If K knots are taken at the same
C   point, then the pp function itself will be discontinuous there.
C
C   Typically, K knots are taken at or to the left of the leftmost data
C   point, K knots at or to the right of the rightmost data point, with
C   the remaining N-K knots in between.  In order for there to be a
C   solution to the interpolation problem the knots T(i), must satisfy
C   certain additional constraints. That is,
C
C                T(i) .LT. X(i) .LT. T(i+K)   i=1,..,N
C
C   Equality is permitted on the left for I=1 and on the right for I=N.
C
C   The two-dimensional interpolant computed by this routine is a tensor
C   product of one-dimensional pp interpolants.  The knots form a grid
C   (TX(i),TY(j)) in the same way that the data points do.  Along lines
C   parallel to the coordinate axes, the interpolant reduces to a one-
C   dimensional pp function with order and knots (KX,TX) or (KY,TY). In
C   this case the appropriate constraints on the knots become
C
C            TX(i) .LT. X(i) .LT. TX(i+KX)   i=1,..,NX
C            TY(i) .LT. Y(i) .LT. TY(i+KY)   i=1,..,NY
C
C   with equality on the left permitted when i=1 and equality on the
C   right permitted when i=NX and i=NY, respectively.
C
C   If these conditions are violated, then DB2INT returns with IFLAG
C   equal to 6 or 10.  The default knot sequence selected by DB2INT
C   always satisfies these conditions.
C
C   When the user sets IFLAG=1 DB2INT selects knots as follows. KX knots
C   are taken at each endpoint in the x direction, not-a-knot end
C   conditions (see references) are used, and the remaining knots are
C   placed at data points if KX is even or at midpoints between data
C   points if KX is odd.  The y direction is treated similarly.  This
C   yields a two-dimensional pp function with KX-2 continuous partial
C   derivatives in x and KY-2 continuous partial derivatives in y.  The
C   interpolant is zero outside the rectangular region defined by the
C   data points, and discontinuous along the boundary of this region.
C
C***SEE ALSO  DB2VAL
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C               Carl de Boor, Efficient computer manipulation of tensor
C                 products, ACM Transactions on Mathematical Software 5,
C                 (1979), pp. 173-182.
C***ROUTINES CALLED  DBKCHK, DBKNOT, DBTPCF, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   820525  DATE WRITTEN
C   900910  Modified to meet SLATEC standards.
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C   951130  Changed library to MSSL. (FNF)
C   951218  Corrected CATEGORY line. (FNF)
C   000503  Changed library to MSSL in error message calls. (ACH)
C   000504  Removed warning comment about SLATEC externals; it no longer
C           applies, as all needed subordinates are now in MSSL. (ACH)
C***END PROLOGUE  DB2INT
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          NX, NY, KX, KY, LDF, IFLAG
      DOUBLE PRECISION X(NX), Y(NY), TX(NX+KX), TY(NY+KY), FCN(LDF,NY),
     +                 WORK(*)
C
C  .. LOCAL VARIABLES
C
      CHARACTER*6      LIBRAR, SUBROU
      PARAMETER        (LIBRAR = 'MSSL', SUBROU = 'DB2INT')
      CHARACTER*50     MESSAG
      INTEGER          I, IW
      LOGICAL          DBKCHK
C
      EXTERNAL         DBKCHK, DBKNOT, DBTPCF, XERMSG
C
C***FIRST EXECUTABLE STATEMENT  DB2INT
C
C  -----------------------
C  CHECK VALIDITY OF INPUT
C  -----------------------
C
      IF ((IFLAG .LT. 0) .OR. (IFLAG .GT. 1))  GO TO 9020
C
      IF ((NX .LT. 2) .OR. (NX .GT. LDF))  GO TO 9030
      IF ((KX .LT. 2) .OR. (KX .GT. NX))  GO TO 9040
      DO 10 I=2,NX
         IF (X(I) .LE. X(I-1))  GO TO 9050
   10 CONTINUE
      IF (IFLAG .EQ. 0) THEN
         IF(.NOT. DBKCHK(X,NX,KX,TX)) GO TO 9060
      ENDIF
C
      IF (NY .LT. 2)  GO TO 9070
      IF ((KY .LT. 2) .OR. (KY .GT. NY))  GO TO 9080
      DO 20 I=2,NY
         IF (Y(I) .LE. Y(I-1))  GO TO 9090
   20 CONTINUE
      IF (IFLAG .EQ. 0) THEN
         IF(.NOT. DBKCHK(Y,NY,KY,TY)) GO TO 9100
      ENDIF
C
C  ------------
C  CHOOSE KNOTS
C  ------------
C
      IF (IFLAG .EQ. 1) THEN
         CALL DBKNOT(X,NX,KX,TX)
         CALL DBKNOT(Y,NY,KY,TY)
      ENDIF
C
C  -------------------------------
C  CONSTRUCT B-SPLINE COEFFICIENTS
C  -------------------------------
C
      IW = NX*NY + 1
      CALL DBTPCF(X,NX,FCN,LDF,NY,TX,KX,WORK,NY,WORK(IW))
      CALL DBTPCF(Y,NY,WORK,NY,NX,TY,KY,FCN,LDF,WORK(IW))
      IFLAG = 0
      GO TO 9999
C
C  -----------
C  ERROR EXITS
C  -----------
C
 9020 CONTINUE
      IFLAG = 2
      MESSAG = 'IFLAG IS OUT OF RANGE'
      GO TO 9900
C
 9030 CONTINUE
      IFLAG = 3
      MESSAG = 'NX OR LDF IS OUT OF RANGE'
      GO TO 9900
C
 9040 CONTINUE
      IFLAG = 4
      MESSAG = 'KX IS OUT OF RANGE'
      GO TO 9900
C
 9050 CONTINUE
      IFLAG = 5
      MESSAG = 'X ARRAY MUST BE STRICTLY INCREASING'
      GO TO 9900
C
 9060 CONTINUE
      IFLAG = 6
      MESSAG = 'TX IS AN ILLEGAL KNOT SEQUENCE'
      GO TO 9900
C
 9070 CONTINUE
      IFLAG = 7
      MESSAG = 'NY IS OUT OF RANGE'
      GO TO 9900
C
 9080 CONTINUE
      IFLAG = 8
      MESSAG = 'KY IS OUT OF RANGE'
      GO TO 9900
C
 9090 CONTINUE
      IFLAG = 9
      MESSAG = 'Y ARRAY MUST BE STRICTLY INCREASING'
      GO TO 9900
C
 9100 CONTINUE
      IFLAG = 10
      MESSAG = 'TY IS AN ILLEGAL KNOT SEQUENCE'
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
*DECK DBKCHK
      LOGICAL FUNCTION DBKCHK (X, N, K, T)
C***BEGIN PROLOGUE  DBKCHK
C***SUBSIDIARY
C***PURPOSE  DBKCHK checks whether piecewise polynomial interpolation
C            of order K at the interpolation points X(i), i=1,..,N
C            is possible using the knot sequence T(i), i=1,..,N+K.
C***LIBRARY   MSSL
C     >> Revised or new BSPLINE routine, provisionally accepted for SLATEC 5.0.
C***TYPE      DOUBLE PRECISION (BKCHK-S, DBKCHK-D)
C***KEYWORDS  KNOT SELECTION, INTERPOLATION, PIECEWISE POLYNOMIALS,
C             SPLINES
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   DBKCHK is a subsidiary routine called by DB2INT and DB3INT.
C
C   DBKCHK checks whether piecewise polynomial interpolation of order K
C   at the interpolation points X(i), i=1,..,N is possible using the
C   knot sequence T(i), i=1,..,N+K.
C
C   Specifically, the following criteria are checked:
C
C     (1) T is non-decreasing
C     (2) T(1) .LE. X(1) .LT. T(1+K)
C     (3) T(i) .LT. X(i) .LT. T(i+K)   i=2,..,N-1
C     (4) T(N) .LT. X(N) .LE. T(N+K)
C
C   DBKCHK returns .TRUE. if the knot sequence is legal and .FALSE. if
C   it is not legal.
C
C
C   I N P U T
C   ---------
C
C   X       Double precision 1D array (size N)
C           Array of data points. Must be strictly increasing.
C
C   N       Integer scalar (N .GE. 2)
C           Number of data points.
C
C   K       Integer scalar (2 .LE. K .LE. N)
C           The order (degree + 1) of polynomial pieces.
C
C   T       Double precision 1D array (size N+K)
C           The selected knot sequence. Non-decreasing.
C
C
C   CAUTION: The constraints on the variables X, N, and K are not
C            checked by this routine.
C
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   900910  DATE WRITTEN
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C   951130  Changed library to MSSL. (FNF)
C   000504  Removed warning comment about SLATEC externals; it no longer
C           applies, as all needed subordinates are now in MSSL. (ACH)
C***END PROLOGUE  DBKCHK
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          N, K
      DOUBLE PRECISION X(N), T(N+K)
C
C  .. LOCAL VARIABLES
C
      INTEGER          I
C
C***FIRST EXECUTABLE STATEMENT  DBKCHK
      DBKCHK = .TRUE.
C
C     ... CHECK WHETHER T IS NON-DECREASING
C
      DO 10 I=2,N+K
         IF (T(I-1) .GT. T(I))  GO TO 100
   10 CONTINUE
C
C     ... CHECK WHETHER THE X ARE PROPERLY DISTRIBUTED
C
      IF ((T(1) .GT. X(1)) .OR. (X(1) .GE. T(1+K)))  GO TO 100
      DO 20 I=2,N-1
         IF ((T(I) .GE. X(I)) .OR. (X(I) .GE. T(I+K)))  GO TO 100
   20 CONTINUE
      IF ((T(N) .GE. X(N)) .OR. (X(N) .GT. T(N+K))) GO TO 100
      GO TO 9999
C
C     ... EXIT THROUGH HERE IF ILLEGAL
C
  100 CONTINUE
      DBKCHK = .FALSE.
C
 9999 CONTINUE
      RETURN
      END
*DECK DBKNOT
      SUBROUTINE DBKNOT (X, N, K, T)
C***BEGIN PROLOGUE  DBKNOT
C***SUBSIDIARY
C***PURPOSE  DBKNOT selects a sequence of N+K knots for use in spline
C            interpolation of order K at a given set of N data points.
C            The selection yields an interpolant with K-2 continuous
C            derivatives, except at the two endpoints.  Not-a-knot end
C            conditions are used.
C***LIBRARY   MSSL
C     >> Revised or new BSPLINE routine, provisionally accepted for SLATEC 5.0.
C***TYPE      DOUBLE PRECISION (BKNOT-S, DBKNOT-D)
C***KEYWORDS  KNOT SELECTION, INTERPOLATION, PIECEWISE POLYNOMIALS,
C             SPLINES
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   DBKNOT is a subsidiary routine called by DB2INT and DB3INT.
C
C   DBKNOT selects a sequence of N+K knots for use in piecewise
C   polynomial interpolation of order K at a given set of N data points.
C
C   The knots T(i), i=1,..,N+K, are selected as follows:
C
C     number                    location
C     ----------------------------------------------------------------
C     K     at X(1)
C     N-K   between X(1+K/2) and X(N-K+K/2); these are placed at data
C           points if N is even and midway between data points if N is
C           odd
C     K     at X(N)
C
C   where X(i), i=1,..,N are the data points.
C
C   This selection yields an interpolant with K-2 continuous deriva-
C   tives, except at X(1) and X(N); the interpolant will be zero outside
C   (X(1),X(N)).  Knot-a-knot end conditions are used, resulting in K-1
C   continuous derivatives at T(K+1) and T(N).
C
C
C   I N P U T
C   ---------
C
C   X       Double precision 1D array (size N)
C           Array of data points. Must be strictly increasing.
C
C   N       Integer scalar (N .GE. 2)
C           Number of data points.
C
C   K       Integer scalar (2 .LE. K .LE. N)
C           The order (degree + 1) of polynomial pieces.
C
C
C   O U T P U T
C   -----------
C
C   T       Double precision 1D array (size N+K)
C           The selected knot sequence. Non-decreasing.
C
C
C   CAUTION: The constraints on the input variables are not checked by
C            this routine.
C
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   820525  DATE WRITTEN
C   900910  Modified to meet SLATEC standards.
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C   951130  Changed library to MSSL. (FNF)
C   000504  Removed warning comment about SLATEC externals; it no longer
C           applies, as all needed subordinates are now in MSSL. (ACH)
C***END PROLOGUE  DBKNOT
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          N, K
      DOUBLE PRECISION X(N), T(N+K)
C
C  .. LOCAL VARIABLES
C
      INTEGER          I, IPJ, J
      DOUBLE PRECISION HALF
C
      PARAMETER ( HALF = 0.50D0 )
C
C  ----------------------------
C  PUT K KNOTS AT EACH ENDPOINT
C  ----------------------------
C
C***FIRST EXECUTABLE STATEMENT  DBKNOT
      DO 100 J=1,K
         T(J)   = X(1)
         T(N+J) = X(N)
  100 CONTINUE
C
C  --------------------------
C  DISTRIBUTE REMAINING KNOTS
C  --------------------------
C
      IF (MOD(K,2) .EQ. 0) THEN
C
C        ... EVEN K --  KNOTS AT DATA POINTS
C
         I = (K/2) - K
         DO 120 J=K+1,N
            T(J) = X(I+J)
  120    CONTINUE
C
      ELSE
C
C        ... ODD K --  KNOTS BETWEEN DATA POINTS
C
         I = (K-1)/2 - K
         DO 160 J=K+1,N
            IPJ = I + J
            T(J) = HALF*( X(IPJ) + X(IPJ+1) )
  160    CONTINUE
C
      ENDIF
C
      RETURN
      END
*DECK DBTPCF
      SUBROUTINE DBTPCF (X, N, FCN, LDF, NF, T, K, BCOEF, LDB, WORK)
C***BEGIN PROLOGUE  DBTPCF
C***SUBSIDIARY
C***PURPOSE  DBTPCF computes NF sets of B-spline coefficients which
C            determine NF distinct piecewise polynomial interpolating
C            functions.  The functions share a common order, set of
C            interpolation points, and knot sequence.
C***LIBRARY   MSSL
C     >> Revised or new BSPLINE routine, provisionally accepted for SLATEC 5.0.
C***TYPE      DOUBLE PRECISION (BTPCF-S, DBTPCF-D)
C***KEYWORDS  INTERPOLATION, PIECEWISE POLYNOMIALS, SPLINES
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   DBTPCF is a subsidiary routine called by DB2INT and DB3INT.
C
C   DBTPCF computes NF sets of B-spline coefficients which determine
C   NF distinct piecewise polynomial interpolating functions of
C   order K.  The functions share a common set of interpolation
C   points X(i), i=1,..,N and knot sequence T(i), i=1,..,N+K.
C
C   The NF sets of N function values to be interpolated are stored in
C   the first NF columns of the array FCN.  On output, the
C   corresponding B-spline coefficients are stored in the ROWS of
C   the array BCOEF.  This is done to facilitate multidimensional
C   interpolation by tensor-products of one-dimensional B-splines.
C   See (de Boor, 1979) for details.
C
C
C   I N P U T
C   ---------
C
C   X       Double precision 1D array (size N)
C           Array of abscissae. Must be strictly increasing.
C
C   N       Integer scalar (N .GE. 2)
C           Number of abscissae.
C
C   FCN     Double precision 2D array (size LDF by NF)
C           Array of function values to interpolate. Each column of
C           FCN contains an independent set of N function values
C           corresponding to the array X of abscissae.  That is, the
C           jth data set is (X(i),FCN(i,j)), i=1,..,N.
C
C   LDF     Integer scalar (LDF .GE. N)
C           The actual leading dimension of FCN used in the calling
C           program.
C
C   NF      Integer scalar (NF .GE. 1)
C           Number of data sets (columns of FCN).
C
C   T       Double precision 1D array (size N+K)
C           The knot sequence.
C
C   K       Integer scalar (2 .LE. K .LE. N)
C           The order (degree + 1) of the piecewise polynomial
C           interpolant.
C
C
C   O U T P U T
C   -----------
C
C   BCOEF   Double precision 2D array (size NF by N)
C           Array of coefficients of the B-spline interpolants.
C           The coefficients of the B-spline interpolating the
C           data in the jth column of FCN are found in the jth ROW
C           of BCOEF.  That is, the B-spline coefficients defining
C           an interpolant to (X(i),FCN(i,j)), i=1,..,N are
C           stored in BCOEF(j,i), i=1,..,N.
C
C   LDB     Integer scalar (LDB .GE. NF)
C           The actual leading dimension of BCOEF used in the calling
C           program.
C
C
C   M I S C E L L A N E O U S
C   -------------------------
C
C   WORK    Double precision 1D array (size 2*K*(N+1))
C           Array of working storage.
C
C
C   CAUTION: The constraints on the input variables are not checked by
C            this routine.
C
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C               Carl de Boor, Efficient computer manipulation of tensor
C                 products, ACM Transactions on Mathematical Software 5,
C                 (1979), pp. 173-182.
C***ROUTINES CALLED  DBINTK, DBNSLV
C***REVISION HISTORY  (YYMMDD)
C   820525  DATE WRITTEN
C   900910  Modified to meet SLATEC standards.
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C   951130  Changed library to MSSL. (FNF)
C   000504  Removed warning comment about SLATEC externals; it no longer
C           applies, as all needed subordinates are now in MSSL. (ACH)
C***END PROLOGUE  DBTPCF
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          N, LDF, NF, K, LDB
      DOUBLE PRECISION X(N), FCN(LDF,NF), T(N+K), BCOEF(LDB,N), WORK(*)
C
C  .. LOCAL VARIABLES
C
      INTEGER          I, IQ, IW, J, K1, K2
C
      EXTERNAL         DBINTK, DBNSLV
C
C  ---------------------------------------------
C  CHECK FOR NULL INPUT AND PARTITION WORK ARRAY
C  ---------------------------------------------
C
C***FIRST EXECUTABLE STATEMENT  DBTPCF
      IF (NF .LE. 0)  GO TO 500
      K1 = K - 1
      K2 = K1 + K
      IQ = 1 + N
      IW = IQ + K2*N + 1
C
C  -----------------------------
C  COMPUTE B-SPLINE COEFFICIENTS
C  -----------------------------
C
C     ... FIRST DATA SET
C
      CALL DBINTK(X,FCN,T,N,K,WORK,WORK(IQ),WORK(IW))
      DO 20 I=1,N
         BCOEF(1,I) = WORK(I)
   20 CONTINUE
C
C     ... ALL REMAINING DATA SETS BY BACK-SUBSTITUTION
C
      IF (NF .GT. 1)  THEN
         DO 100 J=2,NF
            DO 50 I=1,N
               WORK(I) = FCN(I,J)
   50       CONTINUE
            CALL DBNSLV(WORK(IQ),K2,N,K1,K1,WORK)
            DO 60 I=1,N
               BCOEF(J,I) = WORK(I)
   60       CONTINUE
  100    CONTINUE
      ENDIF
C
C  ----
C  EXIT
C  ----
C
  500 CONTINUE
      RETURN
      END
*DECK DBINTK
      SUBROUTINE DBINTK (X, Y, T, N, K, BCOEF, Q, WORK)
C***BEGIN PROLOGUE  DBINTK
C***PURPOSE  Compute the B-representation of a spline which interpolates
C            given data.
C***LIBRARY   MSSL
C     >> Revised or new BSPLINE routine, provisionally accepted for SLATEC 5.0.
C***CATEGORY  E1A
C***TYPE      DOUBLE PRECISION (BINTK-S, DBINTK-D)
C***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C
C         DBINTK is a modification of the SPLINT routine of reference
C         [de Boor 1979].
C
C         DBINTK produces the B-spline coefficients, BCOEF, of the
C         B-spline of order K with knots T(I), I=1,...,N+K, which
C         takes on the value Y(I) at X(I), I=1,...,N.  The spline or
C         any of its derivatives can be evaluated by calls to DB1VAL.
C
C         The I-th equation of the linear system A*BCOEF = B for the
C         coefficients of the interpolant enforces interpolation at
C         X(I), I=1,...,N.  Hence, B(I) = Y(I), for all I, and A is
C         a band matrix with 2K-1 bands if A is invertible.  The matrix
C         A is generated row by row and stored, diagonal by diagonal,
C         in the rows of Q, with the main diagonal going into row K.
C         The banded system is then solved by a call to DBNFAC (which
C         constructs the triangular factorization for A and stores it
C         again in Q), followed by a call to DBNSLV (which then
C         obtains the solution BCOEF by substitution).  DBNFAC does no
C         pivoting, since the total positivity of the matrix A makes
C         this unnecessary.  The linear system to be solved is
C         (theoretically) invertible if and only if
C                 T(I) .LT. X(I) .LT. T(I+K),        for all I.
C         Equality is permitted on the left for I=1 and on the right
C         for I=N when K knots are used at X(1) or X(N).  Otherwise,
C         violation of this condition is certain to lead to an error.
C
C     Description of Arguments
C
C         Input     (X,Y,T are double precision)
C           X      - vector of length N containing data point abscissae
C                    in strictly increasing order.
C           Y      - corresponding vector of length N containing data
C                    point ordinates.
C           T      - knot vector of length N+K
C                    Since T(1),..,T(K) .LE. X(1) and T(N+1),..,T(N+K)
C                    .GE. X(N), this leaves only N-K knots (not nec-
C                    essarily X(I) values) interior to (X(1),X(N))
C           N      - number of data points, N .GE. K
C           K      - order of the spline, K .GE. 1
C
C         Output    (BCOEF,Q,WORK are double precision)
C           BCOEF  - a vector of length N containing the B-spline
C                    coefficients
C           Q      - a work vector of length (2*K-1)*N, containing
C                    the triangular factorization of the coefficient
C                    matrix of the linear system being solved.  The
C                    coefficients for the interpolant of another data
C                    set (X(I),YY(I)), I=1,...,N, with the same
C                    abscissae can be obtained by loading YY into
C                    BCOEF and then executing
C                        CALL DBNSLV(Q,2K-1,N,K-1,K-1,BCOEF)
C           WORK   - work vector of length 2*K
C
C     Error Conditions
C         Improper input is a fatal error, as is a singular system of
C         equations.  Conditions checked and their XERMSG error numbers:
C               1   K is less than 1
C               2   N is less than K
C               3   X values are not distinct or not ordered
C               4   Some abscissa was not in the support of the corres-
C                   ponding basis function, and the system is singular.
C               8   The system solver has detected a singular system,
C                   although the theoretical conditions for a solution
C                   were satisfied.
C         Note: These numbers are printed by XERMSG, but not returned
C         to the calling program.
C
C***REFERENCES  D. E. Amos, Computation with splines and B-splines,
C                 Report SAND78-1968, Sandia Laboratories, March 1979.
C               Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C               Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C***ROUTINES CALLED  DBNFAC, DBNSLV, DBSPVN, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900910  Changed BVALU references to B1VAL in DESCRIPTION.  (RFB)
C   901114  Miscellaneous housekeeping changes.  (FNF)
C           1. Improved readability of prologue and reduced single/
C              double differences.
C           2. Revised declarations.
C           3. Changed XERMSG calls so that different errors have
C              different error numbers, and recorded the numbers in
C              the prologue.
C   920501  Reformatted the REFERENCES section.  (WRB)
C   930527  Modified error return coding to look more like DBINKP. (FNF)
C   951130  Changed library to MSSL. (FNF)
C   000503  Changed library to MSSL in error message calls. (ACH)
C   000504  Removed warning comment about SLATEC externals; it no longer
C           applies, as all needed subordinates are now in MSSL. (ACH)
C***END PROLOGUE  DBINTK
C
C  Implementation note:  The error numbers occur in the order of
C     detection.  The number 8 is historical.  Prior to 901114, the
C     rest all printed number 2.  The errors at statements 80 and 90
C     were detected by de Boor's original SPLINT and the rest were
C     added later by Amos.  (FNF)
C
C  Additional note: This routine is now ready to have IFLAG added to
C     its argument list, as in DBINKP.
C**End
C
      INTEGER  N, K
      DOUBLE PRECISION  X(*), Y(*), T(*), BCOEF(*), Q(*), WORK(*)
C     DIMENSION Q(2*K-1,N), T(N+K)
C
      EXTERNAL  DBNFAC, DBNSLV, DBSPVN, XERMSG
C
      INTEGER  I, IFLAG, ILP1MX, IWORK, J, JJ, KM1, KPKM2, LEFT, LENQ,
     +   NP1
      DOUBLE PRECISION  XI
C
C***FIRST EXECUTABLE STATEMENT  DBINTK
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
        Q(I) = 0
   10 CONTINUE
C
C  ***   LOOP OVER I TO CONSTRUCT THE  N  INTERPOLATION EQUATIONS
      DO 50 I=1,N
        XI = X(I)
        ILP1MX = MIN(I+K,NP1)
C        *** FIND  LEFT  IN THE CLOSED INTERVAL (I,I+K-1) SUCH THAT
C                T(LEFT) .LE. X(I) .LT. T(LEFT+1)
C        MATRIX IS SINGULAR IF THIS IS NOT POSSIBLE
        LEFT = MAX(LEFT,I)
        IF (XI.LT.T(LEFT)) GO TO 80
   20   IF (XI.LT.T(LEFT+1)) GO TO 30
        LEFT = LEFT + 1
        IF (LEFT.LT.ILP1MX) GO TO 20
        LEFT = LEFT - 1
        IF (XI.GT.T(LEFT+1)) GO TO 80
C        *** THE I-TH EQUATION ENFORCES INTERPOLATION AT XI, HENCE
C        A(I,J) = B(J,K,T)(XI), ALL J. ONLY THE  K  ENTRIES WITH  J =
C        LEFT-K+1,...,LEFT ACTUALLY MIGHT BE NONZERO. THESE  K  NUMBERS
C        ARE RETURNED, IN  BCOEF (USED FOR TEMP. STORAGE HERE), BY THE
C        FOLLOWING
   30   CALL DBSPVN(T, K, K, 1, XI, LEFT, BCOEF, WORK, IWORK)
C        WE THEREFORE WANT  BCOEF(J) = B(LEFT-K+J)(XI) TO GO INTO
C        A(I,LEFT-K+J), I.E., INTO  Q(I-(LEFT+J)+2*K,(LEFT+J)-K) SINCE
C        A(I+J,J)  IS TO GO INTO  Q(I+K,J), ALL I,J,  IF WE CONSIDER  Q
C        AS A TWO-DIM. ARRAY , WITH  2*K-1  ROWS (SEE COMMENTS IN
C        DBNFAC). IN THE PRESENT PROGRAM, WE TREAT  Q  AS AN EQUIVALENT
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
      CALL DBNFAC(Q, K+KM1, N, KM1, KM1, IFLAG)
      GO TO (60, 90), IFLAG
C     *** SOLVE  A*BCOEF = Y  BY BACKSUBSTITUTION
   60 DO 70 I=1,N
        BCOEF(I) = Y(I)
   70 CONTINUE
      CALL DBNSLV(Q, K+KM1, N, KM1, KM1, BCOEF)
      IFLAG = 0
      RETURN
C
C     ERROR RETURNS
C
   80 CONTINUE
      IFLAG = 4
      CALL XERMSG ('MSSL', 'DBINTK',
     +   'SOME ABSCISSA WAS NOT IN THE SUPPORT OF THE CORRESPONDING ' //
     +   'BASIS FUNCTION AND THE SYSTEM IS SINGULAR.', IFLAG, 1)
      RETURN
   90 CONTINUE
      IFLAG = 8
      CALL XERMSG ('MSSL', 'DBINTK',
     +   'THE SYSTEM OF SOLVER DETECTS A SINGULAR SYSTEM, ALTHOUGH ' //
     +   'THE THEORETICAL CONDITIONS FOR A SOLUTION WERE SATISFIED.',
     +   IFLAG, 1)
      RETURN
  100 CONTINUE
      IFLAG = 1
      CALL XERMSG ('MSSL', 'DBINTK',
     +   'K DOES NOT SATISFY K.GE.1', IFLAG, 1)
      RETURN
  105 CONTINUE
      IFLAG = 2
      CALL XERMSG ('MSSL', 'DBINTK',
     +   'N DOES NOT SATISFY N.GE.K', IFLAG, 1)
      RETURN
  110 CONTINUE
      IFLAG = 3
      CALL XERMSG ('MSSL', 'DBINTK',
     +   'X(I) DOES NOT SATISFY X(I).LT.X(I+1) FOR SOME I', IFLAG, 1)
      RETURN
      END
*DECK DBNFAC
      SUBROUTINE DBNFAC (W, NROWW, NROW, NBANDL, NBANDU, IFLAG)
C***BEGIN PROLOGUE  DBNFAC
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBINT4 and DBINTK
C***LIBRARY   MSSL
C***TYPE      DOUBLE PRECISION (BNFAC-S, DBNFAC-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C  DBNFAC is the BANFAC routine from
C        * A Practical Guide to Splines *  by C. de Boor
C
C  DBNFAC is a double precision routine
C
C  Returns in  W  the LU-factorization (without pivoting) of the banded
C  matrix  A  of order  NROW  with  (NBANDL + 1 + NBANDU) bands or diag-
C  onals in the work array  W .
C
C *****  I N P U T  ****** W is double precision
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
C *****  O U T P U T  ****** W is double precision
C  IFLAG.....Integer indicating success( = 1) or failure ( = 2) .
C     If  IFLAG = 1, then
C  W.....contains the LU-factorization of  A  into a unit lower triangu-
C        lar matrix  L  and an upper triangular matrix  U (both banded)
C        and stored in customary fashion over the corresponding entries
C        of  A . This makes it possible to solve any particular linear
C        system  A*X = B  for  X  by a
C              CALL DBNSLV ( W, NROWW, NROW, NBANDL, NBANDU, B )
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
C     The routine should NOT be used for an arbitrary banded matrix.
C
C***SEE ALSO  DBINT4, DBINTK
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   000501  Changed library to MSSL (ACH)
C***END PROLOGUE  DBNFAC
C
      INTEGER IFLAG, NBANDL, NBANDU, NROW, NROWW, I, IPK, J, JMAX, K,
     1 KMAX, MIDDLE, MIDMK, NROWM1
      DOUBLE PRECISION W(NROWW,*), FACTOR, PIVOT
C
C***FIRST EXECUTABLE STATEMENT  DBNFAC
      IFLAG = 1
      MIDDLE = NBANDU + 1
C                         W(MIDDLE,.) CONTAINS THE MAIN DIAGONAL OF  A .
      NROWM1 = NROW - 1
C      IF (NROWM1) 120, 110, 10
      IF (NROWM1.LT.0) THEN
         GOTO 120
      ELSE IF (NROWM1.EQ.0) THEN
         GOTO 110
      ELSE
         GOTO 10
      ENDIF
   10 IF (NBANDL.GT.0) GO TO 30
C                A IS UPPER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO .
      DO 20 I=1,NROWM1
        IF (W(MIDDLE,I).EQ.0.0D0) GO TO 120
   20 CONTINUE
      GO TO 110
   30 IF (NBANDU.GT.0) GO TO 60
C              A IS LOWER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO AND
C                 DIVIDE EACH COLUMN BY ITS DIAGONAL .
      DO 50 I=1,NROWM1
        PIVOT = W(MIDDLE,I)
        IF (PIVOT.EQ.0.0D0) GO TO 120
        JMAX = MIN(NBANDL,NROW-I)
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
        IF (PIVOT.EQ.0.0D0) GO TO 120
C                 JMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN COLUMN  I
C                     BELOW THE DIAGONAL .
        JMAX = MIN(NBANDL,NROW-I)
C              DIVIDE EACH ENTRY IN COLUMN  I  BELOW DIAGONAL BY PIVOT .
        DO 70 J=1,JMAX
          W(MIDDLE+J,I) = W(MIDDLE+J,I)/PIVOT
   70   CONTINUE
C                 KMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN ROW  I  TO
C                     THE RIGHT OF THE DIAGONAL .
        KMAX = MIN(NBANDU,NROW-I)
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
  110 IF (W(MIDDLE,NROW).NE.0.0D0) RETURN
  120 IFLAG = 2
      RETURN
      END
*DECK DBNSLV
      SUBROUTINE DBNSLV (W, NROWW, NROW, NBANDL, NBANDU, B)
C***BEGIN PROLOGUE  DBNSLV
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBINT4 and DBINTK
C***LIBRARY   MSSL
C***TYPE      DOUBLE PRECISION (BNSLV-S, DBNSLV-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C  DBNSLV is the BANSLV routine from
C        * A Practical Guide to Splines *  by C. de Boor
C
C  DBNSLV is a double precision routine
C
C  Companion routine to  DBNFAC . It returns the solution  X  of the
C  linear system  A*X = B  in place of  B , given the LU-factorization
C  for  A  in the work array  W from DBNFAC.
C
C *****  I N P U T  ****** W,B are DOUBLE PRECISION
C  W, NROWW,NROW,NBANDL,NBANDU.....Describe the LU-factorization of a
C        banded matrix  A  of order  NROW  as constructed in  DBNFAC .
C        For details, see  DBNFAC .
C  B.....Right side of the system to be solved .
C
C *****  O U T P U T  ****** B is DOUBLE PRECISION
C  B.....Contains the solution  X , of order  NROW .
C
C *****  M E T H O D  ******
C     (With  A = L*U, as stored in  W,) the unit lower triangular system
C  L(U*X) = B  is solved for  Y = U*X, and  Y  stored in  B . Then the
C  upper triangular system  U*X = Y  is solved for  X  . The calcul-
C  ations are so arranged that the innermost loops stay within columns.
C
C***SEE ALSO  DBINT4, DBINTK
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   000501  Changed library to MSSL (ACH)
C***END PROLOGUE  DBNSLV
C
      INTEGER NBANDL, NBANDU, NROW, NROWW, I, J, JMAX, MIDDLE, NROWM1
      DOUBLE PRECISION W(NROWW,*), B(*)
C***FIRST EXECUTABLE STATEMENT  DBNSLV
      MIDDLE = NBANDU + 1
      IF (NROW.EQ.1) GO TO 80
      NROWM1 = NROW - 1
      IF (NBANDL.EQ.0) GO TO 30
C                                 FORWARD PASS
C            FOR I=1,2,...,NROW-1, SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN
C            OF  L )  FROM RIGHT SIDE  (BELOW I-TH ROW) .
      DO 20 I=1,NROWM1
        JMAX = MIN(NBANDL,NROW-I)
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
      JMAX = MIN(NBANDU,I-1)
      DO 70 J=1,JMAX
        B(I-J) = B(I-J) - B(I)*W(MIDDLE-J,I)
   70 CONTINUE
      I = I - 1
      IF (I.GT.1) GO TO 60
   80 B(1) = B(1)/W(MIDDLE,1)
      RETURN
      END
*DECK DBSPVN
      SUBROUTINE DBSPVN (T, JHIGH, K, INDEX, X, ILEFT, VNIKX, WORK,
     +   IWORK)
C***BEGIN PROLOGUE  DBSPVN
C***PURPOSE  Calculate the value of all (possibly) nonzero basis
C            functions at X.
C***LIBRARY   MSSL
C***CATEGORY  E3, K6
C***TYPE      DOUBLE PRECISION (BSPVN-S, DBSPVN-D)
C***KEYWORDS  EVALUATION OF B-SPLINE
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Abstract    **** a double precision routine ****
C         DBSPVN is the BSPLVN routine of the reference.
C
C         DBSPVN calculates the value of all (possibly) nonzero basis
C         functions at X of order MAX(JHIGH,(J+1)*(INDEX-1)), where T(K)
C         .LE. X .LE. T(N+1) and J=IWORK is set inside the routine on
C         the first call when INDEX=1.  ILEFT is such that T(ILEFT) .LE.
C         X .LT. T(ILEFT+1).  A call to DINTRV(T,N+1,X,ILO,ILEFT,MFLAG)
C         produces the proper ILEFT.  DBSPVN calculates using the basic
C         algorithm needed in DBSPVD.  If only basis functions are
C         desired, setting JHIGH=K and INDEX=1 can be faster than
C         calling DBSPVD, but extra coding is required for derivatives
C         (INDEX=2) and DBSPVD is set up for this purpose.
C
C         Left limiting values are set up as described in DBSPVD.
C
C     Description of Arguments
C
C         Input      T,X are double precision
C          T       - knot vector of length N+K, where
C                    N = number of B-spline basis functions
C                    N = sum of knot multiplicities-K
C          JHIGH   - order of B-spline, 1 .LE. JHIGH .LE. K
C          K       - highest possible order
C          INDEX   - INDEX = 1 gives basis functions of order JHIGH
C                          = 2 denotes previous entry with work, IWORK
C                              values saved for subsequent calls to
C                              DBSPVN.
C          X       - argument of basis functions,
C                    T(K) .LE. X .LE. T(N+1)
C          ILEFT   - largest integer such that
C                    T(ILEFT) .LE. X .LT.  T(ILEFT+1)
C
C         Output     VNIKX, WORK are double precision
C          VNIKX   - vector of length K for spline values.
C          WORK    - a work vector of length 2*K
C          IWORK   - a work parameter.  Both WORK and IWORK contain
C                    information necessary to continue for INDEX = 2.
C                    When INDEX = 1 exclusively, these are scratch
C                    variables and can be used for other purposes.
C
C     Error Conditions
C         Improper input is a fatal error.
C
C***REFERENCES  Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  Calls to XERROR changed to calls to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   000501  Changed library to MSSL (ACH)
C***END PROLOGUE  DBSPVN
C
      INTEGER ILEFT, IMJP1, INDEX, IPJ, IWORK, JHIGH, JP1, JP1ML, K, L
      DOUBLE PRECISION T, VM, VMPREV, VNIKX, WORK, X
C     DIMENSION T(ILEFT+JHIGH)
      DIMENSION T(*), VNIKX(*), WORK(*)
C     CONTENT OF J, DELTAM, DELTAP IS EXPECTED UNCHANGED BETWEEN CALLS.
C     WORK(I) = DELTAP(I), WORK(K+I) = DELTAM(I), I = 1,K
C***FIRST EXECUTABLE STATEMENT  DBSPVN
      IF(K.LT.1) GO TO 90
      IF(JHIGH.GT.K .OR. JHIGH.LT.1) GO TO 100
      IF(INDEX.LT.1 .OR. INDEX.GT.2) GO TO 105
      IF(X.LT.T(ILEFT) .OR. X.GT.T(ILEFT+1)) GO TO 110
      GO TO (10, 20), INDEX
   10 IWORK = 1
      VNIKX(1) = 1.0D0
      IF (IWORK.GE.JHIGH) GO TO 40
C
   20 IPJ = ILEFT + IWORK
      WORK(IWORK) = T(IPJ) - X
      IMJP1 = ILEFT - IWORK + 1
      WORK(K+IWORK) = X - T(IMJP1)
      VMPREV = 0.0D0
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
      CALL XERMSG ('MSSL', 'DBSPVN', 'K DOES NOT SATISFY K.GE.1', 2,
     +   1)
      RETURN
  100 CONTINUE
      CALL XERMSG ('MSSL', 'DBSPVN',
     +   'JHIGH DOES NOT SATISFY 1.LE.JHIGH.LE.K', 2, 1)
      RETURN
  105 CONTINUE
      CALL XERMSG ('MSSL', 'DBSPVN', 'INDEX IS NOT 1 OR 2', 2, 1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('MSSL', 'DBSPVN',
     +   'X DOES NOT SATISFY T(ILEFT).LE.X.LE.T(ILEFT+1)', 2, 1)
      RETURN
      END
