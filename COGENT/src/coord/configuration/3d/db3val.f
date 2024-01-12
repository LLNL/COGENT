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
*DECK DB3VAL
      DOUBLE PRECISION FUNCTION DB3VAL (XVAL, YVAL, ZVAL, IDX, IDY, IDZ,
     +   TX, TY, TZ, NX, NY, NZ, KX, KY, KZ, FCN, LDF1, LDF2,
     +   ICONT, IWORK, WORK, IFLAG)
C***BEGIN PROLOGUE  DB3VAL
C***PURPOSE  DB3VAL evaluates the three-dimensional piecewise polynomial
C            interpolating function constructed by the routine DB3INT.
C            Either function values or partial derivative values may be
C            be requested.
C***LIBRARY   MSSL
C     >> Revised or new BSPLINE routine, provisionally accepted for SLATEC 5.0.
C***CATEGORY  E3, K6
C***TYPE      DOUBLE PRECISION (B3VAL-S, DB3VAL-D)
C***KEYWORDS  THREE DIMENSIONS, GRIDDED DATA, SPLINES,
C             PIECEWISE POLYNOMIALS, EVALUATION, DIFFERENTIATION
C***AUTHOR  Boisvert, R. F., (NIST)
C             Computing and Applied Mathematics Laboratory
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   DB3VAL evaluates the tensor product piecewise polynomial interpolant
C   constructed by the routine DB3INT or one of its derivatives at the
C   point (XVAL,YVAL,ZVAL).  The variables IDX, IDY and IDZ indicate
C   which function is to be evaluated.  If B(x,y,z) is the interpolant,
C
C                            IDX+IDY+IDZ
C                           d
C            DB1VAL  =   ----------------- B (XVAL,YVAL,ZVAL)
C                          IDX   IDY   IDZ
C                        dx    dy    dz
C
C   Thus, to evaluate the interpolant itself, set IDX=IDY=IDZ=0. To get
C   the first partial derivative with respect to x, set IDX=1 and
C   IDY=IDZ=0, and so on.

C   Since B is a piecewise polynomial of degree KX-1 in x, KY-1 in y and
C   KZ-1 in z, DB3VAL returns zero whenever IDX.GE.KX or IDY.GE.KY or
C   IDZ.GE.KZ.  DB3VAL also returns zero if (XVAL,YVAL,ZVAL) is out of
C   range, that is, if
C
C            XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR.
C            YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY) .OR.
C            ZVAL.LT.TZ(1) .OR. ZVAL.GT.ZY(NZ+KZ).
C
C   If the knots TX, TY, and TZ were chosen by DB3INT, then this is
C   equivalent to
C
C              XVAL.LT.X(1) .OR. XVAL.GT.X(NX) .OR.
C              YVAL.LT.Y(1) .OR. YVAL.GT.Y(NY) .OR.
C              ZVAL.LT.Z(1) .OR. ZVAL.GT.Z(NZ).
C
C   The input quantities TX, TY, TZ, NX, NY, NZ, KX, KY, KZ, and FCN
C   should be unchanged since the call to DB3INT which produced FCN.
C
C   Note that the derivative functions computed by DB3VAL may be
C   discontinuous when x, y or z correspond to knots.  (When DB3INT
C   selects knots this occurs only when IDX=KX-1, IDY=KY-1 or IDZ=KZ-1.)
C   In these cases DB3VAL returns right limiting values (right
C   derivatives), except at the rightmost knot where left limiting
C   values are returned.
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
C   ZVAL    Double precision scalar
C           Z coordinate of evaluation point.
C
C   IDX     Integer scalar  (IDX .GE. 0)
C           Indicates the x derivative of piecewise polynomial to
C           evaluate: IDX=J for the Jth partial derivative with respect
C           to x.  DB3VAL will return 0 if IDX.GE.KX.
C
C   IDY     Integer scalar  (IDY .GE. 0)
C           Indicates the y derivative of piecewise polynomial to
C           evaluate: IDY=J for the Jth partial derivative with respect
C           to y.  DB3VAL will return 0 if IDY.GE.KY.
C
C   IDZ     Integer scalar  (IDZ .GE. 0)
C           Indicates the z derivative of piecewise polynomial to
C           evaluate: IDZ=J for the Jth partial derivative with respect
C           to z.  DB3VAL will return 0 if IDZ.GE.KZ.
C
C   TX      Double precision 1D array (size NX+KX)
C           Sequence of knots defining the piecewise polynomial in
C           the x direction.  (Unchanged since the call to DB3INT which
C           produced FCN.)
C
C   TY      Double precision 1D array (size NY+KY)
C           Sequence of knots defining the piecewise polynomial in
C           the y direction.  (Unchanged since the call to DB3INT which
C           produced FCN.)
C
C   TZ      Double precision 1D array (size NZ+KZ)
C           Sequence of knots defining the piecewise polynomial in
C           the z direction.  (Unchanged since the call to DB3INT which
C           produced FCN.)
C
C   NX      Integer scalar (NX .GE. KX)
C           The number of interpolation points in x. (Unchanged since
C           the call to DB3INT which produced FCN.)
C
C   NY      Integer scalar (NY .GE. KY)
C           The number of interpolation points in y. (Unchanged since
C           the call to DB3INT which produced FCN.)
C
C   NZ      Integer scalar (NZ .GE. KZ)
C           The number of interpolation points in z. (Unchanged since
C           the call to DB3INT which produced FCN.)
C
C   KX      Integer scalar (KX .GE. 2)
C           Order of polynomial pieces in x. (Unchanged since the call
C           to DB3INT which produced FCN.)
C
C   KY      Integer scalar (KY .GE. 2)
C           Order of polynomial pieces in y. (Unchanged since the call
C           to DB3INT which produced FCN.)
C
C   KZ      Integer scalar (KZ .GE. 2)
C           Order of polynomial pieces in z. (Unchanged since the call
C           to DB3INT which produced FCN.)
C
C   FCN     Double precision 3D array (size LDF1 by LDF2 by NZ)
C           The B-spline coefficients computed by DB3INT.
C
C   LDF1    Integer scalar (LDF1 .GE. NX)
C           The actual leading dimension of FCN used in the calling
C           program.
C
C   LDF2    Integer scalar (LDF2 .GE. NY)
C           The actual second dimension of FCN used in the calling
C           program.
C
C
C   I N P U T / O U T P U T
C   -----------------------
C
C   ICONT   Integer scalar
C           An input flag.  Set ICONT=0 on the first call. Subsequent
C           calls evaluating the same piecewise polynomial function
C           (i.e., the same TX, TY, TZ, NX, NY, NZ, KX, KY, KZ, and
C           FCN) should have ICONT=1.  (See notes on efficiency below.)
C           As a convenience, DB3VAL sets ICONT=1 on output.
C
C   IWORK   Integer 1D array (size 10)
C           A working storage array.  NOTE: This array is used for
C           communication between successive calls of DB3VAL.  It must
C           not be modified between successive calls when ICONT=1.
C           Different interpolants require different work arrays.
C
C   WORK    Double precision 1D array (size KY*KZ + 3*max(KX,KY,KZ) + KZ + 2)
C           A working storage array.  NOTE: This array is used for
C           communication between successive calls of DB3VAL.  It must
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
C           2 : NX or LDF1 out of range
C           3 : KY out of range
C           4 : NY or LDF2 out of range
C           5 : KZ out of range
C           6 : NZ out of range
C           7 : IDX, IDY or IDZ out of range
C
C           Note: Error conditions are checked in the order shown here.
C                 Only the first error detected is reported.
C
C
C   N O T E S   O N   E F F I C I E N C Y
C   -------------------------------------
C
C   DB3VAL is designed so that it will be reasonably efficient when it is
C   being used to evaluate the fitted function on a grid.  The most
C   favorable situation occurs when DB3VAL is called repeatedly with the
C   same XVAL, the same IDX and IDY, and slowly increasing YVAL and ZVAL.
C   When calling DB3VAL in a triple loop over x, y and z :
C
C     (a) vary x in the OUTER loop, and
C     (b) vary z in the INNER loop, and
C     (c) use separate loops or distinct ICONT,IWORK,WORK for
C         each (IDX,IDY) pair
C
C   A typical example of this usage is as follows.
C
C      ICONT = 0
C      DO 30 I=1,NX
C         DO 20 J=1,NY
C            DO 10 K=1,NZ
C               G(I,J) = DB3VAL(X(I),Y(J),Z(K),IDX,IDY,IDZ,TX,TY,TZ,
C     +                        NX,NY,NZ,KX,KY,KZ,FCN,LDF1,LDF2,
C     +                        ICONT,IWORK,WORK,IFLAG)
C               IF (IFLAG .NE. 0) GO TO 9999
C   10       CONTINUE
C   20    CONTINUE
C   30 CONTINUE
C
C   Note that ICONT=0 initially.  This signals DB3VAL that it is starting
C   a new function.  DB3VAL sets ICONT=1 on output, so ICONT=1 on all
C   subsequent calls, which tells DB3VAL to attempt to reuse information
C   from the previous call to speed the evaluation process.
C
C
C   B A C K G R O U N D
C   -------------------
C
C   DB3VAL evaluates the following function or its partial derivatives
C
C                      NX   NY   NZ
C        B(x,y,z)  =  SUM  SUM  SUM  FCN(i,j,k) U (x) V (y) W (z)
C                     i=1  j=1  k=1              i     j     k
C   where
C
C     U (x)  is the ith (one-dimensional) B-spline basis function
C      i     defined by NX, KX, and TX,
C
C     V (y)  is the jth (one-dimensional) B-spline basis function
C      j     defined by NY, KY, and TY, and
C
C     W (z)  is the kth (one-dimensional) B-spline basis function
C      k     defined by NZ, KZ, and TZ.
C
C   See (de Boor, 1978) for a description of the B-spline basis.
C
C  The summation above can be rewritten as
C
C                                  NZ
C  (1)               B(x,y,z)  =  SUM  b  V (z),
C                                 k=1   k  k
C  where
C                         NY
C  (2)            b  =   SUM  c   V (y),   k = 1,..,NZ,
C                  k     j=1   jk  j
C  and
C                     NX
C  (3)        c  =   SUM  FCN(i,j,k) U (x),   j = 1,..,NY,
C              jk    i=1              i       k = 1,..,NZ.
C
C  Note that each summation is the evaluation of a one-dimensional
C  B-spline, which can be done by calls to DB1VAL.  At most KZ basis
C  functions in (1) are nonzero.  The indices of these functions are
C  determined, and only the coefficients b(k) which multiply the
C  nonzero basis functions are be computed using (2).  We also
C  observe that at most KY basis functions in (2) are nonzero.  The
C  indices of these functions also are determined and only those
C  coefficients c(j,k) which multiply nonzero basis functions are
C  computed using (3).
C
C
C   A C K N O W L E D G E M E N T
C   -----------------------------
C
C   Thanks to Fred Fritsch of Lawrence Livermore National Laboratory
C   for his critical review of this code which lead to many
C   improvements.
C
C***SEE ALSO  DB3INT
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
C***END PROLOGUE  DB3VAL
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          IDX, IDY, IDZ, NX, NY, NZ, KX, KY, KZ, LDF1,
     +                 LDF2, ICONT, IWORK(10), IFLAG
      DOUBLE PRECISION XVAL, YVAL, ZVAL, TX(NX+KX), TY(NY+KY),
     +                 TZ(NZ+KZ), FCN(LDF1,LDF2,NZ), WORK(*)
C
C  .. LOCAL VARIABLES
C
      CHARACTER*6      LIBRAR, SUBROU
      PARAMETER        (LIBRAR = 'MSSL', SUBROU = 'DB3VAL')
      CHARACTER*50     MESSAG
      LOGICAL          NEEDCO
      INTEGER          I, IDLAST, IERR, IHAVEX, IHAVEY, ILOY, ILOZ,
     +                 INX, INY, INZ, IPOS, IW, IXY, IZ, J, JDLAST,
     +                 JMAX, JMIN, KMAX, KMIN, K, LEFTY, LEFTZ,
     +                 LYLAST, LZLAST, MFLAG
      DOUBLE PRECISION DB1VAL, XLAST, YLAST
C
      EXTERNAL         DINTRV, DB1VAL, XERMSG
C
C***FIRST EXECUTABLE STATEMENT  DB3VAL
      DB3VAL = 0
C
C  ------------------------------
C  RESET STATE FROM PREVIOUS CALL
C  ------------------------------
C
C  These state variables indicate state of DB3VAL at end of last call.
C  They tell us the range of validity of the coefficients b(k) in (2)
C  and c(j,k) in (3) computed in the previous call and saved in WORK.
C  We want to avoid recomputing these if possible.
C
C  IWORK( 1) == IHAVEX == 1 if coefficients c(j,k) in (3) are available
C  IWORK( 2) == IHAVEY == 1 if coefficients b(k) in (2) are available
C  IWORK( 3) == INX    == knot interval of XVAL (for DB1VAL)
C  IWORK( 4) == ILOY   == knot interval of YVAL (for DINTRV)
C  IWORK( 5) == ILOZ   == knot interval of ZVAL (for DINTRV)
C  IWORK( 6) == LYLAST == value of LEFTY from DINTRV for YVAL
C  IWORK( 7) == LZLAST == value of LEFTZ from DINTRV for ZVAL
C  IWORK( 8) == IDLAST == value of IDX
C  IWORK( 9) == JDLAST == value of IDY
C  IWORK(10) == KMIN   == smallest  index k of nonzero basis b(k) at ZVAL
C   WORK( 1) == XLAST  == value of x variable
C   WORK( 2) == YLAST  == value of y variable
C
      IF ( ICONT .NE. 1 ) THEN
         IHAVEX = 0
         IHAVEY = 0
         INX = 1
         ILOY = 1
         ILOZ = 1
         JMIN = 0
         KMIN = 0
      ELSE
         IHAVEX = IWORK(1)
         IHAVEY = IWORK(2)
         INX    = IWORK(3)
         ILOY   = IWORK(4)
         ILOZ   = IWORK(5)
         LYLAST = IWORK(6)
         LZLAST = IWORK(7)
         IDLAST = IWORK(8)
         JDLAST = IWORK(9)
         KMIN   = IWORK(10)
         XLAST  = WORK(1)
         YLAST  = WORK(2)
      ENDIF
C
C     ... set LEFTY, LEFTZ so that they have values when state is
C         saved after premature termination
C
      LEFTY = 0
      LEFTZ = 0
C
C  -------------
C  SPECIAL CASES
C  -------------
C
C     ... CHECK INPUT FOR ERRORS
C
      IFLAG = 0
      IF (KX. LT. 1) THEN
         IFLAG = 1
         MESSAG = 'KX IS OUT OF RANGE'
      ELSE IF ((NX .LT. KX) .OR. (NX .GT. LDF1)) THEN
         IFLAG = 2
         MESSAG = 'NX OR LDF1 IS OUT OF RANGE'
      ELSE IF (KY. LT. 1) THEN
         IFLAG = 3
         MESSAG = 'KY IS OUT OF RANGE'
      ELSE IF ((NY .LT. KY) .OR. (NY .GT. LDF2)) THEN
         IFLAG = 4
         MESSAG = 'NY OR LDF2 IS OUT OF RANGE'
      ELSE IF (KZ. LT. 1) THEN
         IFLAG = 5
         MESSAG = 'KZ IS OUT OF RANGE'
      ELSE IF (NZ .LT. KZ) THEN
         IFLAG = 6
         MESSAG = 'NZ IS OUT OF RANGE'
      ELSE IF ((IDX .LT. 0) .OR. (IDY .LT. 0) .OR. (IDZ .LT. 0)) THEN
         IFLAG = 7
         MESSAG = 'IDX, IDY OR IDZ IS OUT OF RANGE'
      ELSE
         IFLAG = 0
      ENDIF
C
      IF ( IFLAG .NE. 0 ) THEN
         CALL XERMSG(LIBRAR,SUBROU,MESSAG,IFLAG,1)
         GO TO 9999
      ENDIF
C
C     ... QUICK RETURN WHEN SPLINE EVALUATES TO ZERO
C
      IF ( (IDX  .GE. KX       ) .OR.
     +     (IDY  .GE. KY       ) .OR.
     +     (IDZ  .GT. KZ       ) .OR.
     +     (XVAL .LT. TX(1)    ) .OR.
     +     (XVAL .GT. TX(NX+KX)) .OR.
     +     (YVAL .LT. TY(1)    ) .OR.
     +     (YVAL .GT. TY(NY+KY)) .OR.
     +     (ZVAL .LT. TZ(1)    ) .OR.
     +     (ZVAL .GT. TZ(NZ+KZ))      ) GO TO 9999
C
C  ---------------------
C  EVALUATE THE B-SPLINE
C  ---------------------
C
C     ... PARTITION WORK ARRAY
C
      IXY = 3
      IZ = IXY + KY*KZ
      IW = IZ  + KZ
C
C     ... FIND KNOT INTERVAL CONTAINING ZVAL
C
      CALL DINTRV(TZ,NZ+KZ,ZVAL,ILOZ,LEFTZ,MFLAG)
      IF (MFLAG .NE. 0) THEN
C
C        ... ZVAL .EQ. T(NZ+KZ),  ADJUST TO GET LEFT LIMITING VALUE
C
   10    CONTINUE
         LEFTZ = LEFTZ - 1
         IF (ZVAL. EQ. TZ(LEFTZ)) GO TO 10
      ENDIF
C
C     ... FIND KNOT INTERVAL CONTAINING YVAL
C
      CALL DINTRV(TY,NY+KY,YVAL,ILOY,LEFTY,MFLAG)
      IF (MFLAG .NE. 0) THEN
C
C        ... YVAL .EQ. T(NY+KY),  ADJUST TO GET LEFT LIMITING VALUE
C
   20    CONTINUE
         LEFTY = LEFTY - 1
         IF (YVAL. EQ. TY(LEFTY)) GO TO 20
      ENDIF
C
C     ... CHECK IF WE ALREADY HAVE COEFFICIENTS c(j,k) IN (3)
C
      NEEDCO = (IHAVEX .EQ. 0     ) .OR.
     +         (LEFTZ  .NE. LZLAST) .OR.
     +         (LEFTY  .NE. LYLAST) .OR.
     +         (IDX    .NE. IDLAST) .OR.
     +         (XVAL   .NE. XLAST )
C
      IF ( NEEDCO ) THEN
C
C        ... FIND RANGE OF INDICES OF NONZERO BASIS FUNCTIONS IN (1)
C            (KMIN,KMAX) = (SMALLEST,LARGEST)
C
         IF (LEFTZ .LT. KZ) THEN
            KMIN = 1
            KMAX = KZ
         ELSEIF (LEFTZ .GT. NZ) THEN
            KMIN = NZ - KZ + 1
            KMAX = NZ
         ELSE
            KMIN = LEFTZ-KZ+1
            KMAX = LEFTZ
         ENDIF
C
C        ... FIND RANGE OF INDICES OF NONZERO BASIS FUNCTIONS IN (2)
C            (JMIN,JMAX) = (SMALLEST,LARGEST)
C
         IF (LEFTY .LT. KY) THEN
            JMIN = 1
            JMAX = KY
         ELSEIF (LEFTY .GT. NY) THEN
            JMIN = NY - KY + 1
            JMAX = NY
         ELSE
            JMIN = LEFTY-KY+1
            JMAX = LEFTY
         ENDIF
C
C        ... GET COEFFICIENTS OF NONZERO BASIS FUNCTIONS IN (3)
C
         IPOS = IXY - 1
         DO 50 K=KMIN,KMAX
            DO 51 J=JMIN,JMAX
               IPOS = IPOS + 1
               WORK(IPOS) = DB1VAL(XVAL,IDX,TX,NX,KX,FCN(1,J,K),INX,
     +              WORK(IW),IERR)
   51       CONTINUE     
   50    CONTINUE
         IHAVEX = 1
C
      ENDIF
C
C     ... CHECK IF WE ALREADY HAVE COEFFICIENTS b(j) IN (2)
C
      NEEDCO = NEEDCO               .OR.
     +         (IHAVEY .EQ. 0     ) .OR.
     +         (IDY    .NE. JDLAST) .OR.
     +         (YVAL   .NE. YLAST )
C
      IF ( NEEDCO ) THEN
C
C        ... FIND MIN INDEX OF NONZERO BASIS FUNCTIONS IN (2)
C
         IF (LEFTY .LT. KY) THEN
            JMIN = 1
         ELSEIF (LEFTY .GT. NY) THEN
            JMIN = NY - KY + 1
         ELSE
            JMIN = LEFTY-KY+1
         ENDIF
C
C        ... GET COEFFICIENTS OF NONZERO BASIS FUNCTIONS IN (2)
C
         INY = 1
         I = IXY - KY
         J = IZ - 1
         DO 60 K=1,KZ
            I = I + KY
            J = J + 1
            WORK(J) = DB1VAL(YVAL,IDY,TY(JMIN),KY,KY,WORK(I),INY,
     +                      WORK(IW),IERR)
  60     CONTINUE
         IHAVEY = 1
C
      ENDIF
C
C     ... EVALUATE INTERPOLANT USING (1)
C
      INZ = 1
      DB3VAL = DB1VAL(ZVAL,IDZ,TZ(KMIN),KZ,KZ,WORK(IZ),INZ,WORK(IW),
     +              IERR)
C
C  -------------------
C  SAVE STATE AND EXIT
C  -------------------
C
 9999 CONTINUE
      IWORK(1)  = IHAVEX
      IWORK(2)  = IHAVEY
      IWORK(3)  = INX
      IWORK(4)  = ILOY
      IWORK(5)  = ILOZ
      IWORK(6)  = LEFTY
      IWORK(7)  = LEFTZ
      IWORK(8)  = IDX
      IWORK(9)  = IDY
      IWORK(10) = KMIN
      WORK(1) = XVAL
      WORK(2) = YVAL
      ICONT = 1
C
      RETURN
      END
