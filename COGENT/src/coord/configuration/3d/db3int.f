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
*DECK DB3INT
      SUBROUTINE DB3INT (X, NX, Y, NY, Z, NZ, KX, KY, KZ, TX, TY, TZ,
     +   FCN, LDF1, LDF2, WORK, IFLAG)
C***BEGIN PROLOGUE  DB3INT
C***PURPOSE  DB3INT determines a piecewise polynomial function that
C            interpolates three-dimensional gridded data. Users  specify
C            the polynomial order (degree+1) of the interpolant and
C            (optionally) the knot sequence.
C***LIBRARY   MSSL
C     >> Revised or new BSPLINE routine, provisionally accepted for SLATEC 5.0.
C***CATEGORY  E2A
C***TYPE      DOUBLE PRECISION (B3INT-S, DB3INT-D)
C***KEYWORDS  INTERPOLATION, THREE DIMENSIONS, GRIDDED DATA, SPLINES,
C             PIECEWISE POLYNOMIALS
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   DB3INT determines the parameters of a function that interpolates the
C   three-dimensional  gridded data (X(i),Y(j),Z(k),FCN(i,j,k)) for
C   i=1,..,NX, j=1,..,NY, and k=1,..,NZ. The interpolating function and
C   its derivatives may subsequently be evaluated by function DB3VAL.
C
C   The interpolating function is a piecewise polynomial (pp) function
C   represented as a tensor product of one-dimensional B-splines.  The
C   form of this function is
C
C                        NX   NY   NZ
C          S(x,y,z)  =  SUM  SUM  SUM  a    U (x) V (y) W (z)
C                       i=1  j=1  k=1   ijk  i     j     k
C
C   where the functions U(i), V(j), and W(k) are one-dimensional B-
C   spline basis functions. The coefficients a(i,j,k) are chosen so
C
C   S(X(i),Y(j),Z(k))=FCN(i,j,k) for i=1,..,NX, j=1,..,NY, k=1,..,NZ.
C
C   Note that for fixed y and z S is a pp function of x alone, for
C   fixed x and z S is a pp function of y alone, and for fixed x and y
C   S is a pp function of z alone.  In one dimension a pp function may
C   be created by partitioning a given interval into subintervals and
C   defining a distinct polynomial on each one.  The points where
C   adjacent subintervals meet are called knots.  Each of the functions
C   U(i), V(j), and W(k) above is a piecewise polynomial.
C
C   Users of DB3INT choose the order (degree+1) of the polynomial pieces
C   used to define the interpolant in each of the x, y and z directions
C   (KX, KY, and KZ). Users also may define their own knot sequences in
C   x, y and z separately (TX, TY, and TZ).  If IFLAG=0, however, DB3INT
C   will choose knots that result in a pp interpolant with KX-2, KY-2
C   and KZ-2 continuous partial derivatives in x, y and z respectively.
C   The interpolating function is identically zero outside the rectan-
C   gular region defined by the knots.  See below for more information
C   on knot selection.
C
C   After a call to DB3INT, all information necessary to define the
C   interpolating function is contained in the parameters NX, NY, NZ,
C   KX, KY, KZ, TX, TY, TZ, and FCN. These quantities should not be
C   altered until after the last call of the evaluation routine DB3VAL.
C
C
C   I N P U T
C   ---------
C
C   X       Double precision 1D array (size NX)
C           Array of x abscissae. Must be strictly increasing.
C
C   NX      Integer scalar (2 .LE. NX .LE. LDF1)
C           Number of x abscissae.
C
C   Y       Double precision 1D array (size NY)
C           Array of y abscissae. Must be strictly increasing.
C
C   NY      Integer scalar (2 .LE. NY .LE. LDF2)
C           Number of y abscissae.
C
C   Z       Double precision 1D array (size NZ)
C           Array of z abscissae. Must be strictly increasing.
C
C   NZ      Integer scalar (NZ .GE. 2)
C           Number of z abscissae.
C
C   KX      Integer scalar (2 .LE. KX .LE. NX)
C           The order (degree + 1) of polynomial pieces in x.
C
C   KY      Integer scalar (2 .LE. KY .LE. NY)
C           The order (degree + 1) of polynomial pieces in y.
C
C   KZ      Integer scalar (2 .LE. KZ .LE. NZ)
C           The order (degree + 1) of polynomial pieces in z.
C
C
C   I N P U T   O R   O U T P U T
C   -----------------------------
C
C   TX      Double precision 1D array (size NX+KX)
C           The knots in the x direction. If IFLAG=0 these are chosen
C           by DB3INT. If IFLAG=1 they are specified by the user as
C           input. The knots must be non-decreasing and also satisfy
C                  TX(i).LT.X(i).LT.TX(i+KX),  i=1,..,NX.
C           Equality on the left is permitted for i=1 and equality on
C           the right is permitted for i=NX. (See information below on
C           knot placement.)
C
C   TY      Double precision 1D array (size NY+KY)
C           The knots in the y direction. If IFLAG=0 these are chosen
C           by DB3INT. If IFLAG=1 they are specified by the user as
C           input. The knots must be non-decreasing and also satisfy
C                  TY(i).LT.Y(i).LT.TY(i+KY),  i=1,..,NY.
C           Equality on the left is permitted for i=1 and equality on
C           the right is permitted for i=NY. (See information below on
C           knot placement.)
C
C   TZ      Double precision 1D array (size NZ+KZ)
C           The knots in the y direction. If IFLAG=0 these are chosen
C           by DB3INT. If IFLAG=1 they are specified by the user as
C           input. The knots must be non-decreasing and also satisfy
C                  TZ(i).LT.Z(i).LT.TZ(i+KZ),  i=1,..,NZ.
C           Equality on the left is permitted for i=1 and equality on
C           the right is permitted for i=NZ. (See information below on
C           knot placement.)
C
C
C   I N P U T   A N D   O U T P U T
C   -------------------------------
C
C   FCN     Double precision 3D array (size LDF1 by LDF2 by NZ)
C           Input : Array of function values to interpolate. FCN(I,J,K)
C                   should contain the function value at the point
C                   (X(I),Y(J),Z(K)).
C           Output: Array of coefficients of the B-spline interpolant.
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
C   O T H E R
C   ---------
C
C   WORK    Double precision 1D array
C           Array of working storage. Must be dimensioned of length at
C           least NX*NY*NZ + 2*max( KX*(NX+1), KY*(NY+1), KZ*(NZ+1))
C
C   IFLAG   Integer scalar.
C           Must be set by the user before DB3INT is called.
C           On return IFLAG indicates the status of the output.
C
C           Input :  0 : knot sequence chosen by user
C                    1 : knot sequence chosen by DB3INT
C
C           Output:  0 : successful execution
C                    2 : IFLAG out of range
C                    3 : NX or LDF1 out of range
C                    4 : KX out of range
C                    5 : X not strictly increasing
C                    6 : TX is an illegal knot sequence
C                    7 : NY or LDF2 out of range
C                    8 : KY out of range
C                    9 : Y not strictly increasing
C                   10 : TY is an illegal knot sequence
C                   11 : NZ out of range
C                   12 : KZ out of range
C                   13 : Z not strictly increasing
C                   14 : TZ is an illegal knot sequence
C
C           Note: Error conditions are checked in the order shown here.
C                 Only the first error detected is reported.
C
C           After a successful return, DB3INT may be recalled without
C           resetting IFLAG, provided that only the array FCN has
C           changed.
C
C
C   K N O T   S E L E C T I O N
C   ---------------------------
C
C   In this section we describe the relationship between data points
C   and knots.  Users who choose to let DB3INT select the knot sequences
C   TX, TY and TZ (by setting IFLAG=1) may skip this discussion.
C
C   To describe the relationship between knots and data points we first
C   consider the simpler case of one-dimensional interpolation; in
C   particular, we consider interpolating the data X(i), i=1,..,N with
C   a pp function of order K.
C
C   Knots are the points where the individual polynomial pieces join up,
C   and hence where the pp function may suffer a loss of smoothness.  To
C   define a pp function, one needs N+K knots. If the knots are distinct
C   the interpolant will be as smooth as possible (continuous, with K-2
C   continuous derivatives). If two adjacent knots come together, the
C   smoothness of the function is reduced at that point. In general, if
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
C   Equality is permitted on the left for i=1 and on the right for i=N.
C
C   The three-dimensional interpolant computed by this routine is a
C   tensor product of one-dimensional pp interpolants.  The knots form
C   a grid (TX(i),TY(j),TZ(k)) in the same way that the data points do.
C   Along lines parallel to the coordinate axes, the interpolant reduces
C   to a one-dimensional pp function with order and knots (KX,TX),
C   (KY,TY), or (KZ,TZ). In this case the appropriate constraints on the
C   knots become
C
C            TX(i) .LT. X(i) .LT. TX(i+KX)   i=1,..,NX
C            TY(i) .LT. Y(i) .LT. TY(i+KY)   i=1,..,NY
C            TZ(i) .LT. Z(i) .LT. TZ(i+KZ)   i=1,..,NZ
C
C   with equality on the left permitted when i=1 and equality on the
C   right permitted when i=NX, i=NY and i=NZ respectively.
C
C   If these conditions are violated, then DB3INT returns with IFLAG
C   equal to 6, 10 or 14.  The default knot sequence selected by DB3INT
C   always satisfies these conditions.
C
C   When the user sets IFLAG=1 DB3INT selects knots as follows. KX knots
C   are taken at each endpoint in the x direction, not-a-knot end
C   conditions (see references) are used, and the remaining knots are
C   placed at data points if KX is even or at midpoints between data
C   points if KX is odd.  The y and z directions are treated similarly.
C   This yields a three-dimensional pp function with KX-2, KY-2 and
C   KZ-2 continuous partial derivatives in x, y and z respectively. The
C   interpolant is zero outside the rectangular region defined by the
C   data points, and discontinuous along the boundary of this region.
C
C***SEE ALSO  DB3VAL
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C               Carl de Boor, Efficient computer manipulation of tensor
C                 products, ACM Transactions on Mathematical Software 5,
C                 (1979), pp. 173-182.
C***ROUTINES CALLED  DBKCHK, DBKNOT, DBTPCF, DBUPCK, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   820525  DATE WRITTEN
C   900910  Modified to meet SLATEC standards.
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C   951130  Changed library to MSSL. (FNF)
C   951218  Corrected CATEGORY line. (FNF)
C   000503  Changed library to MSSL in error message calls. (ACH)
C   000504  Removed warning comment about SLATEC externals; it no longer
C           applies, as all needed subordinates are now in MSSL. (ACH)
C***END PROLOGUE  DB3INT
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          NX, NY, NZ, KX, KY, KZ, LDF1, LDF2, IFLAG
      DOUBLE PRECISION X(NX), Y(NY), Z(NZ), TX(NX+KX), TY(NY+KY),
     +                 TZ(NZ+KZ), FCN(LDF1,LDF2,NZ), WORK(*)
C
C  .. LOCAL VARIABLES
C
      CHARACTER*6      LIBRAR, SUBROU
      PARAMETER        (LIBRAR = 'MSSL', SUBROU = 'DB3INT')
      CHARACTER*50     MESSAG
      INTEGER          I, IW, J, K, LOC
      LOGICAL          DBKCHK
C
      EXTERNAL         DBKCHK, DBKNOT, DBTPCF, DBUPCK, XERMSG
C
C***FIRST EXECUTABLE STATEMENT  DB3INT
C
C  -----------------------
C  CHECK VALIDITY OF INPUT
C  -----------------------
C
      IF ((IFLAG .LT. 0) .OR. (IFLAG .GT. 1))  GO TO 9020
C
      IF ((NX .LT. 2) .OR. (NX .GT. LDF1))  GO TO 9030
      IF ((KX .LT. 2) .OR. (KX .GT. NX))  GO TO 9040
      DO 10 I=2,NX
         IF (X(I) .LE. X(I-1))  GO TO 9050
   10 CONTINUE
      IF (IFLAG .EQ. 0) THEN
         IF(.NOT. DBKCHK(X,NX,KX,TX)) GO TO 9060
      ENDIF
C
      IF ((NY .LT. 2) .OR. (NY .GT. LDF2))  GO TO 9070
      IF ((KY .LT. 2) .OR. (KY .GT. NY))  GO TO 9080
      DO 20 I=2,NY
         IF (Y(I) .LE. Y(I-1))  GO TO 9090
   20 CONTINUE
      IF (IFLAG .EQ. 0) THEN
         IF(.NOT. DBKCHK(Y,NY,KY,TY)) GO TO 9100
      ENDIF
C
      IF (NZ .LT. 2)  GO TO 9110
      IF ((KZ .LT. 2) .OR. (KZ .GT. NZ))  GO TO 9120
      DO 30 I=2,NZ
         IF (Z(I) .LE. Z(I-1))  GO TO 9130
   30 CONTINUE
      IF (IFLAG .EQ. 0) THEN
         IF(.NOT. DBKCHK(Z,NZ,KZ,TZ)) GO TO 9140
      ENDIF
C
C  ------------
C  CHOOSE KNOTS
C  ------------
C
      IF (IFLAG .EQ. 1) THEN
         CALL DBKNOT(X,NX,KX,TX)
         CALL DBKNOT(Y,NY,KY,TY)
         CALL DBKNOT(Z,NZ,KZ,TZ)
      ENDIF
C
C  -------------------------------
C  CONSTRUCT B-SPLINE COEFFICIENTS
C  -------------------------------
C
      IW = NX*NY*NZ + 1
C
C     ... COPY FCN TO WORK IN PACKED FORM FOR DBTPCF
C
      LOC = 0
      DO 300 K=1,NZ
         DO 200 J=1,NY
            DO 100 I=1,NX
               WORK(I+LOC) = FCN(I,J,K)
  100       CONTINUE
            LOC = LOC + NX
  200    CONTINUE
  300 CONTINUE
C
C     ... TENSOR-PRODUCT INTERPOLATION
C
      CALL DBTPCF(X,NX,WORK,NX,NY*NZ,TX,KX,FCN, NY*NZ,WORK(IW))
      CALL DBTPCF(Y,NY,FCN, NY,NX*NZ,TY,KY,WORK,NX*NZ,WORK(IW))
      CALL DBTPCF(Z,NZ,WORK,NZ,NX*NY,TZ,KZ,FCN, NX*NY,WORK(IW))
C
C     ... UNPACK FCN
C
      CALL DBUPCK(FCN,NX,NY,NZ,LDF1,LDF2)
C
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
      MESSAG = 'NX OR LDF1 IS OUT OF RANGE'
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
      MESSAG = 'NY OR LDF2 IS OUT OF RANGE'
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
 9110 CONTINUE
      IFLAG = 11
      MESSAG = 'NZ IS OUT OF RANGE'
      GO TO 9900
C
 9120 CONTINUE
      IFLAG = 12
      MESSAG = 'KZ IS OUT OF RANGE'
      GO TO 9900
C
 9130 CONTINUE
      IFLAG = 13
      MESSAG = 'Z ARRAY MUST BE STRICTLY INCREASING'
      GO TO 9900
C
 9140 CONTINUE
      IFLAG = 14
      MESSAG = 'TZ IS AN ILLEGAL KNOT SEQUENCE'
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

*DECK DBUPCK
      SUBROUTINE DBUPCK (A, NX, NY, NZ, LD1, LD2)
C***BEGIN PROLOGUE  DBUPCK
C***SUBSIDIARY
C***PURPOSE  Converts an array packed as though dimensioned A(NX,NY,NZ)
C            to one dimensioned as A(LD1,LD2,NZ).
C***LIBRARY   MSSL
C     >> Revised or new BSPLINE routine, provisionally accepted for SLATEC 5.0.
C***TYPE      DOUBLE PRECISION (BUPCK-S, DBUPCK-D)
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   DBUPCK is a subsidiary routine called by DB3INT.
C
C   DBUPCK converts an array packed as though dimensioned as A(NX,NY,NZ)
C   to one dimensioned as A(LD1,LD2,NZ), where LD1.GE.NX and LD2.GE.NY.
C
C   The routine assumes that the input array is packed, that is, the
C   first NX*NY*NZ locations contain data.  It then unpacks the data
C   so that the calling program can dimension it as A(LD1,LD2,NZ).
C   The operation is done in place.
C
C
C   I N P U T
C   ---------
C
C   A       Double precision 1D array (length at least LD1*LD2*NZ)
C           On input, contains the array to be converted.
C           On output, contains the converted array.
C
C   NX      Integer (NX.GT.0)
C           First dimension of input array.
C
C   NY      Integer (NY.GT.0)
C           Second dimension of input array.
C
C   NZ      Integer (NZ.GT.0)
C           Third dimension of input array.
C
C   LD1     Integer (LD1.GE.NX)
C           First dimension of output array.
C
C   LD2     Integer (LD2.GE.NY)
C           Second dimension of output array.
C
C
C   CAUTION: The constraints on the input variables are not checked by
C            this routine.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   900910  DATE WRITTEN
C   930510  Clarified purpose.  (FNF)
C   930618  Reformatted the AUTHOR section.  (WRB)
C   951130  Changed library to MSSL. (FNF)
C   000504  Removed warning comment about SLATEC externals; it no longer
C           applies, as all needed subordinates are now in MSSL. (ACH)
C***END PROLOGUE  DBUPCK
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          LD1, LD2, NX, NY, NZ
      DOUBLE PRECISION A(*)
C
C  .. LOCAL VARIABLES
C
      INTEGER          I, IF2, IF3, IT2, IT3, J, K
C
C***FIRST EXECUTABLE STATEMENT  DBUPCK
      IF (( NZ .EQ.  1) .AND. (LD1 .EQ. NX))  GO TO 900
      IF ((LD1 .EQ. NX) .AND. (LD2 .EQ. NY))  GO TO 900
C
      DO 300 K=NZ,1,-1
         IF3 = (K-1)*NY - 1
         IT3 = (K-1)*LD2 - 1
         DO 200 J=NY,1,-1
            IF2 = (J + IF3)*NX
            IT2 = (J + IT3)*LD1
            DO 100 I=NX,1,-1
               A(I+IT2) = A(I+IF2)
  100       CONTINUE
  200    CONTINUE
  300 CONTINUE
C
  900 CONTINUE
      RETURN
      END
