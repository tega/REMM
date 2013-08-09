*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  f06.f
*     A subset of the NAG F06 Chapter with some modifications.
*     The routines perform the same function as the NAG F06 routines.
*
*                         Level 0  F06  Scalar routines
*                         -------  ---- ---------------
*     f06aaz+         f06baf/drot3g+  f06bcf/dcsg+    f06blf/ddiv+
*     f06bmf/dnorm+
*
*                         Level 1  BLAS Vector routines
*                         -------  ---- ---------------
*     daxpy           dcopy           ddot            dnrm2    
*     dscal           dswap           idamax          drot
*
*                         Level 1  F06  Vector routines
*                         -------  ---  ---------------
*     f06dbf/iload    f06fbf/dload    f06fcf/ddscl    f06dff/icopy      
*     f06fjf/dssq+    f06fkf          f06flf/dcond    f06klf/idrank+    
*     f06fqf          f06frf/dgrfg+
*
*                         Level 2  BLAS Matrix-vector routines
*                         -------  ---  ----------------------
*     dgemv           dger            dsymv           dsyr     
*     dtrmv           dtrsv
*
*                         Level 2  F06  Matrix routines
*                         -------  ---  ---------------
*     f06qff          f06qgf          f06qhf          f06qjf   
*     f06qkf          f06qnf          f06qrf          f06qsf   
*     f06qtf          f06qvf          f06qwf          f06qxf   
*     f06qzf
*
*    +Differs from the Nag F06 version.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*UPTODATE F06AAZTEXT
      SUBROUTINE F06AAZ( SRNAME, INFO )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 9/25/88.
C     .. Scalar Arguments ..
      INTEGER            INFO
      CHARACTER*13       SRNAME
C     ..
C
C  Purpose
C  =======
C
C  F06AAZ  is an error handler for the Level 2 BLAS routines.
C
C  It is called by the Level 2 BLAS routines if an input parameter is
C  invalid.
C
C  Parameters
C  ==========
C
C  SRNAME - CHARACTER*13.
C           On entry, SRNAME specifies the name of the routine which
C           called F06AAZ.
C
C  INFO   - INTEGER.
C           On entry, INFO specifies the position of the invalid
C           parameter in the parameter-list of the calling routine.
C
C
C  Auxiliary routine for Level 2 Blas.
C
C  Written on 20-July-1986.
C
C     .. Local Scalars ..
      CHARACTER*80       REC (1)
C     ..
C     .. Executable Statements ..
      WRITE (REC (1),99999) SRNAME, INFO
C
      RETURN
C
99999 FORMAT ( ' ** On entry to ', A13, ' parameter number ', I2,
     $         ' had an illegal value' )
C
C     End of F06AAZ.
C
      END
** END OF F06AAZTEXT
      SUBROUTINE F06BAF( X, Y, CS, SN )

      DOUBLE PRECISION   X, Y, CS, SN

      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/

C
C  Note: f06baf/drot3g is different from the Nag routine f06baf.
C
C  DROT3G  generates a plane rotation that reduces the vector (X, Y) to
C  the vector (A, 0),  where A is defined as follows...
C
C     If both X and Y are negligibly small, or
C     if Y is negligible relative to Y,
C     then  A = X,  and the identity rotation is returned.
C
C     If X is negligible relative to Y,
C     then  A = Y,  and the swap rotation is returned.
C
C     Otherwise,  A = sign(X) * sqrt( X**2 + Y**2 ).
C
C  In all cases,  X and Y are overwritten by A and 0,  and CS will lie
C  in the closed interval (0, 1).  Also,  the absolute value of CS and
C  SN (if nonzero) will be no less than the machine precision,  EPS.
C
C  DROT3G  guards against overflow and underflow.
C  It is assumed that  FLMIN .lt. EPS**2  (i.e.  RTMIN .lt. EPS).
C
C  Systems Optimization Laboratory, Stanford University.
C  Original version dated January 1982.
C  F77 version dated 28-June-1986.
C  This version of DROT3G dated 28-June-1986.
C
      DOUBLE PRECISION   A, B, EPS, ONE, RTMIN, ZERO
      LOGICAL            FIRST
      INTRINSIC          ABS, MAX, SQRT
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      SAVE               FIRST , EPS   , RTMIN
      DATA               FIRST / .TRUE. /

      IF( FIRST )THEN
         FIRST = .FALSE.
         EPS    = WMACH(3)
         RTMIN  = WMACH(6)
      END IF

      IF (Y .EQ. ZERO) THEN

         CS = ONE
         SN = ZERO

      ELSE IF (X .EQ. ZERO) THEN

         CS = ZERO
         SN = ONE
         X  = Y

      ELSE

         A      = ABS(X)
         B      = ABS(Y)
         IF (MAX(A,B) .LE. RTMIN) THEN
            CS = ONE
            SN = ZERO
         ELSE
            IF (A .GE. B) THEN
               IF (B .LE. EPS*A) THEN
                  CS = ONE
                  SN = ZERO
                  GO TO 900
               ELSE
                  A  = A * SQRT( ONE + (B/A)**2 )
               END IF
            ELSE
               IF (A .LE. EPS*B) THEN
                  CS = ZERO
                  SN = ONE
                  X  = Y
                  GO TO 900
               ELSE
                  A  = B * SQRT( ONE + (A/B)**2 )
               END IF
            END IF
            IF (X .LT. ZERO) A = - A
            CS = X/A
            SN = Y/A
            X  = A
         END IF
      END IF

  900 Y  = ZERO

      RETURN

*     End of  F06BAF (DROT3G).

      END
*UPTODATE F06BCFTEXT
      SUBROUTINE F06BCF( T, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-602 (MAR 1988).
C     .. Scalar Arguments ..
      DOUBLE PRECISION   C, S, T
C     ..
C
C  F06BCF returns values c and s such that
C
C     c = cos( theta ),   s = sin( theta )
C
C  for a given value of
C
C     t = tan( theta ).
C
C  c is always non-negative and s has the same sign as t, so that
C
C     c = 1.0/sqrt( 1.0 + t**2 ),   s = t/sqrt( 1.0 + t**2 ).
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 28-February-1986.
C     Sven Hammarling, Nag Central Office.
C  -- Modified 19-August-1987.
C     Sven Hammarling and Jeremy Du Croz, Nag Central Office.
C        No longer sets s to zero when t is less than eps.
C  -- Modified 24-July-1991.
C     Philip E. Gill, UCSD.
C        Modified to call mchpar instead of x02ajf
C
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE = 1.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   ABST, EPS, REPS, RRTEPS, RTEPS
      LOGICAL            FIRST
C     .. External Functions ..
C+    DOUBLE PRECISION   X02AJF
C+    EXTERNAL           X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
C     .. Save statement ..
      SAVE               FIRST, EPS, REPS, RTEPS, RRTEPS
C     .. Data statements ..
      DATA               FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( FIRST )THEN
         FIRST  = .FALSE.
         EPS    = WMACH(3)
C+       EPS    = X02AJF( )
         REPS   =  1/EPS
         RTEPS  =  SQRT( EPS )
         RRTEPS =  1/RTEPS
      END IF
C
      ABST = ABS( T )
      IF( ABST.LT.RTEPS )THEN
         C = ONE
         S = T
      ELSE IF( ABST.GT.RRTEPS )THEN
         C = 1/ABST
         S = SIGN( ONE, T )
      ELSE
         C = 1/SQRT( 1 + ABST**2 )
         S = C*T 
      END IF
C
      RETURN
C
C     End of F06BCF. ( SCSG )
C
      END    
** END OF F06BCFTEXT
*UPTODATE F06BLFTEXT
      DOUBLE PRECISION FUNCTION F06BLF( A, B, FAIL )
      DOUBLE PRECISION          DDIV
      ENTRY                     DDIV  ( A, B, FAIL )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 9/25/88.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                  A, B
      LOGICAL                           FAIL
C     ..
C
C  F06BLF returns the value div given by
C
C     div = ( a/b                 if a/b does not overflow,
C           (
C           ( 0.0                 if a .eq. 0.0,
C           (
C           ( sign( a/b )*flmax   if a .ne. 0.0  and a/b would overflow,
C
C  where  flmax  is a large value, via the function name. In addition if
C  a/b would overflow then  fail is returned as true, otherwise  fail is
C  returned as false.
C
C  Note that when  a and b  are both zero, fail is returned as true, but
C  div  is returned as  0.0. In all other cases of overflow  div is such
C  that  abs( div ) = flmax.
C
C  When  b = 0  then  sign( a/b )  is taken as  sign( a ).
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 26-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      ABSB, DIV, FLMAX, FLMIN
      LOGICAL               FIRST
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
C     .. Intrinsic Functions ..
      INTRINSIC             ABS, SIGN
C     .. Save statement ..
      SAVE                  FIRST, FLMIN, FLMAX
C     .. Data statements ..
      DATA                  FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( A.EQ.ZERO )THEN
         DIV = ZERO
         IF( B.EQ.ZERO )THEN
            FAIL = .TRUE.
         ELSE
            FAIL = .FALSE.
         END IF
      ELSE
C
         IF( FIRST )THEN
            FIRST  = .FALSE.
            FLMIN  = WMACH( 5 )
            FLMAX  = WMACH( 7 )
         END IF
C
         IF( B.EQ.ZERO )THEN
            DIV  =  SIGN( FLMAX, A )
            FAIL = .TRUE.
         ELSE
            ABSB = ABS( B )
            IF( ABSB.GE.ONE )THEN
               FAIL = .FALSE.
               IF( ABS( A ).GE.ABSB*FLMIN )THEN
                  DIV = A/B
               ELSE
                  DIV = ZERO
               END IF
            ELSE
               IF( ABS( A ).LE.ABSB*FLMAX )THEN
                  FAIL = .FALSE.
                  DIV  =  A/B
               ELSE
                  FAIL = .TRUE.
                  DIV  = FLMAX
                  IF( ( ( A.LT.ZERO ).AND.( B.GT.ZERO ) ).OR.
     $                ( ( A.GT.ZERO ).AND.( B.LT.ZERO ) )     )
     $               DIV = -DIV
               END IF
            END IF
         END IF
      END IF
C
      F06BLF = DIV
      RETURN
C
C     End of F06BLF. ( DDIV )
C
      END
** END OF F06BLFTEXT
*UPTODATE F06BMFTEXT
      DOUBLE PRECISION FUNCTION F06BMF( SCALE, SSQ )
      DOUBLE PRECISION          DNORM
      ENTRY                     DNORM ( SCALE, SSQ )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 9/25/88.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                  SCALE, SSQ
C     ..
C
C  F06BMF returns the value norm given by
C
C     norm = ( scale*sqrt( ssq ), scale*sqrt( ssq ) .lt. flmax
C            (
C            ( flmax,             scale*sqrt( ssq ) .ge. flmax
C
C  via the function name.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 22-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION      FLMAX, NORM, SQT
      LOGICAL               FIRST
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
C     .. Intrinsic Functions ..
      INTRINSIC             SQRT
C     .. Save statement ..
      SAVE                  FIRST, FLMAX
C     .. Data statements ..
      DATA                  FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( FIRST )THEN
         FIRST = .FALSE.
         FLMAX = WMACH( 7 )
      END IF
C
      SQT = SQRT( SSQ )
      IF( SCALE.LT.FLMAX/SQT )THEN
         NORM = SCALE*SQT
      ELSE
         NORM = FLMAX
      END IF
C
      F06BMF = NORM
      RETURN
C
C     End of F06BMF. ( DNORM )
C
      END
** END OF F06BMFTEXT
*UPTODATE F06ECFTEXT
      SUBROUTINE F06ECF( N, ALPHA, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 9/25/88.
C     .. Entry Points ..
      ENTRY      DAXPY ( N, ALPHA, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
C     ..
C
C  F06ECF performs the operation
C
C     y := alpha*x + y
C
C
C  Nag Fortran 77 version of the Blas routine DAXPY.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 3-September-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.NE.ZERO )THEN
            IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
               DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  Y( IX ) = ALPHA*X( IX ) + Y( IX )
   10          CONTINUE
            ELSE
               IF( INCY.GE.0 )THEN
                  IY = 1
               ELSE
                  IY = 1 - ( N - 1 )*INCY
               END IF
               IF( INCX.GT.0 )THEN
                  DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     Y( IY ) = ALPHA*X( IX ) + Y( IY )
                     IY      = IY            + INCY
   20             CONTINUE
               ELSE
                  IX = 1 - ( N - 1 )*INCX
                  DO 30, I = 1, N
                     Y( IY ) = ALPHA*X( IX ) + Y( IY )
                     IX      = IX            + INCX
                     IY      = IY            + INCY
   30             CONTINUE
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06ECF. ( DAXPY )
C
      END
** END OF F06ECFTEXT
*UPTODATE F06EFFTEXT
      SUBROUTINE F06EFF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DCOPY ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
C     ..
C
C  F06EFF performs the operation
C
C     y := x
C
C
C  Nag Fortran 77 version of the Blas routine DCOPY.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
            DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
               Y( IY ) = X( IY )
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCY.GT.0 )THEN
               DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
                  Y( IY ) = X( IX )
                  IX      = IX      + INCX
   20          CONTINUE
            ELSE
               IY = 1 - ( N - 1 )*INCY
               DO 30, I = 1, N
                  Y( IY ) = X( IX )
                  IY      = IY      + INCY
                  IX      = IX      + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06EFF. ( DCOPY )
C
      END
** END OF F06EFFTEXT
*UPTODATE F06EAFTEXT
      DOUBLE PRECISION FUNCTION F06EAF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      DOUBLE PRECISION          DDOT
      ENTRY                     DDOT  ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER                           INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION                  X( * ), Y( * )
C     ..
C
C  F06EAF returns the value
C
C     F06EAF = x'y
C
C
C  Nag Fortran 77 version of the Blas routine DDOT.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 21-September-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ZERO
      PARAMETER           ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      SUM
      INTEGER               I, IX, IY
C     ..
C     .. Executable Statements ..
      SUM = ZERO
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               SUM = SUM + X( IX )*Y( IX )
   10       CONTINUE
         ELSE
            IF( INCY.GE.0 )THEN
               IY = 1
            ELSE
               IY = 1 - ( N - 1 )*INCY
            END IF
            IF( INCX.GT.0 )THEN
               DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  SUM = SUM + X( IX )*Y( IY )
                  IY  = IY  + INCY
   20          CONTINUE
            ELSE
               IX = 1 - ( N - 1 )*INCX
               DO 30, I = 1, N
                  SUM = SUM + X( IX )*Y( IY )
                  IX  = IX  + INCX
                  IY  = IY  + INCY
   30          CONTINUE
            END IF
         END IF
      END IF
C
      F06EAF = SUM
      RETURN
C
C     End of F06EAF. ( DDOT )
C
      END
** END OF F06EAFTEXT
*UPTODATE F06EJFTEXT
      DOUBLE PRECISION FUNCTION F06EJF( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      DOUBLE PRECISION          DNRM2
      ENTRY                     DNRM2 ( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                           INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION                  X( * )
C     ..
C
C  F06EJF returns the euclidean norm of a vector via the function
C  name, so that
C
C     F06EJF := sqrt( x'*x )
C
C
C  Nag Fortran 77 version of the Blas routine DNRM2.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 25-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      NORM, SCALE, SSQ
C     .. External Functions ..
      DOUBLE PRECISION      F06BMF
      EXTERNAL              F06BMF
C     .. External Subroutines ..
      EXTERNAL              F06FJF
C     .. Intrinsic Functions ..
      INTRINSIC             ABS
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = ABS( X( 1 ) )
      ELSE
         SCALE = ZERO
         SSQ   = ONE
         CALL F06FJF( N, X, INCX, SCALE, SSQ )
         NORM  = F06BMF( SCALE, SSQ )
      END IF
C
      F06EJF = NORM
      RETURN
C
C     End of F06EJF. ( DNRM2 )
C
      END
** END OF F06EJFTEXT
*UPTODATE F06EDFTEXT
      SUBROUTINE F06EDF( N, ALPHA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DSCAL ( N, ALPHA, X, INCX )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06EDF performs the operation
C
C     x := alpha*x
C
C
C  Nag Fortran 77 version of the Blas routine DSCAL.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.EQ.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   10       CONTINUE
         ELSE IF( ALPHA.EQ.( -ONE ) )THEN
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = -X( IX )
   20       CONTINUE
         ELSE IF( ALPHA.NE.ONE )THEN
            DO 30, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ALPHA*X( IX )
   30       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06EDF. ( DSCAL )
C
      END
** END OF F06EDFTEXT
*UPTODATE F06EGFTEXT
      SUBROUTINE F06EGF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DSWAP ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
C     ..
C
C  F06EGF performs the operations
C
C     temp := x,   x := y,   y := temp.
C
C
C  Nag Fortran 77 version of the Blas routine DSWAP.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
            DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
               TEMP    = X( IY )
               X( IY ) = Y( IY )
               Y( IY ) = TEMP
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCY.GT.0 )THEN
               DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
                  TEMP    = X( IX )
                  X( IX ) = Y( IY )
                  Y( IY ) = TEMP
                  IX      = IX      + INCX
   20          CONTINUE
            ELSE
               IY = 1 - ( N - 1 )*INCY
               DO 30, I = 1, N
                  TEMP    = X( IX )
                  X( IX ) = Y( IY )
                  Y( IY ) = TEMP
                  IY      = IY      + INCY
                  IX      = IX      + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06EGF. ( DSWAP )
C
      END
** END OF F06EGFTEXT
*UPTODATE F06JLFTEXT
      INTEGER FUNCTION F06JLF( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      INTEGER          IDAMAX
      ENTRY            IDAMAX( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                  INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION         X( * )
C     ..
C
C  F06JLF returns the smallest value of i such that
C
C     abs( x( i ) ) = max( abs( x( j ) ) )
C                      j
C
C
C  Nag Fortran 77 version of the Blas routine IDAMAX.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 31-May-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION         XMAX
      INTEGER                  I, IMAX, IX
C     .. Intrinsic Functions ..
      INTRINSIC                ABS
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IMAX = 1
         IF( N.GT.1 )THEN
            XMAX = ABS( X( 1 ) )
            IX   = 1
            DO 10, I = 2, N
               IX = IX + INCX
               IF( XMAX.LT.ABS( X( IX ) ) )THEN
                  XMAX = ABS( X( IX ) )
                  IMAX = I
               END IF
   10       CONTINUE
         END IF
      ELSE
         IMAX = 0
      END IF
C
      F06JLF = IMAX
      RETURN
C
C     End of F06JLF. ( IDAMAX )
C
      END
** END OF F06JLFTEXT
*UPTODATE F06EPFTEXT
      SUBROUTINE F06EPF( N, X, INCX, Y, INCY, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DROT  ( N, X, INCX, Y, INCY, C, S )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   C, S
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
C     ..
C
C  F06EPF performs the plane rotation
C
C     ( x  y ) = ( x  y )*( c  -s ).
C                         ( s   c )
C
C
C  Nag Fortran 77 version of the Blas routine DROT.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 23-January-1984.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( S.NE.ZERO ).OR.( C.NE.ONE ) )THEN
            IF( ( C.EQ.ZERO ).AND.( S.EQ.ONE ) )THEN
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP1   = -X( IX )
                     X( IX ) =  Y( IX )
                     Y( IX ) =  TEMP1
   10             CONTINUE
               ELSE
                  IF( INCY.GE.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  IF( INCX.GT.0 )THEN
                     DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                        TEMP1   = -X( IX )
                        X( IX ) =  Y( IY )
                        Y( IY ) =  TEMP1
                        IY      =  IY       + INCY
   20                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 30, I = 1, N
                        TEMP1   = -X( IX )
                        X( IX ) =  Y( IY )
                        Y( IY ) =  TEMP1
                        IX      =  IX      + INCX
                        IY      =  IY      + INCY
   30                CONTINUE
                  END IF
               END IF
            ELSE IF( ( C.EQ.ZERO ).AND.( S.EQ.( -ONE ) ) )THEN
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 40, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP1   =  X( IX )
                     X( IX ) = -Y( IX )
                     Y( IX ) =  TEMP1
   40             CONTINUE
               ELSE
                  IF( INCY.GE.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  IF( INCX.GT.0 )THEN
                     DO 50, IX = 1, 1 + ( N - 1 )*INCX, INCX
                        TEMP1   =  X( IX )
                        X( IX ) = -Y( IY )
                        Y( IY ) =  TEMP1
                        IY      =  IY       + INCY
   50                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 60, I = 1, N
                        TEMP1   =  X( IX )
                        X( IX ) = -Y( IY )
                        Y( IY ) =  TEMP1
                        IX      =  IX      + INCX
                        IY      =  IY      + INCY
   60                CONTINUE
                  END IF
               END IF
            ELSE
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 70, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP1   = X( IX )
                     X( IX ) = S*Y( IX ) + C*TEMP1
                     Y( IX ) = C*Y( IX ) - S*TEMP1
   70             CONTINUE
               ELSE
                  IF( INCY.GE.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  IF( INCX.GT.0 )THEN
                     DO 80, IX = 1, 1 + ( N - 1 )*INCX, INCX
                        TEMP1   = X( IX )
                        X( IX ) = S*Y( IY ) + C*TEMP1
                        Y( IY ) = C*Y( IY ) - S*TEMP1
                        IY      = IY        + INCY
   80                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 90, I = 1, N
                        TEMP1   = X( IX )
                        X( IX ) = S*Y( IY ) + C*TEMP1
                        Y( IY ) = C*Y( IY ) - S*TEMP1
                        IX      = IX        + INCX
                        IY      = IY        + INCY
   90                CONTINUE
                  END IF
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06EPF. ( DROT )
C
      END
** END OF F06EPFTEXT
*UPTODATE F06DBFTEXT
      SUBROUTINE F06DBF( N, CONST, X, INCX )
      ENTRY      ILOAD ( N, CONST, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER            CONST, INCX, N
C     .. Array Arguments ..
      INTEGER            X( * )
C     ..
C                      
C  F06DBF performs the operation
C
C     x = const*e,   e' = ( 1  1 ... 1 ).
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 18-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( CONST.NE.0 )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = CONST
   10       CONTINUE
         ELSE
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = 0
   20       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06DBF. ( ILOAD )
C
      END
** END OF F06DBFTEXT
*UPTODATE F06FBFTEXT
      SUBROUTINE F06FBF( N, CONST, X, INCX )
      ENTRY      DLOAD ( N, CONST, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   CONST
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FBF performs the operation
C
C     x = const*e,   e' = ( 1  1 ... 1 ).
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
C                      
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( CONST.NE.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = CONST
   10       CONTINUE
         ELSE
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   20       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06FBF. ( DLOAD )
C
      END
** END OF F06FBFTEXT
*UPTODATE F06FCFTEXT
      SUBROUTINE F06FCF( N, D, INCD, X, INCX )
      ENTRY      DDSCL ( N, D, INCD, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER            INCD, INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), X( * )
C     ..
C
C  F06FCF performs the operation
C
C     x := diag( d )*x
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, ID, IX
C     .. External Subroutines ..
      EXTERNAL           DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC          ABS
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCD.EQ.0 ).AND.( INCX.NE.0 ) )THEN
            CALL DSCAL( N, D( 1 ), X, ABS( INCX ) )
         ELSE IF( ( INCD.EQ.INCX ).AND.( INCD.GT.0 ) )THEN
            DO 10, ID = 1, 1 + ( N - 1 )*INCD, INCD
               X( ID ) = D( ID )*X( ID )
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCD.GT.0 )THEN
               DO 20, ID = 1, 1 + ( N - 1 )*INCD, INCD
                  X( IX ) = D( ID )*X( IX )
                  IX      = IX              + INCX
   20          CONTINUE
            ELSE
               ID = 1 - ( N - 1 )*INCD
               DO 30, I = 1, N
                  X( IX ) = D( ID )*X( IX )
                  ID      = ID              + INCD
                  IX      = IX              + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06FCF. ( DDSCL )
C            
      END
** END OF F06FCFTEXT
*UPTODATE F06DFFTEXT
      SUBROUTINE F06DFF( N, X, INCX, Y, INCY )
      ENTRY      ICOPY ( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      INTEGER            X( * ), Y( * )
C     ..
C
C  F06DFF performs the operation
C
C     y := x
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 10-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
            DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
               Y( IY ) = X( IY )
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCY.GT.0 )THEN
               DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
                  Y( IY ) = X( IX )
                  IX      = IX      + INCX
   20          CONTINUE
            ELSE
               IY = 1 - ( N - 1 )*INCY
               DO 30, I = 1, N
                  Y( IY ) = X( IX )
                  IY      = IY      + INCY
                  IX      = IX      + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06DFF. ( ICOPY )
C
      END
** END OF F06DFFTEXT
*UPTODATE F06FJFTEXT
      SUBROUTINE F06FJF( N, X, INCX, SCALE, SUMSQ )
      ENTRY      DSSQ  ( N, X, INCX, SCALE, SUMSQ )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   SCALE, SUMSQ
      INTEGER            INCX, N
C     .. Array Arguments ..                   
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FJF returns the values scl and smsq such that
C
C     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
C
C  where x( i ) = X( 1 + ( i - 1 )*INCX ). The value of sumsq is assumed
C  to be at least unity and the value of smsq will then satisfy
C
C     1.0 .le. smsq .le. ( sumsq + n ) .
C
C  scale is assumed to be non-negative and scl returns the value
C
C     scl = max( scale, abs( x( i ) ) ) .
C
C  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
C  scl and smsq are overwritten on SCALE and SUMSQ respectively.
C
C  The routine makes only one pass through the vector X.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   ABSXI
      INTEGER            IX
C     .. Intrinsic Functions ..
      INTRINSIC          ABS
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SUMSQ = 1     + SUMSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ +       ( ABSXI/SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
C
C     End of F06FJF. ( DSSQ )
C
      END
** END OF F06FJFTEXT
*UPTODATE F06FKFTEXT
      DOUBLE PRECISION FUNCTION F06FKF( N, W, INCW, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER                           INCW, INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION                  W( * ), X( * )
C     ..
C
C  F06FKF returns the weighted euclidean norm of a vector via the
C  function name, so that
C
C     F06FKF := sqrt( x'*W*x ),   where   W = diag( w ).
C
C  The elements of w are assumed to be non-negative.
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 25-June-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      ABSYI, NORM, SCALE, SSQ
      INTEGER               I, IW, IX
C     .. External Functions ..
      DOUBLE PRECISION      F06BMF
      EXTERNAL              F06BMF
C     .. Intrinsic Functions ..
      INTRINSIC             ABS, SQRT
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = SQRT( W( 1 ) )*ABS( X( 1 ) )
      ELSE
         IF( INCW.GT.0 )THEN
            IW = 1
         ELSE
            IW = 1 - ( N - 1 )*INCW
         END IF
         IF( INCX.GT.0 )THEN
            IX = 1
         ELSE
            IX = 1 - ( N - 1 )*INCX
         END IF
         SCALE = ZERO
         SSQ   = ONE
         DO 10, I = 1, N
            IF( ( W( IW ).NE.ZERO ).AND.( X( IX ).NE.ZERO ) )THEN
               ABSYI = SQRT( W( IW ) )*ABS( X( IX ) )
               IF( SCALE.LT.ABSYI )THEN
                  SSQ   = 1     + SSQ*( SCALE/ABSYI )**2
                  SCALE = ABSYI
               ELSE
                  SSQ   = SSQ   +     ( ABSYI/SCALE )**2
               END IF
            END IF
            IW = IW + INCW
            IX = IX + INCX
   10    CONTINUE
         NORM = F06BMF( SCALE, SSQ )
      END IF
C
      F06FKF = NORM
      RETURN
C
C     End of F06FKF.
C
      END
** END OF F06FKFTEXT
*UPTODATE F06FLFTEXT
      SUBROUTINE F06FLF( N, X, INCX, XMAX, XMIN )
      ENTRY      DCOND ( N, X, INCX, XMAX, XMIN )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   XMAX, XMIN
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FLF returns the values xmax and xmin given by
C
C     xmax = max( abs( x( i ) ) ),   xmin = min( abs( x( i ) ) ).
C             i                              i
C
C  If n is less than unity then xmax and xmin are returned as zero.
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 27-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         XMAX = ZERO
         XMIN = ZERO
      ELSE
         XMAX = ABS( X( 1 ) )
         XMIN = XMAX
         DO 10 IX = 1 + INCX, 1 + ( N - 1 )*INCX, INCX
            XMAX = MAX( XMAX, ABS( X( IX ) ) )
            XMIN = MIN( XMIN, ABS( X( IX ) ) )
   10    CONTINUE
      END IF
C
      RETURN
C
C     End of F06FLF. ( DCOND )
C
      END
** END OF F06FLFTEXT
*UPTODATE F06KLFTEXT
      INTEGER FUNCTION F06KLF( N, X, INCX, TOL )
      INTEGER          IDRANK                    
      ENTRY            IDRANK( N, X, INCX, TOL )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 9/25/88.
C     .. Scalar Arguments ..
      DOUBLE PRECISION         TOL
      INTEGER                  INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION         X( * )
C     ..
C
C  F06KLF finds the first element of the n element vector x for which
C
C     abs( x( k ) ).le.tol*max( abs( x( 1 ) ), ..., abs( x( k - 1 ) ) )
C
C  and returns the value ( k - 1 ) in the function name F06KLF. If no
C  such k exists then F06KLF is returned as n.
C
C  If tol is supplied as less than zero then the value epsmch, where
C  epsmch is the relative machine precision, is used in place of tol.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 27-February-1986.
C     Sven Hammarling, Nag Central Office.
C                        
C     .. Parameters ..
      DOUBLE PRECISION         ZERO
      PARAMETER              ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION         TL, XMAX
      INTEGER                  IX, K
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
C     .. Intrinsic Functions ..
      INTRINSIC                ABS, MAX
C     ..
C     .. Executable Statements ..
      K = 0
      IF( N.GE.1 )THEN
         IX = 1
         IF( TOL.LT.ZERO )THEN
            TL = WMACH(3)
         ELSE
            TL = TOL
         END IF
         XMAX = ABS( X( IX ) )
C
C+       WHILE( K.LT.N )LOOP
   10    IF   ( K.LT.N )THEN
            IF( ABS( X( IX ) ).LE.TL*XMAX )
     $         GO TO 20
            XMAX = MAX( XMAX, ABS( X( IX ) ) )
            K    = K  + 1
            IX   = IX + INCX
            GO TO 10
         END IF
C+       END WHILE
C
      END IF
C
   20 F06KLF = K
      RETURN
C
C     End of F06KLF. ( IDRANK )
C
      END
** END OF F06KLFTEXT
*UPTODATE F06FQFTEXT
      SUBROUTINE F06FQF( PIVOT, DIRECT, N, ALPHA, X, INCX, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, N
      CHARACTER*1        DIRECT, PIVOT
C     .. Array Arguments ..
      DOUBLE PRECISION   C( * ), S( * ), X( * )
C     ..
C
C  F06FQF generates the parameters of an orthogonal matrix P such that
C
C     when   PIVOT = 'F' or 'f'   and   DIRECT = 'F' or 'f'
C     or     PIVOT = 'V' or 'v'   and   DIRECT = 'B' or 'b'
C
C        P*( alpha ) = ( beta ),
C          (   x   )   (   0  )
C
C     when   PIVOT = 'F' or 'f'   and   DIRECT = 'B' or 'b'
C     or     PIVOT = 'V' or 'v'   and   DIRECT = 'F' or 'f'
C
C        P*(   x   ) = (   0  ),
C          ( alpha ) = ( beta )
C
C  where alpha is a scalar and x is an n element vector.
C
C  When  PIVOT = 'F' or 'f'  ( fixed pivot )
C  and  DIRECT = 'F' or 'f'  ( forward sequence ) then
C        
C     P is given as the sequence of plane rotation matrices
C
C        P = P( n )*P( n - 1 )*...*P( 1 )
C
C     where P( k ) is a plane rotation matrix for the ( 1, k + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'V' or 'v'  ( variable pivot )
C  and  DIRECT = 'B' or 'b'  ( backward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( 1 )*P( 2 )*...*P( n )
C
C     where P( k ) is a plane rotation matrix for the ( k, k + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'F' or 'f'  ( fixed pivot )
C  and  DIRECT = 'B' or 'b'  ( backward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( 1 )*P( 2 )*...*P( n )
C
C     where P( k ) is a plane rotation matrix for the ( k, n + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'V' or 'v'  ( variable pivot )
C  and  DIRECT = 'F' or 'f'  ( forward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( n )*P( n - 1 )*...*P( 1 )
C
C     where P( k ) is a plane rotation matrix for the ( k, k + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  The routine returns the cosine, c( k ), and sine, s( k ) that define
C  the matrix P( k ), such that the two by two rotation part of P( k ),
C  R( k ), has the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  On entry, ALPHA must contain  the scalar alpha and on exit, ALPHA is
C  overwritten by beta. The cosines and sines are returned in the arrays
C  C and S and the vector x is overwritten by the tangents of the plane
C  rotations ( t( k ) = s( k )/c( k ) ).
C
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 19-April-1985.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, IX
C     .. External Subroutines ..
      EXTERNAL           F06BAF
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
            IX = 1 + ( N - 1 )*INCX
            IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
               DO 10, I = N, 2, -1
                  CALL F06BAF( X( IX - INCX ), X( IX ), C( I ), S( I ) )
                  IX = IX - INCX
   10          CONTINUE
               CALL F06BAF( ALPHA, X( IX ), C( 1 ), S( 1 ) )
            ELSE IF( ( PIVOT.EQ.'F' ).OR.( PIVOT.EQ.'f' ) )THEN
C
C              Here we choose c and s so that
C
C                 ( alpha ) := (  c  s )*( alpha  )
C                 (   0   )    ( -s  c ) ( x( i ) )
C
C              which is equivalent to
C
C                 (   0   ) := ( c  -s )*( x( i ) )
C                 ( alpha )    ( s   c ) ( alpha  )
C
C              and so we need to return  s( i ) = -s  in order to make
C              R( i ) look like
C
C                 R( i ) = (  c( i )  s( i ) ).
C                          ( -s( i )  c( i ) )
C
               DO 20, I = N, 1, -1
                  CALL F06BAF( ALPHA, X( IX ), C( I ), S( I ) )
                  S( I )  = -S( I )
                  X( IX ) = -X( IX )
                  IX      =  IX      - INCX
   20          CONTINUE
            END IF
         ELSE IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
            IX = 1
            IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
C
C              Here we choose c and s so that
C
C                 ( x( i + 1 ) ) := (  c  s )*( x( i + 1 ) )
C                 (    0       )    ( -s  c ) ( x( i )     )
C
C              which is equivalent to
C
C                 (    0       ) := ( c  -s )*( x( i )     )
C                 ( x( i + 1 ) )    ( s   c ) ( x( i + 1 ) )
C
C              and so we need to return  s( i ) = -s  in order to make
C              R( i ) look like
C
C                 R( i ) = (  c( i )  s( i ) ).
C                          ( -s( i )  c( i ) )
C
               DO 30, I = 1, N - 1
                  CALL F06BAF( X( IX + INCX ), X( IX ), C( I ), S( I ) )
                  S( I )  = -S( I )
                  X( IX ) = -X( IX )
                  IX      =  IX      + INCX
   30          CONTINUE
               CALL F06BAF( ALPHA, X( IX ), C( N ), S( N ) )
               S( N )  = -S( N )
               X( IX ) = -X( IX )
            ELSE IF( ( PIVOT.EQ.'F' ).OR.( PIVOT.EQ.'f' ) )THEN
               DO 40, I = 1, N
                  CALL F06BAF( ALPHA, X( IX ), C( I ), S( I ) )
                  IX = IX + INCX
   40          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06FQF. ( SSROTG )
C
      END
** END OF F06FQFTEXT
*UPTODATE F06FRFTEXT
      SUBROUTINE F06FRF( N, ALPHA, X, INCX, TOL, ZETA )
      ENTRY      DGRFG ( N, ALPHA, X, INCX, TOL, ZETA )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 9/25/88.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, TOL, ZETA
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FRF generates details of a generalized Householder reflection such
C  that
C
C     P*( alpha ) = ( beta ),   P'*P = I.
C       (   x   )   (   0  )
C
C  P is given in the form
C
C     P = I - ( zeta )*( zeta  z' ),
C             (   z  )
C
C  where z is an n element vector and zeta is a scalar that satisfies
C
C     1.0 .le. zeta .le. sqrt( 2.0 ).
C
C  zeta is returned in ZETA unless x is such that
C
C     max( abs( x( i ) ) ) .le. max( eps*abs( alpha ), tol )
C
C  where eps is the relative machine precision and tol is the user
C  supplied value TOL, in which case ZETA is returned as 0.0 and P can
C  be taken to be the unit matrix.
C
C  beta is overwritten on alpha and z is overwritten on x.
C  the routine may be called with  n = 0  and advantage is taken of the
C  case where  n = 1.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 30-August-1984.
C     Sven Hammarling, Nag Central Office.
C     This version dated 28-September-1984.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   BETA, EPS, SCALE, SSQ
      LOGICAL            FIRST
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
C     .. External Subroutines ..
      EXTERNAL           F06FJF, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
C     .. Save statement ..
      SAVE               EPS, FIRST
C     .. Data statements ..
      DATA               FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         ZETA = ZERO
      ELSE IF( ( N.EQ.1 ).AND.( X( 1 ).EQ.ZERO ) )THEN
         ZETA = ZERO
      ELSE
C
         IF( FIRST )THEN
            FIRST = .FALSE.
            EPS   =  WMACH(3)
         END IF
C
C        Treat case where P is a 2 by 2 matrix specially.
C
         IF( N.EQ.1 )THEN
C
C           Deal with cases where  ALPHA = zero  and
C           abs( X( 1 ) ) .le. max( EPS*abs( ALPHA ), TOL )  first.
C
            IF( ALPHA.EQ.ZERO )THEN
               ZETA   =  ONE
               ALPHA  =  ABS ( X( 1 ) )
               X( 1 ) = -SIGN( ONE, X( 1 ) )
            ELSE IF( ABS( X( 1 ) ).LE.MAX( EPS*ABS( ALPHA ), TOL ) )THEN
               ZETA   =  ZERO
            ELSE
               IF( ABS( ALPHA ).GE.ABS( X( 1 ) ) )THEN
                  BETA = ABS( ALPHA ) *SQRT( 1 + ( X( 1 )/ALPHA )**2 )
               ELSE
                  BETA = ABS( X( 1 ) )*SQRT( 1 + ( ALPHA/X( 1 ) )**2 )
               END IF
               ZETA = SQRT( ( ABS( ALPHA ) + BETA )/BETA )
               IF( ALPHA.GE.ZERO )
     $            BETA = -BETA
               X( 1 ) = -X( 1 )/( ZETA*BETA )
               ALPHA  = BETA
            END IF
         ELSE
C
C           Now P is larger than 2 by 2.
C
            SSQ   = ONE
            SCALE = ZERO
            CALL F06FJF( N, X, INCX, SCALE, SSQ )
C
C           Treat cases where  SCALE = zero,
C           SCALE .le. max( EPS*abs( ALPHA ), TOL )  and
C           ALPHA = zero  specially.
C           Note that  SCALE = max( abs( X( i ) ) ).
C
            IF( ( SCALE.EQ.ZERO ).OR.
     $          ( SCALE.LE.MAX( EPS*ABS( ALPHA ), TOL ) ) )THEN
               ZETA  = ZERO
            ELSE IF( ALPHA.EQ.ZERO )THEN
               ZETA  = ONE
               ALPHA = SCALE*SQRT( SSQ )
               CALL DSCAL( N, -1/ALPHA, X, INCX )
            ELSE
               IF( SCALE.LT.ABS( ALPHA ) )THEN
                  BETA = ABS( ALPHA )*SQRT( 1 + SSQ*( SCALE/ALPHA )**2 )
               ELSE
                  BETA = SCALE       *SQRT( SSQ +   ( ALPHA/SCALE )**2 )
               END IF
               ZETA = SQRT( ( BETA + ABS( ALPHA ) )/BETA )
               IF( ALPHA.GT.ZERO )
     $            BETA = -BETA
               CALL DSCAL( N, -1/( ZETA*BETA ), X, INCX )
               ALPHA = BETA
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06FRF. ( DGRFG )
C
      END
** END OF F06FRFTEXT
*UPTODATE F06PAFTEXT
      SUBROUTINE F06PAF( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DGEMV  performs one of the matrix-vector operations
C
C     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
C
C  where alpha and beta are scalars, x and y are vectors and A is an
C  m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
C
C              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
C
C              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
C           Before entry, the incremented array X must contain the
C           vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
C           Before entry with BETA non-zero, the incremented array Y
C           must contain the vector y. On exit, Y is overwritten by the
C           updated vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PAF/DGEMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
C     up the start points in  X  and  Y.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  y := alpha*A*x + y.
C
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
C
C        Form  y := alpha*A'*x + y.
C
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PAF (DGEMV ).
C
      END
** END OF F06PAFTEXT
*UPTODATE F06PMFTEXT
      SUBROUTINE F06PMF( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DGER   performs the rank 1 operation
C
C     A := alpha*x*y' + A,
C
C  where alpha is a scalar, x is an m element vector, y is an n element
C  vector and A is an m by n matrix.
C
C  Parameters
C  ==========
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( m - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the m
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y.
C           Unchanged on exit.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients. On exit, A is
C           overwritten by the updated matrix.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PMF/DGER  ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06PMF (DGER  ).
C
      END
** END OF F06PMFTEXT
*UPTODATE F06PCFTEXT
      SUBROUTINE F06PCF( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DSYMV  performs the matrix-vector  operation
C
C     y := alpha*A*x + beta*y,
C
C  where alpha and beta are scalars, x and y are n element vectors and
C  A is an n by n symmetric matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the array A is to be referenced as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of A
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of A
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the symmetric matrix and the strictly
C           lower triangular part of A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the symmetric matrix and the strictly
C           upper triangular part of A is not referenced.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y. On exit, Y is overwritten by the updated
C           vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO.EQ.'U' .OR. UPLO.EQ.'u').AND.
     $         .NOT.(UPLO.EQ.'L' .OR. UPLO.EQ.'l')      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 5
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PCF/DSYMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     Set up the start points in  X  and  Y.
C
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
C
C        Form  y  when A is stored in upper triangle.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               DO 50, I = 1, J - 1
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + A( I, J )*X( I )
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*A( J, J ) + ALPHA*TEMP2
   60       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 80, J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
               DO 70, I = 1, J - 1
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + A( I, J )*X( IX )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*A( J, J ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
   80       CONTINUE
         END IF
      ELSE
C
C        Form  y  when A is stored in lower triangle.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP1  = ALPHA*X( J )
               TEMP2  = ZERO
               Y( J ) = Y( J )       + TEMP1*A( J, J )
               DO 90, I = J + 1, N
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + A( I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120, J = 1, N
               TEMP1   = ALPHA*X( JX )
               TEMP2   = ZERO
               Y( JY ) = Y( JY )       + TEMP1*A( J, J )
               IX      = JX
               IY      = JY
               DO 110, I = J + 1, N
                  IX      = IX      + INCX
                  IY      = IY      + INCY
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + A( I, J )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PCF (DSYMV ).
C
      END
** END OF F06PCFTEXT
*UPTODATE F06PPFTEXT
      SUBROUTINE F06PPF( UPLO, N, ALPHA, X, INCX, A, LDA )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DSYR  ( UPLO, N, ALPHA, X, INCX, A, LDA )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, LDA, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
C     ..
C
C  Purpose
C  =======
C
C  DSYR   performs the symmetric rank 1 operation
C
C     A := alpha*x*x' + A,
C
C  where alpha is a real scalar, x is an n element vector and A is an
C  n by n symmetric matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the array A is to be referenced as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of A
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of A
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the symmetric matrix and the strictly
C           lower triangular part of A is not referenced. On exit, the
C           upper triangular part of the array A is overwritten by the
C           upper triangular part of the updated matrix.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the symmetric matrix and the strictly
C           upper triangular part of A is not referenced. On exit, the
C           lower triangular part of the array A is overwritten by the
C           lower triangular part of the updated matrix.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO.EQ.'U' .OR. UPLO.EQ.'u').AND.
     $         .NOT.(UPLO.EQ.'L' .OR. UPLO.EQ.'l')      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PPF/DSYR  ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Set the start point in X if the increment is not unity.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
C
      IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
C
C        Form  A  when A is stored in upper triangle.
C
         IF( INCX.EQ.1 )THEN
            DO 20, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 10, I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE
            JX = KX
            DO 40, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = KX
                  DO 30, I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
   30             CONTINUE
               END IF
               JX = JX + INCX
   40       CONTINUE
         END IF
      ELSE
C
C        Form  A  when A is stored in lower triangle.
C
         IF( INCX.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 50, I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            JX = KX
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = JX
                  DO 70, I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PPF (DSYR  ).
C
      END
** END OF F06PPFTEXT
*UPTODATE F06PFFTEXT
      SUBROUTINE F06PFF( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
C     ..
C
C  Purpose
C  =======
C
C  DTRMV  performs one of the matrix-vector operations
C
C     x := A*x,   or   x := A'*x,
C
C  where x is n element vector and A is an n by n unit, or non-unit,
C  upper or lower triangular matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   x := A*x.
C
C              TRANS = 'T' or 't'   x := A'*x.
C
C              TRANS = 'C' or 'c'   x := A'*x.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit
C           triangular as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular matrix and the strictly lower triangular part of
C           A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular matrix and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C           A are not referenced either, but are assumed to be unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x. On exit, X is overwritten with the
C           tranformed vector x.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO .EQ.'U' .OR. UPLO .EQ.'u').AND.
     $         .NOT.(UPLO .EQ.'L' .OR. UPLO .EQ.'l')      )THEN
         INFO = 1
      ELSE IF( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 2
      ELSE IF( .NOT.(DIAG .EQ.'U' .OR. DIAG .EQ.'u').AND.
     $         .NOT.(DIAG .EQ.'N' .OR. DIAG .EQ.'n')      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PFF/DTRMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  x := A*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
C
C        Form  x := A'*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 90, I = J - 1, 1, -1
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 120, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 110, I = J - 1, 1, -1
                     IX   = IX   - INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  110             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = 1, N
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 130, I = J + 1, N
                     TEMP = TEMP + A( I, J )*X( I )
  130             CONTINUE
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               JX = KX
               DO 160, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = J + 1, N
                     IX   = IX   + INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  150             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  160          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06PFF (DTRMV ).
C
      END
** END OF F06PFFTEXT
*UPTODATE F06PJFTEXT
      SUBROUTINE F06PJF( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
C     ..
C
C  Purpose
C  =======
C
C  DTRSV  solves one of the systems of equations
C
C     A*x = b,   or   A'*x = b,
C
C  where b and x are n element vectors and A is an n by n unit, or
C  non-unit, upper or lower triangular matrix.
C
C  No test for singularity or near-singularity is included in this
C  routine. Such tests must be performed before calling this routine.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the equations to be solved as
C           follows:
C
C              TRANS = 'N' or 'n'   A*x = b.
C
C              TRANS = 'T' or 't'   A'*x = b.
C
C              TRANS = 'C' or 'c'   A'*x = b.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit
C           triangular as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular matrix and the strictly lower triangular part of
C           A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular matrix and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C           A are not referenced either, but are assumed to be unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element right-hand side vector b. On exit, X is overwritten
C           with the solution vector x.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO .EQ.'U' .OR. UPLO .EQ.'u').AND.
     $         .NOT.(UPLO .EQ.'L' .OR. UPLO .EQ.'l')      )THEN
         INFO = 1
      ELSE IF( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 2
      ELSE IF( .NOT.(DIAG .EQ.'U' .OR. DIAG .EQ.'u').AND.
     $         .NOT.(DIAG .EQ.'N' .OR. DIAG .EQ.'n')      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PJF/DTRSV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  x := inv( A )*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO 10, I = J - 1, 1, -1
                        X( I ) = X( I ) - TEMP*A( I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 40, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO 30, I = J - 1, 1, -1
                        IX      = IX      - INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO 50, I = J + 1, N
                        X( I ) = X( I ) - TEMP*A( I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO 70, I = J + 1, N
                        IX      = IX      + INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
C
C        Form  x := inv( A' )*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = 1, N
                  TEMP = X( J )
                  DO 90, I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( I )
   90             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  DO 110, I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   + INCX
  110             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = N, 1, -1
                  TEMP = X( J )
                  DO 130, I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( I )
  130             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  DO 150, I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   - INCX
  150             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  160          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06PJF (DTRSV ).
C
      END
** END OF F06PJFTEXT
*UPTODATE F06QFFTEXT
      SUBROUTINE F06QFF( MATRIX, M, N, A, LDA, B, LDB )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      CHARACTER*1        MATRIX
      INTEGER            M, N, LDA, LDB
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
C     ..
C
C  F06QFF  copies  the  m by n  matrix  A  into  the  m by n  matrix  B.
C
C  If   MATRIX = 'G' or 'g'   then  A  and  B  are  regarded as  general
C                             matrices,
C  if   MATRIX = 'U' or 'u'   then  A  and  B  are  regarded  as   upper
C                             triangular,  and only  elements  for which
C                             i.le.j  are referenced,
C  if   MATRIX = 'L' or 'l'   then  A  and  B  are  regarded  as   lower
C                             triangular,  and only  elements  for which
C                             i.ge.j  are referenced.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 21-November-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
         DO 20 J = 1, N
            DO 10 I = 1, M
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ) )THEN
         DO 40 J = 1, N
            DO 30 I = 1, MIN( M, J )
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ) )THEN
         DO 60 J = 1, MIN( M, N )
            DO 50 I = J, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QFF. ( SMCOPY )
C
      END
** END OF F06QFFTEXT
*UPTODATE F06QGFTEXT
      DOUBLE PRECISION FUNCTION F06QGF( NORM, MATRIX, M, N, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER                           LDA, M, N
      CHARACTER*1                       MATRIX, NORM
C     .. Array Arguments ..
      DOUBLE PRECISION                  A( LDA, * )
C     ..
C
C  Purpose
C  =======
C
C  F06QGF  returns the value of the one norm,  or the Frobenius norm, or
C  the element of  largest absolute value of a real matrix A.  A  may be
C  rectangular,  or square,  or triangular,  or symmetric.
C
C  Description
C  ===========
C
C  F06QGF returns the value
C
C     F06QGF = ( max( abs( a( i, j ) ) ) , NORM = 'M' or 'm'
C              (
C              ( norm1( A )  ,             NORM = '1', 'O' or 'o'
C              (
C              ( normF( A ) ,              NORM = 'F', 'f', 'E' or 'e'
C
C  where norm1 denotes the one norm of a matrix (maximum column sum) and
C  normF denotes the  Frobenius norm of a matrix  (square root of sum of
C  squares).  Note that  max( abs( a( i, j ) ) )  is not a  matrix norm.
C
C  The type of matrix for which  F06QGF is returned is determined by the
C  parameter MATRIX.
C
C  If   MATRIX = 'G' or 'g'   then  A  is regarded as  a general matrix,
C  If   MATRIX = 'U' or 'u'   then  A  is regarded as  upper triangular,
C  If   MATRIX = 'L' or 'l'   then  A  is regarded as  lower triangular,
C  If   MATRIX = 'S' or 's'   then  A  is regarded as symmetric and only
C             or 'H' or 'h'   the  upper triangular part of the array  A
C                             is referenced,
C  If   MATRIX = 'Y' or 'y'   then  A  is regarded as symmetric and only
C             or 'E' or 'e'   the  lower triangular part of the array  A
C                             is referenced.
C
C  Parameters
C  ==========
C
C  NORM  -  CHARACTER*1.
C
C           On entry,  NORM specifies the value to be returned in F06QGF
C           as described above.
C
C           Unchanged on exit.
C
C  MATRIX - CHARACTER*1.
C
C           On entry,  MATRIX  specifies the type of matrix and,  in the
C           case of a  symmetric matrix,  the part of the array in which
C           the matrix is stored as described above.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry,  M  specifies the number of rows of the matrix  A.
C           M  must be at least  zero and when the  matrix is  symmetric
C           then  M must be equal to  N. When  M = 0  then F06QGF is set
C           to zero and an immediate return is effected.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N  must be at least zero. When  N = 0  then F06QGF is set to
C           zero and an immediate return is effected.
C
C           Unchanged on exit.
C
C  A      - REAL array of DIMENSION ( LDA, n ).
C
C           Before entry,  A  must contain the  m by n  matrix for which
C           F06QGF is required.
C
C           If  MATRIX = 'U' or 'u' or 'S' or 's' or 'H' or 'h' then the
C           strictly lower triangular part of A is not referenced.
C
C           If  MATRIX = 'L' or 'l' or 'Y' or 'y' or 'E' or 'e' then the
C           strictly upper triangular part of A is not referenced.
C
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in  the  calling  (sub)  program.  LDA  must be at least  M.
C
C           Unchanged on exit.
C
C  Further comments
C  ================
C
C  If A is part of a matrix B partitioned as
C
C     B = ( B1  B2 ) ,
C         ( B3  A  )
C
C  where  B1 is an l by k matrix  ( l.ge.0, k.ge.0 ),  then this routine
C  may be called with the parameter  A as  b( l + 1, k + 1 ) and  LDA as
C  the first dimension of  B  as declared in the calling  (sub) program.
C
C  This routine  can be inefficient on  paged machines when the one norm
C  is required, the matrix is symmetric and N is large.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION                  ONE, ZERO
      PARAMETER                         ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION                  SCALE, SUM, VALUE
      INTEGER                           I, J
C     .. External Functions ..
      DOUBLE PRECISION                  F06BMF
      EXTERNAL                          F06BMF
C     .. External Subroutines ..
      EXTERNAL                          F06FJF
C     .. Intrinsic Functions ..
      INTRINSIC                         ABS, MAX, MIN
C     ..
C     .. Executable Statements ..
      IF( MIN( M, N ).EQ.0 )THEN
         VALUE = ZERO
      ELSE IF( ( NORM.EQ.'M' ).OR.( NORM.EQ.'m' ) )THEN
C
C        Find  max( abs( a( i, j ) ) ).
C
         VALUE = ZERO
         IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
            DO 20 J = 1, N
               DO 10 I = 1, M
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10          CONTINUE
   20       CONTINUE
         ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ).OR.
     $            ( MATRIX.EQ.'S' ).OR.( MATRIX.EQ.'s' ).OR.
     $            ( MATRIX.EQ.'H' ).OR.( MATRIX.EQ.'h' ) )THEN
            DO 40 J = 1, N
               DO 30 I = 1, MIN( M, J )
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   30          CONTINUE
   40       CONTINUE
         ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ).OR.
     $            ( MATRIX.EQ.'Y' ).OR.( MATRIX.EQ.'y' ).OR.
     $            ( MATRIX.EQ.'E' ).OR.( MATRIX.EQ.'e' ) )THEN
            DO 60 J = 1, MIN( M, N )
               DO 50 I = J, M
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   50          CONTINUE
   60       CONTINUE
         END IF
      ELSE IF( ( NORM.EQ.'1' ).OR.( NORM.EQ.'O' ).OR.
     $         ( NORM.EQ.'o' ) )THEN
C
C        Find  norm1( A ).
C
         VALUE = ZERO
         IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
            DO 80 J = 1, N
               SUM = ZERO
               DO 70 I = 1, M
                  SUM = SUM + ABS( A( I, J ) )
   70          CONTINUE
               VALUE = MAX( VALUE, SUM )
   80       CONTINUE
         ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ) )THEN
            DO 100 J = 1, N
               SUM = ZERO
               DO 90 I = 1, MIN( M, J )
                  SUM = SUM + ABS( A( I, J ) )
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ) )THEN
            DO 120 J = 1, MIN( M, N )
               SUM = ZERO
               DO 110 I = J, M
                  SUM = SUM + ABS( A( I, J ) )
  110          CONTINUE
               VALUE = MAX( VALUE, SUM )
  120       CONTINUE
         ELSE IF( ( MATRIX.EQ.'S' ).OR.( MATRIX.EQ.'s' ).OR.
     $            ( MATRIX.EQ.'H' ).OR.( MATRIX.EQ.'h' ) )THEN
            DO 150 J = 1, N
               SUM = ZERO
               DO 130 I = 1, J
                  SUM = SUM + ABS( A( I, J ) )
  130          CONTINUE
               DO 140 I = J + 1, N
                  SUM = SUM + ABS( A( J, I ) )
  140          CONTINUE
               VALUE = MAX( VALUE, SUM )
  150       CONTINUE
         ELSE IF( ( MATRIX.EQ.'Y' ).OR.( MATRIX.EQ.'y' ).OR.
     $            ( MATRIX.EQ.'E' ).OR.( MATRIX.EQ.'e' ) )THEN
            DO 180 J = 1, N
               SUM = ZERO
               DO 160 I = 1, J - 1
                  SUM = SUM + ABS( A( J, I ) )
  160          CONTINUE
               DO 170 I = J, N
                  SUM = SUM + ABS( A( I, J ) )
  170          CONTINUE
               VALUE = MAX( VALUE, SUM )
  180       CONTINUE
         END IF
      ELSE IF( ( NORM.EQ.'F' ).OR.( NORM.EQ.'f' ).OR.( NORM.EQ.'E' ).OR.
     $         ( NORM.EQ.'e' ) )THEN
C
C        Find  normF( A ).
C
         SCALE = ZERO
         SUM = ONE
         IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
            DO 190 J = 1, N
               CALL F06FJF( M, A( 1, J ), 1, SCALE, SUM )
  190       CONTINUE
         ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ) )THEN
            DO 200 J = 1, N
               CALL F06FJF( MIN( M, J ), A( 1, J ), 1, SCALE, SUM )
  200       CONTINUE
         ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ) )THEN
            DO 210 J = 1, MIN( M, N )
               CALL F06FJF( M - J + 1, A( J, J ), 1, SCALE, SUM )
  210       CONTINUE
         ELSE IF( ( MATRIX.EQ.'S' ).OR.( MATRIX.EQ.'s' ).OR.
     $            ( MATRIX.EQ.'H' ).OR.( MATRIX.EQ.'h' ).OR.
     $            ( MATRIX.EQ.'Y' ).OR.( MATRIX.EQ.'y' ).OR.
     $            ( MATRIX.EQ.'E' ).OR.( MATRIX.EQ.'e' ) )THEN
            IF( ( MATRIX.EQ.'S' ).OR.( MATRIX.EQ.'s' ).OR.
     $          ( MATRIX.EQ.'H' ).OR.( MATRIX.EQ.'h' ) )THEN
               DO 220 J = 2, N
                  CALL F06FJF( J - 1, A( 1, J ), 1, SCALE, SUM )
  220          CONTINUE
            ELSE
               DO 230 J = 1, N - 1
                  CALL F06FJF( N - J, A( J + 1, J ), 1, SCALE, SUM )
  230          CONTINUE
            END IF
            SUM = 2*SUM
            CALL F06FJF( N, A( 1, 1 ), LDA + 1, SCALE, SUM )
         END IF
         VALUE = F06BMF( SCALE, SUM )
      END IF
C
      F06QGF = VALUE
      RETURN
C
C     End of F06QGF. ( SMNRM )
C
      END
** END OF F06QGFTEXT
*UPTODATE F06QHFTEXT
      SUBROUTINE F06QHF( MATRIX, M, N, CONST, DIAG, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      CHARACTER*1        MATRIX
      DOUBLE PRECISION   CONST, DIAG
      INTEGER            LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
C     ..
C
C  F06QHF forms the m by n matrix A given by
C
C     a( i, j ) = (  diag  i.eq.j,
C                 (
C                 ( const  i.ne.j.
C
C  If   MATRIX = 'G' or 'g'   then  A  is regarded  as a general matrix,
C  if   MATRIX = 'U' or 'u'   then  A  is regarded  as upper triangular,
C                             and only  elements  for which  i.le.j  are
C                             referenced,
C  if   MATRIX = 'L' or 'l'   then  A  is regarded  as lower triangular,
C                             and only  elements  for which  i.ge.j  are
C                             referenced.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 21-November-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
         DO 20 J = 1, N
            DO 10 I = 1, M
               A( I, J ) = CONST
   10       CONTINUE
   20    CONTINUE
         IF( CONST.NE.DIAG )THEN
            DO 30 I = 1, MIN( M, N )
               A( I, I ) = DIAG
   30       CONTINUE
         END IF
      ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ) )THEN
         DO 50 J = 1, N
            DO 40 I = 1, MIN( M, J )
               A( I, J ) = CONST
   40       CONTINUE
   50    CONTINUE
         IF( CONST.NE.DIAG )THEN
            DO 60 I = 1, MIN( M, N )
               A( I, I ) = DIAG
   60       CONTINUE
         END IF
      ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ) )THEN
         DO 80 J = 1, MIN( M, N )
            DO 70 I = J, M
               A( I, J ) = CONST
   70       CONTINUE
   80    CONTINUE
         IF( CONST.NE.DIAG )THEN
            DO 90 I = 1, MIN( M, N )
               A( I, I ) = DIAG
   90       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06QHF. ( SMLOAD )
C
      END
** END OF F06QHFTEXT
*UPTODATE F06QJFTEXT
      SUBROUTINE F06QJF( SIDE, TRANS, N, PERM, K, B, LDB )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K, LDB, N
      CHARACTER*1        SIDE, TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * )
      INTEGER            PERM( * )
C     ..
C
C  Purpose
C  =======
C
C  F06QJF performs one of the transformations
C
C     B := P'*B   or   B := P*B,   where B is an m by k matrix,
C
C  or
C
C     B := B*P'   or   B := B*P,   where B is a k by m matrix,
C
C  P being an m by m permutation matrix of the form
C
C     P = P( 1, index( 1 ) )*P( 2, index( 2 ) )*...*P( n, index( n ) ),
C
C  where  P( i, index( i ) ) is the permutation matrix that interchanges
C  items i and index( i ). That is P( i, index( i ) ) is the unit matrix
C  with rows and columns  i and index( i )  interchanged.  Of course, if
C  index( i ) = i  then  P( i, index( i ) ) = I.
C
C  This routine  is intended for use in  conjunction with  Nag auxiliary
C  routines that  perform  interchange  operations,  such  as  pivoting.
C
C  Parameters
C  ==========
C
C  SIDE   - CHARACTER*1.
C  TRANS
C           On entry,  SIDE  ( Left-hand side, or Right-hand side )  and
C           TRANS  ( Transpose, or No transpose )  specify the operation
C           to be performed as follows.
C
C           SIDE = 'L' or 'l'   and   TRANS = 'T' or 't'
C
C              Perform the operation   B := P'*B.
C
C           SIDE = 'L' or 'l'   and   TRANS = 'N' or 'n'
C
C              Perform the operation   B := P*B.
C
C           SIDE = 'R' or 'r'   and   TRANS = 'T' or 't'
C
C              Perform the operation   B := B*P'.
C
C           SIDE = 'R' or 'r'   and   TRANS = 'N' or 'n'
C
C              Perform the operation   B := B*P.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry,  N must specify the value of n. N must be at least
C           zero.  When  N = 0  then an  immediate  return is  effected.
C
C           Unchanged on exit.
C
C  PERM   - INTEGER array of DIMENSION at least ( n ).
C
C           Before  entry,  PERM  must  contain the  n  indices  for the
C           permutation matrices. index( i ) must satisfy
C
C              1 .le. index( i ) .le. m.
C
C           It is usual for index( i ) to be at least i, but this is not
C           necessary for this routine.
C
C           Unchanged on exit.
C
C  K      - INTEGER.
C
C           On entry with  SIDE = 'L' or 'l',  K must specify the number
C           of columns of B and on entry with  SIDE = 'R' or 'r', K must
C           specify the  number of rows of B.  K must be at least  zero.
C           When  K = 0  then an immediate return is effected.
C
C           Unchanged on exit.
C
C  B      - REAL  array of  DIMENSION  ( LDB, ncolb ),  where  ncolb = k
C           when  SIDE = 'L' or 'l'  and  ncolb = m  when  SIDE = 'R' or
C           'r'.
C
C           Before entry  with  SIDE = 'L' or 'l',  the  leading  m by K
C           part  of  the  array   B  must  contain  the  matrix  to  be
C           transformed  and  before entry with  SIDE = 'R' or 'r',  the
C           leading  K by m part of the array  B must contain the matrix
C           to  be  transformed.  On  exit,  B  is  overwritten  by  the
C           transformed matrix.
C
C  LDB    - INTEGER.
C
C           On entry,  LDB  must specify  the  leading dimension  of the
C           array  B  as declared  in the  calling  (sub) program.  When
C           SIDE = 'L' or 'l'   then  LDB  must  be  at  least  m,  when
C           SIDE = 'R' or 'r'   then  LDB  must  be  at  least  k.
C           Unchanged on exit.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, J, L
      LOGICAL            LEFT, NULL, RIGHT, TRNSP
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( MIN( N, K ).EQ.0 )
     $   RETURN
      LEFT = ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' )
      RIGHT = ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' )
      NULL = ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' )
      TRNSP = ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' )
      IF( LEFT )THEN
         IF( TRNSP )THEN
            DO 20 I = 1, N
               IF( PERM( I ).NE.I )THEN
                  L = PERM( I )
                  DO 10 J = 1, K
                     TEMP = B( I, J )
                     B( I, J ) = B( L, J )
                     B( L, J ) = TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE IF( NULL )THEN
            DO 40 I = N, 1, -1
               IF( PERM( I ).NE.I )THEN
                  L = PERM( I )
                  DO 30 J = 1, K
                     TEMP = B( L, J )
                     B( L, J ) = B( I, J )
                     B( I, J ) = TEMP
   30             CONTINUE
               END IF
   40       CONTINUE
         END IF
      ELSE IF( RIGHT )THEN
         IF( TRNSP )THEN
            DO 60 J = N, 1, -1
               IF( PERM( J ).NE.J )THEN
                  L = PERM( J )
                  DO 50 I = 1, K
                     TEMP = B( I, J )
                     B( I, J ) = B( I, L )
                     B( I, L ) = TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE IF( NULL )THEN
            DO 80 J = 1, N
               IF( PERM( J ).NE.J )THEN
                  L = PERM( J )
                  DO 70 I = 1, K
                     TEMP = B( I, L )
                     B( I, L ) = B( I, J )
                     B( I, J ) = TEMP
   70             CONTINUE
               END IF
   80       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06QJF. ( SGEAP )
C
      END
** END OF F06QJFTEXT
*UPTODATE F06QKFTEXT
      SUBROUTINE F06QKF( SIDE, TRANS, N, PERM, K, B, LDB )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K, LDB, N
      CHARACTER*1        SIDE, TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION   PERM( * ), B( LDB, * )
C     ..
C
C  Purpose
C  =======
C
C  F06QKF performs one of the transformations
C
C     B := P'*B   or   B := P*B,   where B is an m by k matrix,
C
C  or
C
C     B := B*P'   or   B := B*P,   where B is a k by m matrix,
C
C  P being an m by m permutation matrix of the form
C
C     P = P( 1, index( 1 ) )*P( 2, index( 2 ) )*...*P( n, index( n ) ),
C
C  where  P( i, index( i ) ) is the permutation matrix that interchanges
C  items i and index( i ). That is P( i, index( i ) ) is the unit matrix
C  with rows and columns  i and  index( i )  interchanged. Of course, if
C  index( i ) = i  then  P( i, index( i ) ) = I.
C
C  This  routine is intended  for use in conjunction with  Nag auxiliary
C  routines  that  perform  interchange  operations,  such  as  sorting.
C
C  Parameters
C  ==========
C
C  SIDE   - CHARACTER*1.
C  TRANS
C           On entry,  SIDE  ( Left-hand side, or Right-hand side )  and
C           TRANS  ( Transpose, or No transpose )  specify the operation
C           to be performed as follows.
C
C           SIDE = 'L' or 'l'   and   TRANS = 'T' or 't'
C
C              Perform the operation   B := P'*B.
C
C           SIDE = 'L' or 'l'   and   TRANS = 'N' or 'n'
C
C              Perform the operation   B := P*B.
C
C           SIDE = 'R' or 'r'   and   TRANS = 'T' or 't'
C
C              Perform the operation   B := B*P'.
C
C           SIDE = 'R' or 'r'   and   TRANS = 'N' or 'n'
C
C              Perform the operation   B := B*P.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N must specify the value of n.  N must be at least
C           zero.  When  N = 0  then an  immediate  return  is effected.
C
C           Unchanged on exit.
C
C  PERM   - REAL             array of DIMENSION at least ( n ).
C
C           Before  entry,  PERM  must  contain  the  n indices  for the
C           permutation matrices. index( i ) must satisfy
C
C              1 .le. index( i ) .le. m.
C
C           It is usual for index( i ) to be at least i, but this is not
C           necessary for this routine. It is assumed that the statement
C           INDEX = PERM( I )  returns the correct integer in  INDEX, so
C           that,  if necessary,  PERM( I )  should contain a real value
C           slightly larger than  INDEX.
C
C           Unchanged on exit.
C
C  K      - INTEGER.
C
C           On entry with  SIDE = 'L' or 'l',  K must specify the number
C           of columns of B and on entry with  SIDE = 'R' or 'r', K must
C           specify the number of rows of  B.  K must be at least  zero.
C           When  K = 0  then an immediate return is effected.
C
C           Unchanged on exit.
C
C  B      - REAL  array  of  DIMENSION ( LDB, ncolb ),  where  ncolb = k
C           when  SIDE = 'L' or 'l'  and  ncolb = m  when  SIDE = 'R' or
C           'r'.
C
C           Before entry  with  SIDE = 'L' or 'l',  the  leading  m by K
C           part  of  the  array   B  must  contain  the  matrix  to  be
C           transformed  and before  entry with  SIDE = 'R' or 'r',  the
C           leading  K by m part of the array  B must contain the matrix
C           to  be  transformed.  On exit,   B  is  overwritten  by  the
C           transformed matrix.
C
C  LDB    - INTEGER.
C
C           On entry,  LDB  must specify  the  leading dimension  of the
C           array  B  as declared  in the  calling  (sub) program.  When
C           SIDE = 'L' or 'l'   then  LDB  must  be  at  least  m,  when
C           SIDE = 'R' or 'r'   then  LDB  must  be  at  least  k.
C           Unchanged on exit.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 11-August-1987.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      LOGICAL            LEFT, NULL, RIGHT, TRNSP
      INTEGER            I, J, L
      DOUBLE PRECISION   TEMP
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( MIN( N, K ).EQ.0 )
     $   RETURN
      LEFT = ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' )
      RIGHT = ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' )
      NULL = ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' )
      TRNSP = ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' )
      IF( LEFT )THEN
         IF( TRNSP )THEN
            DO 20 I = 1, N
               L = PERM( I )
               IF( L.NE.I )THEN
                  DO 10 J = 1, K
                     TEMP = B( I, J )
                     B( I, J ) = B( L, J )
                     B( L, J ) = TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE IF( NULL )THEN
            DO 40 I = N, 1, -1
               L = PERM( I )
               IF( L.NE.I )THEN
                  DO 30 J = 1, K
                     TEMP = B( L, J )
                     B( L, J ) = B( I, J )
                     B( I, J ) = TEMP
   30             CONTINUE
               END IF
   40       CONTINUE
         END IF
      ELSE IF( RIGHT )THEN
         IF( TRNSP )THEN
            DO 60 J = N, 1, -1
               L = PERM( J )
               IF( L.NE.J )THEN
                  DO 50 I = 1, K
                     TEMP = B( I, J )
                     B( I, J ) = B( I, L )
                     B( I, L ) = TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE IF( NULL )THEN
            DO 80 J = 1, N
               L = PERM( J )
               IF( L.NE.J )THEN
                  DO 70 I = 1, K
                     TEMP = B( I, L )
                     B( I, L ) = B( I, J )
                     B( I, J ) = TEMP
   70             CONTINUE
               END IF
   80       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06QKF. ( SGEAPR )
C
      END
** END OF F06QKFTEXT
*UPTODATE F06QNFTEXT
      SUBROUTINE F06QNF( SIDE, N, K1, K2, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), S( * )
C     ..
C
C  F06QNF applies a  sequence  of  pairwise interchanges to either  the
C  left,  or the right,  of the  n by n  upper triangular matrix  U,  to
C  transform U to an  upper Hessenberg matrix. The interchanges are
C  applied in planes k1 up to k2.
C
C  The upper Hessenberg matrix, H, is formed as
C
C     H = P*U,    when   SIDE = 'L' or 'l',  (  Left-hand side )
C
C  where P is a permutation matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 )
C
C  and is formed as
C
C     H = U*P',   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is a permutation matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k ) being a pairwise interchange for the  ( k, k + 1 ) plane.
C  The  two by two
C  interchange part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = ( 0  1 ).
C              ( 1  0 )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of  H.
C
C  The  sub-diagonal elements of  H, h( k + 1, k ),  are returned in the
C  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 16-May-1988.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, TEMP
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the permutations to columns n back to k1.
C
         DO 20 J = N, K1, -1
            IF( J.GE.K2 )THEN
               AIJ = A( K2, J )
            ELSE
C
C              Form  the  additional sub-diagonal element  h( j + 1, j )
C              and store it in s( j ).
C
               AIJ    = ZERO
               S( J ) = A( J, J )
            END IF
            DO 10 I = MIN( K2, J ) - 1, K1, -1
               TEMP          = A( I, J )
               A( I + 1, J ) = TEMP
               AIJ           = AIJ
   10       CONTINUE
            A( K1, J ) = AIJ
   20    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Apply  the  plane interchanges to  columns  k1  up to
C        ( k2 - 1 ) and  form   the   additional  sub-diagonal
C        elements,   storing  h( j + 1, j ) in s( j ).
C
         DO 40 J = K1, K2 - 1
            DO 30 I = 1, J
               TEMP = A( I, J + 1 )
               A( I, J + 1 ) = A( I, J )
               A( I, J )     = TEMP
   30       CONTINUE
            S( J )            = A( J + 1, J + 1 )
            A( J + 1, J + 1 ) = ZERO
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QNF. ( SUTSRH )
C
      END
** END OF F06QNFTEXT
*UPTODATE F06QRFTEXT
      SUBROUTINE F06QRF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QRF restores an upper Hessenberg matrix H to upper triangular form
C  by  applying a sequence of  plane rotations  from either the left, or
C  the right.  The matrix  H  is assumed to have  non-zero  sub-diagonal
C  elements  in  positions  h( k + 1, k ),  k = k1, k1 + 1, ..., k2 - 1,
C  only  and  h( k + 1, k )  must  be  supplied  in  s( k ).
C
C  H is restored to the upper triangular matrix R either as
C
C     R = P*H,   when   SIDE = 'L' or 'l'  (  Left-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  or as
C
C     R = H*P',  when   SIDE = 'R' or 'r'  ( Right-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C  in both cases  P( k )  being a  plane rotation  for the  ( k, k + 1 )
C  plane.  The cosine and sine that define P( k ) are returned in c( k )
C  and  s( k )  respectively.  The two by two  rotation part of  P( k ),
C  Q( k ), is of the form
C
C     Q( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The upper triangular part of the matrix  H  must be supplied in the n
C  by n  leading upper triangular part of  A, and this is overwritten by
C  the upper triangular matrix R.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, STEMP, SUBH, TEMP
      INTEGER            I, J
C     .. External Subroutines ..
      EXTERNAL           F06BAF
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Restore   H  to  upper  triangular  form  by  annihilating  the
C        sub-diagonal elements of H.  The jth rotation is chosen so that
C
C           ( h( j, j ) ) := (  c  s )*( h( j, j )     ).
C           (     0     )    ( -s  c ) ( h( j + 1, j ) )
C
C        Apply the rotations in columns k1 up to n.
C
         DO 20 J = K1, N
            AIJ = A( K1, J )
            DO 10 I = K1, MIN( J, K2 ) - 1
               TEMP = A( I + 1, J )
               A( I, J ) = S( I )*TEMP + C( I )*AIJ
               AIJ = C( I )*TEMP - S( I )*AIJ
   10       CONTINUE
            IF( J.LT.K2 )THEN
C
C              Set up the rotation.
C
               SUBH = S( J )
               CALL F06BAF( AIJ, SUBH, C( J ), S( J ) )
               A( J, J ) = AIJ
            ELSE
               A( K2, J ) = AIJ
            END IF
   20    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Restore   H  to  upper  triangular  form  by  annihilating  the
C        sub-diagonal elements of H.  The jth rotation is chosen so that
C
C           ( h( j + 1, j + 1 ) ) := (  c  s )*( h( j + 1, j + 1 ) ),
C           (         0         )    ( -s  c ) ( h( j + 1, j )     )
C
C        which can be expressed as
C
C           ( 0  h( j + 1, j + 1 ) ) :=
C
C               ( h( j + 1, j )  h( j + 1, j + 1 ) )*(  c  s ).
C                                                    ( -s  c )
C
C        Thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
C        rotation matrix look like
C
C           Q( j ) = (  c( j )  s( j ) ).
C                    ( -s( j )  c( j ) )
C
         DO 40 J = K2 - 1, K1, -1
            SUBH = S( J )
            CALL F06BAF( A( J + 1, J + 1 ), SUBH, CTEMP, STEMP )
            STEMP = -STEMP
            S( J ) = STEMP
            C( J ) = CTEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               DO 30 I = J, 1, -1
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
            END IF
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QRF. ( SUHQR )
C
      END
** END OF F06QRFTEXT
*UPTODATE F06QSFTEXT
      SUBROUTINE F06QSF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QSF restores an upper spiked matrix  H to upper triangular form by
C  applying a sequence of plane rotations, in planes  k1 up to k2,  from
C  either the left, or the right.
C
C  The matrix  H is assumed to have non-zero elements only in the spiked
C  positions, h( k2, k ) for a row spike and h( k + 1, k1 ) for a column
C  spike, k = k1, k1 + 1, ..., k2 - 1, and these must be supplied in the
C  elements s( k ).
C
C  When  SIDE = 'L' or 'l'  (  Left-hand side )
C
C     H  is  assumed  to have a  row spike  and is restored to the upper
C     triangular matrix  R as
C
C        R = P*H,
C
C     where P is an orthogonal matrix of the form
C
C        P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C     P( k )  being a  plane rotation  matrix for the  ( k, k2 )  plane.
C
C  When  SIDE = 'R' or 'r'  ( Right-hand side )
C
C     H  is assumed to have a  column spike and is restored to the upper
C     triangular matrix R as
C
C        R = H*P',
C
C     where P is an orthogonal matrix of the form
C
C        P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C     P( k ) being a plane rotation matrix for the  ( k1, k + 1 ) plane.
C
C  The  two by two  rotation  part of  P( k ),  Q( k ),  is of  the form
C
C     Q( k ) = (  c( k )  s( k ) )
C              ( -s( k )  c( k ) )
C
C  and  c( k ) and s( k ) are returned in the kth elements of the arrays
C  C and S respectively.
C
C  The upper triangular part of the matrix  H must be supplied in the  n
C  by n  leading upper triangular part of  A, and this is overwritten by
C  the upper triangular matrix R.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, SPIKE, STEMP, TEMP
      INTEGER            I, J
C     .. External Subroutines ..
      EXTERNAL           F06BAF
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Restore H to upper triangular form by annihilating the elements
C        in  the  spike  of  H.  The  jth rotation  is  chosen  so  that
C
C        ( h( j, j ) ) := (  c  s )*( h( j , j ) ).
C        (     0     )    ( -s  c ) ( h( k2, j ) )
C
C        Apply the rotations in columns k1 up to ( k2 - 1 ).
C
         DO 20 J = K1, K2 - 1
            SPIKE = S( J )
            DO 10 I = K1, J - 1
               AIJ = A( I, J )
               A( I, J ) = S( I )*SPIKE + C( I )*AIJ
               SPIKE = C( I )*SPIKE - S( I )*AIJ
   10       CONTINUE
C
C           Set up the rotation.
C
            CALL F06BAF( A( J, J ), SPIKE, C( J ), S( J ) )
   20    CONTINUE
C
C        Apply the rotations to columns k2 up to n.
C
         DO 40 J = K2, N
            TEMP = A( K2, J )
            DO 30 I = K1, K2 - 1
               AIJ = A( I, J )
               A( I, J ) = S( I )*TEMP + C( I )*AIJ
               TEMP = C( I )*TEMP - S( I )*AIJ
   30       CONTINUE
            A( K2, J ) = TEMP
   40    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Restore H to upper triangular form by annihilating the spike of
C        H. The jth rotation is chosen so that
C
C           ( h( j, j ) ) := (  c  s )*( h( j, j )  ),
C           (     0     )    ( -s  c ) ( h( j, k1 ) )
C
C        which can be expressed as
C
C           ( 0  h( j, j ) ) := ( h( j, k1 )  h( j, j ) )*(  c  s ).
C                                                         ( -s  c )
C
C        Thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
C        rotation matrix look like
C
C           Q( j ) = (  c( j )  s( j ) ).
C                    ( -s( j )  c( j ) )
C
         DO 70 J = K2, K1 + 1, -1
            CALL F06BAF( A( J, J ), S( J - 1 ), CTEMP, STEMP )
            STEMP = -STEMP
            S( J - 1 ) = STEMP
            C( J - 1 ) = CTEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               DO 50 I = J - 1, K1 + 1, -1
                  SPIKE = S( I - 1 )
                  S( I - 1 ) = STEMP*A( I, J ) + CTEMP*SPIKE
                  A( I, J ) = CTEMP*A( I, J ) - STEMP*SPIKE
   50          CONTINUE
               DO 60 I = K1, 1, -1
                  TEMP = A( I, K1 )
                  A( I, K1 ) = STEMP*A( I, J ) + CTEMP*TEMP
                  A( I, J ) = CTEMP*A( I, J ) - STEMP*TEMP
   60          CONTINUE
            END IF
   70    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QSF. ( SUSQR )
C
      END
** END OF F06QSFTEXT
*UPTODATE F06QTFTEXT
      SUBROUTINE F06QTF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QTF performs the transformation
C
C     R := P*U*Q'  when  SIDE = 'L' or 'l'  (  Left-hand side )
C
C     R := Q*U*P'  when  SIDE = 'R' or 'r'  ( Right-hand side ),
C
C  where  U and R  are  n by n  upper  triangular  matrices,   P  is  an
C  orthogonal matrix,  consisting of a given sequence of plane rotations
C  to be  applied  in  planes  k1 to k2,  and  Q  is  a  unitary  matrix
C  consisting of a sequence of plane rotations, applied in planes  k1 to
C  k2,  chosen to make  R  upper triangular.
C
C  When  SIDE = 'L' or 'l'  then  P  is  given  as a  sequence of  plane
C  rotation matrices
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  where  P( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
C  In this case the matrix Q is given as
C
C     Q = Q( k2 - 1 )*...*Q( k1 + 1 )*Q( k1 ),
C
C  where  Q( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
C
C  When  SIDE = 'R' or 'r'  then  P  is  given  as a  sequence of  plane
C  rotation matrices
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C  where  P( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
C  In this case the matrix Q is given as
C
C     Q = Q( k1 )*Q( k1 + 1 )*...*Q( k2 - 1 ),
C
C  where  Q( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
C
C  The  upper  triangular  matrix  U  must  be  supplied  in the  n by n
C  leading upper triangular part of  A,  and this  is overwritten by the
C  upper triangular matrix  R.  The cosine  and  sine  that  define  the
C  plane rotation matrix  P( k )  must be supplied in  c( k ) and s( k )
C  respectively,  and  the two by two rotation part of  P( k ),  T( k ),
C  is assumed to be of the form
C
C     T( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The cosine  and  sine that define  Q( k )  are overwritten on  c( k )
C  and  s( k )  respectively and the two by two rotation part of  Q( k )
C  will have the form of  T( k )  above.
C
C  If  n or k1  are less  than  unity, or  k1  is not  less than  k2, or
C  k2  is greater than  n  then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 26-November-1987.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, FILL, STEMP, TEMP
      INTEGER            I, I1, J
C     .. External Subroutines ..
      EXTERNAL           F06BAF
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the left-hand transformations,  column by column,  to the
C        triangular part of  U,  but not to  anywhere  that would  cause
C        fill.
C
         DO 20 J = K1 + 1, N
C
C           Apply  P( k1 ) ... P( j - 1 )  to column j.
C
            AIJ = A( K1, J )
            DO 10 I = K1, MIN( J - 1, K2 - 1 )
               A( I, J ) = S( I )*A( I + 1, J ) + C( I )*AIJ
               AIJ = C( I )*A( I + 1, J ) - S( I )*AIJ
   10       CONTINUE
            A( I, J ) = AIJ
   20    CONTINUE
C
C           Now apply each  left-hand tranformation  to form the fill-in
C           elements and apply a  right-hand transformation to eliminate
C           the fill-in element.
C
         DO 40 J = K1, K2 - 1
C
C           Apply  P( j )  to the jth diagonal element  and the  fill-in
C           position.
C
            FILL = -S( J )*A( J, J )
            A( J, J ) = C( J )*A( J, J )
C
C           Now  set up  the rotation  Q( j )  to eliminate the  fill-in
C           element,  and  apply  Q( j )  to  the  jth  and  ( j + 1 )th
C           columns.
C
            CALL F06BAF( A( J + 1, J + 1 ), FILL, CTEMP, STEMP )
            C( J ) = CTEMP
            S( J ) = -STEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               STEMP = -STEMP
               DO 30 I = 1, J
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
            END IF
   40    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        We intermingle the  left and right hand transformations so that
C        at the kth step we form
C
C           A := Q( k )*A*P( k )'.
C
C        First  apply  the  transformations  in  columns  k2 back to k1.
C
         DO 60 J = K2 - 1, K1, -1
C
C           First apply  P( j ).
C
            IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
               CTEMP = C( J )
               STEMP = S( J )
               DO 50 I = 1, J
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   50          CONTINUE
C
C              Next form the fill-in element  a( j + 1, j )  by applying
C              P( j ).
C
               FILL = S( J )*A( J + 1, J + 1 )
               A( J + 1, J + 1 ) = C( J )*A( J + 1, J + 1 )
C
C              Now set up the rotation  Q( j )  to eliminate the fill-in
C              element.
C
               CALL F06BAF( A( J, J ), FILL, C( J ), S( J ) )
            END IF
   60    CONTINUE
C
C        Finally  apply  Q( k2 - 1 ) ... Q( k1 )  to columns  n  back to
C        ( k1 + 1 ).
C
         DO 80 J = N, K1 + 1, -1
            I1 = MIN( K2, J )
            AIJ = A( I1, J )
            DO 70 I = I1 - 1, K1, -1
               TEMP = A( I, J )
               A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
               AIJ = S( I )*AIJ + C( I )*TEMP
   70       CONTINUE
            A( K1, J ) = AIJ
   80    CONTINUE
      END IF
      RETURN
C
C     End of F06QTF. ( SUTSQR )
C
      END
** END OF F06QTFTEXT
*UPTODATE F06QVFTEXT
      SUBROUTINE F06QVF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QVF applies a  given sequence  of  plane rotations  to either  the
C  left,  or the right,  of the  n by n  upper triangular matrix  U,  to
C  transform U to an  upper Hessenberg matrix. The rotations are applied
C  in planes k1 up to k2.
C
C  The upper Hessenberg matrix, H, is formed as
C
C     H = P*U,    when   SIDE = 'L' or 'l',  (  Left-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 )
C
C  and is formed as
C
C     H = U*P',   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k ) being a plane rotation matrix for the  ( k, k + 1 ) plane. The
C  cosine and sine that define P( k ), k = k1, k1 + 1, ..., k2 - 1, must
C  be  supplied  in  c( k )  and  s( k )  respectively.  The  two by two
C  rotation part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of  H.
C
C  The  sub-diagonal elements of  H, h( k + 1, k ),  are returned in the
C  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, STEMP, TEMP
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the plane rotations to columns n back to k1.
C
         DO 20 J = N, K1, -1
            IF( J.GE.K2 )THEN
               AIJ = A( K2, J )
            ELSE
C
C              Form  the  additional sub-diagonal element  h( j + 1, j )
C              and store it in s( j ).
C
               AIJ = C( J )*A( J, J )
               S( J ) = -S( J )*A( J, J )
            END IF
            DO 10 I = MIN( K2, J ) - 1, K1, -1
               TEMP = A( I, J )
               A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
               AIJ = S( I )*AIJ + C( I )*TEMP
   10       CONTINUE
            A( K1, J ) = AIJ
   20    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Apply  the  plane rotations  to  columns  k1  up to  ( k2 - 1 )
C        and  form   the   additional  sub-diagonal  elements,   storing
C        h( j + 1, j ) in s( j ).
C
         DO 40 J = K1, K2 - 1
            IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
               STEMP = S( J )
               CTEMP = C( J )
               DO 30 I = 1, J
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
               S( J ) = STEMP*A( J + 1, J + 1 )
               A( J + 1, J + 1 ) = CTEMP*A( J + 1, J + 1 )
            END IF
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QVF. ( SUTSRH )
C
      END
** END OF F06QVFTEXT
*UPTODATE F06QWFTEXT
      SUBROUTINE F06QWF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QWF applies a  given sequence  of  plane rotations  to either  the
C  left,  or the right,  of the  n by n  upper triangular  matrix  U  to
C  transform  U  to an upper spiked matrix. The rotations are applied in
C  planes k1 up to k2.
C
C  The upper spiked matrix, H, is formed as
C
C     H = P*U,   when   SIDE = 'L' or 'l',  ( Left-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C  P( k ) being a plane rotation matrix for the ( k, k2 ) plane, and is
C  formed as
C
C     H = U*P',   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k )  being a  plane rotation matrix for the  ( k1, k + 1 )  plane.
C
C  The cosine and sine that define  P( k ), k = k1, k1 + 1, ..., k2 - 1,
C  must be  supplied  in  c( k ) and s( k ) respectively. The two by two
C  rotation part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of H.
C
C  When  SIDE = 'L' or 'l'  then a  row spike  is  generated  in  H  and
C  when  SIDE = 'R' or 'r'  then a  column spike is generated. For a row
C  spike the elements  h( k2, k )  and for a  column spike  the elements
C  h( k + 1, k1 ), k = k1, k1 + 1, ..., k2 - 1, are returned in  s( k ).
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, SPIKE, STEMP, TEMP
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the plane rotations to columns n back to k2.
C
         DO 20 J = N, K2, -1
            TEMP = A( K2, J )
            DO 10 I = K2 - 1, K1, -1
               AIJ = A( I, J )
               A( I, J ) = S( I )*TEMP + C( I )*AIJ
               TEMP = C( I )*TEMP - S( I )*AIJ
   10       CONTINUE
            A( K2, J ) = TEMP
   20    CONTINUE
C
C        Form  the spike  and apply the rotations in columns  ( k2 - 1 )
C        back to k1.
C
         DO 40 J = K2 - 1, K1, -1
            SPIKE = -S( J )*A( J, J )
            A( J, J ) = C( J )*A( J, J )
            DO 30 I = J - 1, K1, -1
               AIJ = A( I, J )
               A( I, J ) = S( I )*SPIKE + C( I )*AIJ
               SPIKE = C( I )*SPIKE - S( I )*AIJ
   30       CONTINUE
            S( J ) = SPIKE
   40    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Apply the  plane rotations to columns  ( k1 + 1 ) up to k2  and
C        form the spike.
C
         DO 70 J = K1 + 1, K2
            CTEMP = C( J - 1 )
            STEMP = S( J - 1 )
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               DO 50 I = 1, K1
                  TEMP = A( I, K1 )
                  A( I, K1 ) = STEMP*A( I, J ) + CTEMP*TEMP
                  A( I, J ) = CTEMP*A( I, J ) - STEMP*TEMP
   50          CONTINUE
               DO 60 I = K1 + 1, J - 1
                  SPIKE = S( I - 1 )
                  S( I - 1 ) = STEMP*A( I, J ) + CTEMP*SPIKE
                  A( I, J ) = CTEMP*A( I, J ) - STEMP*SPIKE
   60          CONTINUE
               S( J - 1 ) = STEMP*A( J, J )
               A( J, J ) = CTEMP*A( J, J )
            END IF
   70    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QWF. ( SUTSRS )
C
      END
** END OF F06QWFTEXT
*UPTODATE F06QXFTEXT
      SUBROUTINE F06QXF( SIDE, PIVOT, DIRECT, M, N, K1, K2, C, S, A,
     $                   LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, M, N
      CHARACTER*1        DIRECT, PIVOT, SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QXF  performs the transformation
C
C     A := P*A,   when   SIDE = 'L' or 'l'  (  Left-hand side )
C
C     A := A*P',  when   SIDE = 'R' or 'r'  ( Right-hand side )
C
C  where A is an m by n matrix and P is an orthogonal matrix, consisting
C  of a  sequence  of  plane  rotations,  applied  in  planes  k1 to k2,
C  determined by the parameters PIVOT and DIRECT as follows:
C
C     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
C        plane.
C
C     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
C        plane.
C
C     When  PIVOT  = 'T' or 't'  ( Top pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
C
C        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
C        plane.
C
C     When  PIVOT  = 'T' or 't'  ( Top pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
C        plane.
C
C     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
C
C        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
C        plane.
C
C     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
C        plane.
C
C  c( k ) and s( k )  must contain the  cosine and sine  that define the
C  matrix  P( k ).  The  two by two  plane rotation  part of the  matrix
C  P( k ), R( k ), is assumed to be of the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  If m, n or k1 are less than unity,  or k2 is not greater than k1,  or
C  SIDE = 'L' or 'l'  and  k2  is greater than  m, or  SIDE = 'R' or 'r'
C  and  k2  is greater than  n,  then an  immediate return  is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 20-November-1986.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, STEMP, TEMP
      INTEGER            I, J
      LOGICAL            LEFT, RIGHT
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      LEFT = ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' )
      RIGHT = ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' )
      IF( ( MIN( M, N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $    ( ( LEFT ).AND.( K2.GT.M ) ).OR.
     $    ( ( RIGHT ).AND.( K2.GT.N ) ) )RETURN
      IF( LEFT )THEN
         IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 20 J = 1, N
                  AIJ = A( K1, J )
                  DO 10 I = K1, K2 - 1
                     TEMP = A( I + 1, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     AIJ = C( I )*TEMP - S( I )*AIJ
   10             CONTINUE
                  A( K2, J ) = AIJ
   20          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 40 J = 1, N
                  AIJ = A( K2, J )
                  DO 30 I = K2 - 1, K1, -1
                     TEMP = A( I, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     AIJ = S( I )*AIJ + C( I )*TEMP
   30             CONTINUE
                  A( K1, J ) = AIJ
   40          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'T' ).OR.( PIVOT.EQ.'t' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 60 J = 1, N
                  TEMP = A( K1, J )
                  DO 50 I = K1, K2 - 1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     TEMP = S( I )*AIJ + C( I )*TEMP
   50             CONTINUE
                  A( K1, J ) = TEMP
   60          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 80 J = 1, N
                  TEMP = A( K1, J )
                  DO 70 I = K2 - 1, K1, -1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     TEMP = S( I )*AIJ + C( I )*TEMP
   70             CONTINUE
                  A( K1, J ) = TEMP
   80          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'B' ).OR.( PIVOT.EQ.'b' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 100 J = 1, N
                  TEMP = A( K2, J )
                  DO 90 I = K1, K2 - 1
                     AIJ = A( I, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     TEMP = C( I )*TEMP - S( I )*AIJ
   90             CONTINUE
                  A( K2, J ) = TEMP
  100          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 120 J = 1, N
                  TEMP = A( K2, J )
                  DO 110 I = K2 - 1, K1, -1
                     AIJ = A( I, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     TEMP = C( I )*TEMP - S( I )*AIJ
  110             CONTINUE
                  A( K2, J ) = TEMP
  120          CONTINUE
            END IF
         END IF
      ELSE IF( RIGHT )THEN
         IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 140 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 130 I = 1, M
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 160 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 150 I = M, 1, -1
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'T' ).OR.( PIVOT.EQ.'t' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 180 J = K1 + 1, K2
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 200 J = K2, K1 + 1, -1
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 190 I = M, 1, -1
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'B' ).OR.( PIVOT.EQ.'b' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 220 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, K2 ) + CTEMP*TEMP
                        A( I, K2 ) = CTEMP*A( I, K2 ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 240 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 230 I = M, 1, -1
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, K2 ) + CTEMP*TEMP
                        A( I, K2 ) = CTEMP*A( I, K2 ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06QXF. ( SGESRC )
C
      END
** END OF F06QXFTEXT
      SUBROUTINE F06QZF( HESS, N, K1, K2, C, S, A, LDA )
*     .. Scalar Arguments ..
      CHARACTER*1        HESS
      INTEGER            K1, K2, LDA, N
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
*     ..
*
*  F06QZF  either applies a  given sequence  of  plane rotations  to the
*  right of the n by n reverse lower triangular matrix T, to transform T
*  to a  reverse lower Hessenberg matrix  H, or restores a reverse lower
*  Hessenberg matrix H to reverse lower triangular form T, by applying a
*  sequence of plane rotations from the right.
*
*  The rotations are applied  in planes k1 up to k2.
*
*  When   HESS = 'C' or 'c',   ( Create ),  then   the   reverse   lower
*  Hessenberg matrix, H, is formed as
*
*     H = T*P',
*
*  where P is an orthogonal matrix of the form
*
*     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
*
*  P( k ) being a plane rotation matrix for the  ( k, k + 1 ) plane. The
*  cosine and sine that define P( k ), k = k1, k1 + 1, ..., k2 - 1, must
*  be  supplied  in  c( k )  and  s( k )  respectively.  The  two by two
*  rotation part of P( k ), R( k ), is assumed to have the form
*
*     R( k ) = (  c( k )  s( k ) ).
*              ( -s( k )  c( k ) )
*
*  The matrix  T must be supplied in the n by n reverse lower triangular
*  part  of the array  A,  and this is overwritten by the  reverse lower
*  triangular part of  H.
*
*  The super-diagonal elements of  H, h( n - k, k ), are returned in the
*  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
*
*  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
*  greater than n then an immediate return is effected.
*
*  When   HESS = 'R' or 'r',   ( Remove ),  then   the   reverse   lower
*  Hessenberg matrix  H  is  assumed  to  have  non-zero  super-diagonal
*  elements  in  positions  h( n - k, k ),  k = k1, k1 + 1, ..., k2 - 1,
*  only and  h( n - k, k ) must be supplied in  s( k ). H is restored to
*  the reverse lower triangular matrix T as
*
*     T = H*P',
*
*  where P is an orthogonal matrix of the form
*
*     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
*
*  P( k ) being a plane rotation for the  ( k, k + 1 ) plane. The cosine
*  and  sine  that  define  P( k )  are  returned  in  c( k ) and s( k )
*  respectively.  The  two by two  rotation part of  P( k ),  R( k ), is
*  of the form
*
*     R( k ) = (  c( k )  s( k ) ).
*              ( -s( k )  c( k ) )
*
*  The reverse lower triangular part of the matrix H must be supplied in
*  the  n by n  reverse  lower  triangular  part  of  A,   and  this  is
*  overwritten by the reverse triangular matrix T.
*
*  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
*  greater than n then an immediate return is effected.
*
*  When   n = 7, k1 = 2 and k2 = 5   then  T  and  H  are  of  the  form
*
*     T = ( 0  0  0  0  0  0  X ),   H = ( 0  0  0  0  0  0  X ).
*         ( 0  0  0  0  0  X  X )        ( 0  0  0  0  X  X  X )
*         ( 0  0  0  0  X  X  X )        ( 0  0  0  X  X  X  X )
*         ( 0  0  0  X  X  X  X )        ( 0  0  X  X  X  X  X )
*         ( 0  0  X  X  X  X  X )        ( 0  X  X  X  X  X  X )
*         ( 0  X  X  X  X  X  X )        ( 0  X  X  X  X  X  X )
*         ( X  X  X  X  X  X  X )        ( X  X  X  X  X  X  X )
*
*
*  This routine  is  principally intended  for use  with the  non-linear
*  optimization routines such as E04UCF, in order to help vectorization.
*  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
*
*  -- Written on 10-May-1988.
*     Sven Hammarling, Nag Central Office.
*
*
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     .. External Subroutines ..
      EXTERNAL           F06BAF
*     .. Local Scalars ..
      DOUBLE PRECISION   CTEMP, STEMP, SUPH, TEMP
      INTEGER            I, J
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.( K2.GT.N ) )
     $   RETURN
      IF( ( HESS.EQ.'C' ).OR.( HESS.EQ.'c' ) )THEN
*
*        Apply  the  plane rotations  to  columns  k1  up to  ( k2 - 1 )
*        and  form   the  additional  super-diagonal  elements,  storing
*        h( n - j, j ) in s( j ).
*
         DO 20, J = K1, K2 - 1
            IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
               STEMP             = S( J )
               CTEMP             = C( J )
               S( J )            = STEMP*A( N - J, J + 1 )
               A( N - J, J + 1 ) = CTEMP*A( N - J, J + 1 )
               DO 10, I = N - J + 1, N
                  TEMP          = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J )     = STEMP*TEMP + CTEMP*A( I, J )
   10          CONTINUE
            END IF
   20    CONTINUE
      ELSE IF( ( HESS.EQ.'R' ).OR.( HESS.EQ.'r' ) )THEN
*
*        Restore  H to reverse lower triangular form by annihilating the
*        super-diagonal elements of  H.  The  jth rotation  is chosen so
*        that
*
*          ( h( n - j, n - j ) ) := (  c  s )*( h( n - j, n - j     ) ),
*          (         0         )    ( -s  c ) ( h( n - j, n - j - 1 ) )
*
*        which can be expressed as
*
*           ( 0  h( n - j, n - j ) ) :=
*
*               ( h( n - j, n - j - 1 )  h( n - j, n - j ) )*(  c  s ).
*                                                            ( -s  c )
*
*        Thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
*        rotation matrix look like
*
*           R( j ) = (  c( j )  s( j ) ).
*                    ( -s( j )  c( j ) )
*
         DO 40, J = K2 - 1, K1, -1
            SUPH   =  S( J )
            CALL F06BAF( A( N - J, J + 1 ), SUPH, CTEMP, STEMP )
            STEMP  = -STEMP
            S( J ) =  STEMP
            C( J ) =  CTEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               DO 30, I = N - J + 1, N
                  TEMP          = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J )     = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
            END IF
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of F06QZF.
*
      END
** END OF F06QZFTEXT
