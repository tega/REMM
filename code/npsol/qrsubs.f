*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  QRSUBS FORTRAN
*
*     DGEQR    DGEQRP   DGEAP    DGEAPQ
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DGEQR ( M, N, A, LDA, ZETA, INFORM )
      INTEGER            M, N, LDA, INFORM
      DOUBLE PRECISION   A( LDA, * ), ZETA( * )
C
C  1. Purpose
C     =======
C
C  DGEQR  reduces the  m by n, m.ge.n, matrix A to upper triangular form
C  by means of orthogonal transformations.
C
C  2. Description
C     ===========
C
C  The m by n matrix A is factorized as
C
C     A = Q*( R )   when   m.gt.n,
C           ( 0 )
C
C     A = Q*R       when   m = n,
C
C  where  Q  is an  m by m  orthogonal matrix and  R  is an n by n upper
C  triangular matrix.
C
C  The  factorization  is  obtained  by  Householder's  method. The  kth
C  transformation matrix, Q( k ), which is used to introduce zeros  into
C  the kth column of A is given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',   u( k ) = ( zeta( k ) ),
C                                             (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C  zeta( k )  and  z( k ) are chosen to annhilate the elements below the
C  triangular part of  A.
C
C  The vector  u( k )  is returned in the kth element of ZETA and in the
C  kth column of A, such that zeta( k ) is in ZETA( k ) and the elements
C  of z( k ) are in a( k + 1, k ), ..., a( m, k ). The elements of R are
C  returned in the upper triangular part of  A.
C
C  Q is given by
C
C     Q = ( Q( p )*Q( p - 1 )*...*Q( 1 ) )',
C
C  where p = min( n, m - 1 ).
C
C  3. Parameters
C     ==========
C
C  M      - INTEGER.
C
C           On entry, M must specify the number of rows of  A. M must be
C           at least  n.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N must specify the number of columns of  A. N must
C           be  at  least zero. When  N = 0  then an immediate return is
C           effected.
C
C           Unchanged on exit.
C
C  A      - 'real' array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On exit, the  N by N upper triangular part of A will contain
C           the  upper  triangular  matrix  R  and the  M by N  strictly
C           lower triangular part of  A  will  contain  details  of  the
C           factorization as described above.
C
C  LDA    - INTEGER.
C
C           On entry, LDA  must  specify  the  leading dimension of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least  m.
C
C           Unchanged on exit.
C
C  ZETA   - 'real' array of DIMENSION at least ( n ).
C
C           On  exit, ZETA( k )  contains the scalar  zeta( k )  for the
C           kth  transformation.  If  T( k ) = I  then   ZETA( k ) = 0.0
C           otherwise  ZETA( k )  contains  zeta( k ) as described above
C           and is always in the range ( 1.0, sqrt( 2.0 ) ).
C
C  INFORM - INTEGER.
C
C           On successful  exit  INFORM  will be zero, otherwise  INFORM
C           will  be set to unity indicating that an input parameter has
C           been  incorrectly  set. See  the  next section  for  further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  INFORM = 1
C
C     One or more of the following conditions holds:
C
C        M   .lt. N
C        N   .lt. 0
C        LDA .lt. M
C
C  5. Further information
C     ===================
C
C  Following the use of this routine the operations
C
C     B := Q'*B   and   B := Q*B,
C
C  where  B  is an  m by k  matrix, can  be  performed  by calls to  the
C  auxiliary  linear  algebra routine  DGEAPQ. The  operation  B := Q'*B
C  can be obtained by the call:
C
C     INFORM = 0
C     CALL DGEAPQ( 'Transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, INFORM )
C
C  and  B := Q*B  can be obtained by the call:
C
C     INFORM = 0
C     CALL DGEAPQ( 'No transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, INFORM )
C
C  In  both  cases  WORK  must be a  k  element array  that  is used  as
C  workspace. If  B  is a one-dimensional array (single column) then the
C  parameter  LDB  can be replaced by  M. See routine DGEAPQ for further
C  details.
C
C  Operations involving the matrix  R  are performed by  the
C  Level 2 BLAS  routines  DTRMV  and DTRSV . Note that no test for near
C  singularity of R is incorporated in this routine or in routine  DTRSV
C  and  so it is  strongly recommended that the auxiliary linear algebra
C  routine  DUTCO  be called, prior to solving equations involving R, in
C  order  to determine whether  or not  R  is nearly singular. If  R  is
C  nearly  singular  then  the  auxiliary linear algebra  routine  DUTSV
C  can  be used to  determine  the  singular value decomposition  of  R.
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 13-December-1984.
C     Sven Hammarling, Nag Central Office.
C
      EXTERNAL           DGEMV , DGER  , DGRFG
      INTRINSIC          MIN
      INTEGER            K     , LA
      DOUBLE PRECISION   TEMP
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )

*     Check the input parameters.

      IF( N.EQ.0 )THEN
         INFORM = 0
         RETURN
      END IF
      IF( ( M.LT.N ).OR.( N.LT.0 ).OR.( LDA.LT.M ) )THEN
         INFORM = 1
         RETURN
      END IF

*     Perform the factorization.

      LA = LDA
      DO 20, K = 1, MIN( M - 1, N )

*        Use a Householder reflection to zero the kth column of A.
*        First set up the reflection.

         CALL DGRFG ( M - K, A( K, K ), A( K + 1, K ), 1, ZERO,
     $                ZETA( K ) )
         IF( ( ZETA( K ).GT.ZERO ).AND.( K.LT.N ) )THEN
            IF( ( K + 1 ).EQ.N )
     $         LA = M - K + 1
            TEMP      = A( K, K )
            A( K, K ) = ZETA( K )

*           We now perform the operation  A := Q( k )*A.

*           Let B denote the bottom ( m - k + 1 ) by ( n - k ) part
*           of A.

*           First form  work = B'*u. ( work is stored in the elements
*           ZETA( k + 1 ), ..., ZETA( n ). )

            CALL DGEMV ( 'Transpose', M - K + 1, N - K,
     $                   ONE, A( K, K + 1 ), LA, A( K, K ), 1,
     $                   ZERO, ZETA( K + 1 ), 1 )

*           Now form  B := B - u*work'.

            CALL DGER  ( M - K + 1, N - K, -ONE, A( K, K ), 1,
     $                   ZETA( K + 1 ), 1, A( K, K + 1 ), LA )

*           Restore beta.

            A( K, K ) = TEMP
         END IF
   20 CONTINUE

*     Store the final zeta when m.eq.n.

      IF( M.EQ.N )
     $   ZETA( N ) = ZERO

      INFORM = 0
      RETURN

*     End of DGEQR .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DGEQRP( PIVOT, M, N, A, LDA, ZETA, PERM, WORK, INFORM )
      CHARACTER*1        PIVOT
      INTEGER            M, N, LDA, INFORM
      INTEGER            PERM( * )
      DOUBLE PRECISION   A( LDA, * ), ZETA( * ), WORK( * )

      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/

C  1. Purpose
C     =======
C
C  DGEQRP reduces the  m by n matrix A to upper triangular form by means
C  of orthogonal transformations and column permutations.
C
C  2. Description
C     ===========
C
C  The m by n matrix A is factorized as
C
C     A = Q*( R )*P'      when   m.gt.n,
C           ( 0 )
C
C     A = Q*R*P'          when   m = n,
C
C     A = Q*( R  X )*P'   when   m.lt.n,
C
C  where  Q  is  an  m by m  orthogonal matrix, R  is a  min( m, n )  by
C  min( m, n )  upper triangular matrix and  P is an  n by n permutation
C  matrix.
C
C  The  factorization  is  obtained  by  Householder's  method. The  kth
C  transformation matrix, Q( k ),  which is used to introduce zeros into
C  the kth column of A is given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',   u( k ) = ( zeta( k ) ),
C                                             (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C  zeta( k )  and  z( k ) are chosen to annhilate the elements below the
C  triangular part of  A.
C
C  The vector  u( k )  is returned in the kth element of ZETA and in the
C  kth column of A, such that zeta( k ) is in ZETA( k ) and the elements
C  of z( k ) are in a( k + 1, k ), ..., a( m, k ). The elements of R are
C  returned in the upper triangular part of A.
C
C  Q is given by
C
C     Q = ( Q( p )*Q( p - 1 )*...*Q( 1 ) )',
C
C  where p = min( m - 1, n ).
C
C  Two options are available for the column permutations. In either case
C  the column for which the  sub-diagonal elements are to be annihilated
C  at the  kth step is chosen from the remaining ( n - k + 1 )  columns.
C  The  particular column chosen as the pivot column is either that  for
C  which  the  unreduced  part  ( elements k onwards )  has the  largest
C  Euclidean  length, or  is that for  which the ratio of the  Euclidean
C  length  of the  unreduced part  to the  Euclidean length of the whole
C  column is a maximum.
C
C  3. Parameters
C     ==========
C
C  PIVOT  - CHARACTER*1.
C
C           On  entry, PIVOT  specifies  the  pivoting  strategy  to  be
C           performed as follows.
C
C           PIVOT = 'C' or 'c'
C
C              Column  interchanges  are  to be  incorporated  into  the
C              factorization, such that the  column whose unreduced part
C              has  maximum  Euclidean  length  is chosen  as the  pivot
C              column at each step.
C
C           PIVOT = 'S' or 's'
C
C              Scaled  column interchanges  are to be  incorporated into
C              the  factorization, such  that the  column for which  the
C              ratio  of the  Euclidean  length of the unreduced part of
C              the column to the original Euclidean length of the column
C              is a maximum is chosen as the  pivot column at each step.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry, M  must specify the number of rows of A. M must be
C           at  least  zero. When  M = 0  then  an  immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N  must specify the number of columns of A. N must
C           be  at least zero. When  N = 0  then an immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  A      - 'real' array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On  exit, the  min( M, N ) by min( M, N )  upper  triangular
C           part of A will contain the upper triangular matrix R and the
C           M by min( M, N )  strictly lower triangular part of  A  will
C           contain details  of the  factorization  as  described above.
C           When m.lt.n then the remaining M by ( N - M ) part of A will
C           contain the matrix X.
C
C  LDA    - INTEGER.
C
C           On  entry, LDA  must  specify  the leading dimension of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least  m.
C
C           Unchanged on exit.
C
C  ZETA   - 'real' array of DIMENSION at least ( n ).
C
C           On exit, ZETA( k )  contains the scalar  zeta  for  the  kth
C           transformation. If T( k ) = I then ZETA( k) = 0.0, otherwise
C           ZETA( k )  contains the scalar  zeta( k ) as described above
C           and  is  always  in  the  range  ( 1.0, sqrt( 2.0 ) ).  When
C           n .gt. m  the  elements  ZETA( m + 1 ),  ZETA( m + 2 ), ...,
C           ZETA( n )  are used as internal workspace.
C
C  PERM   - INTEGER array of DIMENSION at least min( m, n ).
C
C           On exit, PERM  contains details of the permutation matrix P,
C           such  that  PERM( k ) = k  if no  column interchange occured
C           at  the  kth  step  and  PERM( k ) = j, ( k .lt. j .le. n ),
C           if  columns  k and j  were  interchanged at  the  kth  step.
C           Note  that, although  there are  min( m - 1, n )  orthogonal
C           transformations, there are min( m, n ) permutations.
C
C  WORK   - 'real' array of DIMENSION at least ( 2*n ).
C
C           Used as internal workspace.
C
C           On exit, WORK( j ), j = 1, 2, ..., n, contains the Euclidean
C           length  of the  jth  column  of the  permuted  matrix  A*P'.
C
C  INFORM - INTEGER.
C
C           On  successful exit, INFORM  will be zero, otherwise  INFORM
C           will  be set to unity indicating that an input parameter has
C           been  incorrectly supplied. See the next section for further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  INFORM = 1
C
C     One or more of the following conditions holds:
C
C        PIVOT .ne. 'C' or 'c' or 'S' or 's'
C        M     .lt. 0
C        N     .lt. 0
C        LDA   .lt. M
C
C  5. Further information
C     ===================
C
C  Following the use of this routine the operations
C
C     B := Q'*B   and   B := Q*B,
C
C  where  B  is an  m by k  matrix, can  be  performed  by calls to  the
C  auxiliary  linear algebra  routine  DGEAPQ. The  operation  B := Q'*B
C  can be obtained by the call:
C
C     INFORM = 0
C     CALL DGEAPQ( 'Transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, INFORM )
C
C  and  B := Q*B  can be obtained by the call:
C
C     INFORM = 0
C     CALL DGEAPQ( 'No transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, INFORM )
C
C  In  both  cases  WORK  must be  a  k  element array  that is used  as
C  workspace. If B is a one-dimensional array ( single column ) then the
C  parameter  LDB  can be replaced by  M. See routine DGEAPQ for further
C  details.
C
C  Also following the use of this routine the operations
C
C     B := P'*B   and   B := P*B,
C
C  where B is an n by k matrix, and the operations
C
C     B := B*P    and   B := B*P',
C
C  where  B is a k by n  matrix, can  be performed by calls to the basic
C  linear  algebra  routine  DGEAP .  The  operation  B := P'*B  can  be
C  obtained by the call:
C
C     CALL DGEAP ( 'Left', 'Transpose', N, MIN( M, N ), PERM,
C    $             K, B, LDB )
C
C  the operation  B := P*B  can be obtained by the call:
C
C     CALL DGEAP ( 'Left', 'No transpose', N, MIN( M, N ), PERM,
C    $             K, B, LDB )
C
C  If  B is a one-dimensional array ( single column ) then the parameter
C  LDB  can be replaced by  N  in the above two calls.
C  The operation  B := B*P  can be obtained by the call:
C
C     CALL DGEAP ( 'Right', 'No transpose', K, MIN( M, N ), PERM,
C    $             M, B, LDB )
C
C  and  B := B*P'  can be obtained by the call:
C
C     CALL DGEAP ( 'Right', 'Transpose', K, MIN( M, N ), PERM,
C    $             M, B, LDB )
C
C  If  B is a one-dimensional array ( single column ) then the parameter
C  LDB  can be replaced by  K  in the above two calls.
C  See routine DGEAP for further details.
C
C  Operations involving  the matrix  R  are performed by  the
C  Level 2 BLAS  routines  DTRSV  and DTRMV.  Note that no test for near
C  singularity of  R is incorporated in this routine or in routine DTRSV
C  and  so it is  strongly recommended that the auxiliary linear algebra
C  routine  DUTCO  be called, prior to solving equations involving R, in
C  order  to determine whether  or not  R  is nearly singular. If  R  is
C  nearly  singular then  the  auxiliary  linear algebra  routine  DUTSV
C  can  be  used  to  determine  the  singular value decomposition of R.
C  Operations  involving  the  matrix  X  can also be  performed  by the
C  Level 2  BLAS  routines.  Matrices  of  the  form   ( R  X )  can  be
C  factorized as
C
C     ( R  X ) = ( T  0 )*S',
C
C  where  T is upper triangular and S is orthogonal, using the auxiliary
C  linear algebra routine  DUTRQ .
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 13-December-1984.
C     Sven Hammarling, Nag Central Office.
C
      EXTERNAL           MCHPAR, DGEMV , DGER  , DGRFG , DNRM2 , DSWAP
      INTRINSIC          ABS   , MAX   , MIN   , SQRT
      INTEGER            J     , JMAX  , K     , LA
      DOUBLE PRECISION   EPS   , MAXNRM, NORM  , DNRM2 , TEMP  , TOL
      DOUBLE PRECISION   LAMDA
      PARAMETER        ( LAMDA = 1.0D-2 )
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )

*     Check the input parameters.

      IF( MIN( M, N ).EQ.0 )THEN
         INFORM = 0
         RETURN
      END IF
      IF( ( ( PIVOT.NE.'C' ).AND.( PIVOT.NE.'c' ).AND.
     $      ( PIVOT.NE.'S' ).AND.( PIVOT.NE.'s' )      ).OR.
     $    ( M.LT.0 ).OR.( N.LT.0 ).OR.( LDA.LT.M )           )THEN
         INFORM = 1
         RETURN
      END IF

*     Compute eps and the initial column norms.

      CALL MCHPAR()
      EPS = WMACH( 3 )
      DO 10, J = 1, N
         WORK( J )     = DNRM2 ( M, A( 1, J ), 1 )
         WORK( J + N ) = WORK( J )
   10 CONTINUE

*     Perform the factorization. TOL is the tolerance for DGRFG .

      LA = LDA
      DO 50, K = 1, MIN( M, N )

*        Find the pivot column.

         MAXNRM = ZERO
         JMAX   = K
         DO 20, J = K, N
            IF( ( PIVOT.EQ.'C' ).OR.( PIVOT.EQ.'c' ) )THEN
               IF( WORK( J + N  ).GT.MAXNRM )THEN
                  MAXNRM = WORK( J + N )
                  JMAX   = J
               END IF
            ELSE IF( WORK( J ).GT.ZERO )THEN
               IF( ( WORK( J + N )/WORK( J ) ).GT.MAXNRM )THEN
                  MAXNRM = WORK( J + N )/WORK( J )
                  JMAX   = J
               END IF
            END IF
   20    CONTINUE
         PERM( K ) = JMAX
         IF( JMAX.GT.K )THEN
            CALL DSWAP ( M, A( 1, K ), 1, A( 1, JMAX ), 1 )
            TEMP             = WORK( K )
            WORK( K )        = WORK( JMAX )
            WORK( JMAX )     = TEMP
            WORK( JMAX + N ) = WORK( K + N )
            PERM( K )        = JMAX
         END IF
         TOL = EPS*WORK( K )
         IF( K.LT.M )THEN

*           Use a Householder reflection to zero the kth column of A.
*           First set up the reflection.

            CALL DGRFG ( M - K, A( K, K ), A( K + 1, K ), 1, TOL,
     $                   ZETA( K ) )
            IF( K.LT.N )THEN
               IF( ZETA( K ).GT.ZERO )THEN
                  IF( ( K + 1 ).EQ.N )
     $               LA = M - K + 1
                  TEMP      = A( K, K )
                  A( K, K ) = ZETA( K )

*                 We now perform the operation  A := Q( k )*A.

*                 Let B denote the bottom ( m - k + 1 ) by ( n - k )
*                 part of A.

*                 First form  work = B'*u. ( work is stored in the
*                 elements ZETA( k + 1 ), ..., ZETA( n ). )

                  CALL DGEMV ( 'Transpose', M - K + 1, N - K,
     $                         ONE, A( K, K + 1 ), LA, A( K, K ), 1,
     $                         ZERO, ZETA( K + 1 ), 1 )

*                 Now form  B := B - u*work'.

                  CALL DGER  ( M - K + 1, N - K, -ONE, A( K, K ), 1,
     $                         ZETA( K + 1 ), 1, A( K, K + 1 ), LA )

*                 Restore beta.

                  A( K, K ) = TEMP
               END IF

*              Update the unreduced column norms. Use the Linpack
*              criterion for when to recompute the norms, except that
*              we retain the original column lengths throughout and use
*              a smaller lamda.

               DO 40, J = K + 1, N
                  IF( WORK( J + N ).GT.ZERO )THEN
                     TEMP = ABS( A( K, J ) )/WORK( J + N )
                     TEMP = MAX( ( ONE + TEMP )*( ONE - TEMP ), ZERO )
                     NORM = TEMP
                     TEMP = ONE +
     $                      LAMDA*TEMP*( WORK( J + N )/WORK( J ) )**2
                     IF( TEMP.GT.ONE )THEN
                        WORK( J + N ) = WORK( J + N )*SQRT( NORM )
                     ELSE
                        WORK( J + N ) = DNRM2 ( M - K,
     $                                          A( K + 1, J ), 1 )
                     END IF
                  END IF
   40          CONTINUE
            END IF
         END IF
   50 CONTINUE

*     Store the final zeta when m.le.n.

      IF( M.LE.N )
     $   ZETA( M ) = ZERO

      INFORM = 0
      RETURN

*     End of DGEQRP.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DGEAP ( SIDE, TRANS, M, N, PERM, K, B, LDB )
*     .. Scalar Arguments ..
      INTEGER            K, LDB, M, N
      CHARACTER*1        SIDE, TRANS
*     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * )
      INTEGER            PERM( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEAP  performs one of the transformations
*
*     B := P'*B   or   B := P*B,   where B is an m by k matrix,
*
*  or
*
*     B := B*P'   or   B := B*P,   where B is a k by m matrix,
*
*  P being an m by m permutation matrix of the form
*
*     P = P( 1, index( 1 ) )*P( 2, index( 2 ) )*...*P( n, index( n ) ),
*
*  where  P( i, index( i ) ) is the permutation matrix that interchanges
*  items i and index( i ). That is P( i, index( i ) ) is the unit matrix
*  with rows and columns  i and index( i )  interchanged.  Of course, if
*  index( i ) = i  then  P( i, index( i ) ) = I.
*
*  This routine  is intended for use in  conjunction with  Nag auxiliary
*  routines that  perform  interchange  operations,  such  as  pivoting.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*  TRANS
*           On entry,  SIDE  ( Left-hand side, or Right-hand side )  and
*           TRANS  ( Transpose, or No transpose )  specify the operation
*           to be performed as follows.
*
*           SIDE = 'L' or 'l'   and   TRANS = 'T' or 't'
*
*              Perform the operation   B := P'*B.
*
*           SIDE = 'L' or 'l'   and   TRANS = 'N' or 'n'
*
*              Perform the operation   B := P*B.
*
*           SIDE = 'R' or 'r'   and   TRANS = 'T' or 't'
*
*              Perform the operation   B := B*P'.
*
*           SIDE = 'R' or 'r'   and   TRANS = 'N' or 'n'
*
*              Perform the operation   B := B*P.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*
*           On entry, M must specify the order of the permutation matrix
*           P.  M must be at least zero.  When  M = 0  then an immediate
*           return is effected.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*
*           On entry,  N must specify the value of n. N must be at least
*           zero.  When  N = 0  then an  immediate  return is  effected.
*
*           Unchanged on exit.
*
*  PERM   - INTEGER array of DIMENSION at least ( n ).
*
*           Before  entry,  PERM  must  contain the  n  indices  for the
*           permutation matrices. index( i ) must satisfy
*
*              1 .le. index( i ) .le. m.
*
*           It is usual for index( i ) to be at least i, but this is not
*           necessary for this routine.
*
*           Unchanged on exit.
*
*  K      - INTEGER.
*
*           On entry with  SIDE = 'L' or 'l',  K must specify the number
*           of columns of B and on entry with  SIDE = 'R' or 'r', K must
*           specify the  number of rows of B.  K must be at least  zero.
*           When  K = 0  then an immediate return is effected.
*
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of  DIMENSION  ( LDB, ncolb ),  where
*           ncolb = k   when   SIDE = 'L' or 'l'  and   ncolb = m   when
*           SIDE = 'R' or 'r'.
*
*           Before entry  with  SIDE = 'L' or 'l',  the  leading  M by K
*           part  of  the  array   B  must  contain  the  matrix  to  be
*           transformed  and  before entry with  SIDE = 'R' or 'r',  the
*           leading  K by M part of the array  B must contain the matrix
*           to  be  transformed.  On  exit,  B  is  overwritten  by  the
*           transformed matrix.
*
*  LDB    - INTEGER.
*
*           On entry,  LDB  must specify  the  leading dimension  of the
*           array  B  as declared  in the  calling  (sub) program.  When
*           SIDE = 'L' or 'l'   then  LDB  must  be  at  least  m,  when
*           SIDE = 'R' or 'r'   then  LDB  must  be  at  least  k.
*           Unchanged on exit.
*
*
*  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
*
*  -- Written on 13-January-1986.
*     Sven Hammarling, Nag Central Office.
*
*
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, J, L
      LOGICAL            LEFT, NULL, RIGHT, TRNSP
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
      IF( MIN( M, N, K ).EQ.0 )
     $   RETURN
      LEFT  = ( SIDE .EQ.'L' ).OR.( SIDE .EQ.'l' )
      RIGHT = ( SIDE .EQ.'R' ).OR.( SIDE .EQ.'r' )
      NULL  = ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' )
      TRNSP = ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' )
      IF( LEFT )THEN
         IF( TRNSP )THEN
            DO 20, I = 1, N
               IF( PERM( I ).NE.I )THEN
                  L = PERM( I )
                  DO 10, J = 1, K
                     TEMP      = B( I, J )
                     B( I, J ) = B( L, J )
                     B( L, J ) = TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE IF( NULL )THEN
            DO 40, I = N, 1, -1
               IF( PERM( I ).NE.I )THEN
                  L = PERM( I )
                  DO 30, J = 1, K
                     TEMP      = B( L, J )
                     B( L, J ) = B( I, J )
                     B( I, J ) = TEMP
   30             CONTINUE
               END IF
   40       CONTINUE
         END IF
      ELSE IF( RIGHT )THEN
         IF( TRNSP )THEN
            DO 60, J = 1, N
               IF( PERM( J ).NE.J )THEN
                  L = PERM( J )
                  DO 50, I = 1, K
                     TEMP      = B( I, J )
                     B( I, J ) = B( L, J )
                     B( L, J ) = TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE IF( NULL )THEN
            DO 80, J = N, 1, -1
               IF( PERM( J ).NE.J )THEN
                  L = PERM( J )
                  DO 70, I = 1, K
                     TEMP      = B( L, J )
                     B( L, J ) = B( I, J )
                     B( I, J ) = TEMP
   70             CONTINUE
               END IF
   80       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEAP . ( F06QJF )
*
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DGEAPQ( TRANS, WHEREZ, M, N, A, LDA, ZETA,
     $                   NCOLB, B, LDB, WORK, INFORM )
      CHARACTER*1        TRANS, WHEREZ
      INTEGER            M, N, LDA, NCOLB, LDB, INFORM
      DOUBLE PRECISION   A( LDA, * ), ZETA( * ), B( LDB, * ), WORK( * )
C
C  1. Purpose
C     =======
C
C  DGEAPQ performs one of the transformations
C
C     B := Q'*B   or   B := Q*B,
C
C  where B is an m by ncolb matrix and Q is an m by m orthogonal matrix,
C  given as the product of  Householder transformation matrices, details
C  of  which are stored in the  m by n ( m.ge.n )  array  A  and, if the
C  parameter  WHEREZ = 'S' or 's', in the array ZETA.
C
C  This  routine is  intended for use following auxiliary linear algebra
C  routines such as  DGEQR , DGEHES and DSLTRI. ( See those routines for
C  example calls. )
C
C  2. Description
C     ===========
C
C  Q is assumed to be given by
C
C     Q = ( Q( p )*Q( p - 1 )*...*Q( 1 ) )',
C
C  Q( k ) being given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',   u( k ) = ( zeta( k ) ),
C                                             (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C
C  z( k )  must  be  supplied  in  the  kth  column  of  A  in  elements
C  a( k + 1, k ), ..., a( m, k )  and  zeta( k ) must be supplied either
C  in  a( k, k )  or in  zeta( k ), depending upon the parameter WHEREZ.
C
C  To obtain Q explicitly B may be set to I and premultiplied by Q. This
C  is more efficient than obtaining Q'.
C
C  3. Parameters
C     ==========
C
C  TRANS  - CHARACTER*1.
C
C           On entry, TRANS  specifies the operation to be performed  as
C           follows.
C
C           TRANS = ' ' or 'N' or 'n'
C
C              Perform the operation  B := Q*B.
C
C           TRANS = 'T' or 't' or 'C' or 'c'
C
C              Perform the operation  B := Q'*B.
C
C           Unchanged on exit.
C
C  WHEREZ - CHARACTER*1.
C
C           On entry, WHEREZ specifies where the elements of zeta are to
C           be found as follows.
C
C           WHEREZ = 'I' or 'i'
C
C              The elements of zeta are in A.
C
C           WHEREZ = 'S' or 's'
C
C              The elements of zeta are separate from A, in ZETA.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry, M  must specify the number of rows of A. M must be
C           at least n.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N  must specify the number of columns of A. N must
C           be  at least zero. When  N = 0  then an immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  A      - 'real' array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  stricly lower  triangular
C           part of the array  A  must contain details of the matrix  Q.
C           In  addition, when  WHEREZ = 'I' or 'i'  then  the  diagonal
C           elements of A must contain the elements of zeta.
C
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C
C           On  entry, LDA  must specify  the leading dimension  of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least m.
C
C           Unchanged on exit.
C
C  ZETA   - 'real' array of DIMENSION at least min( m - 1, n ).
C
C           Before entry with  WHEREZ = 'S' or 's', the array  ZETA must
C           contain the elements of the vector  zeta.
C
C           When  WHEREZ = 'I' or 'i', the array ZETA is not referenced.
C
C           Unchanged on exit.
C
C  NCOLB  - INTEGER.
C
C           On  entry, NCOLB  must specify  the number of columns of  B.
C           NCOLB  must  be  at  least  zero.  When  NCOLB = 0  then  an
C           immediate return is effected.
C
C           Unchanged on exit.
C
C  B      - 'real' array of DIMENSION ( LDB, ncolb ).
C
C           Before entry, the leading  M by NCOLB  part of  the array  B
C           must  contain  the matrix to be  transformed.
C
C           On  exit,  B  is  overwritten  by  the  transformed  matrix.
C
C  LDB    - INTEGER.
C
C           On  entry, LDB  must specify  the  leading dimension of  the
C           array  B as declared in the calling (sub) program. LDB  must
C           be at least m.
C
C           Unchanged on exit.
C
C  WORK   - 'real' array of DIMENSION at least ( ncolb ).
C
C           Used as internal workspace.
C
C  INFORM - INTEGER.
C
C           On  successful exit  INFORM  will be zero, otherwise  INFORM
C           will  be set to unity indicating that an input parameter has
C           been  incorrectly  set. See  the  next  section  for further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  INFORM = 1
C
C     One or more of the following conditions holds:
C
C        TRANS  .ne. ' ' or 'N' or 'n' or 'T' or 't' or 'C' or 'c'
C        WHEREZ .ne. 'I' or 'i' or 'S' or 's'
C        M      .lt. N
C        N      .lt. 0
C        LDA    .lt. M
C        NCOLB  .lt. 0
C        LDB    .lt. M
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 15-November-1984.
C     Sven Hammarling, Nag Central Office.
C
      EXTERNAL           DGEMV , DGER
      INTRINSIC          MIN
      INTEGER            J     , K     , KK    , LB
      DOUBLE PRECISION   TEMP
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )

*     Check the input parameters.

      IF( MIN( N, NCOLB ).EQ.0 )THEN
         INFORM = 0
         RETURN
      END IF
      IF( ( ( TRANS .NE.' ' ).AND.
     $      ( TRANS .NE.'N' ).AND.( TRANS .NE.'n' ).AND.
     $      ( TRANS .NE.'T' ).AND.( TRANS .NE.'t' ).AND.
     $      ( TRANS .NE.'C' ).AND.( TRANS .NE.'c' )      ).OR.
     $    ( ( WHEREZ.NE.'I' ).AND.( WHEREZ.NE.'i' ).AND.
     $      ( WHEREZ.NE.'S' ).AND.( WHEREZ.NE.'s' )      ).OR.
     $    ( M.LT.N ).OR.( N.LT.0 ).OR.( LDA.LT.M ).OR.
     $    ( NCOLB.LT.0 ).OR.( LDB.LT.M )                      )THEN
         INFORM = 1
         RETURN
      END IF

*     Perform the transformation.

      LB = LDB
      DO 20, KK = 1, MIN( M - 1, N )
         IF( ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' ).OR.
     $       ( TRANS.EQ.'C' ).OR.( TRANS.EQ.'c' )     )THEN

*           Q'*B = Q( p )*...*Q( 2 )*Q( 1 )*B,     p = min( m - 1, n ).

            K = KK
         ELSE

*           Q*B  = Q( 1 )'*Q( 2 )'*...*Q( p )'*B,  p = min( m - 1, n ).
*           Note that  Q( k )' = Q( k ).

            K = MIN( N, M - 1 ) + 1 - KK
         END IF
         IF( ( WHEREZ.EQ.'S' ).OR.( WHEREZ.EQ.'s' ) )THEN
            TEMP      = A( K, K )
            A( K, K ) = ZETA( K )
         END IF

*        If ZETA( k ) is zero then Q( k ) = I and we can skip the kth
*        transformation.

         IF( A( K, K ).GT.ZERO )THEN
            IF( NCOLB.EQ.1 )
     $         LB = M - K + 1

*           Let C denote the bottom ( m - k + 1 ) by ncolb part of B.

*           First form  work = C'*u.

            DO 10, J = 1, NCOLB
               WORK( J ) = ZERO
   10       CONTINUE
            CALL DGEMV ( 'Transpose', M - K + 1, NCOLB,
     $                   ONE, B( K, 1 ), LB, A( K, K ), 1,
     $                   ZERO, WORK, 1 )

*           Now form  C := C - u*work'.

            CALL DGER  ( M - K + 1, NCOLB, -ONE, A( K, K ), 1,
     $                   WORK, 1, B( K, 1 ), LB )
         END IF

*        Restore the diagonal element of A.

         IF( ( WHEREZ.EQ.'S' ).OR.( WHEREZ.EQ.'s' ) )
     $      A( K, K ) = TEMP
   20 CONTINUE

      INFORM = 0
      RETURN

*     End of DGEAPQ.

      END
