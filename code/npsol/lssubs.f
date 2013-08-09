*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     file  lssubs.f
*
*     lsadd    lsadds   lsbnds   lschol   lscore   lscrsh   lsdel
*     lsdflt   lsfeas   lsfile   lsgetp   lsgset   lskey    lsloc
*     lsmove   lsmuls   lsoptn   lsprt    lssetx   lssol
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSADD ( UNITQ,
     $                   INFORM, IFIX, IADD, JADD,
     $                   NACTIV, NZ, NFREE, NRANK, NRES, NGQ,
     $                   N, LDA, LDZY, LDR, LDT,
     $                   KX, CONDMX,
     $                   A, R, T, RES, GQM, ZY,
     $                   W, C, S )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            UNITQ
      INTEGER            KX(N)
      DOUBLE PRECISION   A(LDA,*), R(LDR,*), T(LDT,*),
     $                   RES(N,*), GQM(N,*), ZY(LDZY,*)
      DOUBLE PRECISION   W(N), C(N), S(N)

C***********************************************************************
C     LSADD   updates the factorization,  A(free) * (Z Y) = (0 T),  when
C     a constraint is added to the working set.  If  NRANK .gt. 0, the
C     factorization  ( R ) = PCQ  is also updated,  where  C  is the
C                    ( 0 )
C     least squares matrix,  R  is upper-triangular,  and  P  is an
C     orthogonal matrix.  The matrices  C  and  P  are not stored.
C
C     There are three separate cases to consider (although each case
C     shares code with another)...
C
C     (1) A free variable becomes fixed on one of its bounds when there
C         are already some general constraints in the working set.
C
C     (2) A free variable becomes fixed on one of its bounds when there
C         are only bound constraints in the working set.
C
C     (3) A general constraint (corresponding to row  IADD  of  A) is
C         added to the working set.
C
C     In cases (1) and (2), we assume that  KX(IFIX) = JADD.
C     In all cases,  JADD  is the index of the constraint being added.
C
C     If there are no general constraints in the working set,  the
C     matrix  Q = (Z Y)  is the identity and will not be touched.
C
C     If  NRES .GT. 0,  the row transformations are applied to the rows
C     of the  (N by NRES)  matrix  RES.
C     If  NGQ .GT. 0,  the column transformations are applied to the
C     columns of the  (NGQ by N)  matrix  GQM'.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October--1984.
C     Level-2 matrix routines added 25-Apr-1988.
C     This version of  LSADD  dated 28-May-1988.
C***********************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON    /SOL5CM/ ASIZE, DTMAX, DTMIN

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      LOGICAL            BOUND , OVERFL
      EXTERNAL           DDOT  , DDIV  , DNRM2
      INTRINSIC          MAX   , MIN
      PARAMETER         (ZERO = 0.0D+0, ONE = 1.0D+0)

*     If the condition estimator of the updated factors is greater than
*     CONDBD,  a warning message is printed.

      CONDBD = ONE / EPSPT9

      OVERFL = .FALSE.
      BOUND  = JADD .LE. N

      IF (BOUND) THEN
*        ===============================================================
*        A simple bound has entered the working set.  IADD  is not used.
*        ===============================================================
         IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $      WRITE (NOUT, 1010) NACTIV, NZ, NFREE, IFIX, JADD, UNITQ
         NANEW = NACTIV

         IF (UNITQ) THEN

*           Q  is not stored, but KX defines an ordering of the columns
*           of the identity matrix that implicitly define  Q.
*           Define the sequence of pairwise interchanges P that moves
*           the newly-fixed variable to position NFREE.
*           Reorder KX accordingly.

            DO 100, I = 1, NFREE-1
               IF (I .GE. IFIX) THEN
                  W (I) = I + 1
                  KX(I) = KX(I+1)
               ELSE
                  W(I) = I
               END IF
  100       CONTINUE

         ELSE
*           ------------------------------------------------------------
*           Q  is stored explicitly.
*           ------------------------------------------------------------
*           Set  W = the  (IFIX)-th  row of  Q.
*           Move the  (NFREE)-th  row of  Q  to position  IFIX.

            CALL DCOPY ( NFREE, ZY(IFIX,1), LDZY, W, 1 )
            IF (IFIX .LT. NFREE) THEN
               CALL DCOPY ( NFREE, ZY(NFREE,1), LDZY, ZY(IFIX,1), LDZY )
               KX(IFIX) = KX(NFREE)
            END IF
         END IF
         KX(NFREE) = JADD
      ELSE
*        ===============================================================
*        A general constraint has entered the working set.
*        IFIX  is not used.
*        ===============================================================
         IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $      WRITE (NOUT, 1020) NACTIV, NZ, NFREE, IADD, JADD, UNITQ

         NANEW  = NACTIV + 1

*        Transform the incoming row of  A  by  Q'.  Use C as workspace.

         CALL DCOPY ( N, A(IADD,1), LDA, W, 1 )
         CALL CMQMUL( 8, N, NZ, NFREE, LDZY, UNITQ, KX, W, ZY, C )

*        Check that the incoming row is not dependent upon those
*        already in the working set.

         DTNEW  = DNRM2 ( NZ, W, 1 )
         IF (NACTIV .EQ. 0) THEN

*           This is the only general constraint in the working set.

            COND   = DDIV  ( ASIZE, DTNEW, OVERFL )
            TDTMAX = DTNEW
            TDTMIN = DTNEW
         ELSE

*           There are already some general constraints in the working
*           set. Update the estimate of the condition number.

            TDTMAX = MAX( DTNEW, DTMAX )
            TDTMIN = MIN( DTNEW, DTMIN )
            COND   = DDIV  ( TDTMAX, TDTMIN, OVERFL )
         END IF

         IF (COND .GT. CONDMX  .OR.  OVERFL) GO TO 900

         IF (UNITQ) THEN

*           First general constraint added.  Set  Q = I.

            CALL F06QHF( 'General', NFREE, NFREE, ZERO, ONE, ZY, LDZY )
            UNITQ  = .FALSE.
         END IF
      END IF

      IF (BOUND) THEN
         NPIV  = NFREE
      ELSE
         NPIV  = NZ
      END IF

      NT = MIN( NRANK, NPIV )

      IF (UNITQ) THEN
*        ---------------------------------------------------------------
*        Q (i.e., ZY) is not stored explicitly.
*        Apply the sequence of pairwise interchanges P that moves the
*        newly-fixed variable to position NFREE.
*        ---------------------------------------------------------------
         IF (NGQ .GT. 0)
     $      CALL F06QKF( 'Left', 'Transpose', NFREE-1, W, NGQ, GQM, N )
            
         IF (NRANK .GT. 0) THEN

*           Apply the pairwise interchanges to the triangular part of R.
*           The subdiagonal elements generated by this process are
*           stored in  s(1), s(2), ..., s(nt-1).

            CALL F06QNF( 'Right', N, IFIX, NT, S, R, LDR )

            IF (NT .LT. NPIV) THEN

*              R is upper trapezoidal.  Apply the interchanges in
*              columns  nt  thru  npiv.

               DO 200, I = IFIX, NT-1
                  W(I) = I
  200          CONTINUE

               CALL F06QKF( 'Right', 'Normal', NFREE-1, W, NT, R, LDR )
            END IF
            
*           Eliminate the subdiagonal elements of R with a left-hand
*           sweep of rotations P2 in planes (1,2), (2,3), ...,(nt-1,nt).
*           Apply P2 to RES.

            CALL F06QRF( 'Left ', N, IFIX, NT, C, S, R, LDR )
            IF (NRES .GT. 0) 
     $         CALL F06QXF( 'Left', 'Variable', 'Forwards', NT, NRES,
     $                      IFIX, NT, C, S, RES, N )
         END IF
      ELSE
*        ---------------------------------------------------------------
*        Full matrix Q.  Define a sweep of plane rotations P such that
*                           Pw = beta*e(npiv).
*        The rotations are applied in the planes (1,2), (2,3), ...,
*        (npiv-1,npiv).  The rotations must be applied to ZY, R, T
*        and GQM'.
*        ---------------------------------------------------------------
         CALL F06FQF( 'Varble', 'Forwrds', NPIV-1, W(NPIV), W, 1, C, S )

         IF (BOUND  .AND.  NACTIV .GT. 0) THEN

            CALL DCOPY ( NACTIV, S(NZ), 1, W(NZ), 1 )
         
            S(       NZ  ) = S(NZ)*T(NACTIV,NZ+1)
            T(NACTIV,NZ+1) = C(NZ)*T(NACTIV,NZ+1)
         
            CALL F06QZF( 'Create', NACTIV, 1, NACTIV, C(NZ+1), S(NZ+1),
     $                   T(1,NZ+1), LDT )
            CALL DCOPY ( NACTIV, S(NZ), 1, T(NACTIV,NZ), LDT-1 )
         
            CALL DCOPY ( NACTIV, W(NZ), 1, S(NZ), 1 )
         END IF

         IF (NGQ .GT. 0)
     $      CALL F06QXF( 'Left ', 'Variable', 'Forwards', NPIV , NGQ,
     $                   1, NPIV, C, S, GQM, N )
         CALL F06QXF( 'Right', 'Variable', 'Forwards', NFREE, NFREE,
     $                1, NPIV, C, S, ZY, LDZY )

         IF (NRANK .GT. 0) THEN

*           Apply the rotations to the triangular part of R.
*           The subdiagonal elements generated by this process are
*           stored in  s(1),  s(2), ..., s(nt-1).

            NT = MIN( NRANK, NPIV )
            CALL F06QVF( 'Right', N, 1, NT, C, S, R, LDR )

            IF (NT .LT. NPIV) THEN

*              R is upper trapezoidal.  Pretend R is (nt x n) and
*              apply the rotations in columns  nt  thru  npiv.

               CALL F06QXF( 'Right', 'Variable', 'Forwards', NT, N,
     $                      NT, NPIV, C, S, R, LDR )
            END IF

*           Eliminate the subdiagonal elements of R with a left-hand
*           sweep of rotations P2 in planes (1,2), (2,3), ...,(nt-1,nt).
*           Apply P2 to RES.

            CALL F06QRF( 'Left ', N, 1, NT, C, S, R, LDR )
            IF (NRES .GT. 0)
     $         CALL F06QXF( 'Left', 'Variable', 'Forwards', NT, NRES,
     $                      1, NT, C, S, RES, N )
         END IF

         IF (BOUND) THEN

*           The last row and column of ZY has been transformed to plus
*           or minus the unit vector E(NFREE).  We can reconstitute the
*           columns of GQM and R corresponding to the new fixed variable.

            IF (W(NFREE) .LT. ZERO) THEN
               NF = MIN( NRANK, NFREE )
               IF (NF  .GT. 0) CALL DSCAL ( NF , -ONE,   R(1,NFREE), 1 )
               IF (NGQ .GT. 0) CALL DSCAL ( NGQ, -ONE, GQM(NFREE,1), N )
            END IF

*           ------------------------------------------------------------
*           The diagonals of T have been altered.  Recompute the
*           largest and smallest values.
*           ------------------------------------------------------------
            IF (NACTIV .GT. 0) THEN
               CALL DCOND( NACTIV, T(NACTIV,NZ), LDT-1, TDTMAX, TDTMIN )
               COND   = DDIV  ( TDTMAX, TDTMIN, OVERFL )
            END IF
         ELSE
*           ------------------------------------------------------------
*           General constraint.  Install the new row of T.
*           ------------------------------------------------------------
            CALL DCOPY ( NANEW, W(NZ), 1, T(NANEW,NZ), LDT )
         END IF
      END IF

*     ==================================================================
*     Prepare to exit.  Check the magnitude of the condition estimator.
*     ==================================================================
  900 IF (NANEW .GT. 0) THEN
         IF (COND .LT. CONDMX  .AND.  .NOT. OVERFL) THEN

*           The factorization has been successfully updated.

            INFORM = 0
            DTMAX  = TDTMAX
            DTMIN  = TDTMIN
            IF (COND .GE. CONDBD) WRITE (NOUT, 2000) JADD
         ELSE

*           The proposed working set appears to be linearly dependent.

            INFORM = 1
            IF (LSDBG  .AND.  ILSDBG(1) .GT. 0) THEN
               WRITE( NOUT, 3000 )
               IF (BOUND) THEN
                  WRITE (NOUT, 3010) ASIZE, DTMAX, DTMIN
               ELSE
                  IF (NACTIV .GT. 0) THEN
                     WRITE (NOUT, 3020) ASIZE, DTMAX, DTMIN, DTNEW
                  ELSE
                     WRITE (NOUT, 3030) ASIZE, DTNEW
                  END IF
               END IF
            END IF
         END IF
      END IF

      RETURN

 1010 FORMAT(/ ' //LSADD //  Simple bound added.'
     $       / ' //LSADD //  NACTIV    NZ NFREE  IFIX  JADD UNITQ'
     $       / ' //LSADD //  ', 5I6, L6 )
 1020 FORMAT(/ ' //LSADD //  General constraint added.           '
     $       / ' //LSADD //  NACTIV    NZ NFREE  IADD  JADD UNITQ'
     $       / ' //LSADD //  ', 5I6, L6 )
 2000 FORMAT(/ ' XXX  Serious ill-conditioning in the working set',
     $         ' after adding constraint ',  I5
     $       / ' XXX  Overflow may occur in subsequent iterations.'//)
 3000 FORMAT(/ ' //LSADD //  Dependent constraint rejected.' )
 3010 FORMAT(/ ' //LSADD //     ASIZE     DTMAX     DTMIN        '
     $       / ' //LSADD //', 1P, 3E10.2 )
 3020 FORMAT(/ ' //LSADD //     ASIZE     DTMAX     DTMIN     DTNEW'
     $       / ' //LSADD //', 1P, 4E10.2 )
 3030 FORMAT(/ ' //LSADD //     ASIZE     DTNEW'
     $       / ' //LSADD //', 1P, 2E10.2 )

*     End of  LSADD .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSADDS( UNITQ, VERTEX,
     $                   INFORM, K1, K2, NACTIV, NARTIF, NZ, NFREE,
     $                   NRANK, NREJTD, NRES, NGQ,
     $                   N, LDZY, LDA, LDR, LDT,
     $                   ISTATE, KACTIV, KX,
     $                   CONDMX,
     $                   A, R, T, RES, GQM, ZY,
     $                   W, C, S )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            UNITQ, VERTEX
      INTEGER            ISTATE(*), KACTIV(N), KX(N)
      DOUBLE PRECISION   CONDMX
      DOUBLE PRECISION   A(LDA,*), R(LDR,*),
     $                   T(LDT,*), RES(N,*), GQM(N,*), ZY(LDZY,*)
      DOUBLE PRECISION   W(N), C(N), S(N)

C***********************************************************************
C     LSADDS  includes general constraints K1 thru K2 as new rows of
C     the TQ factorization stored in T, ZY.  If NRANK is nonzero, the
C     changes in Q are reflected in NRANK by N triangular factor R such
C     that
C                         C  =  P ( R ) Q,
C                                 ( 0 )
C     where  P  is orthogonal.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written  October-31-1984.
C     This version of LSADDS dated  16-May-1988.
C***********************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
      COMMON    /SOL5CM/ ASIZE, DTMAX, DTMIN

      EXTERNAL           DNRM2
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      RTMAX  = WMACH(8)

*     Estimate the condition number of the constraints that are not
*     to be refactorized.

      IF (NACTIV .EQ. 0) THEN
         DTMAX = ZERO
         DTMIN = ONE
      ELSE
         CALL DCOND ( NACTIV, T(NACTIV,NZ+1), LDT-1, DTMAX, DTMIN )
      END IF

      DO 200, K = K1, K2
         IADD = KACTIV(K)
         JADD = N + IADD
         IF (NACTIV .LT. NFREE) THEN

            CALL LSADD ( UNITQ,
     $                   INFORM, IFIX, IADD, JADD,
     $                   NACTIV, NZ, NFREE, NRANK, NRES, NGQ,
     $                   N, LDA, LDZY, LDR, LDT,
     $                   KX, CONDMX,
     $                   A, R, T, RES, GQM, ZY,
     $                   W, C, S )

            IF (INFORM .EQ. 0) THEN
               NACTIV = NACTIV + 1
               NZ     = NZ     - 1
            ELSE
               ISTATE(JADD) =   0
               KACTIV(K)    = - KACTIV(K)
            END IF
         END IF
  200 CONTINUE

      IF (NACTIV .LT. K2) THEN

*        Some of the constraints were classed as dependent and not
*        included in the factorization.  Re-order the part of KACTIV
*        that holds the indices of the general constraints in the
*        working set.  Move accepted indices to the front and shift
*        rejected indices (with negative values) to the end.

         L      = K1 - 1
         DO 300, K = K1, K2
            I         = KACTIV(K)
            IF (I .GE. 0) THEN
               L      = L + 1
               IF (L .NE. K) THEN
                  ISWAP     = KACTIV(L)
                  KACTIV(L) = I
                  KACTIV(K) = ISWAP
               END IF
            END IF
  300    CONTINUE

*        If a vertex is required, add some temporary bounds.
*        We must accept the resulting condition number of the working
*        set.

         IF (VERTEX) THEN
            CNDMAX = RTMAX
            NZADD  = NZ
            DO 320, IARTIF = 1, NZADD
               IF (UNITQ) THEN
                  IFIX = NFREE
                  JADD = KX(IFIX)
               ELSE
                  ROWMAX = ZERO
                  DO 310, I = 1, NFREE
                     RNORM = DNRM2 ( NZ, ZY(I,1), LDZY )
                     IF (ROWMAX .LT. RNORM) THEN
                        ROWMAX = RNORM
                        IFIX   = I
                     END IF
  310             CONTINUE
                  JADD = KX(IFIX)

                  CALL LSADD ( UNITQ,
     $                         INFORM, IFIX, IADD, JADD,
     $                         NACTIV, NZ, NFREE, NRANK, NRES, NGQ,
     $                         N, LDA, LDZY, LDR, LDT,
     $                         KX, CNDMAX,
     $                         A, R, T, RES, GQM, ZY,
     $                         W, C, S  )
               END IF
               NFREE  = NFREE  - 1
               NZ     = NZ     - 1
               NARTIF = NARTIF + 1
               ISTATE(JADD) = 4
  320       CONTINUE
         END IF
      END IF

      NREJTD = K2 - NACTIV

      RETURN

*     End of  LSADDS.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSBNDS( UNITQ,
     $                   INFORM, NZ, NFREE, NRANK, NRES, NGQ,
     $                   N, LDZY, LDA, LDR, LDT,
     $                   ISTATE, KX, CONDMX,
     $                   A, R, T, RES, GQM, ZY,
     $                   W, C, S )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            UNITQ
      INTEGER            ISTATE(*), KX(N)
      DOUBLE PRECISION   CONDMX
      DOUBLE PRECISION   A(LDA,*), R(LDR,*),
     $                   T(LDT,*), RES(N,*), GQM(N,*), ZY(LDZY,*)
      DOUBLE PRECISION   W(N), C(N), S(N)

************************************************************************
*     LSBNDS updates the factor R as KX is reordered to reflect the
*     status of the bound constraints given by ISTATE.  KX is reordered
*     so that the fixed variables come last.  One of two alternative
*     methods are used to reorder KX. One method needs fewer accesses 
*     to KX, the other gives a matrix  Rz  with more rows and columns.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written  30-December-1985.
*     This version of LSBNDS dated 13-May-88.
************************************************************************

      NFIXED = N - NFREE

      IF (NRANK .LT. N  .AND.  NRANK .GT. 0) THEN
*        ---------------------------------------------------------------
*        R is specified but singular.  Try and keep the dimension of Rz
*        as large as possible.
*        ---------------------------------------------------------------
         NACTV = 0
         NFREE = N
         NZ    = N

         J     = N
*+       WHILE (J .GT. 0  .AND.  N-NFREE .LT. NFIXED) DO
  100    IF    (J .GT. 0  .AND.  N-NFREE .LT. NFIXED) THEN
            IF (ISTATE(J) .GT. 0) THEN
               JADD = J
               DO 110, IFIX = NFREE, 1, -1
                  IF (KX(IFIX) .EQ. JADD) GO TO 120
  110          CONTINUE

*              Add bound JADD.

  120          CALL LSADD ( UNITQ,
     $                      INFORM, IFIX, IADD, JADD,
     $                      NACTV, NZ, NFREE, NRANK, NRES, NGQ,
     $                      N, LDA, LDZY, LDR, LDT,
     $                      KX, CONDMX,
     $                      A, R, T, RES, GQM, ZY,
     $                      W, C, S )

               NFREE = NFREE - 1
               NZ    = NZ    - 1
            END IF
            J = J - 1
            GO TO 100
*+       END WHILE
         END IF
      ELSE     
*        ---------------------------------------------------------------
*        R is of full rank,  or is not specified.
*        ---------------------------------------------------------------
         IF (NFIXED .GT. 0) THEN

*           Order KX so that the free variables come first.

            LSTART = NFREE + 1
            DO 250, K = 1, NFREE
               J = KX(K)
               IF (ISTATE(J) .GT. 0) THEN
                  DO 220, L = LSTART, N
                     J2 = KX(L)
                     IF (ISTATE(J2) .EQ. 0) GO TO 230
  220             CONTINUE

  230             KX(K)  = J2
                  KX(L)  = J
                  LSTART = L + 1

                  IF (NRANK .GT. 0)
     $               CALL CMRSWP( N, NRES, NRANK, LDR, K, L,
     $                            R, RES, C, S )
               END IF
  250       CONTINUE

         END IF
         NZ = NFREE
      END IF

      RETURN

*     End of  LSBNDS.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lschol( ldh, n, nrank, tolrnk, kx, h, inform )

      implicit           double precision (a-h,o-z)
      integer            kx(*)
      double precision   h(ldh,*)

************************************************************************
*     LSCHOL  forms the Cholesky factorization of the positive
*     semi-definite matrix H such that
*                   PHP'  =  R'R
*     where  P  is a permutation matrix and  R  is upper triangular.
*     The permutation P is chosen to maximize the diagonal of R at each
*     stage.  Only the diagonal and super-diagonal elements of H are
*     used.
*
*     Output:
*
*         INFORM = 0   the factorization was computed successfully,
*                      with the Cholesky factor written in the upper
*                      triangular part of H and P stored in KX.
*                  1   the matrix H was indefinite.
*
*     Original version of LSCHOL dated  2-February-1981.
*     Level 2 Blas added 29-June-1986.
*     This version of LSCHOL dated  26-Jun-1989. 
************************************************************************
      common    /sol1cm/ nout
      intrinsic          abs   , max   , sqrt
      external           idamax
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      inform = 0
      nrank  = 0

*     Main loop for computing rows of  R.

      do 200, j = 1, n

*        Find maximum available diagonal.    

         kmax = j - 1 + idamax( n-j+1, h(j,j), ldh+1 )
         dmax = h(kmax,kmax)

         if (dmax .le. tolrnk*abs(h(1,1))) go to 300

*        Perform a symmetric interchange if necessary.

         if (kmax .ne. j) then
            k        = kx(kmax)
            kx(kmax) = kx(j)
            kx(j)    = k

            call dswap ( kmax-j, h(j+1,kmax), 1, h(j,j+1 ), ldh )
            call dswap ( j     , h(1  ,j   ), 1, h(1,kmax), 1   )
            call dswap ( n-kmax+1, h(kmax,kmax), ldh,
     $                             h(j,kmax)   , ldh )

         end if

*        Set the diagonal of  R.

         d      = sqrt( dmax )
         h(j,j) = d
         nrank  = nrank + 1

         if (j .lt. n) then

*           Set the super-diagonal elements of this row of R and update
*           the elements of the block that is yet to be factorized.
                                                          
            call dscal ( n-j,   (one/d), h(j  ,j+1), ldh )
            call dsyr  ( 'u', n-j, -one, h(j  ,j+1), ldh,
     $                                   h(j+1,j+1), ldh )
         end if

  200 continue
*     ------------------------------------------------------------------
*     Check for the semi-definite case.
*     ------------------------------------------------------------------
  300 if (nrank .lt. n) then

*        Find the largest element in the unfactorized block.

         supmax = zero
         do 310, i = j, n-1
            k      = i + idamax( n-i, h(i,i+1), ldh )
            supmax = max( supmax, abs(h(i,k)) )
  310    continue

         if (supmax .gt. tolrnk*abs(h(1,1))) then
            write (nout, 1000) dmax, supmax
            inform = 1
         end if
      end if

      return

 1000 format(' XXX  Hessian appears to be indefinite.'
     $      /' XXX  Maximum diagonal and off-diagonal ignored',
     $             ' in the Cholesky factorization:', 1p, 2e22.14 )

*     End of LSCHOL.

      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lscore( prbtyp, named, names, linobj, unitQ,
     $                   inform, iter, jinf, nclin, nctotl,
     $                   nactiv, nfree, nrank, nz, nrz,
     $                   n, ldA, ldR,
     $                   istate, kactiv, kx,
     $                   ctx, ssq, ssq1, suminf, numinf, xnorm,
     $                   bl, bu, A, clamda, ax,
     $                   featol, R, x, iw, w )

      implicit           double precision(a-h,o-z)
      character*2        prbtyp
      character*8        names(*)
      integer            istate(nctotl), kactiv(n), kx(n)
      integer            iw(*)
      double precision   bl(nctotl), bu(nctotl), a(ldA,*),
     $                   clamda(nctotl), ax(*),
     $                   featol(nctotl), r(ldR,*), x(n)
      double precision   w(*)
      logical            named, linobj, unitQ
************************************************************************
*     LSCORE  is a subroutine for linearly constrained linear-least
*     squares.  On entry, it is assumed that an initial working set of
*     linear constraints and bounds is available.
*     The arrays ISTATE, KACTIV and KX will have been set accordingly
*     and the arrays T and ZY will contain the TQ factorization of
*     the matrix whose rows are the gradients of the active linear
*     constraints with the columns corresponding to the active bounds
*     removed.  the TQ factorization of the resulting (NACTIV by NFREE)
*     matrix is  A(free)*Q = (0 T),  where Q is (NFREE by NFREE) and T
*     is reverse-triangular.
*
*     Values of ISTATE(J) for the linear constraints.......
*
*     ISTATE(J)
*     ---------
*          0    constraint J is not in the working set.
*          1    constraint J is in the working set at its lower bound.
*          2    constraint J is in the working set at its upper bound.
*          3    constraint J is in the working set as an equality.
*
*     Constraint J may be violated by as much as FEATOL(J).
*
*     Systems Optimization Laboratory, Stanford University.
*     This version of  LSCORE  dated  25-Aug-1991.
*
*     Copyright  1984  Stanford University.
*
*  This material may be reproduced by or for the U.S. Government pursu-
*  ant to the copyright license under DAR clause 7-104.9(a) (1979 Mar).
*
*  This material is based upon work partially supported by the National
*  Science Foundation under grants MCS-7926009 and ECS-8012974; the
*  Department of Energy Contract AM03-76SF00326, PA No. DE-AT03-
*  76ER72018; and the Army Research Office Contract DAA29-79-C-0110.
************************************************************************
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
      common    /sol1cm/ nout
      common    /sol3cm/ lennam, ldT   , ncolt , ldzy
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5cm/ asize, dtmax  , dtmin

      integer            locls
      parameter         (lenls = 20)
      common    /sol1ls/ locls(lenls)

      logical            cmdbg, lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg
      common    /cmdebg/ icmdbg(ldbg), cmdbg
*-----------------------------------------------------------------------
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   idbgls, itmax1, itmax2, lcrash, ldbgls, lprob ,
     $                   msgls , nn    , nnclin, nprob , ipadls(20)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolrnk, rpadls(23)

      equivalence       (iprmls(1), idbgls), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*-----------------------------------------------------------------------
      equivalence   (msgls , msglvl), (idbgls, idbg), (ldbgls, msgdbg)

      external           ddiv  , ddot  , dnrm2
      intrinsic          abs   , max   , sqrt
      logical            convrg, cyclin, error , firstv, hitcon,
     $                   hitlow, needfg, overfl, prnt  , prnt1 , rowerr
      logical            singlr, stall , statpt, unbndd, uncon , unitgz,
     $                   weak
      parameter        ( zero   =0.0d+0, half   =0.5d+0, one   =1.0d+0 )
      parameter        ( mrefn  =1     , mstall =50                    )

*     Specify the machine-dependent parameters.

      epsmch = wmach(3)
      flmax  = wmach(7)
      rtmax  = wmach(8)

      lanorm = locls( 2)
      lap    = locls( 3)
      lpx    = locls( 4)
      lres   = locls( 5)
      lres0  = locls( 6)
      lhz    = locls( 7)
      lgq    = locls( 8)
      lcq    = locls( 9)
      lrlam  = locls(10)
      lt     = locls(11)
      lzy    = locls(12)
      lwtinf = locls(13)
      lwrk   = locls(14)

*     Set up the adresses of the contiguous arrays  ( res0, res )
*     and  ( gq, cq ).

      nres   = 0
      if (nrank .gt. 0) nres = 2
      ngq    = 1
      if (linobj) ngq = 2

*     Initialize.

      irefn  =   0
      iter   =   0
      itmax  =   itmax1
      jadd   =   0
      jdel   =   0
      ncnln  =   0
      nphase =   1
      nstall =   0
      numinf = - 1
      nrz    =   0

      alfa   = zero
      condmx = flmax
      drzmax = one
      drzmin = one
      ssq    = zero

      cyclin = .false.
      error  = .false.
      firstv = .false.
      prnt   = .true.
      prnt1  = .true.
      needfg = .true.
      stall  = .true.
      uncon  = .false.
      unbndd = .false.

*     If debug output is required,  print nothing until iteration IDBG.

      msgsvd = msglvl
      if (idbg .gt. 0  .and.  idbg .le. itmax) then
         msglvl = 0
      end if

*======================== start of the main loop =======================
*
*      cyclin = false
*      unbndd = false
*      error  = false
*      k      = 0
*
*      repeat
*            repeat
*                  compute Z'g,  print details of this iteration
*                  stat pt = (Z'g .eq. 0)
*                  if (not stat pt) then
*                     error =  k .ge. itmax
*                     if (not error) then
*                        compute p, alfa
*                        error = unbndd  or  cyclin
*                        if (not error) then
*                           k = k + 1
*                           x = x + alfa p
*                           if (feasible) update Z'g
*                           if necessary, add a constraint
*                        end if
*                     end if
*                  end if
*            until  stat pt  or  error
*
*            compute lam1, lam2, smllst
*            optmul =  smllst .gt. 0
*            if ( not (optmul .or. error) ) then
*                  delete an artificial or regular constraint
*            end if
*      until optmul  or  error
*
*=======================================================================

*     repeat
*        repeat
  100       if (needfg) then
               if (nrank .gt. 0) then
                  resnrm = dnrm2 ( nrank, w(lres), 1 )
                  ssq    = half*(ssq1**2 + resnrm**2 )
               end if

               if (numinf .ne. 0) then

*                 Compute the transformed gradient of either the sum of
*                 of infeasibilities or the objective.  Initialize
*                 singlr and unitgz.

                  call lsgset( prbtyp, linobj, singlr, unitgz, unitQ,
     $                         n, nclin, nfree,
     $                         ldA, ldzy, ldR, nrank, nz, nrz,
     $                         istate, kx,
     $                         bigbnd, tolrnk, numinf, suminf,
     $                         bl, bu, A, w(lres), featol,
     $                         w(lgq), w(lcq), R, x, w(lwtinf),
     $                         w(lzy), w(lwrk) )

                  if (prbtyp .ne. 'FP'  .and.  numinf .eq. 0
     $                                  .and.  nphase .eq. 1) then
                     itmax  = iter + itmax2
                     nphase = 2
                  end if
               end if
            end if

            gznorm = zero
            if (nz  .gt. 0 ) gznorm = dnrm2 ( nz, w(lgq), 1 )

            if (nrz .eq. nz) then
               gz1nrm = gznorm
            else
               gz1nrm = zero
               if (nrz .gt. 0) gz1nrm = dnrm2 ( nrz, w(lgq), 1 )
            end if

            gfnorm = gznorm
            if (nfree .gt. 0  .and.  nactiv .gt. 0)
     $         gfnorm = dnrm2 ( nfree, w(lgq), 1 )

*           ------------------------------------------------------------
*           Print the details of this iteration.
*           ------------------------------------------------------------
*           Define small quantities that reflect the size of x, R and
*           the constraints in the working set.  If feasible,  estimate
*           the rank and condition number of Rz1.
*           Note that nrz .le. nrank + 1.

            if (nrz .eq. 0) then
               singlr = .false.
            else
               if (numinf .gt. 0  .or.  nrz .gt. nrank) then
                  absrzz = zero
                  singlr = .true.
               else
                  call dcond ( nrz, R, ldR+1, drzmax, drzmin )
                  absrzz = abs( r(nrz,nrz) )
                  rownrm = dnrm2 ( n, r(1,1), ldR )
                  singlr =       absrzz      .le. drzmax*tolrnk 
     $                     .or.  rownrm      .le.        tolrnk 
     $                     .or.  abs(r(1,1)) .le. rownrm*tolrnk
               end if

               if (lsdbg  .and.  ilsdbg(1) .gt. 0)
     $            write (nout, 9100) singlr, absrzz, drzmax, drzmin
            end if

            condrz = ddiv  ( drzmax, drzmin, overfl )
            condt  = one
            if (nactiv .gt. 0)
     $         condt  = ddiv  ( dtmax , dtmin , overfl )

            if (prnt) then
               call lsprt ( prbtyp, prnt1, isdel, iter, jadd, jdel,
     $                      msglvl, nactiv, nfree, n, nclin,
     $                      nrank, ldR, ldT, nz, nrz,
     $                      istate,
     $                      alfa, condrz, condt, gfnorm, gznorm, gz1nrm,
     $                      numinf, suminf, ctx, ssq,
     $                      ax, R, w(lt), x, w(lwrk) )
               jdel  = 0
               jadd  = 0
               alfa  = zero
            end if

            if (numinf .gt. 0) then
               dinky  = zero
            else
               objsiz = one  + abs( ssq + ctx )
               wssize = zero
               if (nactiv .gt. 0) wssize = dtmax
               dinky  = epspt8 * max( wssize, objsiz, gfnorm )
               if (uncon) then
                  unitgz = gz1nrm .le. dinky
               end if
            end if

            if (lsdbg  .and.  ilsdbg(1) .gt. 0)
     $         write (nout, 9000) unitgz, irefn, gz1nrm, dinky

*           If the projected gradient  Z'g  is small and Rz is of full
*           rank, X is a minimum on the working set.  An additional
*           refinement step is allowed to take care of an inaccurate
*           value of dinky.

            statpt = .not. singlr  .and.  gz1nrm .le. dinky
     $                             .or.   irefn  .gt. mrefn

            if (.not. statpt) then
*              ---------------------------------------------------------
*              Compute a search direction.
*              ---------------------------------------------------------
               prnt  = .true.

               error = iter .ge. itmax
               if (.not. error) then

                  irefn = irefn + 1
                  iter  = iter  + 1

                  if (iter .eq. idbg) then
                     lsdbg  = .true.
                     cmdbg  =  lsdbg
                     msglvl =  msgsvd
                  end if

                  call lsgetp( linobj, singlr, unitgz, unitQ,
     $                         n, nclin, nfree,
     $                         ldA, ldzy, ldR, nrank, numinf, nrz,
     $                         istate, kx, ctp, pnorm,
     $                         A, w(lap), w(lres), w(lhz), w(lpx),
     $                         w(lgq), w(lcq), R, w(lzy), w(lwrk) )

*                 ------------------------------------------------------
*                 Find the constraint we bump into along p.
*                 Update x and Ax if the step alfa is nonzero.
*                 ------------------------------------------------------
*                 alfhit is initialized to bigalf.  If it remains
*                 that way after the call to cmalf, it will be
*                 regarded as infinite.

                  bigalf = ddiv  ( bigdx, pnorm, overfl )

                  call cmalf ( firstv, hitlow,
     $                         istate, inform, jadd, n, ldA,
     $                         nclin, nctotl, numinf,
     $                         alfhit, palfa, atphit,
     $                         bigalf, bigbnd, pnorm,
     $                         w(lanorm), w(lap), ax,
     $                         bl, bu, featol, w(lpx), x )

*                 If  Rz1  is nonsingular,  alfa = 1.0  will be the
*                 step to the least-squares minimizer on the
*                 current subspace. If the unit step does not violate
*                 the nearest constraint by more than featol,  the
*                 constraint is not added to the working set.

                  hitcon = singlr  .or.  palfa  .le. one
                  uncon  = .not. hitcon

                  if (hitcon) then
                     alfa = alfhit
                  else
                     jadd   = 0
                     alfa   = one
                  end if

*                 Check for an unbounded solution or negligible step.

                  unbndd =  alfa .ge. bigalf
                  stall  = abs( alfa*pnorm ) .le. epspt9*xnorm
                  if (stall) then
                     nstall = nstall + 1
                     cyclin = nstall .gt. mstall
                  else
                     nstall = 0
                  end if

                  error = unbndd  .or.  cyclin
                  if (.not.  error) then
*                    ---------------------------------------------------
*                    Set x = x + alfa*p.  Update Ax, gq, res and ctx.
*                    ---------------------------------------------------
                     if (alfa .ne. zero)
     $                  call lsmove( hitcon, hitlow, linobj, unitgz,
     $                               nclin, nrank, nrz,
     $                               n, ldR, jadd, numinf,
     $                               alfa, ctp, ctx, xnorm,
     $                               w(lap), ax, bl, bu, w(lgq),
     $                               w(lhz), w(lpx), w(lres),
     $                               R, x, w(lwrk) )

                     if (hitcon) then
*                       ------------------------------------------------
*                       Add a constraint to the working set.
*                       Update the TQ factors of the working set.
*                       Use p as temporary work space.
*                       ------------------------------------------------
*                       Update  istate.

                        if (bl(jadd) .eq. bu(jadd)) then
                           istate(jadd) = 3
                        else if (hitlow) then
                           istate(jadd) = 1
                        else
                           istate(jadd) = 2
                        end if
                        iadd = jadd - n
                        if (jadd .le. n) then

                           do 510, ifix = 1, nfree
                              if (kx(ifix) .eq. jadd) go to 520
  510                      continue
  520                   end if

                        call lsadd ( unitQ,
     $                               inform, ifix, iadd, jadd,
     $                               nactiv, nz, nfree, nrank, nres,ngq,
     $                               n, ldA, ldzy, ldR, ldT,
     $                               kx, condmx,
     $                               A, R, w(lt), w(lres),w(lgq),w(lzy),
     $                               w(lwrk), w(lrlam), w(lpx) )

                        nrz    = nrz - 1
                        nz     = nz  - 1

                        if (jadd .le. n) then

*                          A simple bound has been added.

                           nfree  = nfree  - 1
                        else

*                          A general constraint has been added.

                           nactiv = nactiv + 1
                           kactiv(nactiv) = iadd
                        end if
                        irefn  = 0
                     end if
*                    ---------------------------------------------------
*                    Check the feasibility of constraints with non-
*                    negative istate values.  If some violations have
*                    occurred.  Refine the current x and set inform so
*                    that feasibility is checked in lsgset.
*                    ---------------------------------------------------
                     call lsfeas( n, nclin, istate,
     $                            bigbnd, cnorm, err1, jmax1, nviol,
     $                            ax, bl, bu, featol, x, w(lwrk) )

                     if (err1 .gt. featol(jmax1)) then
                        call lssetx( linobj, rowerr, unitQ,
     $                               nclin, nactiv, nfree, nrank, nz,
     $                               n, nctotl, ldzy, ldA, ldR, ldT,
     $                               istate, kactiv, kx,
     $                               jmax1, err2, ctx, xnorm,
     $                               A, ax, bl, bu, w(lcq),
     $                               w(lres), w(lres0), featol, r,
     $                               w(lt), x, w(lzy), w(lpx), w(lwrk) )

                        if (lsdbg  .and.  ilsdbg(1) .gt. 0)
     $                     write (nout, 2100) err1, err2
                        if (rowerr)       write (nout, 2200)
                        uncon  =   .false.
                        irefn  =   0
                        numinf = - 1
                     end if
                     needfg = alfa .ne. zero
                  end if
               end if
            end if
*        until      statpt  .or.  error
         if (.not. (statpt  .or.  error) ) go to 100

*        ===============================================================
*        Try and find the index jdel of a constraint to drop from
*        the working set.
*        ===============================================================
         jdel   = 0

         if (numinf .eq. 0  .and.  prbtyp .eq. 'FP') then
            if (n .gt. nz)
     $         call dload ( n-nz, (zero), w(lrlam), 1 )
            jtiny  = 0
            jsmlst = 0
            jbigst = 0
         else

            call lsmuls( prbtyp,
     $                   msglvl, n, nactiv, nfree,
     $                   ldA, ldT, numinf, nz, nrz,
     $                   istate, kactiv, kx, dinky,
     $                   jsmlst, ksmlst, jinf, jtiny,
     $                   jbigst, kbigst, trulam,
     $                   A, w(lanorm), w(lgq), w(lrlam),
     $                   w(lt), w(lwtinf) )
         end if

         if (.not. error) then
            if (     jsmlst .gt. 0) then

*              LSMULS found a regular constraint with multiplier less
*              than (-dinky).

               jdel   = jsmlst
               kdel   = ksmlst
               isdel  = istate(jdel)
               istate(jdel) = 0

            else if (jsmlst .lt. 0) then

               jdel   = jsmlst

            else if (numinf .gt. 0  .and.  jbigst .gt. 0) then

*              No feasible point exists for the constraints but the
*              sum of the constraint violations may be reduced by
*              moving off constraints with multipliers greater than 1.

               jdel   = jbigst
               kdel   = kbigst
               isdel  = istate(jdel)
               if (trulam .le. zero) is = - 1
               if (trulam .gt. zero) is = - 2
               istate(jdel) = is
               firstv = .true.
               numinf = numinf + 1
            end if

            if      (jdel .ne. 0  .and.  singlr) then

*              Cannot delete a constraint when Rz is singular.
*              Probably a weak minimum.

               jdel = 0
            else if (jdel .ne. 0               ) then

*              Constraint jdel has been deleted.
*              Update the matrix factorizations.

               call lsdel ( unitQ,
     $                      n, nactiv, nfree, nres, ngq, nz, nrz,
     $                      ldA, ldzy, ldR, ldT, nrank,
     $                      jdel, kdel, kactiv, kx,
     $                      A, w(lres), R, w(lt), w(lgq), w(lzy),
     $                      w(lwrk), w(lpx), w(lrlam) )
            end if
         end if

         irefn  =  0
         convrg =  jdel .eq. 0
         prnt   = .false.
         uncon  = .false.
         needfg = .false.

*     until       convrg  .or.  error
      if (.not.  (convrg  .or.  error)) go to 100

*  .........................End of main loop............................
      weak = jtiny .gt. 0  .or.  singlr

      if (error) then
         if (unbndd) then
            inform = 2
            if (numinf .gt. 0) inform = 3
         else if (iter .ge. itmax) then
            inform = 4
         else if (cyclin) then
            inform = 5
         end if
      else if (convrg) then
         inform = 0
         if (numinf .gt. 0) then
            inform = 3
         else if (prbtyp .ne. 'FP'  .and.  weak) then
            inform = 1
         end if
      end if

*     ------------------------------------------------------------------
*     Set   clamda.  Print the full solution.
*     ------------------------------------------------------------------
      msglvl = msgsvd
      if (msglvl .gt. 0) write (nout, 2000) prbtyp, iter, inform

      call cmprt ( msglvl, nfree, ldA,
     $             n, nclin, ncnln, nctotl, bigbnd,
     $             named, names, lennam,
     $             nactiv, istate, kactiv, kx,
     $             A, bl, bu, x, clamda, w(lrlam), x )

      return

 2000 format(/ ' Exit from ', a2, ' problem after ', i4, ' iterations.',
     $         '  inform =', i3 )
 2100 format(  ' XXX  Iterative refinement.  Maximum errors before and',
     $         ' after refinement are ',  1p, 2e14.2 )
 2200 format(  ' XXX  Warning.  Cannot satisfy the constraints to the',
     $         ' accuracy requested.')
 9000 format(/ ' //lscore//  unitgz irefn     gz1nrm      dinky'
     $       / ' //lscore//  ', l6, i6, 1p, 2e11.2 )
 9100 format(/ ' //lscore//  singlr   abs(rzz1)      drzmax      drzmin'
     $       / ' //lscore//  ', l6,     1p, 3e12.4 )

*     End of  LSCORE.

      end                         
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSCRSH( COLD, VERTEX,
     $                   NCLIN, NCTOTL, NACTIV, NARTIF,
     $                   NFREE, N, LDA,
     $                   ISTATE, KACTIV,
     $                   BIGBND, TOLACT,
     $                   A, AX, BL, BU, X, WX, WORK )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            COLD, VERTEX
      INTEGER            ISTATE(NCTOTL), KACTIV(N)
      DOUBLE PRECISION   A(LDA,*), AX(*), BL(NCTOTL), BU(NCTOTL),
     $                   X(N), WX(N), WORK(N)

************************************************************************
*     LSCRSH  computes the quantities  ISTATE (optionally), KACTIV,
*     NACTIV, NZ and NFREE  associated with the working set at X.
*     The computation depends upon the value of the input parameter
*     COLD,  as follows...
*
*     COLD = TRUE.  An initial working set will be selected. First,
*                   nearly-satisfied or violated bounds are added.
*                   Next,  general linear constraints are added that
*                   have small residuals.
*
*     COLD = FALSE. The quantities KACTIV, NACTIV, NZ and NFREE are
*                   computed from ISTATE,  specified by the user.
*
*     Values of ISTATE(j)....
*
*        - 2         - 1         0           1          2         3
*     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of LSCRSH dated 21-Nov-1990.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
      COMMON    /SOL1CM/ NOUT

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      EXTERNAL           DDOT
      INTRINSIC          ABS, MIN
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      FLMAX  =   WMACH(7)
      BIGLOW = - BIGBND
      BIGUPP =   BIGBND
      CALL DCOPY ( N, X, 1, WX, 1 )

      IF (LSDBG) THEN
         IF (ILSDBG(1) .GT. 0)
     $      WRITE (NOUT, 1000) COLD, NCLIN, NCTOTL
         IF (ILSDBG(2) .GT. 0)
     $      WRITE (NOUT, 1100) (WX(J), J = 1, N)
      END IF

      NFIXED = 0
      NACTIV = 0
      NARTIF = 0

*     If a cold start is being made, initialize  ISTATE.
*     If  BL(j) = BU(j),  set  ISTATE(j)=3  for all variables and linear
*     constraints.

      IF (COLD) THEN
         DO 100, J = 1, NCTOTL
            ISTATE(J) = 0
            IF (BL(J) .EQ. BU(J)) ISTATE(J) = 3
  100    CONTINUE
      ELSE
         DO 110, J = 1, NCTOTL
            IF (ISTATE(J) .GT. 3  .OR.  ISTATE(J) .LT. 0) ISTATE(J) = 0
  110    CONTINUE
      END IF

*     Initialize NFIXED, NFREE and KACTIV.
*     Ensure that the number of bounds and general constraints in the
*     working set does not exceed N.

      DO 200, J = 1, NCTOTL
         IF (NFIXED + NACTIV .EQ. N) ISTATE(J) = 0
         IF (ISTATE(J) .GT. 0) THEN
            IF (J .LE. N) THEN
               NFIXED = NFIXED + 1
               IF (ISTATE(J) .EQ. 1) WX(J) = BL(J)
               IF (ISTATE(J) .GE. 2) WX(J) = BU(J)
            ELSE
               NACTIV = NACTIV + 1
               KACTIV(NACTIV) = J - N
            END IF
         END IF
  200 CONTINUE

*     ------------------------------------------------------------------
*     If a cold start is required,  attempt to add as many
*     constraints as possible to the working set.
*     ------------------------------------------------------------------
      IF (COLD) THEN

*        See if any bounds are violated or nearly satisfied.
*        If so,  add these bounds to the working set and set the
*        variables exactly on their bounds.

         J = N
*+       WHILE (J .GE. 1  .AND.  NFIXED + NACTIV .LT. N) DO
  300    IF    (J .GE. 1  .AND.  NFIXED + NACTIV .LT. N) THEN
            IF (ISTATE(J) .EQ. 0) THEN
               B1     = BL(J)
               B2     = BU(J)
               IS     = 0
               IF (B1 .GT. BIGLOW) THEN
                  IF (WX(J) - B1 .LE. (ONE + ABS( B1 ))*TOLACT) IS = 1
               END IF
               IF (B2 .LT. BIGUPP) THEN
                  IF (B2 - WX(J) .LE. (ONE + ABS( B2 ))*TOLACT) IS = 2
               END IF
               IF (IS .GT. 0) THEN
                  ISTATE(J) = IS
                  IF (IS .EQ. 1) WX(J) = B1
                  IF (IS .EQ. 2) WX(J) = B2
                  NFIXED = NFIXED + 1
               END IF
            END IF
            J = J - 1
            GO TO 300
*+       END WHILE
         END IF

*        ---------------------------------------------------------------
*        The following loop finds the linear constraint (if any) with
*        smallest residual less than or equal to TOLACT  and adds it
*        to the working set.  This is repeated until the working set
*        is complete or all the remaining residuals are too large.
*        ---------------------------------------------------------------
*        First, compute the residuals for all the constraints not in the
*        working set.

         IF (NCLIN .GT. 0  .AND.  NACTIV+NFIXED .LT. N) THEN
            DO 410, I = 1, NCLIN
               IF (ISTATE(N+I) .LE. 0)
     $         AX(I) = DDOT  (N, A(I,1), LDA, WX, 1 )
  410       CONTINUE

            IS     = 1
            TOOBIG = TOLACT + TOLACT

*+          WHILE (IS .GT. 0  .AND.  NFIXED + NACTIV .LT. N) DO
  500       IF    (IS .GT. 0  .AND.  NFIXED + NACTIV .LT. N) THEN
               IS     = 0
               RESMIN = TOLACT

               DO 520, I = 1, NCLIN
                  J      = N + I
                  IF (ISTATE(J) .EQ. 0) THEN
                     B1     = BL(J)
                     B2     = BU(J)
                     RESL   = TOOBIG
                     RESU   = TOOBIG
                     IF (B1 .GT. BIGLOW)
     $                  RESL  = ABS( AX(I) - B1 ) / (ONE + ABS( B1 ))
                     IF (B2 .LT. BIGUPP)
     $                  RESU  = ABS( AX(I) - B2 ) / (ONE + ABS( B2 ))
                     RESIDL   = MIN( RESL, RESU )
                     IF(RESIDL .LT. RESMIN) THEN
                        RESMIN = RESIDL
                        IMIN   = I
                        IS     = 1
                        IF (RESL .GT. RESU) IS = 2
                     END IF
                  END IF
  520          CONTINUE

               IF (IS .GT. 0) THEN
                  NACTIV = NACTIV + 1
                  KACTIV(NACTIV) = IMIN
                  J         = N + IMIN
                  ISTATE(J) = IS
               END IF
               GO TO 500
*+          END WHILE
            END IF
         END IF
      END IF
            
      IF (VERTEX  .AND.  NACTIV+NFIXED .LT. N) THEN
*        ---------------------------------------------------------------
*        Find an initial vertex by temporarily fixing some variables.
*        Infeasible variables are moved inside their bounds.
*        ---------------------------------------------------------------
*        Compute lengths of columns of selected linear constraints
*        (just the ones corresponding to variables eligible to be
*        temporarily fixed).        

         DO 630, J = 1, N
            IF (ISTATE(J) .EQ. 0) THEN
               COLSIZ = ZERO
               DO 620, K = 1, NCLIN
                  IF (ISTATE(N+K) .GT. 0)
     $            COLSIZ = COLSIZ + ABS( A(K,J) )
  620          CONTINUE
               WORK(J) = COLSIZ
            END IF
  630    CONTINUE
         
*        Find the  NARTIF  smallest such columns.
*        This is an expensive loop.  Later we can replace it by a
*        4-pass process (say), accepting the first col that is within
*        T  of  COLMIN, where  T = 0.0, 0.001, 0.01, 0.1 (say).
*        (This comment written in 1980).
         
*+       WHILE (NFIXED + NACTIV .LT. N) DO
  640    IF    (NFIXED + NACTIV .LT. N) THEN
            COLMIN = FLMAX
            DO 650, J = 1, N
               IF (ISTATE(J) .EQ. 0) THEN
                  IF (NCLIN .EQ. 0) GO TO 660
                  COLSIZ = WORK(J)
                  IF (COLMIN .GT. COLSIZ) THEN
                     COLMIN = COLSIZ
                     JMIN   = J
                  END IF
               END IF
  650       CONTINUE
            J      = JMIN
  660       ISTATE(J) = 4
            NARTIF = NARTIF + 1
            NFIXED = NFIXED + 1
            B1     = BL(J)
            B2     = BU(J)
         
            IF (B1 .GT. BIGLOW) THEN
               IF (X(J) .LT. B1) X(J) = B1
            END IF
      
            IF (B2 .LT. BIGUPP) THEN
               IF (X(J) .GT. B2) X(J) = B2
            END IF
            GO TO 640
*+       END WHILE
         END IF
      END IF
      
      NFREE = N - NFIXED

      IF (LSDBG) THEN
         IF (ILSDBG(1) .GT. 0)
     $       WRITE (NOUT, 1300) NFIXED, NACTIV, NARTIF
         IF (ILSDBG(2) .GT. 0)
     $       WRITE (NOUT, 1200) (WX(J), J = 1, N)
      END IF

      RETURN

 1000 FORMAT(/ ' //LSCRSH// COLD NCLIN NCTOTL'
     $       / ' //LSCRSH// ', L4, I6, I7 )
 1100 FORMAT(/ ' //LSCRSH// Variables before crash... '/ (5G12.3))
 1200 FORMAT(/ ' //LSCRSH// Variables after  crash... '/ (5G12.3))
 1300 FORMAT(/ ' //LSCRSH// Working set selected ...             '
     $       / ' //LSCRSH// NFIXED NACTIV NARTIF      '
     $       / ' //LSCRSH// ', I6, 2I7 )

*     End of  LSCRSH.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsdel ( unitQ,
     $                   n, nactiv, nfree, nres, ngq, nz, nrz,
     $                   ldA, ldzy, ldR, ldT, nrank,
     $                   jdel, kdel, kactiv, kx,
     $                   A, res, R, t, gq, zy,
     $                   work, c, s )

      implicit           double precision(a-h,o-z)
      logical            unitQ
      integer            kactiv(n), kx(n)
      double precision   a(ldA,*), res(n,*), r(ldR,*), T(ldT,*),
     $                   gq(n,*), zy(ldzy,*)
      double precision   work(n), c(n), s(n)

C***********************************************************************
C     LSDEL   updates the least-squares factor R and the factorization
C     A(free) (Z Y) = (0 T) when a regular, temporary or artificial
C     constraint is deleted from the working set.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     Level-2 matrix routines added 25-Apr-1988.
C     This version of LSDEL dated  25-Aug-1991.
C***********************************************************************
      common    /sol1cm/ nout
      common    /sol5cm/ asize, dtmax, dtmin

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg
      
      intrinsic          max   , min
      external           idamax
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      if (jdel .gt. 0) then
*        ---------------------------------------------------------------
*        Regular constraint or temporary bound deleted.
*        ---------------------------------------------------------------

         if (jdel .le. n) then

*           Case 1.  A simple bound has been deleted.
*           =======  Columns nfree+1 and ir of r must be swapped.

            ir     = nz    + kdel
            if (lsdbg  .and.  ilsdbg(1) .gt. 0)
     $         write (nout, 1100) nactiv, nz, nfree, ir, jdel, unitQ

            itdel  = 1
            nfree  = nfree + 1

            if (nfree .lt. ir) then
               kx(ir)    = kx(nfree)
               kx(nfree) = jdel
               if (nrank .gt. 0)
     $            call cmrswp( n, nres, nrank, ldR, nfree, ir,
     $                         R, res, c, s )
               call dswap ( ngq, gq(nfree,1), n, gq(ir,1), n )
            end if

            if (.not. unitQ) then

*              Copy the incoming column of  A(free)  into the end of T.

               do 130, ka = 1, nactiv
                  i = kactiv(ka)
                  T(ka,nfree) = a(i,jdel)
  130          continue

*              Expand Q by adding a unit row and column.

               if (nfree .gt. 1) then
                  call dload ( nfree-1, zero, zy(nfree,1), ldzy )
                  call dload ( nfree-1, zero, zy(1,nfree), 1  )
               end if
               zy(nfree,nfree) = one
            end if
         else        

*           Case 2.  A general constraint has been deleted.
*           =======

            if (lsdbg  .and.  ilsdbg(1) .gt. 0)
     $         write (nout, 1200) nactiv, nz, nfree, kdel, jdel, unitQ

            itdel  = kdel
            nactiv = nactiv - 1

*           Delete row  kdel  of T and move up the ones below it.
*           T becomes reverse lower Hessenberg.

            do 220, i = kdel, nactiv
               kactiv(i) = kactiv(i+1)
               ld        = nfree - i
               call dcopy ( i+1, T(i+1,ld), ldT, T(i,ld), ldT )
  220       continue
         end if

         nz    = nz     + 1

         if (nactiv .eq. 0) then
            dtmax = one
            dtmin = one
         else
*           ------------------------------------------------------------
*           Restore the nactiv by (nactiv+1) reverse-Hessenberg matrix 
*           T  to reverse-triangular form.  The last  nactiv  super-
*           diagonal elements are removed using a backward sweep of
*           plane rotations.  The rotation for the singleton in the 
*           first column is generated separately.
*           ------------------------------------------------------------
            nsup   = nactiv - itdel + 1
                                              
            if (nsup .gt. 0) then
               npiv   = nfree  - itdel + 1

               if (nsup .gt. 1) then
                  call dcopy ( nsup-1, T(nactiv-1,nz+1), ldT-1, 
     $                         s(nz+1), 1)
                  call f06qzf( 'Remove', nactiv, 1, nsup, 
     $                         c(nz+1), s(nz+1), T(1,nz+1), ldT )
               end if

               call f06baf( T(nactiv,nz+1), T(nactiv,nz), cs, sn )
               T(nactiv,nz) = zero
               s(nz)   = - sn
               c(nz)   =   cs
         
               call f06qxf( 'Right', 'Variable', 'Backwards', 
     $                      nfree, nfree, nz, npiv, c, s, zy, ldzy )
               call f06qxf( 'Left ', 'Variable', 'Backwards', 
     $                      npiv , ngq  , nz, npiv, c, s, gq, n    )
            
               nt = min( nrank, npiv )
               
               if (nt .lt. npiv  .and.  nt .gt. 0) then
               
*                 R is upper trapezoidal, pretend R is (nt x n) and 
*                 apply the rotations in columns  max(nt,nz)  thru npiv.
               
                  call f06qxf( 'Right', 'Variable', 'Backwards', nt, n,
     $                         max(nt,nz), npiv, c, s, R, ldR )
               end if
               
*              Apply the column transformations to the triangular part 
*              of  R.  The arrays  c  and  s  containing the column
*              rotations are overwritten by the row rotations that
*              restore  R  to upper-triangular form.
               
               if (nz .lt. nt) then
                  call f06qtf( 'Right', nt, nz, nt, c, s, R, ldR )
               end if
               
*              Apply the row rotations to the remaining rows of R.
               
               call f06qxf( 'Left', 'Variable', 'Backwards', nt, n-nt,
     $                      nz, nt, c, s, r(1,nt+1), ldR )
               
               if (nres .gt. 0)
     $            call f06qxf( 'Left', 'Variable', 'Backwards', 
     $                         nt, nres, nz, nt, c, s, res, n )
            end if
            
            call dcond ( nactiv, T(nactiv,nz+1), ldT-1, dtmax, dtmin )
         end if
      end if

      nrz1 = nrz + 1

      if (nz .gt. nrz) then
         if (jdel .gt. 0) then
            jart =   nrz1 - 1 + idamax( nz-nrz1+1, gq(nrz1,1), 1 )
         else
            jart = - jdel
         end if

         if (lsdbg  .and.  ilsdbg(1) .gt. 0)
     $      write( nout, 1000 ) nz, nrz1, jart

         if (jart .gt. nrz1) then

*           Swap columns NRZ1 and JART of R.

            if (unitQ) then
               k        = kx(nrz1)
               kx(nrz1)  = kx(jart)
               kx(jart) = k
            else
               call dswap ( nfree, zy(1,nrz1), 1, zy(1,jart), 1 )
            end if

            call dswap ( ngq, gq(nrz1,1), n, gq(jart,1), n )
            if (nrank .gt. 0)
     $         call cmrswp( n, nres, nrank, ldR, nrz1, jart,
     $                      R, res, c, s )
         end if
      end if

      nrz = nrz1

      return

 1000 format(/ ' //lsdel //  Artificial constraint deleted.      '
     $       / ' //lsdel //      nz   nrz   jart                 '
     $       / ' //lsdel //  ', 3i6 )
 1100 format(/ ' //lsdel //  Simple bound deleted.               '
     $       / ' //lsdel //  nactiv    nz nfree    ir  jdel unitQ'
     $       / ' //lsdel //  ', 5i6, l6 )
 1200 format(/ ' //lsdel //  General constraint deleted.         '
     $       / ' //lsdel //  nactiv    nz nfree  kdel  jdel unitQ'
     $       / ' //lsdel //  ', 5i6, l6 )

*     End of  LSDEL .

      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSDFLT( M, N, NCLIN, TITLE )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)

      CHARACTER*(*)      TITLE

************************************************************************
*  LSDFLT  loads the default values of parameters not set by the user.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original Fortran 77 version written 17-September-1985.
*  This version of LSDFLT dated   9-September-1986.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9

      LOGICAL            CMDBG, LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
      COMMON    /CMDEBG/ ICMDBG(LDBG), CMDBG

      LOGICAL            NEWOPT
      COMMON    /SOL3LS/ NEWOPT
      SAVE      /SOL3LS/

*-----------------------------------------------------------------------
      PARAMETER         (MXPARM = 30)
      INTEGER            IPRMLS(MXPARM), IPSVLS
      DOUBLE PRECISION   RPRMLS(MXPARM), RPSVLS

      COMMON    /LSPAR1/ IPSVLS(MXPARM),
     $                   IDBGLS, ITMAX1, ITMAX2, LCRASH, LDBGLS, LPROB ,
     $                   MSGLS , NN    , NNCLIN, NPROB , IPADLS(20)

      COMMON    /LSPAR2/ RPSVLS(MXPARM),
     $                   BIGBND, BIGDX , BNDLOW, BNDUPP, TOLACT, TOLFEA,
     $                   TOLRNK, RPADLS(23)

      EQUIVALENCE       (IPRMLS(1), IDBGLS), (RPRMLS(1), BIGBND)

      SAVE      /LSPAR1/, /LSPAR2/
*-----------------------------------------------------------------------
      EQUIVALENCE   (MSGLS , MSGLVL), (IDBGLS, IDBG), (LDBGLS, MSGDBG)

      LOGICAL            CDEFND
      CHARACTER*4        ICRSH(0:2)
      CHARACTER*3        LSTYPE(1:10)
      CHARACTER*16       KEY
      INTRINSIC          LEN    ,  MAX   , MOD
      PARAMETER        ( ZERO   =  0.0D+0, TEN    = 10.0D+0)
      PARAMETER        ( RDUMMY = -11111., IDUMMY = -11111 )
      PARAMETER        ( GIGANT = 1.0D+10*.99999           )
      PARAMETER        ( WRKTOL = 1.0D-2                   )
      DATA               ICRSH(0), ICRSH(1), ICRSH(2)
     $                 /'COLD'   ,'WARM'   ,'HOT '   /
      DATA               LSTYPE(1), LSTYPE(2)
     $                 /' FP'     ,' LP'     /
      DATA               LSTYPE(3), LSTYPE(4), LSTYPE(5), LSTYPE(6)
     $                 /'QP1'     ,'QP2'     ,'QP3'     ,'QP4'     /
      DATA               LSTYPE(7), LSTYPE(8), LSTYPE(9), LSTYPE(10)
     $                 /'LS1'     ,'LS2'     ,'LS3'     ,'LS4'     /

      EPSMCH = WMACH( 3)

*     Make a dummy call to LSKEY to ensure that the defaults are set.

      CALL LSKEY ( NOUT, '*', KEY )
      NEWOPT = .TRUE.

*     Save the optional parameters set by the user.  The values in
*     RPRMLS and IPRMLS may be changed to their default values.

      CALL ICOPY ( MXPARM, IPRMLS, 1, IPSVLS, 1 )
      CALL DCOPY ( MXPARM, RPRMLS, 1, RPSVLS, 1 )

      IF (       LPROB  .LT. 0      )  LPROB   = 7
                                       CDEFND  = LPROB .EQ. 2*(LPROB/2)
      IF (       LCRASH .LT. 0
     $    .OR.   LCRASH .GT. 2      )  LCRASH  = 0
      IF (       ITMAX1 .LT. 0      )  ITMAX1  = MAX(50, 5*(N+NCLIN))
      IF (       ITMAX2 .LT. 0      )  ITMAX2  = MAX(50, 5*(N+NCLIN))
      IF (       MSGLVL .EQ. IDUMMY )  MSGLVL  = 10
      IF (       IDBG   .LT. 0
     $    .OR.   IDBG   .GT. ITMAX1 + ITMAX2
     $                              )  IDBG    = 0
      IF (       MSGDBG .LT. 0      )  MSGDBG  = 0
      IF (       MSGDBG .EQ. 0      )  IDBG    = ITMAX1 + ITMAX2 + 1
      IF (       TOLACT .LT. ZERO   )  TOLACT  = WRKTOL
      IF (       TOLFEA .EQ. RDUMMY
     $    .OR.  (TOLFEA .GE. ZERO
     $    .AND.  TOLFEA .LT. EPSMCH))  TOLFEA  = EPSPT5
      IF (       TOLRNK .LE. ZERO
     $    .AND.  CDEFND             )  TOLRNK  = EPSPT5
      IF (       TOLRNK .LE. ZERO   )  TOLRNK  = TEN*EPSMCH
      IF (       BIGBND .LE. ZERO   )  BIGBND  = GIGANT
      IF (       BIGDX  .LE. ZERO   )  BIGDX   = MAX(GIGANT, BIGBND)

      LSDBG = IDBG .EQ. 0
      CMDBG = LSDBG
      K     = 1
      MSG   = MSGDBG
      DO 200, I = 1, LDBG
         ILSDBG(I) = MOD( MSG/K, 10 )
         ICMDBG(I) = ILSDBG(I)
         K = K*10
  200 CONTINUE

      IF (MSGLVL .GT. 0) THEN

*        Print the title.

         LENT = LEN( TITLE )
         IF (LENT .GT. 0) THEN
            NSPACE = (81 - LENT)/2 + 1
            WRITE (NOUT, '(///// (80A1) )')
     $         (' ', J=1, NSPACE), (TITLE(J:J), J=1,LENT)
            WRITE (NOUT, '(80A1 //)')
     $         (' ', J=1, NSPACE), ('='       , J=1,LENT)
         END IF

         WRITE (NOUT, 2000)
         WRITE (NOUT, 2100) LSTYPE(LPROB),
     $                      NCLIN , TOLFEA, ICRSH(LCRASH),
     $                      N     , BIGBND, TOLACT,
     $                      M     , BIGDX , TOLRNK
         WRITE (NOUT, 2200) EPSMCH, ITMAX1, MSGLVL,
     $                              ITMAX2
      END IF

      RETURN

 2000 FORMAT(
     $//' Parameters'
     $/ ' ----------' )
 2100 FORMAT(
     $/ ' Problem type...........', 7X, A3
     $/ ' Linear constraints.....', I10,     6X,
     $  ' Feasibility tolerance..', 1PE10.2, 6X,
     $  1X, A4, ' start.............'
     $/ ' Variables..............', I10,     6X,
     $  ' Infinite bound size....', 1PE10.2, 6X,
     $  ' Crash tolerance........', 1PE10.2
     $/ ' Objective matrix rows..', I10,     6X,
     $  ' Infinite step size.....', 1PE10.2, 6X,
     $  ' Rank tolerance.........', 1PE10.2 )
 2200 FORMAT(
     $/ ' EPS (machine precision)', 1PE10.2, 6X,
     $  ' Feasibility phase itns.', I10, 6X,
     $  ' Print level............', I10
     $/ 40X,
     $  ' Optimality  phase itns.', I10 )

*     End of  LSDFLT.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSFEAS( N, NCLIN, ISTATE,
     $                   BIGBND, CVNORM, ERRMAX, JMAX, NVIOL,
     $                   AX, BL, BU, FEATOL, X, WORK )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER            ISTATE(N+NCLIN)
      DOUBLE PRECISION   AX(*), BL(N+NCLIN), BU(N+NCLIN)
      DOUBLE PRECISION   FEATOL(N+NCLIN), X(N)
      DOUBLE PRECISION   WORK(N+NCLIN)

************************************************************************
*  LSFEAS  computes the following...
*  (1)  The number of constraints that are violated by more
*       than  FEATOL  and the 2-norm of the constraint violations.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original version      April    1984.
*  This version of  LSFEAS  dated  17-October-1985.
************************************************************************
      COMMON    /SOL1CM/ NOUT

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      EXTERNAL           IDAMAX, DNRM2
      INTRINSIC          ABS
      PARAMETER        ( ZERO = 0.0D+0 )

      BIGLOW = - BIGBND
      BIGUPP =   BIGBND

*     ==================================================================
*     Compute NVIOL,  the number of constraints violated by more than
*     FEATOL,  and CVNORM,  the 2-norm of the constraint violations and
*     residuals of the constraints in the working set.
*     ==================================================================
      NVIOL  = 0

      DO 200, J = 1, N+NCLIN
         FEASJ  = FEATOL(J)
         IS     = ISTATE(J)
         RES    = ZERO

         IF (IS .GE. 0  .AND.  IS .LT. 4) THEN
            IF (J .LE. N) THEN
               CON =  X(J)
            ELSE
               I   = J - N
               CON = AX(I)
            END IF

            TOLJ   = FEASJ

*           Check for constraint violations.

            IF (BL(J) .GT. BIGLOW) THEN
               RES    = BL(J) - CON
               IF (RES .GT.   FEASJ ) NVIOL = NVIOL + 1
               IF (RES .GT.    TOLJ ) GO TO 190
            END IF

            IF (BU(J) .LT. BIGUPP) THEN
               RES    = BU(J) - CON
               IF (RES .LT. (-FEASJ)) NVIOL = NVIOL + 1
               IF (RES .LT.  (-TOLJ)) GO TO 190
            END IF

*           This constraint is satisfied,  but count the residual as a
*           violation if the constraint is in the working set.

            IF (IS .LE. 0) RES = ZERO
            IF (IS .EQ. 1) RES = BL(J) - CON
            IF (IS .GE. 2) RES = BU(J) - CON
            IF (ABS( RES ) .GT. FEASJ) NVIOL = NVIOL + 1
         END IF
  190    WORK(J) = RES
  200 CONTINUE

      JMAX   = IDAMAX( N+NCLIN, WORK, 1 )
      ERRMAX = ABS ( WORK(JMAX) )

      IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $   WRITE (NOUT, 1000) ERRMAX, JMAX

      CVNORM  = DNRM2 ( N+NCLIN, WORK, 1 )

      RETURN

 1000 FORMAT(/ ' //LSFEAS//  The maximum violation is ', 1PE14.2,
     $                     ' in constraint', I5 )

*     End of  LSFEAS.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSFILE( IOPTNS, INFORM )
      INTEGER            IOPTNS, INFORM

************************************************************************
*     LSFILE  reads the options file from unit  IOPTNS  and loads the
*     options into the relevant elements of  IPRMLS  and  RPRMLS.
*
*     If  IOPTNS .lt. 0  or  IOPTNS .gt. 99  then no file is read,
*     otherwise the file associated with unit  IOPTNS  is read.
*
*     Output:
*
*         INFORM = 0  if a complete  OPTIONS  file was found
*                     (starting with  BEGIN  and ending with  END);
*                  1  if  IOPTNS .lt. 0  or  IOPTNS .gt. 99;
*                  2  if  BEGIN  was found, but end-of-file
*                     occurred before  END  was found;
*                  3  if end-of-file occurred before  BEGIN  or
*                     ENDRUN  were found;
*                  4  if  ENDRUN  was found before  BEGIN.
************************************************************************
      LOGICAL             NEWOPT
      COMMON     /SOL3LS/ NEWOPT
      SAVE       /SOL3LS/

      DOUBLE PRECISION    WMACH(15)
      COMMON     /SOLMCH/ WMACH
      SAVE       /SOLMCH/

      EXTERNAL            MCHPAR, LSKEY
      LOGICAL             FIRST
      SAVE                FIRST , NOUT
      DATA                FIRST /.TRUE./

*     If first time in, set NOUT.
*     NEWOPT is true first time into LSFILE or LSOPTN
*     and just after a call to LSSOL.

      IF (FIRST) THEN
         FIRST  = .FALSE.
         NEWOPT = .TRUE.
         CALL MCHPAR()
         NOUT = WMACH(11)
      END IF

      CALL OPFILE( IOPTNS, NOUT, INFORM, LSKEY )

      RETURN

*     End of  LSFILE.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSGETP( LINOBJ, SINGLR, UNITGZ, UNITQ,
     $                   N, NCLIN, NFREE,
     $                   LDA, LDZY, LDR, NRANK, NUMINF, NRZ,
     $                   ISTATE, KX, CTP, PNORM,
     $                   A, AP, RES, HZ, P,
     $                   GQ, CQ, R, ZY, WORK )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            LINOBJ, SINGLR, UNITGZ, UNITQ
      INTEGER            ISTATE(N+NCLIN), KX(N)
      DOUBLE PRECISION   A(LDA,*), AP(*), RES(*), HZ(*), P(N),
     $                   GQ(N), CQ(*), R(LDR,*), ZY(LDZY,*)
      DOUBLE PRECISION   WORK(N)

************************************************************************
*     LSGETP  computes the following quantities for  LSCORE.
*     (1) The vector  (hz1) = (Rz1)(pz1).
*         If X is not yet feasible,  the product is computed directly.
*         If  Rz1 is singular,  hz1  is zero.  Otherwise  hz1  satisfies
*         the equations
*                        Rz1'hz1 = -gz1,
*         where  g  is the total gradient.  If there is no linear term
*         in the objective,  hz1  is set to  dz1  directly.
*     (2) The search direction P (and its 2-norm).  The vector P is
*         defined as  Z*(pz1), where  (pz1)  depends upon whether or
*         not X is feasible and the nonsingularity of  (Rz1).
*         If  NUMINF .GT. 0,  (pz1)  is the steepest-descent direction.
*         Otherwise,  x  is the solution of the  NRZ*NRZ  triangular
*         system   (Rz1)*(pz1) = (hz1).
*     (3) The vector Ap,  where A is the matrix of linear constraints.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     Level 2 Blas added 11-June-1986.
*     This version of LSGETP dated  5-Aug-1987.
************************************************************************
      COMMON    /SOL1CM/ NOUT

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      EXTERNAL           DDOT  , DNRM2
      INTRINSIC          MIN
      PARAMETER        ( ZERO = 0.0D+0, ONE  = 1.0D+0 )

      IF (SINGLR) THEN
*        ---------------------------------------------------------------
*        The triangular factor for the current objective function is
*        singular,  i.e., the objective is linear along the last column
*        of Z1.  This can only occur when UNITGZ is TRUE.
*        ---------------------------------------------------------------
         IF (NRZ .GT. 1) THEN
            CALL DCOPY ( NRZ-1, R(1,NRZ), 1, P, 1 )
            CALL DTRSV ( 'U', 'N', 'N', NRZ-1, R, LDR, P, 1 )
         END IF
         P(NRZ) = - ONE

         GTP = DDOT  ( NRZ, GQ, 1, P, 1 )
         IF (GTP .GT. ZERO) CALL DSCAL ( NRZ, (-ONE), P, 1 )

         IF (NRZ .LE. NRANK) THEN
            IF (NUMINF .EQ. 0) THEN
               IF (UNITGZ) THEN
                  HZ(NRZ) = R(NRZ,NRZ)*P(NRZ)
               ELSE
                  CALL DLOAD ( NRZ, (ZERO), HZ, 1 )
               END IF
            ELSE
               HZ(1)   = R(1,1)*P(1)
            END IF
         END IF
      ELSE
*        ---------------------------------------------------------------
*        The objective is quadratic in the space spanned by Z1.
*        ---------------------------------------------------------------
         IF (LINOBJ) THEN
            IF (UNITGZ) THEN
               IF (NRZ .GT. 1)
     $            CALL DLOAD ( NRZ-1, (ZERO), HZ, 1 )
               HZ(NRZ) = - GQ(NRZ)/R(NRZ,NRZ)
            ELSE
               CALL DCOPY ( NRZ, GQ  , 1, HZ, 1 )
               CALL DSCAL ( NRZ, (-ONE), HZ, 1 )
               CALL DTRSV ( 'U', 'T', 'N', NRZ, R, LDR, HZ, 1 )
            END IF
         ELSE
            CALL DCOPY ( NRZ, RES, 1, HZ, 1 )
         END IF

*        Solve  Rz1*pz1 = hz1.

         CALL DCOPY ( NRZ, HZ, 1, P, 1 )
         CALL DTRSV ( 'U', 'N', 'N', NRZ, R, LDR, P, 1 )
      END IF

*     Compute  p = Z1*pz1  and its norm.

      IF (LINOBJ)
     $   CTP = DDOT  ( NRZ, CQ, 1, P, 1 )
      PNORM  = DNRM2 ( NRZ, P, 1 )

      CALL CMQMUL( 1, N, NRZ, NFREE, LDZY, UNITQ, KX, P, ZY, WORK )

      IF (LSDBG  .AND.  ILSDBG(2) .GT. 0)
     $   WRITE (NOUT, 1000) (P(J), J = 1, N)

*     Compute  Ap.

      IF (NCLIN .GT. 0) THEN
         CALL DGEMV ( 'No transpose', NCLIN, N, ONE, A, LDA,
     $                P, 1, ZERO, AP, 1 )

         IF (LSDBG  .AND.  ILSDBG(2) .GT. 0)
     $      WRITE (NOUT, 1100) (AP(I), I = 1, NCLIN)
      END IF

      RETURN

 1000 FORMAT(/ ' //LSGETP//   P ... ' / (1P, 5E15.5))
 1100 FORMAT(/ ' //LSGETP//  AP ... ' / (1P, 5E15.5))

*     End of  LSGETP.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSGSET( PRBTYP, LINOBJ, SINGLR, UNITGZ, UNITQ,
     $                   N, NCLIN, NFREE,
     $                   LDA, LDZY, LDR, NRANK, NZ, NRZ,
     $                   ISTATE, KX,
     $                   BIGBND, TOLRNK, NUMINF, SUMINF,
     $                   BL, BU, A, RES, FEATOL,
     $                   GQ, CQ, R, X, WTINF, ZY, WRK )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*2        PRBTYP
      LOGICAL            LINOBJ, SINGLR, UNITGZ, UNITQ
      INTEGER            ISTATE(*), KX(N)
      DOUBLE PRECISION   BL(*), BU(*), A(LDA,*),
     $                   RES(*), FEATOL(*)
      DOUBLE PRECISION   GQ(N), CQ(*), R(LDR,*), X(N), WTINF(*),
     $                   ZY(LDZY,*)
      DOUBLE PRECISION   WRK(N)

************************************************************************
*     LSGSET  finds the number and weighted sum of infeasibilities for
*     the bounds and linear constraints.   An appropriate transformed
*     gradient vector is returned in  GQ.
*
*     Positive values of  ISTATE(j)  will not be altered.  These mean
*     the following...
*
*               1             2           3
*           a'x = bl      a'x = bu     bl = bu
*
*     Other values of  ISTATE(j)  will be reset as follows...
*           a'x lt bl     a'x gt bu     a'x free
*              - 2           - 1           0
*
*     If  x  is feasible,  LSGSET computes the vector Q(free)'g(free),
*     where  g  is the gradient of the the sum of squares plus the
*     linear term.  The matrix Q is of the form
*                    ( Q(free)  0       ),
*                    (   0      I(fixed))
*     where  Q(free)  is the orthogonal factor of  A(free)  and  A  is
*     the matrix of constraints in the working set.  The transformed
*     gradients are stored in GQ.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     Level 2 Blas added 11-June-1986.
*     This version of LSGSET dated 28-May-1988.
************************************************************************
      EXTERNAL           DDOT  , IDRANK
      INTRINSIC          ABS   , MAX   , MIN
      PARAMETER        ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
                     
      BIGUPP =   BIGBND
      BIGLOW = - BIGBND

      NUMINF =   0
      SUMINF =   ZERO
      CALL DLOAD ( N, ZERO, GQ, 1 )

      DO 200, J = 1, N+NCLIN
         IF (ISTATE(J) .LE. 0) THEN
            FEASJ  = FEATOL(J)
            IF (J .LE. N) THEN
               CTX = X(J)
            ELSE
               K   = J - N
               CTX = DDOT  ( N, A(K,1), LDA, X, 1 )
            END IF
            ISTATE(J) = 0

*           See if the lower bound is violated.

            IF (BL(J) .GT. BIGLOW) THEN
               S = BL(J) - CTX
               IF (S     .GT. FEASJ ) THEN
                  ISTATE(J) = - 2
                  WEIGHT    = - WTINF(J)
                  GO TO 160
               END IF
            END IF

*           See if the upper bound is violated.

            IF (BU(J) .GE. BIGUPP) GO TO 200
            S = CTX - BU(J)
            IF (S     .LE. FEASJ ) GO TO 200
            ISTATE(J) = - 1
            WEIGHT    =   WTINF(J)

*           Add the infeasibility.

  160       NUMINF = NUMINF + 1
            SUMINF = SUMINF + ABS( WEIGHT ) * S
            IF (J .LE. N) THEN
               GQ(J) = WEIGHT
            ELSE
               CALL DAXPY ( N, WEIGHT, A(K,1), LDA, GQ, 1 )
            END IF
         END IF
  200 CONTINUE

*     ------------------------------------------------------------------
*     Install  GQ,  the transformed gradient.
*     ------------------------------------------------------------------
      SINGLR = .FALSE.
      UNITGZ = .TRUE.

      IF (NUMINF .GT. 0) THEN
         CALL CMQMUL( 6, N, NZ, NFREE, LDZY, UNITQ, KX, GQ, ZY, WRK )
      ELSE IF (NUMINF .EQ. 0  .AND.  PRBTYP .EQ. 'FP') THEN
         CALL DLOAD ( N, ZERO, GQ, 1 )
      ELSE

*        Ready for the Optimality Phase.
*        Set NRZ so that Rz1 is nonsingular.

         IF (NRANK .EQ. 0) THEN
            IF (LINOBJ) THEN
               CALL DCOPY ( N, CQ, 1, GQ, 1 )
            ELSE
               CALL DLOAD ( N, ZERO, GQ, 1 )
            END IF
            NRZ    = 0
         ELSE

*           Compute  GQ = - R' * (transformed residual)

            CALL DCOPY ( NRANK, RES, 1, GQ, 1 )
            CALL DSCAL ( NRANK, (-ONE), GQ, 1 )
            CALL DTRMV ( 'U', 'T', 'N', NRANK, R, LDR, GQ, 1 )
            IF (NRANK .LT. N)
     $         CALL DGEMV( 'T', NRANK, N-NRANK, -ONE,R(1,NRANK+1),LDR,
     $                      RES, 1, ZERO, GQ(NRANK+1), 1 )

            IF (LINOBJ) CALL DAXPY ( N, ONE, CQ, 1, GQ, 1 )
            UNITGZ = .FALSE.

            ROWNRM = DNRM2 ( N, R(1,1), LDR )
            IF (         ROWNRM  .LE.        TOLRNK
     $          .OR. ABS(R(1,1)) .LE. ROWNRM*TOLRNK) THEN
               NRZ = 0
            ELSE
               NRZ = IDRANK( MIN(NRANK, NZ), R, LDR+1, TOLRNK )
            END IF
         END IF
      END IF

      RETURN

*     End of  LSGSET.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSKEY ( NOUT, BUFFER, KEY )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*(*)      BUFFER

************************************************************************
*     LSKEY   decodes the option contained in  BUFFER  in order to set
*     a parameter value in the relevant element of  IPRMLS  or  RPRMLS.
*
*
*     Input:
*
*     NOUT   A unit number for printing error messages.
*            NOUT  must be a valid unit.
*
*     Output:
*
*     KEY    The first keyword contained in BUFFER.
*
*
*     LSKEY  calls OPNUMB and the subprograms
*                 LOOKUP, SCANNR, TOKENS, UPCASE
*     (now called OPLOOK, OPSCAN, OPTOKN, OPUPPR)
*     supplied by Informatics General, Inc., Palo Alto, California.
*
*     Systems Optimization Laboratory, Stanford University.
*     This version dated Jan 22, 1986.
************************************************************************
*-----------------------------------------------------------------------
      PARAMETER         (MXPARM = 30)
      INTEGER            IPRMLS(MXPARM), IPSVLS
      DOUBLE PRECISION   RPRMLS(MXPARM), RPSVLS

      COMMON    /LSPAR1/ IPSVLS(MXPARM),
     $                   IDBGLS, ITMAX1, ITMAX2, LCRASH, LDBGLS, LPROB ,
     $                   MSGLS , NN    , NNCLIN, NPROB , IPADLS(20)

      COMMON    /LSPAR2/ RPSVLS(MXPARM),
     $                   BIGBND, BIGDX , BNDLOW, BNDUPP, TOLACT, TOLFEA,
     $                   TOLRNK, RPADLS(23)

      EQUIVALENCE       (IPRMLS(1), IDBGLS), (RPRMLS(1), BIGBND)

      SAVE      /LSPAR1/, /LSPAR2/
*-----------------------------------------------------------------------

      EXTERNAL           OPNUMB
      LOGICAL            FIRST , MORE  , NUMBER, OPNUMB, SORTED
      SAVE               FIRST

      PARAMETER         (     MAXKEY = 27,  MAXTIE = 10,   MAXTOK = 10,
     $                        MAXTYP = 16)
      CHARACTER*16       KEYS(MAXKEY), TIES(MAXTIE), TOKEN(MAXTOK),
     $                   TYPE(MAXTYP)
      CHARACTER*16       KEY, KEY2, KEY3, VALUE

      PARAMETER         (IDUMMY = -11111,  RDUMMY = -11111.0,
     $                   SORTED = .TRUE.,  ZERO   =  0.0     )

      DATA                FIRST
     $                  /.TRUE./
      DATA   KEYS
     $ / 'BEGIN           ',
     $   'COLD            ', 'CONSTRAINTS     ', 'CRASH           ',
     $   'DEBUG           ', 'DEFAULTS        ', 'END             ',
     $   'FEASIBILITY     ', 'HOT             ', 'INFINITE        ',
     $   'IPRMLS          ', 'ITERATIONS      ', 'ITERS:ITERATIONS',
     $   'ITNS :ITERATIONS', 'LINEAR          ', 'LIST            ',
     $   'LOWER           ', 'NOLIST          ', 'OPTIMALITY      ',
     $   'PRINT           ', 'PROBLEM         ', 'RANK            ',
     $   'RPRMLS          ', 'START           ', 'UPPER           ',
     $   'VARIABLES       ', 'WARM            '/

      DATA   TIES
     $ / 'BOUND           ', 'CONSTRAINTS     ',
     $   'NO              ', 'NO.      :NUMBER', 'NUMBER          ',
     $   'PHASE           ', 'STEP            ',
     $   'TOLERANCE       ', 'TYPE            ', 'YES             '/

      DATA   TYPE
     $ / 'FP              ',
     $   'LEAST       :LS1', 'LINEAR       :LP', 'LP              ',
     $   'LS          :LS1', 'LS1             ', 'LS2             ',
     $   'LS3             ', 'LS4             ', 'LSQ         :LS1',
     $   'QP          :QP2', 'QP1             ', 'QP2             ',
     $   'QP3             ', 'QP4             ', 'QUADRATIC   :QP2'/
*-----------------------------------------------------------------------

      IF (FIRST) THEN
         FIRST  = .FALSE.
         DO 10, I = 1, MXPARM
            IPRMLS(I) = IDUMMY
            RPRMLS(I) = RDUMMY
   10    CONTINUE
      END IF

*     Eliminate comments and empty lines.
*     A '*' appearing anywhere in BUFFER terminates the string.

      I      = INDEX( BUFFER, '*' )
      IF (I .EQ. 0) THEN
         LENBUF = LEN( BUFFER )
      ELSE
         LENBUF = I - 1
      END IF
      IF (LENBUF .LE. 0) THEN
         KEY = '*'
         GO TO 900
      END IF

*     ------------------------------------------------------------------
*     Extract up to MAXTOK tokens from the record.
*     NTOKEN returns how many were actually found.
*     KEY, KEY2, KEY3 are the first tokens if any, otherwise blank.
*     ------------------------------------------------------------------
      NTOKEN = MAXTOK
      CALL OPTOKN( BUFFER(1:LENBUF), NTOKEN, TOKEN )
      KEY    = TOKEN(1)
      KEY2   = TOKEN(2)
      KEY3   = TOKEN(3)

*     Certain keywords require no action.

      IF (KEY .EQ. ' '     .OR.  KEY .EQ. 'BEGIN' ) GO TO 900
      IF (KEY .EQ. 'LIST'  .OR.  KEY .EQ. 'NOLIST') GO TO 900
      IF (KEY .EQ. 'END'                          ) GO TO 900

*     Most keywords will have an associated integer or real value,
*     so look for it no matter what the keyword.

      I      = 1
      NUMBER = .FALSE.

   50 IF (I .LT. NTOKEN  .AND.  .NOT. NUMBER) THEN
         I      = I + 1
         VALUE  = TOKEN(I)
         NUMBER = OPNUMB( VALUE )
         GO TO 50
      END IF

      IF (NUMBER) THEN
         READ (VALUE, '(BN, E16.0)') RVALUE
      ELSE
         RVALUE = ZERO
      END IF

*     Convert the keywords to their most fundamental form
*     (upper case, no abbreviations).
*     SORTED says whether the dictionaries are in alphabetic order.
*     LOCi   says where the keywords are in the dictionaries.
*     LOCi = 0 signals that the keyword wasn't there.

      CALL OPLOOK( MAXKEY, KEYS, SORTED, KEY , LOC1 )
      CALL OPLOOK( MAXTIE, TIES, SORTED, KEY2, LOC2 )

*     ------------------------------------------------------------------
*     Decide what to do about each keyword.
*     The second keyword (if any) might be needed to break ties.
*     Some seemingly redundant testing of MORE is used
*     to avoid compiler limits on the number of consecutive ELSE IFs.
*     ------------------------------------------------------------------
      MORE   = .TRUE.
      IF (MORE) THEN
         MORE   = .FALSE.
         IF (KEY .EQ. 'COLD        ') THEN
            LCRASH = 0
         ELSE IF (KEY .EQ. 'CONSTRAINTS ') THEN
            NNCLIN = RVALUE
         ELSE IF (KEY .EQ. 'CRASH       ') THEN
            TOLACT = RVALUE
         ELSE IF (KEY .EQ. 'DEBUG       ') THEN
            LDBGLS = RVALUE
         ELSE IF (KEY .EQ. 'DEFAULTS    ') THEN
            DO 20, I = 1, MXPARM
               IPRMLS(I) = IDUMMY
               RPRMLS(I) = RDUMMY
   20       CONTINUE
         ELSE IF (KEY .EQ. 'FEASIBILITY ') THEN
              IF (KEY2.EQ. 'PHASE       ') ITMAX1 = RVALUE
              IF (KEY2.EQ. 'TOLERANCE   ') TOLFEA = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE
            MORE   = .TRUE.
         END IF
      END IF

      IF (MORE) THEN
         MORE   = .FALSE.
         IF (KEY .EQ. 'HOT         ') THEN
            LCRASH = 2
         ELSE IF (KEY .EQ. 'INFINITE    ') THEN
              IF (KEY2.EQ. 'BOUND       ') BIGBND = RVALUE * 0.99999
              IF (KEY2.EQ. 'STEP        ') BIGDX  = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'IPRMLS      ') THEN
*           Allow things like  IPRMLS 21 = 100  to set IPRMLS(21) = 100
            IVALUE = RVALUE
            IF (IVALUE .GE. 1  .AND. IVALUE .LE. MXPARM) THEN
               READ (KEY3, '(BN, I16)') IPRMLS(IVALUE)
            ELSE
               WRITE(NOUT, 2400) IVALUE
            END IF
         ELSE IF (KEY .EQ. 'ITERATIONS  ') THEN
            ITMAX2 = RVALUE
         ELSE IF (KEY .EQ. 'LINEAR      ') THEN
            NNCLIN = RVALUE
         ELSE IF (KEY .EQ. 'LOWER       ') THEN
            BNDLOW = RVALUE
         ELSE
            MORE   = .TRUE.
         END IF
      END IF

      IF (MORE) THEN
         MORE   = .FALSE.
         IF      (KEY .EQ. 'OPTIMALITY  ') THEN
            ITMAX2 = RVALUE
         ELSE IF (KEY .EQ. 'PROBLEM     ') THEN
            IF      (KEY2 .EQ. 'NUMBER') THEN
               NPROB  = RVALUE
            ELSE IF (KEY2 .EQ. 'TYPE  ') THEN

*              Recognize     Problem type = LP     etc.

               CALL OPLOOK( MAXTYP, TYPE, SORTED, KEY3, LOC3 )
               IF (KEY3 .EQ. 'FP' ) LPROB = 1
               IF (KEY3 .EQ. 'LP' ) LPROB = 2
               IF (KEY3 .EQ. 'QP1') LPROB = 3
               IF (KEY3 .EQ. 'QP2') LPROB = 4
               IF (KEY3 .EQ. 'QP3') LPROB = 5
               IF (KEY3 .EQ. 'QP4') LPROB = 6
               IF (KEY3 .EQ. 'LS1') LPROB = 7
               IF (KEY3 .EQ. 'LS2') LPROB = 8
               IF (KEY3 .EQ. 'LS3') LPROB = 9
               IF (KEY3 .EQ. 'LS4') LPROB = 10
               IF (LOC3 .EQ.  0  ) WRITE(NOUT, 2330) KEY3
            ELSE
               WRITE(NOUT, 2320) KEY2
            END IF
         ELSE
            MORE   = .TRUE.
         END IF
      END IF

      IF (MORE) THEN
         MORE   = .FALSE.
         IF      (KEY .EQ. 'PRINT       ') THEN
            MSGLS  = RVALUE
         ELSE IF (KEY .EQ. 'RANK        ') THEN
            TOLRNK = RVALUE
         ELSE IF (KEY .EQ. 'RPRMLS      ') THEN
*           Allow things like  RPRMLS 21 = 2  to set RPRMLS(21) = 2.0
            IVALUE = RVALUE
            IF (IVALUE .GE. 1  .AND. IVALUE .LE. MXPARM) THEN
               READ (KEY3, '(BN, E16.0)') RPRMLS(IVALUE)
            ELSE
               WRITE(NOUT, 2400) IVALUE
            END IF
         ELSE IF (KEY .EQ. 'START       ') THEN
            IDBGLS = RVALUE
         ELSE IF (KEY .EQ. 'UPPER       ') THEN
            BNDUPP = RVALUE
         ELSE IF (KEY .EQ. 'VARIABLES   ') THEN
            NN     = RVALUE
         ELSE IF (KEY .EQ. 'WARM        ') THEN
            LCRASH = 1
         ELSE
            WRITE(NOUT, 2300) KEY
         END IF
      END IF

  900 RETURN

 2300 FORMAT(' XXX  Keyword not recognized:         ', A)
 2320 FORMAT(' XXX  Second keyword not recognized:  ', A)
 2330 FORMAT(' XXX  Third  keyword not recognized:  ', A)
 2400 FORMAT(' XXX  The PARM subscript is out of range:', I10)

*     End of LSKEY

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSLOC ( LPROB, N, NCLIN, LITOTL, LWTOTL )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)

************************************************************************
*     LSLOC   allocates the addresses of the work arrays for  LSCORE.
*
*     Note that the arrays  ( GQ, CQ )  and  ( RES, RES0, HZ )  lie in
*     contiguous areas of workspace.
*     RES, RES0 and HZ are not needed for LP.
*     CQ is defined when the objective has an explicit linear term.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written  29-October-1984.
*     This version of LSLOC dated 16-February-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL3CM/ LENNAM, LDT   , NCOLT, LDZY

      PARAMETER        ( LENLS = 20 )
      COMMON    /SOL1LS/ LOCLS(LENLS)

      LOGICAL            LSDBG
      PARAMETER        ( LDBG = 5 )
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      MINIW     = LITOTL + 1
      MINW      = LWTOTL + 1


*     Assign array lengths that depend upon the problem dimensions.

      IF (NCLIN .EQ. 0) THEN
         LENT  = 0
         LENZY = 0
      ELSE
         LENT  = LDT *NCOLT
         LENZY = LDZY*LDZY
      END IF

      LENCQ  = 0
      IF (LPROB .EQ. 2*(LPROB/2)) LENCQ  = N
      LENRES = 0
      IF (LPROB .GT. 2          ) LENRES = N

      LKACTV    = MINIW
      MINIW     = LKACTV + N

      LANORM    = MINW
      LAP       = LANORM + NCLIN
      LPX       = LAP    + NCLIN
      LGQ       = LPX    + N
      LCQ       = LGQ    + N
      LRES      = LCQ    + LENCQ
      LRES0     = LRES   + LENRES
      LHZ       = LRES0  + LENRES
      LRLAM     = LHZ    + LENRES
      LT        = LRLAM  + N
      LZY       = LT     + LENT
      LWTINF    = LZY    + LENZY
      LWRK      = LWTINF + N  + NCLIN
      LFEATL    = LWRK   + N  + NCLIN
      MINW      = LFEATL + N  + NCLIN

      LOCLS( 1) = LKACTV
      LOCLS( 2) = LANORM
      LOCLS( 3) = LAP
      LOCLS( 4) = LPX
      LOCLS( 5) = LRES
      LOCLS( 6) = LRES0
      LOCLS( 7) = LHZ
      LOCLS( 8) = LGQ
      LOCLS( 9) = LCQ
      LOCLS(10) = LRLAM
      LOCLS(11) = LT
      LOCLS(12) = LZY
      LOCLS(13) = LWTINF
      LOCLS(14) = LWRK
      LOCLS(15) = LFEATL

      LITOTL    = MINIW - 1
      LWTOTL    = MINW  - 1

      RETURN

*     End of  LSLOC .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSMOVE( HITCON, HITLOW, LINOBJ, UNITGZ,
     $                   NCLIN, NRANK, NRZ,
     $                   N, LDR, JADD, NUMINF,
     $                   ALFA, CTP, CTX, XNORM,
     $                   AP, AX, BL, BU, GQ, HZ, P, RES,
     $                   R, X, WORK )

      IMPLICIT           DOUBLE PRECISION (A-H,O-Z)
      LOGICAL            HITCON, HITLOW, LINOBJ, UNITGZ
      DOUBLE PRECISION   AP(*), AX(*), BL(*), BU(*), GQ(*), HZ(*),
     $                   P(N), RES(*), R(LDR,*), X(N)
      DOUBLE PRECISION   WORK(*)

************************************************************************
*     LSMOVE  changes X to X + ALFA*P and updates CTX, AX, RES and GQ
*     accordingly.
*
*     If a bound was added to the working set,  move X exactly on to it,
*     except when a negative step was taken (CMALF may have had to move
*     to some other closer constraint.)
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 27-December-1985.
*     Level 2 BLAS added 11-June-1986.
*     This version of LSMOVE dated 11-June-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      EXTERNAL           DDOT  , DNRM2
      INTRINSIC          ABS   , MIN
      PARAMETER        ( ZERO  = 0.0D+0, ONE = 1.0D+0 )

      CALL DAXPY ( N, ALFA, P, 1, X, 1 )
      IF (LINOBJ) CTX = CTX + ALFA*CTP

      IF (HITCON  .AND.  JADD .LE. N) THEN
         BND = BU(JADD)
         IF (HITLOW) BND = BL(JADD)
         IF (ALFA .GE. ZERO) X(JADD) = BND
      END IF
      XNORM  = DNRM2 ( N, X, 1 )

      IF (NCLIN .GT. 0)
     $   CALL DAXPY ( NCLIN, ALFA, AP, 1, AX, 1 )

      IF (NRZ .LE. NRANK) THEN
         IF (UNITGZ) THEN
            RES(NRZ) = RES(NRZ) - ALFA*HZ(NRZ)
         ELSE
            CALL DAXPY ( NRZ, (-ALFA), HZ, 1, RES, 1  )
         END IF

         IF (NUMINF .EQ. 0) THEN

*           Update the transformed gradient GQ so that
*           GQ = GQ + ALFA*R'( HZ ).
*                            ( 0  )

            IF (UNITGZ) THEN
               CALL DAXPY ( N-NRZ+1, ALFA*HZ(NRZ), R(NRZ,NRZ), LDR,
     $                                             GQ(NRZ)   , 1      )
            ELSE
               CALL DCOPY ( NRZ, HZ, 1, WORK, 1 )
               CALL DTRMV ( 'U', 'T', 'N', NRZ, R, LDR, WORK, 1 )
               IF (NRZ .LT. N)
     $            CALL DGEMV ( 'T', NRZ, N-NRZ, ONE, R(1,NRZ+1), LDR,
     $                         HZ, 1, ZERO, WORK(NRZ+1), 1 )
               CALL DAXPY ( N, ALFA, WORK, 1, GQ, 1 )
            END IF
         END IF
      END IF

      RETURN

*     End of  LSMOVE.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSMULS( PRBTYP,
     $                   MSGLVL, N, NACTIV, NFREE,
     $                   LDA, LDT, NUMINF, NZ, NRZ,
     $                   ISTATE, KACTIV, KX, DINKY,
     $                   JSMLST, KSMLST, JINF, JTINY,
     $                   JBIGST, KBIGST, TRULAM,
     $                   A, ANORMS, GQ, RLAMDA, T, WTINF )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*2        PRBTYP
      INTEGER            ISTATE(*), KACTIV(N), KX(N)
      DOUBLE PRECISION   A(LDA,*), ANORMS(*),
     $                   GQ(N), RLAMDA(N), T(LDT,*), WTINF(*)

************************************************************************
*     LSMULS  first computes the Lagrange multiplier estimates for the
*     given working set.  It then determines the values and indices of
*     certain significant multipliers.  In this process, the multipliers
*     for inequalities at their upper bounds are adjusted so that a
*     negative multiplier for an inequality constraint indicates non-
*     optimality.  All adjusted multipliers are scaled by the 2-norm
*     of the associated constraint row.  In the following, the term
*     minimum refers to the ordering of numbers on the real line,  and
*     not to their magnitude.
*
*     JSMLST  is the index of the minimum of the set of adjusted
*             multipliers with values less than  - DINKY.  A negative
*             JSMLST defines the index in Q'g of the artificial
*             constraint to be deleted.
*     KSMLST  marks the position of general constraint JSMLST in KACTIV.
*
*     JBIGST  is the index of the largest of the set of adjusted
*             multipliers with values greater than (1 + DINKY).
*     KBIGST  marks its position in KACTIV.
*
*     On exit,  elements 1 thru NACTIV of RLAMDA contain the unadjusted
*     multipliers for the general constraints.  Elements NACTIV onwards
*     of RLAMDA contain the unadjusted multipliers for the bounds.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of LSMULS dated  30-June-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      INTRINSIC          ABS, MIN
      PARAMETER        ( ZERO   =0.0D+0,ONE    =1.0D+0 )

      NFIXED =   N - NFREE

      JSMLST =   0
      KSMLST =   0
      SMLLST = - DINKY

      TINYLM =   DINKY
      JTINY  =   0

      JBIGST =   0
      KBIGST =   0
      BIGGST =   ONE + DINKY

      IF (NRZ .LT. NZ) THEN
*        ---------------------------------------------------------------
*        Compute JSMLST for the artificial constraints.
*        ---------------------------------------------------------------
         DO 100, J = NRZ+1, NZ
            RLAM = - ABS( GQ(J) )
            IF (RLAM .LT. SMLLST) THEN
               SMLLST =   RLAM
               JSMLST = - J
            ELSE IF (RLAM .LT. TINYLM) THEN
               TINYLM =   RLAM
               JTINY  =   J
            END IF
  100    CONTINUE

         IF (MSGLVL .GE. 20)
     $      WRITE (NOUT, 1000) (GQ(K), K=NRZ+1,NZ)

      END IF

*     ------------------------------------------------------------------
*     Compute JSMLST for regular constraints and temporary bounds.
*     ------------------------------------------------------------------
*     First, compute the Lagrange multipliers for the general
*     constraints in the working set, by solving  T'*lamda = Y'g.

      IF (N .GT. NZ)
     $   CALL DCOPY ( N-NZ, GQ(NZ+1), 1, RLAMDA, 1 )
      IF (NACTIV .GT. 0)
     $   CALL CMTSOL( 2, LDT, NACTIV, T(1,NZ+1), RLAMDA )

*     -----------------------------------------------------------------
*     Now set elements NACTIV, NACTIV+1,... of  RLAMDA  equal to
*     the multipliers for the bound constraints.
*     -----------------------------------------------------------------
      DO 190, L = 1, NFIXED
         J     = KX(NFREE+L)
         BLAM  = RLAMDA(NACTIV+L)
         DO 170, K = 1, NACTIV
            I    = KACTIV(K)
            BLAM = BLAM - A(I,J)*RLAMDA(K)
  170    CONTINUE
         RLAMDA(NACTIV+L) = BLAM
  190 CONTINUE

*     -----------------------------------------------------------------
*     Find JSMLST and KSMLST.
*     -----------------------------------------------------------------
      DO 330, K = 1, N - NZ
         IF (K .GT. NACTIV) THEN
            J = KX(NZ+K)
         ELSE
            J = KACTIV(K) + N
         END IF

         IS   = ISTATE(J)

         I    = J - N
         IF (J .LE. N) ANORMJ = ONE
         IF (J .GT. N) ANORMJ = ANORMS(I)

         RLAM = RLAMDA(K)

*        Change the sign of the estimate if the constraint is in
*        the working set at its upper bound.

         IF (IS .EQ. 2) RLAM =      - RLAM
         IF (IS .EQ. 3) RLAM =   ABS( RLAM )
         IF (IS .EQ. 4) RLAM = - ABS( RLAM )

         IF (IS .NE. 3) THEN
            SCDLAM = RLAM * ANORMJ
            IF      (SCDLAM .LT. SMLLST) THEN
               SMLLST = SCDLAM
               JSMLST = J
               KSMLST = K
            ELSE IF (SCDLAM .LT. TINYLM) THEN
               TINYLM = SCDLAM
               JTINY  = J
            END IF
         END IF

         IF (NUMINF .GT. 0  .AND.  J .GT. JINF) THEN
            SCDLAM = RLAM/WTINF(J)
            IF (SCDLAM .GT. BIGGST) THEN
               BIGGST = SCDLAM
               TRULAM = RLAMDA(K)
               JBIGST = J
               KBIGST = K
            END IF
         END IF
  330 CONTINUE

*     -----------------------------------------------------------------
*     If required, print the multipliers.
*     -----------------------------------------------------------------
      IF (MSGLVL .GE. 20) THEN
         IF (NFIXED .GT. 0)
     $      WRITE (NOUT, 1100) PRBTYP, (KX(NFREE+K),
     $                         RLAMDA(NACTIV+K), K=1,NFIXED)
         IF (NACTIV .GT. 0)
     $      WRITE (NOUT, 1200) PRBTYP, (KACTIV(K),
     $                         RLAMDA(K), K=1,NACTIV)
      END IF

      IF (LSDBG  .AND.  ILSDBG(1) .GT. 0) THEN
         WRITE (NOUT, 9000) JSMLST, SMLLST, KSMLST
         WRITE (NOUT, 9100) JBIGST, BIGGST, KBIGST
         WRITE (NOUT, 9200) JTINY , TINYLM
      END IF

      RETURN

 1000 FORMAT(/ ' Multipliers for the artificial constraints        '
     $       / 4(5X, 1PE11.2))
 1100 FORMAT(/ ' Multipliers for the ', A2, ' bound  constraints   '
     $       / 4(I5, 1PE11.2))
 1200 FORMAT(/ ' Multipliers for the ', A2, ' linear constraints   '
     $       / 4(I5, 1PE11.2))
 9000 FORMAT(/ ' //LSMULS//  JSMLST     SMLLST     KSMLST (Scaled) '
     $       / ' //LSMULS//  ', I6, 1PE11.2, 5X, I6 )
 9100 FORMAT(  ' //LSMULS//  JBIGST     BIGGST     KBIGST (Scaled) '
     $       / ' //LSMULS//  ', I6, 1PE11.2, 5X, I6 )
 9200 FORMAT(  ' //LSMULS//   JTINY     TINYLM                     '
     $       / ' //LSMULS//  ', I6, 1PE11.2)

*     End of  LSMULS.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSOPTN( STRING )
      CHARACTER*(*)      STRING

************************************************************************
*     LSOPTN  loads the option supplied in  STRING  into the relevant
*     element of  IPRMLS  or  RPRMLS.
************************************************************************

      LOGICAL             NEWOPT
      COMMON     /SOL3LS/ NEWOPT
      SAVE       /SOL3LS/

      DOUBLE PRECISION    WMACH(15)
      COMMON     /SOLMCH/ WMACH
      SAVE       /SOLMCH/

      EXTERNAL            MCHPAR
      CHARACTER*16        KEY
      CHARACTER*72        BUFFER
      LOGICAL             FIRST , PRNT
      SAVE                FIRST , NOUT  , PRNT
      DATA                FIRST /.TRUE./

*     If first time in, set  NOUT.
*     NEWOPT  is true first time into  LSFILE  or  LSOPTN
*     and just after a call to  LSSOL.
*     PRNT    is set to true whenever  NEWOPT  is true.

      IF (FIRST) THEN
         FIRST  = .FALSE.
         NEWOPT = .TRUE.
         CALL MCHPAR()
         NOUT   =  WMACH(11)
      END IF
      BUFFER = STRING

*     Call  LSKEY   to decode the option and set the parameter value.
*     If NEWOPT is true, reset PRNT and test specially for NOLIST.

      IF (NEWOPT) THEN
         NEWOPT = .FALSE.
         PRNT   = .TRUE.
         CALL LSKEY ( NOUT, BUFFER, KEY )

         IF (KEY .EQ. 'NOLIST') THEN
            PRNT   = .FALSE.
         ELSE
            WRITE (NOUT, '(// A / A /)')
     $         ' Calls to LSOPTN',
     $         ' ---------------'
            WRITE (NOUT, '( 6X, A )') BUFFER
         END IF
      ELSE
         IF (PRNT)
     $      WRITE (NOUT, '( 6X, A )') BUFFER
         CALL LSKEY ( NOUT, BUFFER, KEY )

         IF (KEY .EQ.   'LIST') PRNT = .TRUE.
         IF (KEY .EQ. 'NOLIST') PRNT = .FALSE.
      END IF

      RETURN

*     End of  LSOPTN.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSPRT ( PRBTYP, PRNT1, ISDEL, ITER, JADD, JDEL,
     $                   MSGLVL, NACTIV, NFREE, N, NCLIN,
     $                   NRANK, LDR, LDT, NZ, NRZ, ISTATE,
     $                   ALFA, CONDRZ, CONDT, GFNORM, GZNORM, GZ1NRM,
     $                   NUMINF, SUMINF, CTX, SSQ,
     $                   AX, R, T, X, WORK )

      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*2        PRBTYP
      LOGICAL            PRNT1
      INTEGER            ISTATE(*)
      DOUBLE PRECISION   AX(*), R(LDR,*), T(LDT,*), X(N)
      DOUBLE PRECISION   WORK(N)

************************************************************************
*  LSPRT  prints various levels of output for  LSCORE.
*
*           Msg        Cumulative result
*           ---        -----------------
*
*        le   0        no output.
*
*        eq   1        nothing now (but full output later).
*
*        eq   5        one terse line of output.
*
*        ge  10        same as 5 (but full output later).
*
*        ge  20        constraint status,  x  and  Ax.
*
*        ge  30        diagonals of  T  and  R.
*
*
*  Debug printing is performed depending on the logical variable  LSDBG.
*  LSDBG  is set true when  IDBG  major iterations have been performed.
*  At this point,  printing is done according to a string of binary
*  digits of the form  SVT  (stored in the integer array  ILSDBG).
*
*  S  set 'on'  gives information from the maximum step routine  CMALF.
*  V  set 'on'  gives various vectors in  LSCORE  and its auxiliaries.
*  T  set 'on'  gives a trace of which routine was called and an
*               indication of the progress of the run.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original version written 31-October-1984.
*  This version of LSPRT dated 14-January-1985.
************************************************************************
      COMMON    /SOL1CM/ NOUT

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      CHARACTER*2        LADD, LDEL
      CHARACTER*2        LSTATE(0:5)
      DATA               LSTATE(0), LSTATE(1), LSTATE(2)
     $                  /'  '     , 'L '     , 'U '     /
      DATA               LSTATE(3), LSTATE(4), LSTATE(5)
     $                  /'E '     , 'T '     , 'Z '     /

      IF (MSGLVL .GE. 15) WRITE (NOUT, 1000) PRBTYP, ITER

      IF (MSGLVL .GE. 5) THEN
         IF      (JDEL .GT. 0) THEN
            KDEL =   ISDEL
         ELSE IF (JDEL .LT. 0) THEN
            JDEL = - JDEL
            KDEL =   5
         ELSE
            KDEL =   0
         END IF

         IF (JADD .GT. 0) THEN
            KADD = ISTATE(JADD)
         ELSE
            KADD = 0
         END IF

         LDEL   = LSTATE(KDEL)
         LADD   = LSTATE(KADD)

         IF (NUMINF .GT. 0) THEN
            OBJ    = SUMINF
         ELSE
            OBJ    = SSQ + CTX
         END IF

*        ---------------------------------------------------------------
*        Print the terse line.
*        ---------------------------------------------------------------
         IF (NRANK .EQ. 0) THEN
            IF (PRNT1  .OR.  MSGLVL .GE. 15) WRITE (NOUT, 1100)
            WRITE (NOUT, 1200) ITER, JDEL, LDEL, JADD, LADD,
     $                         ALFA, NUMINF, OBJ, N-NFREE, NACTIV,
     $                         NZ, NRZ, GFNORM, GZ1NRM, CONDT
         ELSE
            IF (PRNT1  .OR.  MSGLVL .GE. 15) WRITE (NOUT, 1110)
            WRITE (NOUT, 1200) ITER, JDEL, LDEL, JADD, LADD,
     $                         ALFA, NUMINF, OBJ, N-NFREE, NACTIV,
     $                         NZ, NRZ, GFNORM, GZ1NRM, CONDT, CONDRZ
         END IF

         IF (MSGLVL .GE. 20) THEN
            WRITE (NOUT, 2000) PRBTYP
            WRITE (NOUT, 2100) (X(J) , ISTATE(J)  ,  J=1,N)
            IF (NCLIN .GT. 0)
     $      WRITE (NOUT, 2200) (AX(K), ISTATE(N+K), K=1,NCLIN )

            IF (MSGLVL .GE. 30) THEN
*              ---------------------------------------------------------
*              Print the diagonals of  T  and  R.
*              ---------------------------------------------------------
               IF (NACTIV .GT. 0) THEN
                  CALL DCOPY ( NACTIV, T(NACTIV,NZ+1), LDT-1, WORK,1 )
                  WRITE (NOUT, 3000) PRBTYP, (WORK(J), J=1,NACTIV)
               END IF
               IF (NRANK  .GT. 0)
     $            WRITE (NOUT, 3100) PRBTYP, (R(J,J) , J=1,NRANK )
            END IF
            WRITE (NOUT, 5000)
         END IF
      END IF

      PRNT1 = .FALSE.

      RETURN

 1000 FORMAT(/// ' ', A2, ' iteration', I5
     $         / ' =================' )
 1100 FORMAT(// '  Itn Jdel  Jadd      Step',
     $          ' Ninf  Sinf/Objective', '  Bnd', '  Lin', '    Nz',
     $          '   Nz1   Norm Gf  Norm Gz1   Cond T' )
 1110 FORMAT(// '  Itn Jdel  Jadd      Step',
     $          ' Ninf  Sinf/Objective', '  Bnd', '  Lin', '    Nz',
     $          '   Nz1   Norm Gf  Norm Gz1   Cond T Cond Rz1' )
 1200 FORMAT(I5, I5, A1, I5, A1, 1PE9.1, I5, 1X, 1PE15.6, 2I5,
     $       2I6, 1P, 2E10.2, 1P, 2E9.1 )
 2000 FORMAT(/ ' Values and status of the ', A2, ' constraints'
     $       / ' ---------------------------------------' )
 2100 FORMAT(/ ' Variables...'                 /   (1X, 5(1PE15.6, I5)))
 2200 FORMAT(/ ' General linear constraints...'/   (1X, 5(1PE15.6, I5)))
 3000 FORMAT(/ ' Diagonals of ' , A2,' working set factor T'
     $       /   (1P, 5E15.6))
 3100 FORMAT(/ ' Diagonals of ' , A2, ' triangle R         '
     $       /   (1P, 5E15.6))
 5000 FORMAT(/// ' ---------------------------------------------------',
     $           '--------------------------------------------' )

*     End of  LSPRT .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lssetx( linobj, rowerr, unitQ,
     $                   nclin, nactiv, nfree, nrank, nz,
     $                   n, nctotl, ldzy, ldA, ldR, ldT,
     $                   istate, kactiv, kx,
     $                   jmax, errmax, ctx, xnorm,
     $                   A, ax, bl, bu, cq, res, res0, featol,
     $                   R, t, x, zy, p, work )

      implicit           double precision (a-h,o-z)
      logical            linobj, rowerr, unitQ
      integer            istate(nctotl), kactiv(n), kx(n)
      double precision   a(ldA,*), ax(*), bl(nctotl), bu(nctotl),
     $                   cq(*), res(*), res0(*), featol(nctotl), p(n),
     $                   r(ldR,*), T(ldT,*), zy(ldzy,*), x(n)
      double precision   work(nctotl)

************************************************************************
*  LSSETX  computes the point on a working set that is closest to the
*  input vector  x  (in the least-squares sense).  The norm of  x, the
*  transformed residual vector  Pr - RQ'x,  and the constraint values
*  Ax  are also initialized.
*
*  If the computed point gives a row error of more than the feasibility
*  tolerance, an extra step of iterative refinement is used.  If  x  is
*  still infeasible,  the logical variable  ROWERR  is set.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original version written 31-October-1984.
*  This version dated 29-December-1985.
************************************************************************
      common    /sol1cm/ nout

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      external           idamax, ddot
      intrinsic          abs, min
      parameter        ( ntry  = 2 )
      parameter        ( zero  = 0.0d+0, one = 1.0d+0 )

*     ------------------------------------------------------------------
*     Move  x  onto the simple bounds in the working set.
*     ------------------------------------------------------------------
      do 100, k = nfree+1, n
          j   = kx(k)
          is  = istate(j)
          bnd = bl(j)
          if (is .ge. 2) bnd  = bu(j)
          if (is .ne. 4) x(j) = bnd
  100 continue

*     ------------------------------------------------------------------
*     Move  x  onto the general constraints in the working set.
*     We shall make  ntry  tries at getting acceptable row errors.
*     ------------------------------------------------------------------
      ktry   = 1
      jmax   = 1
      errmax = zero

*     repeat
  200    if (nactiv .gt. 0) then

*           Set  work = residuals for constraints in the working set.
*           Solve for p, the smallest correction to x that gives a point
*           on the constraints in the working set.  Define  p = Y*(py),
*           where  py  solves the triangular system  T*(py) = residuals.

            do 220, i = 1, nactiv
               k   = kactiv(i)
               j   = n + k
               bnd = bl(j)
               if (istate(j) .eq. 2) bnd = bu(j)
               work(i) = bnd - ddot  ( n, a(k,1), ldA, x, 1 )
  220       continue

            call cmtsol( 1, ldT, nactiv, T(1,nz+1), work )
            call dload ( n, zero, p, 1 )
            call dcopy ( nactiv, work, 1, p(nz+1), 1 )

            call cmqmul( 2, n, nz, nfree, ldzy, unitQ, kx, p, zy, work )
            call daxpy ( n, one, p, 1, x, 1 )
         end if

*        ---------------------------------------------------------------
*        Compute the 2-norm of  x.
*        Initialize  Ax  for all the general constraints.
*        ---------------------------------------------------------------
         xnorm  = dnrm2 ( n, x, 1 )
         if (nclin .gt. 0)
     $      call dgemv ( 'N', nclin, n, one, A, ldA,
     $                   x, 1, zero, ax, 1 )

*        ---------------------------------------------------------------
*        Check the row residuals.
*        ---------------------------------------------------------------
         if (nactiv .gt. 0) then
            do 300, k = 1, nactiv
               i   = kactiv(k)
               j   = n + i
               is  = istate(j)
               if (is .eq. 1) work(k) = bl(j) - ax(i)
               if (is .ge. 2) work(k) = bu(j) - ax(i)
  300       continue

            jmax   = idamax( nactiv, work, 1 )
            errmax = abs( work(jmax) )
         end if

         ktry = ktry + 1
*     until    (errmax .le. featol(jmax) .or. ktry .gt. ntry
      if (.not.(errmax .le. featol(jmax) .or. ktry .gt. ntry)) go to 200

      rowerr = errmax .gt. featol(jmax)

*     ==================================================================
*     Compute the linear objective value  c'x  and the transformed
*     residual  Pr  -  RQ'x = RES0  -  RQ'x.
*     ==================================================================
      if (nrank .gt. 0  .or.  linobj) then
         call dcopy ( n, x, 1, p, 1 )
         call cmqmul( 6, n, nz, nfree, ldzy, unitQ, kx, p, zy, work )
      end if

      ctx = zero
      if (linobj)
     $   ctx = ddot  ( n, cq, 1, p, 1 )

      if (nrank .gt. 0) then

         call dtrmv ( 'U', 'N', 'N', nrank, R, ldR, p, 1 )
         if (nrank .lt. n)
     $      call dgemv ( 'N', nrank, n-nrank, one, r(1,nrank+1), ldR,
     $                   p(nrank+1), 1, one, p, 1 )

         call dcopy ( nrank,       res0, 1, res, 1 )
         call daxpy ( nrank, -one, p   , 1, res, 1 )

      end if

      if (lsdbg  .and.  ilsdbg(2) .gt. 0)
     $   write (nout, 2200) (x(j), j = 1, n)

      return

 2200 format(/ ' //LSSETX// Variables after refinement ... '/ (5G12.3))

*     End of  lssetx.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lssol ( mm, n,
     $                   nclin, ldA, ldR,
     $                   A, bl, bu, cvec,
     $                   istate, kx, x, R, b,
     $                   inform, iter, obj, clamda,
     $                   iw, leniw, w, lenw )

      implicit           double precision(a-h,o-z)
      integer            leniw, lenw
      integer            istate(n+nclin), kx(n)
      integer            iw(leniw)
      double precision   bl(n+nclin), bu(n+nclin), a(ldA,*)
      double precision   clamda(n+nclin), cvec(*)
      double precision   r(ldR,*), x(n), b(*)
      double precision   w(lenw)

************************************************************************
*  LSSOL  solves problems of the form
*
*           Minimize               F(x)
*              x
*                                 (  x )
*           subject to    bl  .le.(    ).ge.  bu,
*                                 ( Ax )
*
*  where  '  denotes the transpose of a column vector,  x  denotes the
*  n-vector of parameters and  F(x) is one of the following functions..
*
*  FP =  None                         (find a feasible point).
*  LP =  c'x
*  QP1=        1/2 x'Rx                R  n times n, symmetric pos. def.
*  QP2=  c'x + 1/2 x'Rx                .  .   ..        ..       ..  ..
*  QP3=        1/2 x'R'Rx              R  m times n, upper triangular.
*  QP4=  c'x + 1/2 x'R'Rx              .  .   ..  .   ..      ...
*  LS1=        1/2 (b - Rx)'(b - Rx)   R  m times n, rectangular.
*  LS2=  c'x + 1/2 (b - Rx)'(b - Rx)   .  .   ..  .     ...
*  LS3=        1/2 (b - Rx)'(b - Rx)   R  m times n, upper triangular.
*  LS4=  c'x + 1/2 (b - Rx)'(b - Rx)   .  .   ..  .   ..      ...
*
*  The matrix  R  is entered as the two-dimensional array  R  (of row
*  dimension  LDR).  If  LDR = 0,  R  is not accessed.
*
*  The vector  c  is entered in the one-dimensional array  CVEC.
*
*  NCLIN  is the number of general linear constraints (rows of  A).
*  (NCLIN may be zero.)
*
*  The first  N  components of  BL  and   BU  are lower and upper
*  bounds on the variables.  The next  NCLIN  components are
*  lower and upper bounds on the general linear constraints.
*
*  The matrix  A  of coefficients in the general linear constraints
*  is entered as the two-dimensional array  A  (of dimension
*  LDA by N).  If NCLIN = 0, A is not accessed.
*
*  The vector  x  must contain an initial estimate of the solution,
*  and will contain the computed solution on output.
*
*
*  Complete documentation for  LSSOL  is contained in Report SOL 86-1,
*  Users Guide for LSSOL (Version 1.0), by P.E. Gill, S. J. Hammarling,
*  W. Murray, M.A. Saunders and M.H. Wright, Department of
*  Operations Research, Stanford University, Stanford, California 94305.
*
*  Systems Optimization Laboratory, Stanford University.
*  Version 1.00 Dated  30-Jan-1986.
*  Version 1.01 Dated  30-Jun-1986.   Level-2 BLAS added
*  Version 1.02 Dated  13-May-1988.   Level-2 matrix routines added.
*  Version 1.03 Dated  19-Jun-1989.   Some obscure bugs fixed.      
*  Version 1.04 Dated  26-Aug-1991.   nrank bug fixed.      
*
*  Copyright  1984  Stanford University.
*
*  This material may be reproduced by or for the U.S. Government pursu-
*  ant to the copyright license under DAR clause 7-104.9(a) (1979 Mar).
*
*  This material is based upon work partially supported by the National
*  Science Foundation under Grants MCS-7926009 and ECS-8312142; the
*  Department of Energy Contract AM03-76SF00326, PA No. DE-AT03-
*  76ER72018; the Army Research Office Contract DAA29-84-K-0156;
*  and the Office of Naval Research Grant N00014-75-C-0267.
************************************************************************
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
      common    /sol1cm/ nout
      common    /sol3cm/ lennam, ldT, ncolt, ldzy
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5cm/ asize, dtmax, dtmin

      parameter         (lenls = 20)
      common    /sol1ls/ locls(lenls)

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg
*-----------------------------------------------------------------------
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   idbgls, itmax1, itmax2, lcrash, ldbgls, lprob ,
     $                   msgls , nn    , nnclin, nprob , ipadls(20)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolrnk, rpadls(23)

      equivalence       (iprmls(1), idbgls), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*-----------------------------------------------------------------------
      equivalence   (msgls , msglvl), (idbgls, idbg), (ldbgls, msgdbg)

      intrinsic          max, min

*     Local variables.

      logical            cold  , factrz, linobj, named , rowerr,
     $                   unitQ , vertex
      character*2        prbtyp
      character*8        names(1)
      parameter        ( zero   =0.0d+0, point1 =0.1d+0, point3 =3.3d-1)
      parameter        ( point8 =0.8d+0, point9 =0.9d+0, one    =1.0d+0)

      character*40       title
      data               title
     $                 / 'SOL/LSSOL  ---  Version 1.04   Aug  1991' /

*     Set the machine-dependent constants.

      call mchpar()

      epsmch = wmach( 3)
      rteps  = wmach( 4)
      nout   = wmach(11)

      epspt3 = epsmch**point3
      epspt5 = rteps
      epspt8 = epsmch**point8
      epspt9 = epsmch**point9

      named  = .false.

      inform = 0
      iter   = 0

      condmx = one / epspt5

      nctotl = n + nclin

*     Set the default values of the parameters.

      call lsdflt( mm, n, nclin, title )

*     Set all parameters determined by the problem type.

      if      (lprob .eq. 1 ) then
         prbtyp    = 'FP'
         m      = 0
         linobj = .false.
         factrz = .true.
      else if (lprob .eq. 2 ) then
         prbtyp    = 'LP'
         m      = 0
         linobj = .true.
         factrz = .true.
      else if (lprob .eq. 3 ) then
         prbtyp    = 'QP'
         m      = mm
         linobj = .false.
         factrz = .true.
      else if (lprob .eq. 4 ) then
         prbtyp    = 'QP'
         m      = mm
         linobj = .true.
         factrz = .true.
      else if (lprob .eq. 5 ) then
         prbtyp    = 'QP'
         m      = mm
         linobj = .false.
         factrz = .false.
      else if (lprob .eq. 6 ) then
         prbtyp    = 'QP'
         m      = mm
         linobj = .true.
         factrz = .false.
      else if (lprob .eq. 7 ) then
         prbtyp    = 'LS'
         m      = mm
         linobj = .false.
         factrz = .true.
      else if (lprob .eq. 8 ) then
         prbtyp    = 'LS'
         m      = mm
         linobj = .true.
         factrz = .true.
      else if (lprob .eq. 9 ) then
         prbtyp    = 'LS'
         m      = mm
         linobj = .false.
         factrz = .false.
      else if (lprob .eq. 10) then
         prbtyp    = 'LS'
         m      = mm
         linobj = .true.
         factrz = .false.
      end if

*     Assign the dimensions of arrays in the parameter list of LSCORE.
*     Economies of storage are possible if the minimum number of active
*     constraints and the minimum number of fixed variables are known in
*     advance.  The expert user should alter MINACT and MINFXD
*     accordingly.
*     If a linear program is being solved and the matrix of general
*     constraints is fat,  i.e.,  NCLIN .LT. N,  a non-zero value is
*     known for MINFXD.  Note that in this case, VERTEX must be
*     set  .TRUE..

      minact = 0
      minfxd = 0

      vertex = .false.
      if (      (prbtyp .eq. 'LP'  .or.  prbtyp .eq. 'FP')
     $    .and.  nclin  .lt. n   ) then
         minfxd = n - nclin - 1
         vertex = .true.
      end if

      mxfree = n - minfxd
      maxact = max( 1, min( n, nclin ) )
      maxnz  = n - ( minfxd + minact )

      if (nclin .eq. 0) then
         ldzy   = 1
         ldT    = 1
         ncolt  = 1
         vertex = .false.
      else
         ldzy   = max( 1, mxfree )
         ldT    = max( maxnz, maxact )
         ncolt  = mxfree
      end if

      ncnln  = 0
      lennam = 1

*     Allocate certain arrays that are not done in LSLOC.

      litotl = 0

      lax    = 1
      lwtotl = lax + nclin  - 1

*     Allocate remaining work arrays.

      call lsloc ( lprob, n, nclin, litotl, lwtotl )

      cold  = lcrash .eq. 0

*     Check input parameters and storage limits.

      call cmchk ( nerror, msglvl, lcrash, (.not.factrz),
     $             leniw, lenw, litotl, lwtotl,
     $             n, nclin, ncnln,
     $             istate, kx, named, names, lennam,
     $             bl, bu, x )

      if (nerror .gt. 0) then
         inform = 6
         go to 800
      end if

      lkactv = locls( 1)

      lanorm = locls( 2)
      lpx    = locls( 4)
      lres   = locls( 5)
      lres0  = locls( 6)
      lgq    = locls( 8)
      lcq    = locls( 9)
      lrlam  = locls(10)
      lt     = locls(11)
      lzy    = locls(12)
      lwtinf = locls(13)
      lwrk   = locls(14)
      lfeatl = locls(15)

      if (tolfea .gt. zero)
     $   call dload ( n+nclin, (tolfea), w(lfeatl), 1 )

      ianrmj = lanorm
      do 200, j = 1, nclin
         w(ianrmj) = dnrm2 ( n, a(j,1), ldA )
         ianrmj    = ianrmj + 1
  200 continue
      if (nclin .gt. 0)
     $   call dcond ( nclin, w(lanorm), 1, asize, amin )

      call dcond ( nctotl, w(lfeatl), 1, feamax, feamin )
      call dcopy ( nctotl, w(lfeatl), 1, w(lwtinf), 1 )
      call dscal ( nctotl, (one/feamin), w(lwtinf), 1 )

      ssq1   = zero

      if (factrz) then
*        ===============================================================
*        Factorize R using QR or Cholesky.  kx must be initialized.
*        ===============================================================
         do 210, i = 1, n
            kx(i) = i
  210    continue

         if      (prbtyp .eq. 'LP'  .or.  prbtyp .eq. 'FP') then
            nrank = 0
         else if (prbtyp .eq. 'QP') then
*           ------------------------------------------------------------
*           Compute the Cholesky factorization of R.  The Hessian is
*           M times M and resides in the upper left-hand corner of R.
*           ------------------------------------------------------------
            do 220, j = m+1, n
               call dload ( m, (zero), r(1,j), 1 )
  220       continue

            call lschol( ldR, m, nrank, tolrnk, kx, R, info )

            if (nrank .gt. 0)
     $         call dload ( nrank, (zero), w(lres0), 1 )
         else if (prbtyp .eq. 'LS') then
*           ------------------------------------------------------------
*           Compute the orthogonal factorization PRQ = ( U ),  where P
*                                                      ( 0 )
*           is an orthogonal matrix and Q is a permutation matrix.
*           Overwrite R with the upper-triangle U.  The orthogonal
*           matrix P is applied to the residual and discarded.  The
*           permutation is stored in the array KX.  Once U has been
*           computed we need only work with vectors of length N within
*           LSCORE.  However, it is necessary to store the sum of
*           squares of the terms  B(NRANK+1),...,B(M),  where B = Pr.
*           ------------------------------------------------------------
            call dgeqrp( 'Column iterchanges', m, n, R, ldR,
     $                   w(lwrk), iw(lkactv), w(lgq), info )

            lj  = lkactv
            do 230, j = 1, n
               jmax = iw(lj)
               if (jmax .gt. j) then
                  jsave    = kx(jmax)
                  kx(jmax) = kx(j)
                  kx(j)    = jsave
               end if
               lj = lj + 1
  230       continue

            call dgeapq( 'Transpose', 'Separate', m, min( n,m-1 ), 
     $                   R, ldR, w(lwrk), 1, b, m, w(lgq), info )

            rownrm = dnrm2 ( n, R(1,1), ldR )
            if (          rownrm  .le.        tolrnk
     $          .or.  abs(R(1,1)) .le. rownrm*tolrnk) then
               nrank = 0
            else
               nrank = idrank( min(n, m), R, ldR+1, tolrnk )
            end if

            if (m .gt. nrank) ssq1 = dnrm2 ( m-nrank, b(nrank+1), 1 )

            if (nrank .gt. 0)
     $         call dcopy ( nrank, b, 1, w(lres0), 1 )
         end if
      else
*        ===============================================================
*        R is input as an upper-triangular matrix with m rows.
*        ===============================================================
         nrank = m
         if (nrank .gt. 0) then
            if      (prbtyp .eq. 'QP') then
               call dload ( nrank, (zero), w(lres0), 1 )
            else if (prbtyp .eq. 'LS') then
               call dcopy ( nrank, b, 1, w(lres0), 1 )
            end if
         end if
      end if

      if (       msglvl .gt. 0     .and.  nrank  .lt. n
     $    .and.  prbtyp .ne. 'LP'  .and.  prbtyp .ne. 'FP')
     $   write (nout, 9000) nrank

*     ------------------------------------------------------------------
*     Find an initial working set.
*     ------------------------------------------------------------------
      call lscrsh( cold, vertex,
     $             nclin, nctotl, nactiv, nartif,
     $             nfree, n, ldA,
     $             istate, iw(lkactv),
     $             bigbnd, tolact,
     $             A, w(lax), bl, bu, x, w(lgq), w(lwrk) )

*     ------------------------------------------------------------------
*     Compute the TQ factorization of the constraints while keeping R in
*     upper-triangular form.  Transformations associated with Q are
*     applied to CQ.  Transformations associated with P are applied to
*     RES0.  If some simple bounds are in the working set,  KX is
*     re-ordered so that the free variables come first.
*     ------------------------------------------------------------------
*     First, add the bounds. To save a bit of work, CQ is not loaded
*     until after KX has been re-ordered.

      ngq   = 0
      nres  = 0
      if (nrank .gt. 0) nres = 1
      unitQ = .true.

      call lsbnds( unitQ,
     $             inform, nz, nfree, nrank, nres, ngq,
     $             n, ldzy, ldA, ldR, ldT,
     $             istate, kx, condmx,
     $             A, R, w(lt), w(lres0), w(lcq), w(lzy),
     $             w(lwrk), w(lpx), w(lrlam) )

      if (linobj) then

*        Install the transformed linear term in CQ.
*        CMQMUL applies the permutations in KX to CVEC.

         ngq = 1
         call dcopy ( n, cvec, 1, w(lcq), 1 )
         call cmqmul( 6, n, nz, nfree, ldzy, unitQ,
     $                kx, w(lcq), w(lzy), w(lwrk) )
      end if

      if (nactiv .gt. 0) then
         nact1  = nactiv
         nactiv = 0

         call lsadds( unitQ, vertex,
     $                inform, 1, nact1, nactiv, nartif, nz, nfree,
     $                nrank, nrejtd, nres, ngq,
     $                n, ldzy, ldA, ldR, ldT,
     $                istate, iw(lkactv), kx, condmx,
     $                A, R, w(lt), w(lres0), w(lcq), w(lzy),
     $                w(lwrk), w(lpx), w(lrlam) )
      end if

*     ------------------------------------------------------------------
*     Move the initial  x  onto the constraints in the working set.
*     Compute the transformed residual vector  Pr = Pb - RQ'x.
*     ------------------------------------------------------------------
      call lssetx( linobj, rowerr, unitQ,
     $             nclin, nactiv, nfree, nrank, nz,
     $             n, nctotl, ldzy, ldA, ldR, ldT,
     $             istate, iw(lkactv), kx,
     $             jmax, errmax, ctx, xnorm,
     $             A, w(lax), bl, bu, w(lcq), w(lres), w(lres0),
     $             w(lfeatl), R, w(lt), x, w(lzy), w(lpx), w(lwrk) )

      jinf = 0

      call lscore( prbtyp, named, names, linobj, unitQ,
     $             inform, iter, jinf, nclin, nctotl,
     $             nactiv, nfree, nrank, nz, nrz,
     $             n, ldA, ldR,
     $             istate, iw(lkactv), kx,
     $             ctx, obj, ssq1,
     $             suminf, numinf, xnorm,
     $             bl, bu, A, clamda, w(lax),
     $             w(lfeatl), R, x, iw, w )

      obj    = obj    + ctx
      if (prbtyp .eq. 'LS'  .and.  nrank .gt. 0)
     $   call dcopy ( nrank, w(lres), 1, b, 1 )

*     ==================================================================
*     Print messages if required.
*     ==================================================================
  800 if (msglvl .gt.   0) then
         if (inform .eq.   0) then
            if (prbtyp .eq. 'FP') then
               write (nout, 2001)
            else
               write (nout, 2002) prbtyp
            end if
         end if
         if (inform .eq.   1) write (nout, 2010) prbtyp
         if (inform .eq.   2) write (nout, 2020) prbtyp
         if (inform .eq.   3) write (nout, 2030)
         if (inform .eq.   4) write (nout, 2040)
         if (inform .eq.   5) write (nout, 2050)
         if (inform .eq.   6) write (nout, 2060) nerror

         if (inform .lt.   6) then
            if      (numinf .eq. 0) then
                if (prbtyp .ne. 'FP') write (nout, 3000) prbtyp, obj
            else if (inform .eq. 3) then
               write (nout, 3010) suminf
            else
               write (nout, 3020) suminf
            end if
            if (numinf .gt. 0) obj = suminf
         end if
      end if

*     Recover the optional parameters set by the user.

      call icopy ( mxparm, ipsvls, 1, iprmls, 1 )
      call dcopy ( mxparm, rpsvls, 1, rprmls, 1 )

      return

 2001 format(/ ' Exit LSSOL - Feasible point found.     ')
 2002 format(/ ' Exit LSSOL - Optimal ', A2, ' solution.')
 2010 format(/ ' Exit LSSOL - Weak ',    A2, ' solution.')
 2020 format(/ ' Exit LSSOL - ', A2,         ' solution is unbounded.' )
 2030 format(/ ' Exit LSSOL - Cannot satisfy the linear constraints. ' )
 2040 format(/ ' Exit LSSOL - Too many iterations.')
 2050 format(/ ' Exit LSSOL - Too many iterations without changing X.' )
 2060 format(/ ' Exit LSSOL - ', I10, ' errors found in the input',
     $         ' parameters.  Problem abandoned.'         )
 3000 format(/ ' Final ', A2, ' objective value =', G16.7 )
 3010 format(/ ' Minimum sum of infeasibilities =', G16.7 )
 3020 format(/ ' Final sum of infeasibilities =',   G16.7 )

 9000 format(/ ' Rank of the objective function data matrix = ', I5 )

*     End of  LSSOL .

      end
