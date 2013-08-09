*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     file  npsubs.f
*
*     npsol    npalf    npchkd   npcore   npcrsh   npdflt   npfd     
*     npfeas   npfile   npgetr   npiqp    npkey    nploc    npmrt    
*     npoptn   npprt    nprset   npsavr   npsetx   npsrch   npupdt
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npsol ( n, nclin, ncnln, lda, ldcju, ldr,
     $                   a, bl, bu,
     $                   confun, objfun,
     $                   inform, iter, istate,
     $                   c, cjacu, clamda, objf, gradu, r, x,
     $                   iw, leniw, w, lenw )

      implicit           double precision (a-h,o-z)
      external           confun, objfun
      integer            istate(n+nclin+ncnln)
      integer            iw(leniw)
      double precision   a(lda,*), bl(n+nclin+ncnln),
     $                   bu(n+nclin+ncnln)
      double precision   c(*), cjacu(ldcju,*), clamda(n+nclin+ncnln)
      double precision   gradu(n), r(ldr,*), x(n)
      double precision   w(lenw)

*-----------------------------------------------------------------------
*
*  NPSOL   solves the nonlinear programming problem
*
*            minimize                   F(x)
*
*                                    (    x  )
*            subject to    bl  .le.  (  A*x  )  .le.  bu
*                                    (  c(x) )
*
*  where  F(x)  is a smooth scalar function,  A  is a constant matrix
*  and  c(x)  is a vector of smooth nonlinear functions.  The feasible
*  region is defined by a mixture of linear and nonlinear equality or
*  inequality constraints on  x.
*
*  The dimensions of the problem are...
*
*  N        the number of variables (dimension of  x),
*
*  NCLIN    the number of linear constraints (rows of the matrix  A),
*
*  NCNLN    the number of nonlinear constraints (dimension of  c(x)),
*
*
*  NPSOL   uses a sequential quadratic programming algorithm, with a
*  positive-definite quasi-Newton approximation to the transformed
*  Hessian  Q'HQ  of the Lagrangian function (which will be stored in
*  the array  R).
*
*
*  Complete documentation for  NPSOL  is contained in Report
*  SOL 86-2, Users guide for NPSOL (Version 4.0), by P.E. Gill,
*  W. Murray, M.A. Saunders and M.H. Wright, Department of Operations
*  Research,  Stanford University, Stanford, California 94305.
*
*  Systems Optimization Laboratory, Stanford University.
*  Version 1.1,  April     12, 1983. (The less said about this one.....)
*  Version 2.0,  April     30, 1984.
*  Version 3.0,  March     20, 1985. (First Fortran 77 version).
*  Version 3.2,  August    20, 1985.
*  Version 4.0,  April     16, 1986. (First version with differences).
*  Version 4.01, June      30, 1986. (Level 2 BLAS + F77 linesearch).
*  Version 4.02, August     5, 1986. (Reset SSBFGS. One call to CHFD).
*  Version 4.03, June      14, 1987. (Step limit).
*  Version 4.04, June      28, 1989. (Vectorizable BLAS).
*  Version 4.05, November  28, 1989. (Load and save files added).
*  Version 4.06, November   5, 1991. (New versions of srchq and srchc).
*
*  Copyright  1983  Stanford University.
*
*  This material may be reproduced by or for the U.S. Government pursu-
*  ant to the copyright license under DAR Clause 7-104.9(a) (1979 Mar).
*
*  This material is based upon work partially supported by the National
*  Science Foundation under Grants MCS-7926009 and ECS-8312142; the
*  Department of Energy Contract AM03-76SF00326, PA No. DE-AT03-
*  76ER72018; the Army Research Office Contract DAA29-84-K-0156;
*  and the Office of Naval Research Grant N00014-75-C-0267.
*  ---------------------------------------------------------------------

*  Common blocks.

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
*-----include npparm----------------------------------------------------
      integer            iprmnp(mxparm), ipsvnp
      double precision   rprmnp(mxparm), rpsvnp

      common    /nppar1/ ipsvnp(mxparm),
     $                   idbgnp, itmxnp, jvrfy1, jvrfy2, jvrfy3, jvrfy4,
     $                   ldbgnp, lformh, lvlder, lverfy, msgnp , nlnf  ,
     $                   nlnj  , nlnx  , nncnln, nsave , nload , ksave ,
     $                   ipadnp(12)

      common    /nppar2/ rpsvnp(mxparm),
     $                   cdint , ctol  , dxlim , epsrf , eta   , fdint ,
     $                   ftol  , hcndbd, rpadnp(22)

      equivalence       (iprmnp(1), idbgnp), (rprmnp(1), cdint)

      save      /nppar1/, /nppar2/
*-----------------------------------------------------------------------
      equivalence  (idbgnp, idbg  ), (itmxnp, nmajor), (itmax2, nminor)
      equivalence  (ldbgls, mnrdbg), (ldbgnp, mjrdbg), (msgls , msgqp )

      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      common    /sol1cm/ nout
      common    /sol3cm/ lennam, ldt   , ncolt , ldQ
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5cm/ asize , dtmax , dtmin
      common    /sol6cm/ rcndbd, rfrobn, drmax , drmin

      logical            unitQ
      common    /sol1sv/ nactiv, nfree , nz   , unitQ
      save      /sol1sv/

      parameter         (lenls = 20)
      common    /sol1ls/ locls(lenls)

      parameter         (lennp = 35)
      common    /sol1np/ locnp(lennp)
      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset
      common    /sol5np/ lvrfyc, jverfy(4)
      logical            incrun
      common    /sol6np/ rhomax, rhonrm, rhodmp, scale, incrun

      logical            cmdbg, lsdbg, npdbg
      parameter         (ldbg = 5)
      common    /npdebg/ inpdbg(ldbg), npdbg
      common    /lsdebg/ ilsdbg(ldbg), lsdbg
      common    /cmdebg/ icmdbg(ldbg), cmdbg

      intrinsic          abs   , max   , min   , mod   , sqrt  , real

*     Local variables.

      external           ddiv  , ddot  , dnorm , dnrm2
      character*8        names(1)
      logical            cold  , linobj, named , overfl, rowerr, vertex
      logical            needfd
      parameter         (zero   =0.0d+0, point1 =0.1d+0, point3 =3.3d-1)
      parameter         (point8 =0.8d+0, point9 =0.9d+0, one    =1.0d+0)
      parameter         (growth =1.0d+2                                )

      character*36       title
      data               title
     $                 / 'NPSOL  ---  Version 4.06   Nov  1991' /

*     Set the machine-dependent constants.

      call mchpar()

      epsmch = wmach( 3)
      rteps  = wmach( 4)
      nout   = wmach(11)

      epspt3 = epsmch**point3
      epspt5 = rteps
      epspt8 = epsmch**point8
      epspt9 = epsmch**point9

      rhomax = one/epsmch
      rootn  = sqrt(real(n))

*     Default names will be provided for variables during printing.

      named  = .false.
      inform = 0

*     Set the default values for the parameters.

      call npdflt( n, nclin, ncnln, leniw, lenw, title )

      needfd = lvlder .eq. 0  .or.  lvlder .eq. 2
     $                        .or. (lvlder .eq. 1  .and.  ncnln .gt. 0)
      cold   = lcrash .eq. 0
      lvldif = 0
      if (needfd) lvldif = 1

      nplin  = n     + nclin
      nctotl = nplin + ncnln

*     Assign the dimensions of arrays in the parameter list of NPCORE.
*     Economies of storage are possible if the minimum number of active
*     constraints and the minimum number of fixed variables are known in
*     advance.  The expert user should alter MINACT and MINFXD
*     accordingly.

      minact = 0
      minfxd = 0

      mxfree = n - minfxd
      maxact = max( 1, min( n, nclin ) )
      maxnz  = n - ( minfxd + minact )

      if (nclin + ncnln .eq. 0) then
         ldQ   = 1
         ldt   = 1
         ncolt = 1
      else
         ldQ   = max( 1, mxfree )
         ldt   = max( maxnz, maxact )
         ncolt = mxfree
      end if

      lennam = 1

      ldaqp = max( nclin+ncnln, 1 )
      if (ncnln .eq. 0  .and.  nclin .gt. 0) ldaqp = lda

*     nploc  defines the arrays that contain the locations of various
*     work arrays within  w  and  iw.

      litotl = 0
      lwtotl = 0
      call nploc( n, nclin, ncnln, nctotl, litotl, lwtotl)

*     Allocate certain addresses that are not allocated in NPLOC.

      lax    = lwtotl + 1
      lwtotl = lax    + nclin - 1
      lax    = min( lax, lwtotl )

*     Check input parameters and storage limits.

      call cmchk ( nerror, msgnp, lcrash, .false.,
     $             leniw, lenw, litotl, lwtotl,
     $             n, nclin, ncnln,
     $             istate, iw, named, names, lennam,
     $             bl, bu, x )

      if (nerror .gt. 0) then
         inform = 9
         go to 800
      end if

      lkactv = locls( 1)
      lanorm = locls( 2)
      lcjdx  = locls( 3)
      lres   = locls( 5)
      lres0  = locls( 6)
      lgq    = locls( 9)
      lrlam  = locls(10)
      lt     = locls(11)
      lQ     = locls(12)
      lwtinf = locls(13)
      lwrk1  = locls(14)

      lkx    = locnp( 1)
      liperm = locnp( 2)
      laqp   = locnp( 3)
      ldx    = locnp( 7)
      lfeatl = locnp(10)
      lwrk2  = locnp(12)

      lcmul  = locnp(16)
      lrho   = locnp(20)
      lwrk3  = locnp(21)
      lneedc = locnp(24)
      lhfrwd = locnp(25)
      lhctrl = locnp(26)
      lcjac  = locnp(27)
      lgrad  = locnp(28)

      ldcj  = max ( ncnln, 1 )

      tolrnk = zero
      rcndbd = sqrt( hcndbd )

*     ==================================================================
*     If a unit number for a load file has been set, read initial values
*     from an old run.  These values override existing settings. 
*     ==================================================================
      if (nload .gt. 0) then
         call npgetr( nerror, unitQ, n, nclin, ncnln, ldr, ldQ,
     $                nfree0, iter, istate, iw(lkx),
     $                w(lhfrwd), w(lhctrl),
     $                w(lcmul), r, w(lrho), x, w(lQ) )

         if (nerror .gt. 0) then
            inform = 9
            go to 800
         end if
      end if

*     ==================================================================
*     Load the arrays of feasibility tolerances.
*     ==================================================================
      if (tolfea .gt. zero)
     $   call dload ( nplin, tolfea, w(lfeatl), 1 )

      if (ncnln .gt. 0  .and.  ctol .gt. zero)
     $   call dload ( ncnln, ctol, w(lfeatl+nplin), 1 )

      if (lfdset .eq. 0) then
         fdchk = sqrt( epsrf )
      else if (lfdset .eq. 1) then
         fdchk = fdint
      else
         fdchk = w(lhfrwd)
      end if

      nfun   = 0
      ngrad  = 0
      nstate = 1

*     ------------------------------------------------------------------
*     If required,  compute the problem functions.
*     If the constraints are nonlinear,  the first call of CONFUN
*     sets up any constant elements in the Jacobian matrix.  A copy of
*     the Jacobian (with constant elements set) is placed in  CJACU.
*     ------------------------------------------------------------------
      if (lverfy .ge. 10) then
         xnorm  = dnrm2 ( n, x, 1 )
         lvrfyc = lverfy - 10

         call npchkd( info, msgnp, nstate, lvlder, nfun, ngrad,
     $                ldcj, ldcju, n, ncnln,
     $                confun, objfun, iw(lneedc),
     $                bigbnd, epsrf, cdint, fdint,
     $                fdchk, fdnorm, objf, xnorm,
     $                bl, bu, c, w(lwrk3), w(lcjac), cjacu, w(lcjdx),
     $                w(ldx), w(lgrad), gradu, w(lhfrwd), w(lhctrl),
     $                x, w(lwrk1), w(lwrk2), w, lenw )

         if (info .ne. 0) then
            if (info .gt. 0) inform = 7
            if (info .lt. 0) inform = info
            go to 800
         end if
         nstate = 0
      end if

      call icopy ( ldbg, ilsdbg, 1, icmdbg, 1 )

      if (nclin .gt. 0) then
         ianrmj = lanorm
         do 110, j = 1, nclin
            w(ianrmj) = dnrm2 ( n, a(j,1), lda )
            ianrmj    = ianrmj + 1
  110    continue
         call dcond ( nclin, w(lanorm), 1, asize, amin )
      end if

      call dcond ( nplin, w(lfeatl), 1, feamax, feamin )
      call dcopy ( nplin, w(lfeatl), 1, w(lwtinf), 1 )
      call dscal ( nplin, (one/feamin), w(lwtinf), 1 )

*     ==================================================================
*     The input values of x and (optionally)  istate are used by
*     lscrsh  to define an initial working set.
*     ==================================================================
      vertex = .false.
      call lscrsh( cold, vertex,
     $             nclin, nplin, nactiv, nartif,
     $             nfree, n, lda,
     $             istate, iw(lkactv),
     $             bigbnd, tolact,
     $             a, w(lax), bl, bu, x, w(lwrk1), w(lwrk2) )

      nres   = 0
      ngq    = 0
      condmx = one / epspt5

      if (lcrash .le. 1) then
*        ===============================================================
*        Cold or warm start. The upper-triangular matrix R is the factor
*        of an approximate Lagrangian Hessian.
*        ===============================================================
         unitQ  = .true.
         iter   = 0

         ikx    = lkx
         do 120, i = 1, n
            iw(ikx) = i
            ikx     = ikx + 1
  120    continue

         if (cold) then
            call f06qhf( 'Upper-triangular', n, n, zero, one, r, ldr )
            rfrobn = rootn

            nrank  = 0
            if (ncnln .gt. 0) call dload ( ncnln, (zero), w(lcmul), 1 )
         else

*           R will be updated while finding a feasible x.

            nrank  = nlnx
            call dload ( nlnx, (zero), w(lres0), 1 )
            if (ncnln .gt. 0)
     $         call dcopy ( ncnln, clamda(nplin+1), 1, w(lcmul), 1 )

         end if

         incrun = .true.
         rhonrm =  zero
         rhodmp =  one
         scale  =  one
         call dload ( ncnln, (zero), w(lrho), 1 )

*        ---------------------------------------------------------------
*        Re-order kx so that the free variables come first.
*        If a warm start is required, nrank will be nonzero and the
*        factor R will be updated.
*        ---------------------------------------------------------------
         call lsbnds( unitQ,
     $                inform, nz, nfree, nrank, nres, ngq,
     $                n, ldQ, lda, ldr, ldt,
     $                istate, iw(lkx), condmx,
     $                a, r, w(lt), w(lres0), w(lgq), w(lQ),
     $                w(lwrk1), w(lwrk2), w(lrlam) )             

      else
*        ===============================================================
*        Hot start.
*        Stop if the computed and input values of nfree don't match.
*        ===============================================================
         if (nfree0 .ne. nfree) then
            nerror = 1
            inform = 9
            go to 800
         end if
      end if

*     ------------------------------------------------------------------
*     Factorize the linear constraints in the initial working set.
*     ------------------------------------------------------------------
      if (nactiv .gt. 0) then
         nact1  = nactiv
         nactiv = 0

         call lsadds( unitQ, vertex,
     $                inform, 1, nact1, nactiv, nartif, nz, nfree,
     $                nrank, nrejtd, nres, ngq,
     $                n, ldQ, lda, ldr, ldt,
     $                istate, iw(lkactv), iw(lkx), condmx,
     $                a, r, w(lt), w(lres0), w(lgq), w(lQ),
     $                w(lwrk1), w(lwrk2), w(lrlam) )
      end if

      if (lcrash .le. 1) then
*        ===============================================================
*        Cold or warm start.  Move  x  on to the linear constraints and
*        find a feasible point.
*        ===============================================================
         ssq1   =  zero
         linobj = .false.
         call lssetx( linobj, rowerr, unitQ,
     $                nclin, nactiv, nfree, nrank, nz,
     $                n, nplin, ldQ, lda, ldr, ldt,
     $                istate, iw(lkactv), iw(lkx),
     $                jmax, errmax, ctx, xnorm,
     $                a, w(lax), bl, bu, w(lgq), w(lres), w(lres0),
     $                w(lfeatl), r, w(lt), x, w(lQ),w(lwrk1),w(lwrk2) )

*        ---------------------------------------------------------------
*        Call  lscore  to find a feasible  x.
*        ---------------------------------------------------------------
*        Use  work2  as the multiplier vector.

         jinf   = 0
         lclam  = lwrk2

         idbgsv = idbg
         if (idbg .gt. 0) then
            idbg = nminor + 1
         end if

         call lscore( 'FP problem', named, names, linobj, unitQ,
     $                nlperr, itns, jinf, nclin, nplin,
     $                nactiv, nfree, nrank, nz, nz1,
     $                n, lda, ldr,
     $                istate, iw(lkactv), iw(lkx),
     $                ctx, obj, ssq1, suminf, numinf, xnorm,
     $                bl, bu, a, w(lclam), w(lax),
     $                w(lfeatl), r, x, iw, w )

         if (nlperr .gt. 0) then
            inform = 2
            go to 800
         end if

         idbg  = idbgsv
         call icopy ( ldbg, inpdbg, 1, icmdbg, 1 )
      else
*        ---------------------------------------------------------------
*        Hot start.
*        The point  x  is preassigned.  Compute the 2-norm of  x.
*        Initialize  Ax  for the linear constraints.
*        ---------------------------------------------------------------
         nrank  = nlnx
         xnorm  = dnrm2 ( n, x, 1 )
         if (nclin .gt. 0)
     $      call dgemv ( 'N', nclin, n, one, a, lda,
     $                   x, 1, zero, w(lax), 1 )
      end if

      if (lcrash .gt. 0) then

*        Check for a bad R.

         rfrobn = f06qgf( 'Frobenius norm', 'Upper', n, n, r, ldr )
         call dcond ( n, r, ldr+1, drmax, drmin )
         cond   = ddiv  ( drmax, drmin, overfl )

         if (      cond   .gt. rcndbd
     $       .or.  rfrobn .gt. rootn*growth*drmax) then
*           ------------------------------------------------------------
*           Refactorize the Hessian and bound the condition estimator.
*           ------------------------------------------------------------
            write( nout, 9000 )
            call nprset( unitQ,
     $                   n, nfree, nz, ldQ, ldr,
     $                   iw(liperm), iw(lkx),
     $                   w(lgq), r, w(lQ), w(lwrk1), w(lres0) )
         end if
      end if

*     ==================================================================
*     Now we can check the gradients at a feasible x.
*     ==================================================================
      lvrfyc = lverfy
      if (lverfy .ge. 10) lvrfyc = -1

      call npchkd( info, msgnp, nstate, lvlder, nfun, ngrad,
     $             ldcj, ldcju, n, ncnln,
     $             confun, objfun, iw(lneedc),
     $             bigbnd, epsrf, cdint, fdint,
     $             fdchk, fdnorm, objf, xnorm,
     $             bl, bu, c, w(lwrk3), w(lcjac), cjacu, w(lcjdx),
     $             w(ldx), w(lgrad), gradu, w(lhfrwd), w(lhctrl),
     $             x, w(lwrk1), w(lwrk2), w, lenw )

      if (info .ne. 0) then
         if (info .gt. 0) inform = 7
         if (info .lt. 0) inform = info
         go to 800
      end if

      call dcopy ( n, w(lgrad), 1, w(lgq), 1 )
      call cmqmul( 6, n, nz, nfree, ldQ, unitQ,
     $             iw(lkx), w(lgq), w(lQ), w(lwrk1) )

*     ==================================================================
*     Solve the problem.
*     ==================================================================
      if (ncnln .eq. 0) then
*        ---------------------------------------------------------------
*        The problem has only linear constraints and bounds.
*        ---------------------------------------------------------------
         call npcore( named, names, unitQ, inform, iter,
     $                n, nclin, ncnln, nctotl, nactiv, nfree, nz,
     $                lda, ldcj, ldcju, ldaqp, ldr,
     $                nfun, ngrad, istate, iw(lkactv), iw(lkx),
     $                objf, fdnorm, xnorm, confun, objfun,
     $                a, w(lax), bl, bu, c, w(lcjac), cjacu, clamda,
     $                w(lfeatl), w(lgrad), gradu, r, x, iw, w, lenw )
      else
*        ---------------------------------------------------------------
*        The problem has some nonlinear constraints.
*        ---------------------------------------------------------------
         if (nclin .gt. 0)
     $      call f06qff( 'General', nclin, n, a, lda, w(laqp), ldaqp )

*        Try and add some nonlinear constraint indices to KACTIV.
*
         call npcrsh( cold, n, nclin, ncnln,
     $                nctotl, nactiv, nfree, nz,
     $                istate, iw(lkactv), bigbnd, tolact,
     $                bl, bu, c )

         call npcore( named, names, unitQ, inform, iter,
     $                n, nclin, ncnln, nctotl, nactiv, nfree, nz,
     $                lda, ldcj, ldcju, ldaqp, ldr,
     $                nfun, ngrad, istate, iw(lkactv),iw(lkx),
     $                objf, fdnorm, xnorm, confun, objfun,
     $                w(laqp), w(lax), bl, bu, c, w(lcjac),cjacu,clamda,
     $                w(lfeatl), w(lgrad), gradu, r, x, iw, w, lenw )

      end if

*     ==================================================================
*     If a unit number for a save file has been set, save the details of
*     this run.
*     ==================================================================
      if (nsave .gt. 0  .and.  ksave .gt. nmajor) then
         call npsavr( unitQ, n, nclin, ncnln, ldr, ldQ,
     $                nfree, nsave, iter, istate, iw(lkx),
     $                w(lhfrwd), w(lhctrl),
     $                w(lcmul), r, w(lrho), x, w(lQ) )
      end if

*     ------------------------------------------------------------------
*     If required, form the triangular factor of the Hessian.
*     ------------------------------------------------------------------
*     First,  form the square matrix  R  such that  H = R'R.
*     Compute the  QR  factorization of  R.

      if (lformh .gt. 0) then
         lv     = lwrk2
         do 400 j = 1, n
            if (j .gt. 1)
     $         call dload ( j-1, zero, w(lv), 1 )

            lvj = lv + j - 1
            call dcopy ( n-j+1, r(j,j), ldr, w(lvj), 1     )
            call cmqmul( 3, n, nz, nfree, ldQ, unitQ,
     $                   iw(lkx), w(lv), w(lQ), w(lwrk1) )
            call dcopy ( n    , w(lv) , 1    , r(j,1), ldr )
  400    continue

         call dgeqr ( n, n, r, ldr, w(lwrk1), info )
      end if

*     ==================================================================
*     Print messages if required.
*     ==================================================================
  800 if (msgnp .gt.   0) then
         if (inform .lt.   0) write (nout, 3000)
         if (inform .eq.   0) write (nout, 4000)
         if (inform .eq.   1) write (nout, 4100)
         if (inform .eq.   2) write (nout, 4200)
         if (inform .eq.   3) write (nout, 4300)
         if (inform .eq.   4) write (nout, 4400)
         if (inform .eq.   5) write (nout, 4500)
         if (inform .eq.   6) write (nout, 4600)
         if (inform .eq.   7) write (nout, 4700)
         if (inform .eq.   9) write (nout, 4900) nerror

         if (inform .ge. 0  .and.  inform .lt. 7) then
            if (nlperr .eq. 0) then
               write (nout, 5000) objf
            else
               if (nlperr .eq. 3) then
                  write (nout, 5010) suminf
               else
                  write (nout, 5020) suminf
               end if
            end if
         end if
      end if

*     Recover the optional parameters set by the user.
*     Load the final constraint and objective gradients.

      call icopy ( mxparm, ipsvls, 1, iprmls, 1 )
      call dcopy ( mxparm, rpsvls, 1, rprmls, 1 )
      call icopy ( mxparm, ipsvnp, 1, iprmnp, 1 )
      call dcopy ( mxparm, rpsvnp, 1, rprmnp, 1 )

      if (ncnln .gt. 0) then
         call f06qff( 'General', ncnln, n, w(lcjac), ldcj, cjacu, ldcju)
      end if
      call dcopy ( n, w(lgrad), 1, gradu, 1 )

      return

 3000 format(/ ' Exit NPSOL - User requested termination.'          )
 4000 format(/ ' Exit NPSOL - Optimal solution found.'              )
 4100 format(/ ' Exit NPSOL - Optimal solution found, ',
     $         ' but the requested accuracy could not be achieved.' )
 4200 format(/ ' Exit NPSOL - No feasible point for the linear',
     $         ' constraints.')
 4300 format(/ ' Exit NPSOL - No feasible point for the nonlinear',
     $         ' constraints.')
 4400 format(/ ' Exit NPSOL - Too many major iterations.             ')
 4500 format(/ ' Exit NPSOL - Problem is unbounded (or badly scaled).')
 4600 format(/ ' Exit NPSOL - Current point cannot be improved upon. ')
 4700 format(/ ' Exit NPSOL - Large errors found in the derivatives. ')

 4900 format(/ ' Exit NPSOL - ', I10, ' errors found in the input',
     $         ' parameters.  Problem abandoned.')
 5000 format(/ ' Final nonlinear objective value =', G16.7 )
 5010 format(/ ' Minimum sum of infeasibilities =',  G16.7 )
 5020 format(/ ' Final sum of infeasibilities =',    G16.7 )

 9000 format(  ' XXX  Bad initial Hessian,   R  refactorized.' )

*     End of  npsol .

      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npalf ( inform, n, nclin, ncnln,
     $                   alfa, alfmin, alfmax, bigbnd, dxnorm,
     $                   anorm, adx, ax, bl, bu,
     $                   dslk, dx, slk, x )

      implicit           double precision (a-h,o-z)
      double precision   anorm(*), adx(*), ax(*), bl(*), bu(*),
     $                   dslk(*), dx(n), slk(*), x(n)

************************************************************************
*     npalf   finds a step alfa such that the point x + alfa*p reaches
*     one of the slacks or linear constraints.  The step alfa is the
*     maximum step that can be taken without violating one of the slacks
*     or linear constraints that is currently satisfied.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 77 version written  June 1986.
*     This version of npalf dated  13-Jun-1987.
************************************************************************
      common    /sol1cm/ nout
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      logical            cmdbg
      integer            lcmdbg
      parameter         (lcmdbg = 5)
      common    /cmdebg/ icmdbg(lcmdbg), cmdbg

      intrinsic          abs, max, min
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      if (cmdbg  .and.  icmdbg(3) .gt. 0) write (nout, 1000)

      alfa   = alfmax
      j      = 1

*+    while (j .le. n+nclin+ncnln .and. alfa .gt. alfmin) do
  100 if    (j .le. n+nclin+ncnln .and. alfa .gt. alfmin) then

         if      (j .le. n      ) then
            axi    =  x(j)
            adxi   = dx(j)
            rownrm = one
         else if (j .le. n+nclin) then

            i      = j - n
            axi    = ax(i)
            adxi   = adx(i)
            rownrm = anorm(i) + one
         else

            i      = j - n - nclin
            axi    = slk(i)
            adxi   = dslk(i)
            rownrm = one
         end if

         res = - one
         if (     adxi .le. - epspt9*rownrm*dxnorm) then

*           Constraint decreasing.

            adxi = - adxi
            if (bl(j) .gt. -bigbnd) res = axi   - bl(j)
         else if (adxi .gt.   epspt9*rownrm*dxnorm) then

*           Constraint increasing.

            if (bu(j) .lt.  bigbnd) res = bu(j) - axi
         end if

         if (res .gt. zero  .and.  alfa*adxi .gt. res)
     $      alfa  = res / adxi

         if (cmdbg  .and.  icmdbg(3) .gt. 0)
     $      write (nout, 1200) j, res, adxi, alfa

         j = j + 1
         go to 100
*+    end while
      end if

*     ==================================================================
*     Determine alfa, the bound on the step to be taken.
*     ==================================================================
      alfa   = max( alfa, alfmin )

      inform = 0
      if (alfa .ge. alfmax) inform = 1

      if (cmdbg  .and.  icmdbg(1) .gt. 0  .and.  inform .gt. 0)
     $   write (nout, 9010) alfa

      return

 1000 format(/ ' npalf  entered'
     $       / '    j            res             Ap           alfa '/)
 1200 format( i5, 3g15.5 )
 9010 format(/ ' //npalf //  No finite step.'
     $       / ' //npalf //             alfa'
     $       / ' //npalf //  ', g15.4 )

*     End of  npalf .
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npchkd( inform, msgnp, nstate, lvlder, nfun, ngrad,
     $                   ldcj, ldcju, n, ncnln,
     $                   confun, objfun, needc,
     $                   bigbnd, epsrf, cdint, fdint,
     $                   fdchk, fdnorm, objf, xnorm,
     $                   bl, bu, c, c1, cjac, cjacu, cjdx,
     $                   dx, grad, gradu, hforwd, hcntrl,
     $                   x, wrk1, wrk2, w, lenw )

      implicit           double precision (a-h,o-z)
      integer            needc(*)
      double precision   c(*), c1(*), cjac(ldcj,*), cjacu(ldcju,*),
     $                   cjdx(*)
      double precision   bl(n), bu(n), dx(n), grad(n), gradu(n), x(n)
      double precision   hforwd(*), hcntrl(*)
      double precision   wrk1(n), wrk2(*), w(lenw)
      external           confun, objfun

C***********************************************************************
C     npchkd  performs the following...
C     (1)  Computes the objective and constraint values objf and c.
C     (2)  Evaluates the user-provided gradients in cjacu and gradu.
C     (3)  Counts the missing gradients.
C     (4)  Loads the known gradients into grad and cjac.
C     (5)  Checks that the known gradients are programmed correctly.
C     (6)  Computes the missing gradient elements.
C
C     Systems Optimization Laboratory, Stanford University, California.
C     Original version written 4-September-1985.
C     This version of npchkd dated 26-Nov-1989.
C***********************************************************************
      common    /sol1cm/ nout
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset
      common    /sol5np/ lvrfyc, jverfy(4)

      logical            centrl, needfd
      parameter        ( rdummy =-11111.0d+0)

      infog  = 0
      infocj = 0
      nfdiff = 0
      ncdiff = 0
      ncset  = n*ncnln

      if (ncnln .gt. 0) then
*        ===============================================================
*        Compute the constraints and Jacobian matrix.
*        ===============================================================
*        If some derivatives are missing, load the Jacobian with dummy
*        values.  Any elements left unaltered after the call to confun
*        must be estimated.  A record of the missing jacobian elements
*        is stored in  cjacu.

         needfd = lvlder .eq. 0  .or.  lvlder .eq. 1

         if (needfd)
     $      call f06qhf( 'General', ncnln, n, rdummy, rdummy,
     $                   cjacu, ldcju )

         call iload ( ncnln, (1), needc, 1 )

         mode   = 2
         call confun( mode, ncnln, n, ldcju,
     $                needc, x, c, cjacu, nstate )
         if (mode .lt. 0) go to 9999

         call f06qff( 'General', ncnln, n, cjacu, ldcju, cjac, ldcj )

         if (needfd) then

*           Count the number of missing Jacobian elements.

            do 220, j = 1, n
               do 210, i = 1, ncnln
                  if (cjacu(i,j) .eq. rdummy) ncdiff = ncdiff + 1
  210          continue
  220       continue

            ncset = ncset - ncdiff
            if (nstate .eq. 1) then
               if (ncdiff .eq. 0) then
                  if (lvlder .eq. 0) lvlder = 2
                  if (lvlder .eq. 1) lvlder = 3
                  if (msgnp .gt. 0)
     $               write (nout, 1000) lvlder
               else
                  if (msgnp .gt. 0)
     $               write (nout, 1100) ncset, n*ncnln, ncdiff
               end if
            end if
         end if
      end if

*     ==================================================================
*     Repeat the procedure above for the objective function.
*     ==================================================================
      needfd = lvlder .eq. 0  .or.  lvlder .eq. 2

      if (needfd)
     $   call dload ( n, rdummy, gradu, 1 )

      mode  = 2
      call objfun( mode, n, x, objf, gradu, nstate )
      if (mode .lt. 0) go to 9999

      call dcopy ( n, gradu, 1, grad, 1 )

      if (needfd) then

*        Count the number of missing gradient elements.

         do 300, j = 1, n
            if (gradu(j) .eq. rdummy) nfdiff = nfdiff + 1
  300    continue

         if (nstate .eq. 1) then
            if (nfdiff .eq. 0) then
               if (lvlder .eq. 0) lvlder = 1
               if (lvlder .eq. 2) lvlder = 3
               if (msgnp .gt. 0)
     $            write (nout, 2000) lvlder
            else
               if (msgnp .gt. 0)
     $            write (nout, 2100) n - nfdiff, n, nfdiff
            end if
         end if
      end if

      nfun  = nfun  + 1
      ngrad = ngrad + 1

*     ==================================================================
*     Check whatever gradient elements have been provided.
*     ==================================================================
      if (lvrfyc .ge. 0) then
         if (ncset .gt. 0) then
            call chcjac( mode, lvlder, msgnp,
     $                   ncset, n, ncnln, ldcj, ldcju,
     $                   bigbnd, epsrf, epspt3, fdchk, xnorm,
     $                   confun, needc,
     $                   bl, bu, c, c1, cjac, cjacu, cjdx,
     $                   dx, wrk2, x, wrk1, w, lenw )
            if (mode .lt. 0) go to 9999
            infocj = mode
         end if

         if (nfdiff .lt. n) then
            call chfgrd( mode, msgnp, n,
     $                   bigbnd, epsrf, epspt3, fdchk, objf, xnorm,
     $                   objfun,
     $                   bl, bu, grad, gradu, dx, x, wrk1, w, lenw )
            if (mode .lt. 0) go to 9999
            infog   = mode
         end if
      end if

      needfd = ncdiff .gt. 0  .or.  nfdiff .gt. 0
      if (needfd) then
*        ===============================================================
*        Compute the missing gradient elements.
*        ===============================================================
*        chfd computes the finite-difference intervals or checks
*        preassigned (scalar) values.

         call chfd  ( mode, msgnp, lvlder,
     $                n, ncnln, ldcj, ldcju,
     $                bigbnd, epsrf, fdnorm, objf,
     $                objfun, confun, needc,
     $                bl, bu, c, c1, cjdx, cjac, cjacu,
     $                grad, gradu, hforwd, hcntrl, x,
     $                dx, w, lenw )

         if (mode .lt. 0) go to 9999

         if (lfdset .gt. 0) then
            centrl = lvldif .eq. 2
            call npfd  ( centrl, mode,
     $                   ldcj, ldcju, n, ncnln,
     $                   bigbnd, cdint, fdint, fdnorm, objf,
     $                   confun, objfun, needc,
     $                   bl, bu, c, c1, cjdx, cjac, cjacu,
     $                   grad, gradu, hforwd, hcntrl, x,
     $                   w, lenw )

            if (mode .lt. 0) go to 9999
         end if
      end if

      inform = infocj + infog

      return

*     The user requested termination.

 9999 inform = mode
      return

 1000 format(//' All Jacobian elements have been set.  ',
     $         ' Derivative level increased to ', i4 )
 1100 format(//' The user sets ', i6, '   out of', i6,
     $         '   Jacobian elements.'
     $       / ' Each iteration, ', i6,
     $         '   Jacobian elements will be estimated numerically.' )
 2000 format(//' All objective gradient elements have been set.  ',
     $         ' Derivative level increased to ', i4 )
 2100 format(//' The user sets ', i6, '   out of', i6,
     $         '   objective gradient elements.'
     $       / ' Each iteration, ', i6,
     $         '   gradient elements will be estimated numerically.' )

*     End of  npchkd.
      end
