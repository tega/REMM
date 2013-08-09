*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npcore( named, names, unitq, inform, majits,
     $                   n, nclin, ncnln, nctotl, nactiv, nfree, nz,
     $                   lda, ldcj, ldcju, ldaqp, ldr,
     $                   nfun, ngrad, istate, kactiv, kx,
     $                   objf, fdnorm, xnorm, confun, objfun,
     $                   aqp, ax, bl, bu, c, cjac, cjacu, clamda,
     $                   featol, grad, gradu, r, x, iw, w, lenw )

      implicit           double precision (a-h,o-z)
      logical            named
      integer            istate(*), kactiv(n), kx(n)
      integer            iw(*)
      double precision   aqp(ldaqp,*), ax(*), bl(nctotl), bu(nctotl),
     $                   c(*), cjac(ldcj,*), cjacu(ldcju,*)
      double precision   clamda(nctotl), featol(nctotl), grad(n),
     $                   gradu(n), r(ldr,*), x(n)
      double precision   w(lenw)
      external           confun, objfun

      double precision   asize, dtmax, dtmin
      character*8        names(*)
               
************************************************************************
*     npcore  is the core routine for  npsol,  a sequential quadratic
*     programming (SQP) method for nonlinearly constrained optimization.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version      February-1982.
*     This version of NPCORE dated 22-Oct-91.
************************************************************************
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      common    /sol1cm/ nout
      common    /sol3cm/ lennam, ldt   , ncolt , ldq
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5cm/ asize , dtmax , dtmin
      common    /sol6cm/ rcndbd, rfrobn, drmax, drmin

      parameter         (lenls = 20)
      common    /sol1ls/ locls(lenls)

      parameter         (lennp = 35)
      common    /sol1np/ locnp(lennp)
      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset
      common    /sol5np/ lvrfyc, jverfy(4)
      logical            incrun
      common    /sol6np/ rhomax, rhonrm, rhodmp, scale, incrun

      parameter         (ldbg = 5)
      logical            cmdbg, npdbg
      common    /npdebg/ inpdbg(ldbg), npdbg
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

      logical            goodgq, newgq
      logical            centrl, convrg, convpt, done  , error , feasqp
      logical            infeas, needfd, optiml, overfl, unitq
      logical            ktcond(2)
      intrinsic          abs   , max   , min   , mod   , real  , sqrt
      external           ddiv  , ddot  , dnrm2

      character*25       lsumry
      character*2        job
      parameter        ( job  = 'NP' )
      parameter        ( zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0 )
      parameter        ( growth=1.0d+2                              )

*     specify machine-dependent parameters.

      epsmch = wmach(3)
      flmax  = wmach(7)
      rtmax  = wmach(8)

      lanorm = locls( 2)
      lrpq   = locls( 5)
      lqrwrk = locls( 6)
      lhpq   = locls( 8)
      lgq    = locls( 9)
      lrlam  = locls(10)
      lt     = locls(11)
      lq     = locls(12)
      lwtinf = locls(13)
      lwrk1  = locls(14)
      lqptol = locls(15)

      liperm = locnp( 2)
      laqp   = locnp( 3)
      ladx   = locnp( 4)
      lbl    = locnp( 5)
      lbu    = locnp( 6)
      ldx    = locnp( 7)
      lgq1   = locnp( 8)
      lx1    = locnp(11)
      lwrk2  = locnp(12)
      lcs1   = locnp(13)
      lcs2   = locnp(14)
      lc1mul = locnp(15)
      lcmul  = locnp(16)
      lcjdx1 = locnp(17)
      ldlam  = locnp(18)
      ldslk  = locnp(19)
      lrho   = locnp(20)
      lwrk3  = locnp(21)
      lslk1  = locnp(22)
      lslk   = locnp(23)
      lneedc = locnp(24)
      lhfrwd = locnp(25)
      lhctrl = locnp(26)

      lcjac1 = laqp   + nclin
      lcjdx  = ladx   + nclin
      lvioln = lwrk3

*     Initialize

      lsumry = '                         '
      nqpinf = 0

      majit0 = majits
      nplin  = n     + nclin
      ncqp   = nclin + ncnln
      nl     = min( nplin + 1, nctotl )

      ldcj1 = max( ncqp , 1 )

      needfd = lvlder .eq. 0  .or.  lvlder .eq. 2
     $                        .or. (lvlder .eq. 1  .and.  ncnln .gt. 0)

      alfa   = zero
      alfdx  = zero
      rtftol = sqrt( ftol )
      rootn  = sqrt( real(n) )

*     If debug printing is required,  turn off any extensive printing
*     until iteration  IDBG.

      msgsv1 = msgnp
      msgsv2 = msgqp
      if (idbg .le. nmajor  .and.  idbg .gt. 0) then
         msgnp = 0
         if (msgsv1 .ge. 5) msgnp = 5
         msgqp = 0
         if (msgsv2 .ge. 5) msgqp = 5
      end if

*     ------------------------------------------------------------------
*     Information from the feasibility phase will be used to generate a
*     hot start for the first QP subproblem.
*     ------------------------------------------------------------------
      call dcopy ( nctotl, featol, 1, w(lqptol), 1 )

      nstate = 0

      objalf = objf
      if (ncnln .gt. 0) then
         objalf = objalf - ddot  ( ncnln, w(lcmul), 1, c, 1 )
      end if

      newgq  = .false.

**    ==================================================================
*+    repeat                             (until converged or error exit)

*        ===============================================================
*        See if we want to save the details of this iteration.
*        ===============================================================
  100    if (mod(majits,ksave) .eq. 0 .and. majits .ne. majit0) then
            call npsavr( unitq, n, nclin, ncnln, ldr, ldq,
     $                   nfree, nsave, majits, istate, kx,
     $                   w(lhfrwd), w(lhctrl),
     $                   w(lcmul), r, w(lrho), x, w(lq) )
         end if
         
**       ===============================================================
*+       repeat                         (Until a good gradient is found)

  110       centrl = lvldif .eq. 2

            if (newgq) then
               if (needfd) then
*                 ------------------------------------------------------
*                 compute any missing gradient elements and the
*                 transformed gradient of the objective.
*                 ------------------------------------------------------
                  call npfd  ( centrl, mode,
     $                         ldcj, ldcju, n, ncnln,
     $                         bigbnd, cdint, fdint, fdnorm, objf,
     $                         confun, objfun, iw(lneedc),
     $                         bl, bu, c, w(lwrk2), w(lwrk3),cjac,cjacu,
     $                         grad, gradu, w(lhfrwd), w(lhctrl), x,
     $                         w, lenw )
                  inform = mode
                  if (mode .lt. 0) go to 800

               end if

               call dcopy ( n, grad, 1, w(lgq), 1 )
               call cmqmul( 6, n, nz, nfree, ldq, unitq,
     $                      kx, w(lgq), w(lq), w(lwrk1) )

               newgq  = .false.
            end if

*           ============================================================
*           (1) Solve an inequality quadratic program (IQP) for the
*               search direction and multiplier estimates.
*           (2) For each nonlinear inequality constraint,  compute
*               the slack variable for which the merit function is
*               minimized.
*           (3) Compute the search direction for the slack variables
*               and multipliers.
*
*           Note that the array VIOLN is WRK3.
*           ============================================================
            call npiqp ( feasqp, unitq, nqperr, minits,
     $                   n, nclin, ncnln, lda, ldcj, ldaqp, ldr,
     $                   linact, nlnact, nactiv, nfree, nz, numinf,
     $                   istate, kactiv, kx,
     $                   dxnorm, gdx, qpcurv,
     $                   aqp, w(ladx), w(lanorm), ax, bl, bu,
     $                   c, cjac, clamda, w(lcmul), w(lcs1),
     $                   w(ldlam), w(ldslk), w(ldx), w(lbl), w(lbu),
     $                   w(lqptol), r, w(lrho), w(lslk), w(lvioln), x,
     $                   w(lwtinf), iw, w )

            if (feasqp) then
               nqpinf = 0
            else
               nqpinf = nqpinf + 1
               lsumry(2:22) = 'Infeasible subproblem'
            end if

*           ============================================================
*           Compute quantities needed for the convergence test.
*           ============================================================
*           Compute the norms of the projected gradient and the
*           gradient with respect to the free variables.

            gznorm = zero
            if (nz .gt. 0)
     $         gznorm = dnrm2 ( nz   , w(lgq), 1 )
            gfnorm = gznorm
            if (nfree .gt. 0  .and.  nactiv .gt. 0)
     $         gfnorm = dnrm2 ( nfree, w(lgq), 1 )

*           If the forward-difference estimate of the transformed
*           gradient of the Lagrangian function is small,  switch to
*           central differences, recompute the derivatives and re-solve
*           the QP.

            goodgq = .true.
            if (needfd  .and.  .not. centrl) then

               glnorm = dnrm2 ( n, w(lhpq), 1 )
               if (ncnln .eq. 0) then
                  cnorm = zero
               else
                  cnorm = dnrm2 ( ncnln, c, 1 )
               end if

               gltest = (one + abs(objf) + abs(cnorm))*epsrf/fdnorm
               if (glnorm .le. gltest) then
                  goodgq       = .false.
                  lsumry(3:21) = 'Central differences'
                  lvldif       = 2
                  newgq        = .true.
               end if

            end if

*+       until     (goodgq)
         if (.not.  goodgq ) go to 110

*        ===============================================================
*        (1) Compute the number of constraints that are violated by more
*            than FEATOL.
*        (2) Compute the 2-norm of the residuals of the constraints in
*            the QP working set.
*        ===============================================================
         call npfeas( n, nclin, ncnln, istate,
     $                bigbnd, cvnorm, errmax, jmax, nviol,
     $                ax, bl, bu, c, featol, x, w(lwrk2) )

*        Define small quantities that reflect the magnitude of OBJF and
*        the norm of GRAD(free).

         objsiz = one + abs( objf )
         xsize  = one +  xnorm
         gtest  = max( objsiz, gfnorm )
         dinky  = rtftol * gtest

         if (nactiv .eq. 0) then
            condt = zero
         else if (nactiv .eq. 1) then
            condt = dtmin
         else
            condt = ddiv  ( dtmax, dtmin, overfl )
         end if

         call dcond ( n, r, ldr+1, drmax, drmin )

         condh = ddiv  ( drmax, drmin, overfl )
         if (condh .lt. rtmax) then
            condh = condh*condh
         else
            condh = flmax
         end if

         if (nz .eq. 0) then
            condhz = one
         else if (nz .eq. n) then
            condhz = condh
         else
            call dcond ( nz, r, ldr+1, drzmax, drzmin )
            condhz = ddiv  ( drzmax, drzmin, overfl )
            if (condhz .lt. rtmax) then
               condhz = condhz*condhz
            else
               condhz = flmax
            end if
         end if

*        ---------------------------------------------------------------
*        Test for convergence.
*        The point test CONVPT checks for a K-T point at the initial
*        point or after a large change in X.
*        ---------------------------------------------------------------
         convpt    = gznorm .le. epspt8*gtest  .and.  nviol  .eq. 0

         ktcond(1) = gznorm .lt. dinky
         ktcond(2) = nviol  .eq. 0
         optiml    = ktcond(1)  .and.  ktcond(2)

         convrg    = majits .gt. 0  .and.  alfdx .le. rtftol*xsize

         infeas    =       convrg         .and.  .not. feasqp
     $               .or.  nqpinf .gt. 7

         done      = convpt  .or.  (convrg  .and. optiml)
     $                       .or.   infeas

         objalf = objf
         grdalf = gdx
         glf1   = gdx
         if (ncnln .gt. 0) then
            glf1   = glf1
     $                 - ddot( ncnln, w(lcjdx), 1, clamda(nl), 1 )

*           Compute the value and directional derivative of the
*           augmented Lagrangian merit function.
*           The penalty parameters may be increased or decreased.

            call npmrt ( feasqp, n, nclin, ncnln,
     $                   objalf, grdalf, qpcurv,
     $                   istate,
     $                   w(lcjdx), w(lcmul), w(lcs1),
     $                   w(ldlam), w(lrho), w(lvioln),
     $                   w(lwrk1), w(lwrk2) )
         end if

*        ===============================================================
*        Print the details of this iteration.
*        ===============================================================
         call npprt ( ktcond, convrg, lsumry, msgnp, msgqp,
     $                ldr, ldt, n, nclin, ncnln,
     $                nctotl, nactiv, linact, nlnact, nz, nfree,
     $                majit0, majits, minits, istate, alfa, nfun,
     $                condhz, condh, condt, objalf, objf,
     $                gfnorm, gznorm, cvnorm,
     $                ax, c, r, w(lt), w(lvioln), x, w(lwrk1) )

         alfa  = zero
         error = majits .ge. nmajor

         if (.not. (done  .or.  error)) then
            majits = majits + 1

            if (majits .eq. idbg) then
               npdbg = .true.
               cmdbg =  npdbg
               msgnp =  msgsv1
               msgqp =  msgsv2
            end if

*           Make copies of information needed for the BFGS update.

            call dcopy ( n, x     , 1, w(lx1) , 1 )
            call dcopy ( n, w(lgq), 1, w(lgq1), 1 )

            if (ncnln .gt. 0) then
               call dcopy ( ncnln, w(lcjdx), 1, w(lcjdx1), 1 )
               call dcopy ( ncnln, w(lcmul), 1, w(lc1mul), 1 )
               call dcopy ( ncnln, w(lslk) , 1, w(lslk1) , 1 )
            end if

*           ============================================================
*           Compute the parameters for the linesearch.
*           ============================================================
*           alfmin is the smallest allowable step predicted by the QP
*           subproblem.

            alfmin = one
            if (.not. feasqp) alfmin = zero

*           ------------------------------------------------------------
*           alfmax is the largest feasible steplength subject to a user-
*           defined limit alflim on the change in x.
*           ------------------------------------------------------------
            if (ncnln .gt. 0  .and.  needfd) then
               alfmax = one
            else
               alfmax = ddiv  ( bigdx, dxnorm, overfl )
               call npalf ( info, n, nclin, ncnln,
     $                      alfa, alfmin, alfmax, bigbnd, dxnorm,
     $                      w(lanorm), w(ladx), ax, bl, bu,
     $                      w(ldslk), w(ldx), w(lslk), x )
               alfmax = alfa
               if (alfmax .lt. one + epspt3  .and.  feasqp)
     $            alfmax = one
            end if

*           ------------------------------------------------------------
*           alfbnd is a tentative upper bound on the steplength.  If the
*           merit function is decreasing at ALFBND and certain
*           conditions hold,  ALFBND will be increased in multiples of
*           two (subject to not being greater than ALFMAX).
*           ------------------------------------------------------------
            if (ncnln .eq. 0) then
               alfbnd = alfmax
            else
               alfbnd = min( one, alfmax )
            end if

*           ------------------------------------------------------------
*           alfsml trips the computation of central differences.  If a
*           trial steplength falls below ALFSML, the linesearch is
*           terminated.
*           ------------------------------------------------------------
            alfsml = zero
            if (needfd  .and. .not. centrl) then
               alfsml = ddiv  ( fdnorm, dxnorm, overfl )
               alfsml = min   ( alfsml, alfmax )
            end if

*           ============================================================
*           Compute the steplength using safeguarded interpolation.
*           ============================================================
            alflim = ddiv ( (one+xnorm)*dxlim, dxnorm, overfl )
            alfa   = min  ( alflim, one )

            call npsrch( needfd, nlserr, n, ncnln,
     $                   ldcj, ldcju, nfun, ngrad,
     $                   iw(lneedc), confun, objfun,
     $                   alfa, alfbnd, alfmax, alfsml, dxnorm,
     $                   epsrf, eta, gdx, grdalf, glf1, glf2,
     $                   objf, objalf, qpcurv, xnorm,
     $                   c, w(lwrk1), cjac, cjacu, w(lcjdx), w(lwrk3),
     $                   w(lc1mul), w(lcmul), w(lcs1),
     $                   w(lcs2), w(ldx), w(ldlam), w(ldslk), grad,
     $                   gradu, clamda(nl), w(lrho),
     $                   w(lslk1), w(lslk), w(lx1), x,
     $                   w(lwrk2), w, lenw )

*           ------------------------------------------------------------
*           npsrch  sets nlserr to the following values...
*
*           < 0  if the user wants to stop.
*             1  if the search is successful and alfa < alfmax.
*             2  if the search is successful and alfa = alfmax.
*             3  if a better point was found but too many functions
*                were needed (not sufficient decrease).
*
*           Values of nlserr occurring with a nonzero value of alfa.
*             4  if alfmax < tolabs (too small to do a search).
*             5  if alfa  < alfsml (srchq only -- maybe want to switch
*                to central differences to get a better direction).
*             6  if the search found that there is no useful step.
*                The interval of uncertainty is less than 2*tolabs.
*                The minimizer is very close to alfa = zero
*                or the gradients are not sufficiently accurate.
*             7  if there were too many function calls.
*             8  if the input parameters were bad
*                (alfmax le toltny  or  uphill).
*           ------------------------------------------------------------
            if (nlserr .lt. 0) then
               inform = nlserr
               go to 800
            end if

            if (alfa .gt. alflim) lsumry(4:4) = 'L'

            error  = nlserr .ge. 4
            if (error) then
*              ---------------------------------------------------------
*              The linesearch failed to find a better point.
*              If exact gradients or central differences are being used,
*              or the KT conditions are satisfied, stop.  Otherwise,
*              switch to central differences and solve the QP again.
*              ---------------------------------------------------------
               if (needfd  .and.  .not. centrl) then
                  if (.not. optiml) then
                     error        = .false.
                     lsumry(3:21) = 'Central differences'
                     lvldif       = 2
                     newgq        = .true.
                  end if
               end if
            else
               if (needfd) then
*                 ======================================================
*                 Compute the missing gradients.
*                 ======================================================
                  mode  = 1
                  ngrad = ngrad + 1

                  if (ncnln .gt. 0) then
                     call iload ( ncnln, (1), iw(lneedc), 1 )

                     call confun( mode, ncnln, n, ldcju, iw(lneedc),
     $                            x, w(lwrk1), cjacu, nstate )
                     inform = mode
                     if (mode .lt. 0) go to 800

                     call f06qff( 'General', ncnln, n, cjacu, ldcju,
     $                            cjac, ldcj )
                  end if

                  call objfun( mode, n, x, obj, gradu, nstate )
                  inform = mode
                  if (mode .lt. 0) go to 800

                  call dcopy ( n, gradu, 1, grad, 1 )

                  call npfd  ( centrl, mode,
     $                         ldcj, ldcju, n, ncnln,
     $                         bigbnd, cdint, fdint, fdnorm, objf,
     $                         confun, objfun, iw(lneedc),
     $                         bl, bu, c, w(lwrk2), w(lwrk3),cjac,cjacu,
     $                         grad, gradu, w(lhfrwd), w(lhctrl), x,
     $                         w, lenw )

                  inform = mode
                  if (mode .lt. 0) go to 800

                  gdx  =  ddot( n, grad, 1, w(ldx), 1 )
                  glf2 =  gdx
                  if (ncnln .gt. 0) then
                     call dgemv ( 'N', ncnln, n, one, cjac, ldcj,
     $                            w(ldx), 1, zero, w(lcjdx), 1 )
                     glf2 = glf2 -
     $                      ddot( ncnln, w(lcjdx), 1, clamda(nl), 1 )
                  end if
               end if

               call dcopy ( n, grad, 1, w(lgq), 1 )
               call cmqmul( 6, n, nz, nfree, ldq, unitq,
     $                      kx, w(lgq), w(lq), w(lwrk1) )

               xnorm  = dnrm2 ( n, x, 1 )

               if (ncnln .gt. 0  .and.  alfa .ge. one)
     $            call dcopy ( ncnln, clamda(nl), 1, w(lcmul), 1 )

               if (nclin .gt. 0)
     $            call daxpy ( nclin, alfa, w(ladx), 1, ax, 1 )
               alfdx   = alfa * dxnorm

*              =========================================================
*              Update the factors of the approximate Hessian of the
*              Lagrangian function.
*              =========================================================
               call npupdt( lsumry, unitq,
     $                      n, ncnln, nfree, nz,
     $                      ldcj1, ldcj, ldq, ldr, kx,
     $                      alfa, glf1, glf2, qpcurv,
     $                      w(lcjac1), cjac, w(lcjdx1), w(lcjdx),
     $                      w(lcs1), w(lcs2), w(lgq1), w(lgq),
     $                      w(lhpq), w(lrpq), clamda(nl), r,
     $                      w(lwrk3), w(lq), w(lwrk2), w(lwrk1) )

               call dcond ( n, r, ldr+1, drmax, drmin )
               cond   = ddiv  ( drmax, drmin, overfl )

               if (      cond   .gt. rcndbd
     $             .or.  rfrobn .gt. rootn*growth*drmax) then
*                 ------------------------------------------------------
*                 reset the condition estimator and range-space
*                 partition of Q'HQ.
*                 ------------------------------------------------------
                  if (npdbg  .and.  inpdbg(1) .gt. 0)
     $               write (nout, 9000) rfrobn, drmax, drmin,cond,rcndbd

                  LSUMRY(5:23) = 'Refactorize Hessian'

                  call nprset( unitq,
     $                         n, nfree, nz, ldq, ldr,
     $                         iw(liperm), kx,
     $                         w(lgq), r, w(lq), w(lwrk1), w(lqrwrk) )
               end if
            end if
         end if

*+    until     (done  .or.  error)
      if (.not. (done  .or.  error) ) go to 100

*     ======================end of main loop============================

      if (done) then
         if (convpt  .or.  optiml) then
            inform = 0
         else if (infeas) then
            inform = 3
         end if
      else if (error) then
         if (majits .ge. nmajor) then
            inform = 4
         else if (optiml) then
            inform = 1
         else
            inform = 6
         end if
      end if

*     ------------------------------------------------------------------
*     Set  clamda.  Print the full solution.
*     ------------------------------------------------------------------
  800 msgnp = msgsv1
      msgqp = msgsv2
      if (msgnp .gt. 0)
     $   write (nout, 2100) inform, majits, nfun, ngrad

      call cmprt ( msgnp, nfree, ldaqp,
     $             n, nclin, ncnln, nctotl, bigbnd,
     $             named, names, lennam,
     $             nactiv, istate, kactiv, kx,
     $             aqp, bl, bu, c, clamda, w(lrlam), x )
      if (ncnln .gt. 0)
     $call dcopy ( ncnln, w(lcmul), 1, clamda(n+nclin+1), 1 )

      return

 2100 format(/ ' Exit  NP phase.  inform = ', i2, ' majits = ', i5,
     $         '   nfun = ', i5, '   ngrad = ', i5 )

 9000 format(/ ' //npcore//        rfrobn         drmax         drmin'
     $       / ' //npcore//', 1p, 3e14.2
     $       / ' //npcore//          cond        rcndbd'
     $       / ' //npcore//', 1p, 2e14.2 )

*     end of  npcore.
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npcrsh( cold, n, nclin, ncnln,
     $                   nctotl, nactiv, nfree, nz,
     $                   istate, kactiv, bigbnd, tolact,
     $                   bl, bu, c )

      implicit           double precision (a-h,o-z)
      logical            cold
      integer            istate(nctotl), kactiv(n)
      double precision   c(*), bl(nctotl), bu(nctotl)

************************************************************************
*     npcrsh  adds indices of nonlinear constraints to the initial
*     working set.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version   14-February 1985.
*     This version of  NPCRSH  dated 14-November-1985.
************************************************************************
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      common    /sol1cm/ nout

      logical            npdbg
      parameter        ( ldbg = 5 )
      common    /npdebg/ inpdbg(ldbg), npdbg

      intrinsic          abs, min
      parameter        ( one = 1.0d+0 )

      nfixed = n      - nfree
      linact = nactiv
      nplin  = n      + nclin

*     If a cold start is being made, initialize the status of the QP
*     working set.  First,  if  BL(j) = BU(j),  set ISTATE(j)=3.

      if (cold) then
         do 130, j = nplin+1, nctotl
            istate(j) = 0
            if (bl(j) .eq. bu(j)) istate(j) = 3
  130    continue
      end if

*     Increment NACTIV and KACTIV.
*     Ensure that the number of bounds and general constraints in the
*     QP  working set does not exceed N.

      do 200, j = nplin+1, nctotl
         if (nfixed + nactiv .eq. n) istate(j) = 0
         if (istate(j) .gt. 0) then
            nactiv = nactiv + 1
            kactiv(nactiv) = j - n
         end if
  200 continue

      if (cold) then

*        ---------------------------------------------------------------
*        If a cold start is required, an attempt is made to add as many
*        nonlinear constraints as possible to the working set.
*        ---------------------------------------------------------------
*        The following loop finds the most violated constraint.  If
*        there is room in kactiv, it will be added to the working set
*        and the process will be repeated.

         is     =   1
         biglow = - bigbnd
         bigupp =   bigbnd
         toobig =   tolact + tolact

*        while (is .gt. 0  .and.  nfixed + nactiv .lt. n) do
  500    if    (is .gt. 0  .and.  nfixed + nactiv .lt. n) then
            is   = 0
            cmin = tolact

            do 520, i = 1, ncnln
               j      = nplin + i
               if (istate(j) .eq. 0) then
                  b1     = bl(j)
                  b2     = bu(j)
                  resl   = toobig
                  resu   = toobig
                  if (b1 .gt. biglow)
     $            resl   = abs( c(i) - b1 ) / (one + abs( b1 ))
                  if (b2 .lt. bigupp)
     $            resu   = abs( c(i) - b2 ) / (one + abs( b2 ))
                  res    = min( resl, resu )
                  if (res .lt. cmin) then
                     cmin = res
                     imin = i
                     is   = 1
                     if (resl .gt. resu) is = 2
                  end if
               end if
  520       continue

            if (is .gt. 0) then
               nactiv         = nactiv + 1
               kactiv(nactiv) = nclin  + imin
               j              = nplin  + imin
               istate(j)      = is
            end if
            go to 500
*        end while
         end if
      end if

*     ------------------------------------------------------------------
*     An initial working set has now been selected.
*     ------------------------------------------------------------------
      nlnact = nactiv - linact
      nz     = nfree  - nactiv
      if (npdbg  .and.  inpdbg(1) .gt. 0)
     $   write (nout, 1000) nfixed, linact, nlnact

      return

 1000 format(/ ' //npcrsh//  Working set selected....'
     $       / ' //npcrsh// nfixed linact nlnact     '
     $       / ' //npcrsh//', 3i7 )

*     end of  npcrsh
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE NPDFLT( N, NCLIN, NCNLN, LENIW, LENW, TITLE )

      IMPLICIT           DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*(*)      TITLE

************************************************************************
*     NPDFLT  loads the default values of parameters not set in the
*     options file.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 77 version written 10-September-1985.
*     This version of NPDFLT dated 26-Nov-89.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/

      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
                 
      COMMON    /SOL4NP/ LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON    /SOL5NP/ LVRFYC, JVERFY(4)

      LOGICAL            CMDBG, LSDBG, NPDBG
      PARAMETER        ( LDBG = 5 )
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
      COMMON    /CMDEBG/ ICMDBG(LDBG), CMDBG

      LOGICAL            NEWOPT
      COMMON    /SOL7NP/ NEWOPT
      SAVE      /SOL7NP/

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
*-----INCLUDE NPPARM----------------------------------------------------
      INTEGER            IPRMNP(MXPARM), IPSVNP
      DOUBLE PRECISION   RPRMNP(MXPARM), RPSVNP

      COMMON    /NPPAR1/ IPSVNP(MXPARM),
     $                   IDBGNP, ITMXNP, JVRFY1, JVRFY2, JVRFY3, JVRFY4,
     $                   LDBGNP, LFORMH, LVLDER, LVERFY, MSGNP , NLNF  ,
     $                   NLNJ  , NLNX  , NNCNLN, NSAVE , NLOAD , KSAVE ,
     $                   IPADNP(12)

      COMMON    /NPPAR2/ RPSVNP(MXPARM),
     $                   CDINT , CTOL  , DXLIM , EPSRF , ETA   , FDINT ,
     $                   FTOL  , HCNDBD, RPADNP(22)

      EQUIVALENCE       (IPRMNP(1), IDBGNP), (RPRMNP(1), CDINT)

      SAVE      /NPPAR1/, /NPPAR2/
*-----------------------------------------------------------------------
      EQUIVALENCE  (IDBGNP, IDBG  ), (ITMXNP, NMAJOR), (ITMAX2, NMINOR)
      EQUIVALENCE  (LDBGLS, MNRDBG), (LDBGNP, MJRDBG), (MSGLS , MSGQP )

      INTRINSIC          ABS    , LEN    , MAX   , MOD   , REAL
      PARAMETER        ( ZERO   =  0.0D+0, ONE    =  1.0D+0 )
      PARAMETER        ( POINT3 =  3.3D-1, POINT8 =  0.8D+0 )
      PARAMETER        ( POINT9 =  0.9D+0                   )
      PARAMETER        ( TWO    =  2.0D+0, TEN    = 10.0D+0 )
      PARAMETER        ( TENP6  =  1.0D+6, HUNDRD = 10.0D+1 )
      PARAMETER        ( RDUMMY = -11111., IDUMMY = -11111  )
      PARAMETER        ( GIGANT =  1.0D+20*.99999           )
      PARAMETER        ( WRKTOL =  1.0D-2                   )

      CHARACTER*4        ICRSH(0:2)
      CHARACTER*16       KEY
      DATA                ICRSH(0),  ICRSH(1),  ICRSH(2)
     $                 / 'COLD'   , 'WARM'   , 'HOT '    /

      EPSMCH = WMACH( 3)
      NOUT   = WMACH(11)

      CONDBD = MAX ( ONE/(HUNDRD*EPSMCH*REAL(N)), TENP6 )

      NCQP   = NCLIN + NCNLN
      NPLIN  = N     + NCLIN
      NCTOTL = NPLIN + NCNLN

*     Make a dummy call NPKEY to ensure that the defaults are set.

      CALL NPKEY ( NOUT, '*', KEY )
      NEWOPT = .TRUE.

*     Save the optional parameters set by the user.  The values in
*     IPRMLS, RPRMLS, IPRMNP and RPRMNP may be changed to their
*     default values.

      CALL ICOPY ( MXPARM, IPRMLS, 1, IPSVLS, 1 )
      CALL DCOPY ( MXPARM, RPRMLS, 1, RPSVLS, 1 )
      CALL ICOPY ( MXPARM, IPRMNP, 1, IPSVNP, 1 )
      CALL DCOPY ( MXPARM, RPRMNP, 1, RPSVNP, 1 )

      IF (          LCRASH .LT. 0
     $    .OR.      LCRASH .GT. 2     )   LCRASH  =  0
      IF (          LVLDER .LT. 0
     $    .OR.      LVLDER .GT. 3     )   LVLDER  =  3
      IF (          LFORMH .LT. 0
     $    .OR.      LFORMH .GT. 1     )   LFORMH  =  0

      IF (          NMAJOR .LT. 0     )   NMAJOR  = MAX(50, 3*NPLIN+
     $                                                     10*NCNLN )
      IF (          NMINOR .LT. 1     )   NMINOR  = MAX(50, 3*NCTOTL)
      IF (          MJRDBG .LT. 0     )   MJRDBG  =  0
      IF (          MNRDBG .LT. 0     )   MNRDBG  =  0
      IF (          IDBG   .LT. 0
     $    .OR.      IDBG   .GT. NMAJOR)   IDBG    =  0
      IF (          MJRDBG .EQ. 0
     $    .AND.     MNRDBG .EQ. 0     )   IDBG    = NMAJOR + 1
      IF (          MSGNP  .EQ. IDUMMY)   MSGNP   = 10
      IF (          MSGQP  .EQ. IDUMMY)   MSGQP   =  0
                                          NLNF    =  N
                                          NLNJ    =  N
                                          NLNX    =  N
      IF (          JVRFY2 .LT. 0
     $    .OR.      JVRFY2 .GT. N     )   JVRFY2  =  N
      IF (          JVRFY1 .LT. 0
     $    .OR.      JVRFY1 .GT. JVRFY2)   JVRFY1  =  1
      IF (          JVRFY4 .LT. 0
     $    .OR.      JVRFY4 .GT. N     )   JVRFY4  =  N
      IF (          JVRFY3 .LT. 0
     $    .OR.      JVRFY3 .GT. JVRFY4)   JVRFY3  =  1
      IF (          LVERFY .EQ. IDUMMY
     $    .OR.      LVERFY .GT. 13    )   LVERFY  =  0

      IF (          KSAVE  .LE. 0     )   KSAVE   =  NMAJOR + 1
      IF (          NSAVE  .LT. 0     )   NSAVE   =  0
      IF (          NSAVE  .EQ. 0     )   KSAVE   =  NMAJOR + 1
      IF (          NLOAD  .LT. 0     )   NLOAD   =  0
      IF (          LCRASH .LE. 1     )   NLOAD   =  0
      IF (          NLOAD  .EQ. 0 
     $    .AND.     LCRASH .EQ. 2     )   LCRASH  =  0

      IF (          TOLACT .LT. ZERO
     $    .OR.      TOLACT .GE. ONE   )   TOLACT  =  WRKTOL
      IF (          TOLFEA .LT. EPSMCH
     $    .OR.      TOLFEA .GE. ONE   )   TOLFEA  =  EPSPT5
      IF (          EPSRF  .LT. EPSMCH
     $    .OR.      EPSRF  .GE. ONE   )   EPSRF   =  EPSPT9
                                          LFDSET  =  0
      IF (          FDINT  .LT. ZERO  )   LFDSET  =  2
      IF (          FDINT  .EQ. RDUMMY)   LFDSET  =  0
      IF (          FDINT  .GE. EPSMCH
     $    .AND.     FDINT  .LT. ONE   )   LFDSET  =  1
      IF (          LFDSET .EQ. 1
     $    .AND.    (CDINT  .LT. EPSMCH
     $    .OR.      CDINT  .GE. ONE  ))   CDINT   = EPSRF**POINT3
      IF (          BIGBND .LE. ZERO  )   BIGBND  = GIGANT
      IF (          BIGDX  .LE. ZERO  )   BIGDX   = MAX( GIGANT,BIGBND )
      IF (          DXLIM  .LE. ZERO  )   DXLIM   = TWO
      IF (          ETA    .LT. ZERO
     $    .OR.      ETA    .GE. ONE   )   ETA     = POINT9
      IF (          FTOL   .LT. EPSRF
     $    .OR.      FTOL   .GE. ONE   )   FTOL    = EPSRF**POINT8

      IF (          HCNDBD .LT. ONE   )   HCNDBD  = CONDBD

                                          DCTOL   = EPSPT5
      IF (          LVLDER .LT. 2     )   DCTOL   = EPSPT3
      IF (          CTOL   .LT. EPSMCH
     $    .OR.      CTOL   .GE. ONE   )   CTOL    = DCTOL

      ITMAX1    = MAX( 50, 3*(N + NCLIN + NCNLN) )
      JVERFY(1) = JVRFY1
      JVERFY(2) = JVRFY2
      JVERFY(3) = JVRFY3
      JVERFY(4) = JVRFY4

      NPDBG = IDBG .EQ. 0
      CMDBG = NPDBG

      K     = 1
      MSG1  = MJRDBG
      MSG2  = MNRDBG
      DO 200 I = 1, LDBG
         INPDBG(I) = MOD( MSG1/K, 10 )
         ICMDBG(I) = INPDBG(I)
         ILSDBG(I) = MOD( MSG2/K, 10 )
         K = K*10
  200 CONTINUE

      IF (MSGNP .GT. 0) THEN

*        Print the title.  If no hot start is specified,  the parameters
*        are final and can be printed.

         LENT = LEN( TITLE )
         IF (LENT .GT. 0) THEN
            NSPACE = (81 - LENT)/2 + 1
            WRITE (NOUT, '(///// (80A1) )')
     $         (' ', J=1, NSPACE), (TITLE(J:J), J=1,LENT)
            WRITE (NOUT, '(80A1 //)')
     $         (' ', J=1, NSPACE), ('='       , J=1,LENT)
         END IF


         IF (LCRASH .LE. 1) THEN
            WRITE (NOUT, 2000)
            WRITE (NOUT, 2100) NCLIN , TOLFEA, ICRSH(LCRASH) ,
     $                         N     , BIGBND, TOLACT,
     $                         DXLIM , BIGDX
            WRITE (NOUT, 2200) NCNLN , FTOL  , EPSRF ,
     $                         NLNJ  , CTOL  ,
     $                         NLNF  , ETA   ,
     $                         EPSMCH,
     $                         LVLDER, LVERFY
            WRITE (NOUT, 2300) NMAJOR, MSGNP ,
     $                         NMINOR, MSGQP ,
     $                         NLOAD , NSAVE , KSAVE
            
            IF (LVLDER .LT. 3) THEN
               IF      (LFDSET .EQ. 0) THEN
                  WRITE (NOUT, 2400)
               ELSE IF (LFDSET .EQ. 1) THEN
                  WRITE (NOUT, 2401) FDINT, CDINT
               ELSE IF (LFDSET .EQ. 2) THEN
                  WRITE (NOUT, 2402)
               END IF
            END IF
         END IF
      END IF

      RETURN

 2000 FORMAT(
     $//' Parameters'
     $/ ' ----------' )
 2100 FORMAT(
     $/ ' Linear constraints.....', I10,     6X,
     $  ' Linear feasibility.....', 1PE10.2, 6X,
     $  1X, A4, ' start.............'
     $/ ' Variables..............', I10,     6X,
     $  ' Infinite bound size....', 1PE10.2, 6X,
     $  ' Crash tolerance........', 1PE10.2
     $/ ' Step limit.............', 1PE10.2, 6X,
     $  ' Infinite step size.....', 1PE10.2  )
 2200 FORMAT(
     $/ ' Nonlinear constraints..', I10,     6X,
     $  ' Optimality tolerance...', 1PE10.2, 6X,
     $  ' Function precision.....', 1PE10.2
     $/ ' Nonlinear Jacobian vars', I10,     6X,
     $  ' Nonlinear feasibility..', 1PE10.2
     $/ ' Nonlinear objectiv vars', I10,     6X,
     $  ' Linesearch tolerance...', 1PE10.2
     $/ ' EPS (machine precision)', 1PE10.2, 6X,
     $  ' Derivative level.......', I10,     6X,
     $  ' Verify level...........', I10)
 2300 FORMAT(
     $/ ' Major iterations limit.', I10, 6X,
     $  ' Major print level......', I10
     $/ ' Minor iterations limit.', I10, 6X,
     $  ' Minor print level......', I10
     $/ ' RUN loaded from file...', I10, 6X,
     $  ' RUN to be saved on file', I10, 6X,
     $  ' Save frequency.........', I10)

 2400 FORMAT(/ ' Difference intervals to be computed.' )
 2401 FORMAT(/ ' Difference interval....', 1PE10.2, 6X,
     $         ' Central diffce interval', 1PE10.2 )
 2402 FORMAT(/ ' User-supplied difference intervals.' )

*     End of  NPDFLT.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npfd  ( centrl, inform,
     $                   ldcj, ldcju, n, ncnln,
     $                   bigbnd, cdint, fdint, fdnorm, objf,
     $                   confun, objfun, needc,
     $                   bl, bu, c, c1, c2, cjac, cjacu,
     $                   grad, gradu, hforwd, hcntrl, x,
     $                   w, lenw )

      implicit           double precision (a-h,o-z)
      logical            centrl
      integer            needc(*)

      double precision   bl(n), bu(n), c(*), c1(*), c2(*),
     $                   cjac(ldcj,*), cjacu(ldcju,*)
      double precision   grad(n), gradu(n), hforwd(n), hcntrl(n), x(n)
      double precision   w(lenw)
      external           confun, objfun

************************************************************************
*     npfd   evaluates any missing gradients.
*
*     Systems Optimization Laboratory, Stanford University, California.
*     Original version written 3-July-1986.
*     This version of npfd   dated 14-July-1986.
************************************************************************

      common    /sol1cm/ nout
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset

      logical            npdbg
      parameter         (ldbg = 5)
      common    /npdebg/ inpdbg(ldbg), npdbg

      intrinsic          abs   , max

      parameter         (rdummy=-11111.0d+0)
      parameter         (zero  = 0.0d+0, half  = 0.5d+0, one   = 1.0d+0)
      parameter         (two   = 2.0d+0, three = 3.0d+0, four  = 4.0d+0)

      inform = 0

*     ==================================================================
*     Use the pre-assigned difference intervals to approximate the
*     derivatives.
*     ==================================================================
*     Use either the same interval for each component (lfdset = 1),
*     or the intervals already in hforwd or hcntrl (lfdset = 0 or 2).

      nstate =   0
      mode   =   0

      biglow = - bigbnd
      bigupp =   bigbnd

      fdnorm =   zero

      do 340, j = 1, n

         xj     = x(j)
         nfound = 0
         if (ncdiff .gt. 0) then
            do 310, i = 1, ncnln
               if (cjacu(i,j) .eq. rdummy) then
                  needc(i) = 1
                  nfound   = nfound + 1
               else
                  needc(i) = 0
               end if
  310       continue
         end if

         if (nfound .gt. 0  .or.  gradu(j) .eq. rdummy) then
            stepbl = biglow
            stepbu = bigupp
            if (bl(j) .gt. biglow) stepbl = bl(j) - xj
            if (bu(j) .lt. bigupp) stepbu = bu(j) - xj

            if (centrl) then
               if (lfdset .eq. 1) then
                  delta = cdint
               else
                  delta = hcntrl(j)
               end if
            else
               if (lfdset .eq. 1) then
                  delta = fdint
               else
                  delta = hforwd(j)
               end if
            end if

            delta  = delta*(one + abs(xj))
            fdnorm = max (fdnorm, delta)
            if (half*(stepbl + stepbu) .lt. zero) delta =  - delta

            x(j) = xj + delta
            if (nfound .gt. 0) then
               call confun( mode, ncnln, n, ldcju,
     $                      needc, x, c1, cjacu, nstate )
               if (mode .lt. 0) go to 999
            end if

            if (gradu(j) .eq. rdummy) then
               call objfun( mode, n, x, objf1, gradu, nstate )
               if (mode .lt. 0) go to 999
            end if

            if (centrl) then
*              ---------------------------------------------------------
*              Central differences.
*              ---------------------------------------------------------
               x(j)  = xj + delta + delta

               if (nfound .gt. 0) then
                  call confun( mode, ncnln, n, ldcju,
     $                         needc, x, c2, cjacu, nstate )
                  if (mode .lt. 0) go to 999

                  do 320, i = 1, ncnln
                     if (needc(i) .eq. 1)
     $                  cjac(i,j) = (four*c1(i) - three*c(i) - c2(i))
     $                                  / (delta + delta)
  320             continue
               end if

               if (gradu(j) .eq. rdummy) then
                  call objfun( mode, n, x, objf2, gradu, nstate )
                  if (mode .lt. 0) go to 999

                  grad(j) = (four*objf1 - three*objf - objf2)
     $                                  / (delta + delta)

               end if
            else
*              ---------------------------------------------------------
*              Forward Differences.
*              ---------------------------------------------------------
               if (nfound .gt. 0) then
                  do 330, i = 1, ncnln
                     if (needc(i) .eq. 1)
     $                  cjac(i,j) = (c1(i) -  c(i))/  delta
  330             continue
               end if

               if (gradu(j) .eq. rdummy)
     $            grad(j) = (objf1 - objf) /  delta

            end if
         end if
         x(j) = xj

  340 continue

      return

  999 inform = mode
      return

*     end of  npfd
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npfeas( n, nclin, ncnln, istate,
     $                   bigbnd, cvnorm, errmax, jmax, nviol,
     $                   ax, bl, bu, c, featol, x, work )

      implicit           double precision (a-h,o-z)
      integer            istate(n+nclin+ncnln)
      double precision   ax(*), bl(n+nclin+ncnln), bu(n+nclin+ncnln)
      double precision   c(*), featol(n+nclin+ncnln), x(n)
      double precision   work(n+nclin+ncnln)

************************************************************************
*     npfeas  computes the following...
*     (1)  The number of constraints that are violated by more
*          than  featol  and the 2-norm of the constraint violations.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version      April    1984.
*     This version of  npfeas  dated  16-October-1985.
************************************************************************
      common    /sol1cm/ nout

      logical            npdbg
      parameter        ( ldbg = 5 )
      common    /npdebg/ inpdbg(ldbg), npdbg

      external           idamax, dnrm2
      intrinsic          abs
      parameter        ( zero = 0.0d+0 )

      biglow = - bigbnd
      bigupp =   bigbnd

*     ==================================================================
*     Compute nviol, the number of constraints violated by more than
*     featol,  and cvnorm,  the 2-norm of the constraint
*     violations and residuals of the constraints in the qp working set.
*     ==================================================================
      nviol  = 0

      do 200, j = 1, n+nclin+ncnln
         feasj  = featol(j)
         res    = zero

         if (j .le. n + nclin) then

*           Bound or general linear constraint.

            if (j .le. n) then
               con =  x(j)
            else
               con = ax(j-n)
            end if

            tolj   = feasj
         else

*           Nonlinear constraint.

            con    = c(j-n-nclin)
            tolj   = zero
         end if

*        Check for constraint violations.

         if (bl(j) .gt. biglow) then
            res    = bl(j) - con
            if (res .gt.   feasj ) nviol = nviol + 1
            if (res .gt.    tolj ) go to 190
         end if

         if (bu(j) .lt. bigupp) then
            res    = bu(j) - con
            if (res .lt. (-feasj)) nviol = nviol + 1
            if (res .lt.  (-tolj)) go to 190
         end if

*        This constraint is satisfied,  but count the residual as a
*        violation if the constraint is in the working set.

         is     = istate(j)

         if (is .eq. 0) then
            res = zero
         else if (is .eq. 1  .or.  is .le. -2) then
            res = bl(j) - con
         else if (is .ge. 2  .or.  is .eq. -1) then
            res = bu(j) - con
         end if

         if (abs( res ) .gt. feasj) nviol = nviol + 1

*        Set the array of violations.

  190    work(j) = res
  200 continue

      jmax   = idamax( n+nclin+ncnln, work, 1 )
      errmax = abs ( work(jmax) )

      if (npdbg  .and.  inpdbg(1) .gt. 0)
     $   write (nout, 1000) errmax, jmax

      cvnorm = dnrm2 ( n+nclin+ncnln, work, 1 )

      return

 1000 format(/ ' //npfeas//  The maximum violation is ', 1p, e14.2,
     $                     ' in constraint', i5 )

*     end of  npfeas
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npfile( ioptns, inform )
      integer            ioptns, inform

************************************************************************
*     npfile  reads the options file from unit  ioptns  and loads the
*     options into the relevant elements of  iprmnp  and  rprmnp.
*
*     If  ioptns .lt. 0  or  ioptns .gt. 99  then no file is read,
*     otherwise the file associated with unit  ioptns  is read.
*
*     Output:
*
*         inform = 0  if a complete  options  file was found
*                     (starting with  begin  and ending with  end);
*                  1  if  ioptns .lt. 0  or  ioptns .gt. 99;
*                  2  if  begin  was found, but end-of-file
*                     occurred before  end  was found;
*                  3  if end-of-file occurred before  begin  or
*                     endrun  were found;
*                  4  if  endrun  was found before  begin.
************************************************************************
      logical             newopt
      common     /sol7np/ newopt
      save       /sol7np/

      double precision    wmach(15)
      common     /solmch/ wmach
      save       /solmch/

      external            mchpar, npkey
      logical             first
      save                first , nout
      data                first /.true./

*     If first time in, set  nout.
*     newopt is true first time into npfile or npoptn
*     and just after a call to npsol.

      if (first) then
         first  = .false.
         newopt = .true.
         call mchpar()
         nout = wmach(11)
      end if

      call opfile( ioptns, nout, inform, npkey )

      return

*     end of  npfile
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npgetr( inform, unitq, n, nclin, ncnln, ldr, ldq,
     $                   nfree, iter, istate, kx,
     $                   hforwd, hcntrl,
     $                   cmul, r, rho, x, zy )

      implicit           double precision (a-h,o-z)
      logical            unitq
      integer            istate(n+nclin+ncnln), kx(n)
      double precision   r(ldr,*), x(n), zy(ldq,*)
      double precision   hforwd(*), hcntrl(*)
      double precision   cmul(*), rho(*)

C***********************************************************************
C     npgetr  loads details of a previous run from unit nload.
C
C     Original version   24-Nov-89.
C     This version of  npgetr  dated 26-Nov-89.
C***********************************************************************
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
      common    /sol1cm/ nout
      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset
      logical            incrun
      common    /sol6np/ rhomax, rhonrm, rhodmp, scale, incrun

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

      character*4        icrsh(0:2)
      data                icrsh(0),  icrsh(1),  icrsh(2)
     $                 / 'COLD'   , 'WARM'   , 'HOT '    /

      if (nload .le. 0) return

      write( nout , 4000 ) nload

      read ( nload, 1000, end=999 ) iter, nfree, lfdset, lvldif, unitq
      do 110, j = 1, n
         read ( nload, 1010, end=999 ) jold, kx(j), istate(j), x(j)
  110 continue

      if (jold .ne. n) then
         write( nout, 9000 )
         inform = 1
         return
      end if

      do 120, j = n+1, n+nclin
         read ( nload, 1020, end=999 ) jold, istate(j)
  120 continue

      if (jold .ne. n+nclin) then
         write( nout, 9000 )
         inform = 1
         return
      end if

      if (ncnln .gt. 0) then
         k = 1
         do 130, j = n+nclin+1, n+nclin+ncnln
            read( nload, 1030, end=999 ) jold, istate(j), cmul(k),
     $                                   rho(k)
            k = k + 1
  130    continue

         if (jold .ne. n+nclin+ncnln) then
            write( nout, 9000 )
            inform = 1
            return
         end if
     
         read( nload, 1040, end=999 ) rhomax, rhonrm, rhodmp, scale,
     $                                incrun
      end if

*     ------------------------------------------------------------------
*     Read   Q(free)  and the factor of  Q'HQ.
*     ------------------------------------------------------------------
      if (.not. unitq) then
         do 160, j = 1, nfree
            do 150, i = 1, nfree
              read( nload, 1050, end=999 ) iold, jold, zy(i,j)
  150       continue     
  160    continue
      end if 

      do 180, j = 1, n
         do 170, i = 1, j
           read ( nload, 1050, end=999 ) iold, jold, r(i,j)
  170    continue
  180 continue

      if (jold .ne. n  .or. iold .ne. n) then
         write( nout, 9000 )
         inform = 1
         return
      end if

*     ------------------------------------------------------------------
*     Read the finite-difference intervals.  
*     ------------------------------------------------------------------
      if (lvldif .gt. 0) then
         if (lfdset .eq. 0  .or.  lfdset .eq. 2) then
            do 190, j = 1, n
               read( nload, 1090 ) jold, hforwd(j), hcntrl(j)
  190       continue
         end if
      end if                                            

      write ( nout , 4010 ) n, iter

      call mcclos( nload )

*     ------------------------------------------------------------------
*     Now that all values have been set, we can print the parameters.
*     ------------------------------------------------------------------
      if (msgnp .gt. 0) then
         epsmch = wmach( 3)
         write (nout, 2000)
         write (nout, 2100) nclin , tolfea, icrsh(lcrash) ,
     $                      n     , bigbnd, tolact,
     $                      dxlim , bigdx
         write (nout, 2200) ncnln , ftol  , epsrf ,
     $                      nlnj  , ctol  ,
     $                      nlnf  , eta   ,
     $                      epsmch,
     $                      lvlder, lverfy
         write (nout, 2300) nmajor, msgnp ,
     $                      nminor, msgqp ,
     $                      nload , nsave , ksave

         if (lvlder .lt. 3) then
            if      (lfdset .eq. 0) then
               write (nout, 2400)
            else if (lfdset .eq. 1) then
               write (nout, 2401) fdint, cdint
            else if (lfdset .eq. 2) then
               write (nout, 2402)
            end if
         end if
      end if

      return
                       
 999  write( nout, 9000 )
      inform = 1 
      return

 1000 format(2i8, 1x, 2i2, 1x, l1 )
 1010 format(2i8, 1x,  i2, 1p, 2e24.14 )
 1020 format( i8, 1x,  i2 )
 1030 format( i8, 1x,  i2, 1p, 2e24.14 )
 1040 format( 1p, 4e24.14, 1x, l1 )
 1050 format(2i8, 1x, 1p,  e24.14 )
 1090 format( i8, 1x, 1p, 2e24.14 )

 2000 format(
     $//' Parameters'
     $/ ' ----------' )
 2100 format(
     $/ ' Linear constraints.....',     i10  , 6x,
     $  ' Linear feasibility.....', 1p, e10.2, 6x,
     $  1x, a4, ' start.............'
     $/ ' Variables..............',     i10,   6x,
     $  ' Infinite bound size....',     e10.2, 6x,
     $  ' Crash tolerance........',     e10.2
     $/ ' Step limit.............',     e10.2, 6x,
     $  ' Infinite step size.....',     e10.2  )
 2200 format(
     $/ ' Nonlinear constraints..',     i10,   6x,
     $  ' Optimality tolerance...', 1p, e10.2, 6x,
     $  ' Function precision.....',     e10.2
     $/ ' Nonlinear Jacobian vars',     i10,   6x,
     $  ' Nonlinear feasibility..',     e10.2
     $/ ' Nonlinear objectiv vars',     i10,   6x,
     $  ' Linesearch tolerance...',     e10.2
     $/ ' EPS (machine precision)',     e10.2, 6x,
     $  ' Derivative level.......',     i10,   6x,
     $  ' Verify level...........',     i10)
 2300 format(
     $/ ' Major iterations limit.',     i10, 6x,
     $  ' Major print level......',     i10
     $/ ' Minor iterations limit.',     i10, 6x,
     $  ' Minor print level......',     i10
     $/ ' RUN loaded from file...',     i10, 6x,
     $  ' RUN to be saved on file',     i10, 6x,
     $  ' Save frequency.........',     i10)

 2400 format(/ ' Difference intervals to be computed.' )
 2401 format(/ ' Difference interval....', 1p, e10.2, 6x,
     $         ' Central diffce interval',     e10.2 )
 2402 format(/ ' User-supplied difference intervals.' )
 4000 format(/ ' OLD RUN to be loaded from file', i4)
 4010 format(/ ' Number of variables loaded = ', i5, '    Itn = ', i8)
 9000 format(/ ' Exit NPSOL - The OLD RUN file dimensions do not match'
     $       , '  this problem.' )

*     end of  npgetr
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npiqp ( feasqp, unitq, nqperr, minits,
     $                   n, nclin, ncnln, lda, ldcj, ldaqp, ldr,
     $                   linact, nlnact, nactiv, nfree, nz, numinf,
     $                   istate, kactiv, kx,
     $                   dxnorm, gdx, qpcurv,
     $                   aqp, adx, anorm, ax, bl, bu,
     $                   c, cjac, clamda, cmul, cs,
     $                   dlam, dslk, dx, qpbl, qpbu, qptol,
     $                   r, rho, slk, violn, x,
     $                   wtinf, iw, w )

      implicit           double precision (a-h,o-z)
      logical            feasqp, unitq
      integer            istate(*), kactiv(n), kx(n)
      integer            iw(*)
      double precision   aqp(ldaqp,*), adx(*), anorm(*), ax(*),
     $                   bl(*), bu(*),
     $                   c(*), cjac(ldcj,*), clamda(*), cmul(*), cs(*)
      double precision   dlam(*), dslk(*), dx(n)
      double precision   qpbl(*), qpbu(*),
     $                   qptol(*), r(ldr,*), rho(*), slk(*),
     $                   violn(*), x(n), wtinf(*)
      double precision   w(*)

************************************************************************
*     NPIQP   does the following:
*
*     (1)  Generate the upper and lower bounds for the QP  subproblem.
*
*     (2)  Compute the  TQ  factors of the rows of  AQP  specified by
*          the array  ISTATE.  The part of the factorization defined by
*          the first contiguous group of linear constraints does not
*          need to be recomputed.  The remaining rows (which could be
*          comprised of both linear and nonlinear constraints) are
*          included as new rows of the  TQ  factorization stored in
*          T and ZY.  Note that if there are no nonlinear constraints,
*          no factorization is required.
*
*     (3)  Solve the  QP  subproblem.
*                 minimize     1/2 (W p - d)'(Wp - d) + g'p
*
*                 subject to   qpbl .le. (  p ) .le. qpbu,
*                                        ( Ap )
*
*          where  W  is a matrix (not stored) such that  W'W = H  and
*          WQ = R,  d  is the zero vector,  and  g  is the gradient.
*          If the subproblem is infeasible, compute the point which
*          minimizes the sum of infeasibilities.
*
*    (4)   Find the value of each slack variable for which the merit
*          function is minimized.
*
*    (5)   Compute  DSLK,  DLAM  and  DX,  the search directions for
*          the slack variables, the multipliers and the variables.
*
*     Systems Optimization Laboratory, Stanford University.
*     Fortran 66 version written 10-January-1983.
*     This version of NPIQP dated 10-May-90.
************************************************************************
      common    /sol1cm/ nout
      common    /sol3cm/ lennam, ldt   , ncolt , ldq
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5cm/ asize , dtmax , dtmin
      common    /sol6cm/ rcndbd, rfrobn, drmax , drmin
      
      integer            locls
      parameter         (lenls = 20)
      common    /sol1ls/ locls(lenls)

      logical            incrun
      common    /sol6np/ rhomax, rhonrm, rhodmp, scale, incrun

      logical            cmdbg, lsdbg, npdbg
      parameter        ( ldbg = 5 )
      common    /npdebg/ inpdbg(ldbg), npdbg
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

      character*8        names(1)
      logical            linobj, overfl, qpnamd, vertex
      intrinsic          abs   , min   , max
      external           ddiv  , ddot  , dnrm2
      parameter        ( qpnamd =.false.,vertex =.false. )
      parameter        ( zero   =0.0d+0, one    =1.0d+0, two   =2.0d+0 )
      parameter        ( hundrd =1.0d+2                                )

      idbgsv = idbg
      if (npdbg) then
         idbg   = 0
      else
         idbg = nminor + 1
      end if
      lsdbg  = npdbg
      cmdbg  = npdbg
      call icopy ( ldbg, ilsdbg, 1, icmdbg, 1 )

      lrpq   = locls( 5)
      lrpq0  = locls( 6)
      lhpq   = locls( 8)
      lgq    = locls( 9)
      lrlam  = locls(10)
      lt     = locls(11)
      lq     = locls(12)
      lwrk1  = locls(14)

      nrpq   = 0
      ngq    = 1

      feasqp =  .true.
      linobj =  .true.

      biglow = - bigbnd
      bigupp =   bigbnd
      ssq1   =   zero

      nplin  = n     + nclin
      nctotl = nplin + ncnln
      ncqp   = nclin + ncnln
      nrank  = n
      nrejtd = 0                           

*     ==================================================================
*     Generate the upper and lower bounds upon the search direction, the
*     weights on the sum of infeasibilities and the nonlinear constraint
*     violations.
*     ==================================================================
      wscale = - one
      do 170 j = 1, nctotl

         if (j .le. n) then                
            con = x(j)
         else if (j .le. nplin) then
            con = ax(j-n)
         else
            con = c(j-nplin)
         end if

         blj = bl(j)
         buj = bu(j)
         if (blj .gt. biglow) blj = blj - con
         if (buj .lt. bigupp) buj = buj - con

         weight = one
         if (j .le. nplin) then
            if (abs(blj) .le. qptol(j)) blj = zero
            if (abs(buj) .le. qptol(j)) buj = zero
         else
            i    = j - nplin
            viol = zero
            if (bl(j) .gt. biglow) then
               if (blj .gt. zero) then
                  viol   = blj
                  if (rho(i) .gt. zero) then
                     weight =   viol*rho(i)
                  else
                     weight =   viol
                  end if
                  wscale = max( wscale,   weight )
                  go to 160
               end if
            end if

            if (bu(j) .lt. bigupp) then
               if (buj .lt. zero) then
                  viol   =   buj
                  if (rho(i) .gt. zero) then
                     weight = - viol*rho(i)
                  else
                     weight = - viol
                  end if
                  wscale = max( wscale, weight )
               end if
            end if

*           Set the vector of nonlinear constraint violations.

  160       violn(i) = viol
         end if

         wtinf(j) = weight
         qpbl(j)  = blj
         qpbu(j)  = buj
  170 continue

      if (wscale .gt. zero) then
         wscale = one/wscale
         call dscal ( nctotl, (wscale), wtinf, 1 )
      end if

      call dcond ( nctotl, wtinf, 1, wtmax, wtmin )
      wtmin  = epspt9*wtmax
      do 180, j = 1, nctotl
         wtinf(j) = max( wtinf(j), wtmin )
  180 continue

*     Set the maximum allowable condition estimator of the constraints
*     in the working set.  Note that a relatively well-conditioned
*     working set is used to start the QP iterations.

      condmx = max( one/epspt3, hundrd )

      if (ncnln .gt. 0) then
*        ===============================================================
*        Refactorize part of the  QP  constraint matrix.
*        ===============================================================
*        Load the new Jacobian into the  QP  matrix  A.  Compute the
*        2-norms of the rows of the Jacobian.

         call f06qff( 'General', ncnln, n, cjac, ldcj,
     $                aqp(nclin+1,1), ldaqp )

         do 190 j = nclin+1, ncqp
            anorm(j) = dnrm2 ( n, aqp(j,1), ldaqp )
  190    continue

*        Count the number of linear constraints in the working set and
*        move them to the front of kactiv.  Compute the norm of the
*        matrix of constraints in the working set.
*        Let k1  point to the first nonlinear constraint.  Constraints
*        with indices kactiv(k1),..., kactiv(nactiv)  must be
*        refactorized.

         asize  = zero
         linact = 0
         k1     = nactiv + 1
         do 200 k = 1, nactiv
            i     = kactiv(k)
            asize = max( asize, anorm(i) )

            if (i .le. nclin) then
               linact = linact + 1
               if (linact .ne. k) then
                  iswap  = kactiv(linact)
                  kactiv(linact) = i
                  kactiv(k)      = iswap
               end if
            else

*              Record the old position of the 1st. nonlinear constraint.

               if (k1 .gt. nactiv) k1 = k
            end if
  200    continue

         if (nactiv .le. 1 )
     $      call dcond ( ncqp, anorm, 1, asize, amin )

*        Compute the absolute values of the nonlinear constraints in
*        the working set.  Use DX as workspace.

         do 210 k = linact+1, nactiv
            j = n + kactiv(k)
            if (istate(j) .eq. 1) dx(k) = abs( qpbl(j) )
            if (istate(j) .ge. 2) dx(k) = abs( qpbu(j) )
  210    continue

*        Sort the elements of KACTIV corresponding to nonlinear
*        constraints in descending order of violation (i.e.,
*        the first element of KACTIV for a nonlinear constraint
*        is associated with the most violated constraint.)
*        In this way, the rows of the Jacobian corresponding
*        to the more violated constraints tend to be included
*        in the  TQ  factorization.

*        The sorting procedure is taken from the simple insertion
*        sort in D. Knuth, ACP Volume 3, Sorting and Searching,
*        Page 81.  It should be replaced by a faster sort if the
*        number of active nonlinear constraints becomes large.

         do 230 k = linact+2, nactiv
            l     = k
            viol  = dx(l)
            kviol = kactiv(l)
*           while (l .gt. linact+1  .and.  dx(l-1) .lt. viol) do
  220       if    (l .gt. linact+1                          ) then
               if (                        dx(l-1) .lt. viol) then
                  dx(l)     = dx(l-1)
                  kactiv(l) = kactiv(l-1)
                  l         = l - 1
                  go to 220
               end if
*           end while
            end if
            dx(l)     = viol
            kactiv(l) = kviol
  230    continue

         k2     = nactiv
         nactiv = k1     - 1
         nz     = nfree  - nactiv

*        Update the factors  R,  T  and  Q  to include constraints
*        K1  through  K2.

         if (k1 .le. k2)
     $      call lsadds( unitq, vertex,
     $                   inform, k1, k2, nactiv, nartif, nz, nfree,
     $                   nrank, nrejtd, nrpq, ngq,
     $                   n, ldq, ldaqp, ldr, ldt,
     $                   istate, kactiv, kx, condmx,
     $                   aqp, r, w(lt), w(lrpq), w(lgq), w(lq),
     $                   w(lwrk1), dx, w(lrlam) )
      end if

*     ==================================================================
*     Solve for dx, the vector of minimum two-norm that satisfies the
*     constraints in the working set.
*     ==================================================================
      call npsetx( unitq,
     $             ncqp, nactiv, nfree, nz,
     $             n, nlnx, nctotl, ldq, ldaqp, ldr, ldt,
     $             istate, kactiv, kx,
     $             dxnorm, gdx,
     $             aqp, adx, qpbl, qpbu, w(lrpq), w(lrpq0), dx, w(lgq),
     $             r, w(lt), w(lq), w(lwrk1) )

*     ==================================================================
*     Solve a quadratic program for the search direction  DX  and
*     multiplier estimates  clamda.
*     ==================================================================
*     If there is no feasible point for the subproblem,  the sum of
*     infeasibilities is minimized subject to the linear constraints
*     (1  thru  jinf)  being satisfied.

      jinf  = n + nclin

      ntry  = 1
*+    repeat
  450    call lscore( 'QP subproblem', qpnamd, names, linobj, unitq,
     $                nqperr, minits, jinf, ncqp, nctotl,
     $                nactiv, nfree, nrank, nz, nz1,
     $                n, ldaqp, ldr,
     $                istate, kactiv, kx,
     $                gdx, ssq, ssq1, suminf, numinf, dxnorm,
     $                qpbl, qpbu, aqp, clamda, adx,
     $                qptol, r, dx, iw, w )

         if (npdbg  .and.  inpdbg(1) .gt. 0)
     $      write (nout, 8000) nqperr

         nviol = 0
         if (numinf .gt. 0) then

*           Count the violated linear constraints.

            do 460 j = 1, nplin
               if (istate(j) .lt. 0) nviol = nviol + 1
  460       continue

            if (nviol .gt. 0) then
               ntry   = ntry + 1
               unitq  = .true.
               nactiv = 0
               nfree  = n
               nz     = n
               call iload ( nctotl, (0), istate, 1 )

               call npsetx( unitq,
     $                      ncqp, nactiv, nfree, nz,
     $                      n, nlnx, nctotl, ldq, ldaqp, ldr, ldt,
     $                      istate, kactiv, kx,
     $                      dxnorm, gdx,
     $                      aqp, adx, qpbl, qpbu, w(lrpq), w(lrpq0),
     $                      dx, w(lgq), r, w(lt), w(lq), w(lwrk1) )
            end if
         end if
      if (.not. (nviol .eq. 0  .or.  ntry .gt. 2)) go to 450
*+    until (    nviol .eq. 0  .or.  ntry .gt. 2)

*     ==================================================================
*     Count the number of nonlinear constraint gradients in the  QP
*     working set.  Make sure that all small  QP  multipliers associated
*     with nonlinear inequality constraints have the correct sign.
*     ==================================================================
      nlnact  = 0
      if (nactiv .gt. 0  .and.  ncnln .gt. 0) then
         do 500 k = 1, nactiv
            l     = kactiv(k)
            if (l .gt. nclin) then
               nlnact = nlnact + 1
               j      = n      + l
               if (istate(j) .eq. 1) clamda(j) = max( zero, clamda(j) )
               if (istate(j) .eq. 2) clamda(j) = min( zero, clamda(j) )
            end if
  500    continue
      end if

      linact = nactiv - nlnact

*     ------------------------------------------------------------------
*     Extract various useful quantities from the QP solution.
*     ------------------------------------------------------------------
*     Compute  HPQ = R'R(pq)  from the transformed gradient of the QP
*     objective function and  R(pq)  from the transformed residual.

      call dscal ( n, (-one), w(lrpq), 1 )
      call daxpy ( n, (-one), w(lgq) , 1, w(lhpq), 1 )
      qpcurv = two*ssq

      if (ncnln .gt. 0) then
         if (numinf .gt. 0) then
            feasqp = .false.
            call dload ( nctotl, (zero), clamda, 1 )

            if (nz .gt. 0) then
*              ---------------------------------------------------------
*              Compute a null space component for the search direction
*              as the solution of  Z'HZ(pz) = -Z'g - Z'HY(py).
*              ---------------------------------------------------------
*              Overwrite DX with the transformed search direction
*              Q'(dx).  The first NZ components of DX are zero.

               call cmqmul( 6, n, nz, nfree, ldq, unitq,
     $                      kx, dx, w(lq), w(lwrk1) )

*              Overwrite the first NZ components of DX with the solution
*              of  (Rz)u = -(v + w),  where  (Rz)'w = Z'g  and  v  is
*              vector of first NZ components of  R(pq).

               call dcopy ( nz, w(lgq), 1, dx, 1 )
               call dtrsv ( 'U', 'T', 'N', nz, r, ldr, dx, 1 )

               call daxpy ( nz, (one), w(lrpq), 1, dx, 1 )

               call dtrsv ( 'U', 'N', 'N', nz, r, ldr, dx, 1 )
               call dscal ( nz, (-one), dx, 1 )

*              Recompute rpq, hpq, gdx and qpcurv.

               call dcopy ( nlnx, dx, 1, w(lrpq), 1 )
               call dtrmv ( 'U', 'N', 'N', nlnx, r, ldr, w(lrpq), 1 )
               if (nlnx .lt. n)
     $            call dgemv( 'N', nlnx, n-nlnx, one, r(1,nlnx+1), ldr,
     $                        dx(nlnx+1), 1, one, w(lrpq), 1 )

               gdx    = ddot  ( n, w(lgq) , 1, dx     , 1 )
               qpcurv = ddot  ( n, w(lrpq), 1, w(lrpq), 1 )

               call cmqmul( 3, n, nz, nfree, ldq, unitq,
     $                      kx, dx, w(lq), w(lwrk1) )

*              ---------------------------------------------------------
*              Recompute ADX and the 2-norm of DX.
*              ---------------------------------------------------------
               dxnorm  = dnrm2 ( n, dx, 1 )
               if (ncqp .gt. 0)
     $            call dgemv ( 'N', ncqp, n, one, aqp, ldaqp,
     $                         dx, 1, zero, adx, 1 )

               if (npdbg  .and.  inpdbg(2) .gt. 0)
     $            write (nout, 8100) (dx(j), j = 1, n)
            end if

            call dcopy ( nlnx, w(lrpq), 1, w(lhpq), 1 )
            call dtrmv ( 'U', 'T', 'N', nlnx, r, ldr, w(lhpq), 1 )
            if (nlnx .lt. n)
     $         call dgemv ( 'T', nlnx, n-nlnx, one, r(1,nlnx+1), ldr,
     $                      w(lrpq), 1, zero, w(lhpq+nlnx), 1 )
         end if

*        ===============================================================
*        For given values of the objective function and constraints,
*        attempt to minimize the merit function with respect to each
*        slack variable.
*        ===============================================================
         do 600 i = 1, ncnln
            j      = nplin + i
            con    = c(i)

            if (      .not. feasqp  .and.
     $          violn(i) .ne. zero  .and.  rho(i) .le. zero )
     $         rho(i) = one

            quotnt = ddiv  ( cmul(i), scale*rho(i), overfl )

*           Define the slack variable to be  CON - MULT / RHO.
*           Force each slack to lie within its upper and lower bounds.

            if (bl(j) .gt. biglow) then
               if (qpbl(j) .ge. - quotnt) then
                  slk(i) = bl(j)
                  go to 550
               end if
            end if

            if (bu(j) .lt. bigupp) then
               if (qpbu(j) .le. - quotnt) then
                  slk(i) = bu(j)
                  go to 550
               end if
            end if

            slk(i) = con - quotnt

*           The slack has been set within its bounds.

  550       cs(i)  = con - slk(i)

*           ------------------------------------------------------------
*           Compute the search direction for the slacks and multipliers.
*           ------------------------------------------------------------
            dslk(i) = adx(nclin+i) + cs(i)

            if (feasqp) then
*
*              If any constraint is such that  (DLAM)*(C - S)  is
*              positive,  the merit function may be reduced immediately
*              by substituting the QP multiplier.
*
               dlam(i)  = clamda(j) - cmul(i)
               if (dlam(i) * cs(i) .ge. zero) then
                  cmul(i) = clamda(j)
                  dlam(i) = zero
               end if
            else

*              The  QP  subproblem was infeasible.

               dlam(i) = zero

               if (istate(j) .lt. 0  .or.  violn(i) .ne. zero)
     $            dslk(i)  = zero

            end if
  600    continue

         if (.not. feasqp)
     $      rhonrm = dnrm2 ( ncnln, rho, 1 )

         if (npdbg  .and.  inpdbg(2) .gt. 0) then
            write (nout, 8200) (violn(i), i=1,ncnln)
            write (nout, 8300) (slk(i)  , i=1,ncnln)
         end if
      end if

      call icopy ( ldbg, inpdbg, 1, icmdbg, 1 )
      idbg   = idbgsv

      return

 8000 format(/ ' //npiqp // nqperr'
     $       / ' //npiqp // ',  i6 )
 8100 format(/ ' //npiqp // dx recomputed with null space portion...'
     $       / (5g12.3))
 8200 format(/ ' //npiqp // Violations = '/ (1p, 5e15.6))
 8300 format(/ ' //npiqp // Slacks     = '/ (1p, 5e15.6))

*     end of  npiqp
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE NPKEY ( NOUT, BUFFER, KEY )

      IMPLICIT           DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*)      BUFFER

************************************************************************
*     NPKEY   decodes the option contained in  BUFFER  in order to set
*     a parameter value in the relevant element of the parameter arrays.
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
*     NPKEY  calls OPNUMB and the subprograms
*                 LOOKUP, SCANNR, TOKENS, UPCASE
*     (now called OPLOOK, OPSCAN, OPTOKN, OPUPPR)
*     supplied by Informatics General, Inc., Palo Alto, California.
*
*     Systems Optimization Laboratory, Stanford University.
*     This version of NPKEY  dated 24-Nov-89.
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
*-----INCLUDE NPPARM----------------------------------------------------
      INTEGER            IPRMNP(MXPARM), IPSVNP
      DOUBLE PRECISION   RPRMNP(MXPARM), RPSVNP

      COMMON    /NPPAR1/ IPSVNP(MXPARM),
     $                   IDBGNP, ITMXNP, JVRFY1, JVRFY2, JVRFY3, JVRFY4,
     $                   LDBGNP, LFORMH, LVLDER, LVERFY, MSGNP , NLNF  ,
     $                   NLNJ  , NLNX  , NNCNLN, NSAVE , NLOAD , KSAVE ,
     $                   IPADNP(12)

      COMMON    /NPPAR2/ RPSVNP(MXPARM),
     $                   CDINT , CTOL  , DXLIM , EPSRF , ETA   , FDINT ,
     $                   FTOL  , HCNDBD, RPADNP(22)

      EQUIVALENCE       (IPRMNP(1), IDBGNP), (RPRMNP(1), CDINT)

      SAVE      /NPPAR1/, /NPPAR2/
*-----------------------------------------------------------------------
      EQUIVALENCE  (IDBGNP, IDBG  ), (ITMXNP, NMAJOR), (ITMAX2, NMINOR)
      EQUIVALENCE  (LDBGLS, MNRDBG), (LDBGNP, MJRDBG), (MSGLS , MSGQP )

      EXTERNAL           OPNUMB
      LOGICAL            FIRST , MORE  , NUMBER, OPNUMB, SORTED
      SAVE               FIRST

      PARAMETER         (     MAXKEY = 42,  MAXTIE = 21,   MAXTOK = 10)
      CHARACTER*16       KEYS(MAXKEY), TIES(MAXTIE), TOKEN(MAXTOK)
      CHARACTER*16       KEY, KEY2, KEY3, VALUE

      PARAMETER         (IDUMMY = -11111,  RDUMMY = -11111.0d+0,
     $                   SORTED = .TRUE.,  ZERO   =  0.0     )

      DATA                FIRST
     $                  /.TRUE./
      DATA   KEYS
     $ / 'BEGIN           ',
     $   'CENTRAL         ', 'COLD            ', 'CONDITION       ',
     $   'CONSTRAINTS     ',
     $   'CRASH           ', 'DEBUG           ', 'DEFAULTS        ',
     $   'DERIVATIVE      ', 'DIFFERENCE      ', 'END             ',
     $   'FEASIBILITY     ', 'FUNCTION        ', 'HESSIAN         ',
     $   'HOT             ', 'INFINITE        ', 'IPRMLS          ',
     $   'ITERATIONS      ', 'ITERS:ITERATIONS', 'ITNS :ITERATIONS',
     $   'LINEAR          ', 'LINESEARCH      ', 'LIST            ',
     $   'LOAD            ',
     $   'LOWER           ', 'MAJOR           ', 'MINOR           ',
     $   'NOLIST          ', 'NONLINEAR       ', 'OPTIMALITY      ',
     $   'PRINT           ', 'PROBLEM         ', 'ROW             ',
     $   'RPRMLS          ', 'SAVE            ', 'START           ',
     $   'STEP            ',
     $   'STOP            ', 'UPPER           ', 'VARIABLES       ',
     $   'VERIFY          ', 'WARM            '/

      DATA   TIES
     $ / 'BOUND           ', 'CONSTRAINTS     ', 'DEBUG           ',
     $   'FEASIBILITY     ', 'FREQUENCY       ', 'GRADIENTS       ',
     $   'ITERATIONS      ', 'ITERS:ITERATIONS',
     $   'ITNS :ITERATIONS', 'JACOBIAN        ', 'LEVEL           ',
     $   'NO              ',
     $   'NO.      :NUMBER',
     $   'NUMBER          ', 'OBJECTIVE       ', 'PRINT           ',
     $   'RUN             ', 'STEP            ', 'TOLERANCE       ',
     $   'VARIABLES       ', 'YES             '/
*-----------------------------------------------------------------------

      IF (FIRST) THEN
         FIRST  = .FALSE.
         DO 10 I = 1, MXPARM
            RPRMLS(I) = RDUMMY
            IPRMLS(I) = IDUMMY
            RPRMNP(I) = RDUMMY
            IPRMNP(I) = IDUMMY
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
         IF      (KEY .EQ. 'CENTRAL     ') THEN
            CDINT  = RVALUE
         ELSE IF (KEY .EQ. 'COLD        ') THEN
            LCRASH = 0
         ELSE IF (KEY .EQ. 'CONDITION   ') THEN
            HCNDBD = RVALUE
         ELSE IF (KEY .EQ. 'CONSTRAINTS ') THEN
            NNCLIN = RVALUE
         ELSE IF (KEY .EQ. 'CRASH       ') THEN
            TOLACT = RVALUE
         ELSE IF (KEY .EQ. 'DEBUG       ') THEN
            IDBG   = RVALUE
         ELSE IF (KEY .EQ. 'DEFAULTS    ') THEN
            DO 20 I = 1, MXPARM
               IPRMLS(I) = IDUMMY
               RPRMLS(I) = RDUMMY
               IPRMNP(I) = IDUMMY
               RPRMNP(I) = RDUMMY
   20       CONTINUE
         ELSE IF (KEY .EQ. 'DERIVATIVE  ') THEN
            LVLDER = RVALUE
         ELSE IF (KEY .EQ. 'DIFFERENCE  ') THEN
            FDINT  = RVALUE
         ELSE IF (KEY .EQ. 'FEASIBILITY ') THEN
            TOLFEA = RVALUE
            CTOL   = RVALUE
         ELSE IF (KEY .EQ. 'FUNCTION    ') THEN
            EPSRF  = RVALUE
         ELSE
            MORE   = .TRUE.
         END IF
      END IF

      IF (MORE) THEN
         MORE   = .FALSE.
         IF      (KEY .EQ. 'HESSIAN     ') THEN
            LFORMH = 1
            IF   (KEY2.EQ. 'NO          ') LFORMH = 0
         ELSE IF (KEY .EQ. 'HOT         ') THEN
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
            NMAJOR = RVALUE
         ELSE IF (KEY .EQ. 'LINEAR      ') THEN
            IF (KEY2  .EQ. 'CONSTRAINTS ') NNCLIN = RVALUE
            IF (KEY2  .EQ. 'FEASIBILITY ') TOLFEA = RVALUE
            IF (LOC2 .EQ.  0             ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'LINESEARCH  ') THEN
            ETA    = RVALUE
         ELSE IF (KEY .EQ. 'LOWER       ') THEN
            BNDLOW = RVALUE
         ELSE IF (KEY .EQ. 'LOAD        ') THEN
            NLOAD  = RVALUE
         ELSE
            MORE   = .TRUE.
         END IF
      END IF

      IF (MORE) THEN
         MORE   = .FALSE.
         IF      (KEY .EQ. 'MAJOR       ') THEN
              IF (KEY2.EQ. 'DEBUG       ') MJRDBG = RVALUE
              IF (KEY2.EQ. 'ITERATIONS  ') NMAJOR = RVALUE
              IF (KEY2.EQ. 'PRINT       ') MSGNP  = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'MINOR       ') THEN
              IF (KEY2.EQ. 'DEBUG       ') MNRDBG = RVALUE
              IF (KEY2.EQ. 'ITERATIONS  ') NMINOR = RVALUE
              IF (KEY2.EQ. 'PRINT       ') MSGQP  = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'NONLINEAR   ') THEN
              IF (KEY2.EQ. 'CONSTRAINTS ') NNCNLN = RVALUE
              IF (KEY2.EQ. 'FEASIBILITY ') CTOL   = RVALUE
              IF (KEY2.EQ. 'JACOBIAN    ') NLNJ   = RVALUE
              IF (KEY2.EQ. 'OBJECTIVE   ') NLNF   = RVALUE
              IF (KEY2.EQ. 'VARIABLES   ') NLNX   = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'OPTIMALITY  ') THEN
            FTOL   = RVALUE
         ELSE
            MORE   = .TRUE.
         END IF
      END IF

      IF (MORE) THEN
         MORE   = .FALSE.
         IF      (KEY .EQ. 'PRINT       ') THEN
            MSGNP  = RVALUE
         ELSE IF (KEY .EQ. 'PROBLEM     ') THEN
              IF (KEY2.EQ. 'NUMBER      ') NPROB  = RVALUE
         ELSE IF (KEY .EQ. 'ROW         ') THEN
              IF (KEY2.EQ. 'TOLERANCE   ') CTOL   = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'RPRMLS      ') THEN
*           Allow things like  RPRMLS 21 = 2  to set RPRMLS(21) = 2.0
            IVALUE = RVALUE
            IF (IVALUE .GE. 1  .AND. IVALUE .LE. MXPARM) THEN
               READ (KEY3, '(BN, E16.0)') RPRMLS(IVALUE)
            ELSE
               WRITE(NOUT, 2400) IVALUE
            END IF
         ELSE IF (KEY .EQ. 'SAVE        ') THEN
              IF (KEY2.EQ. 'RUN         ') NSAVE  = RVALUE
              IF (KEY2.EQ. 'FREQUENCY   ') KSAVE  = RVALUE
         ELSE IF (KEY .EQ. 'START       ') THEN
              IF (KEY2.EQ. 'CONSTRAINTS ') JVRFY3 = RVALUE
              IF (KEY2.EQ. 'OBJECTIVE   ') JVRFY1 = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'STEP        ') THEN
            DXLIM  = RVALUE
         ELSE IF (KEY .EQ. 'STOP        ') THEN
              IF (KEY2.EQ. 'CONSTRAINTS ') JVRFY4 = RVALUE
              IF (KEY2.EQ. 'OBJECTIVE   ') JVRFY2 = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'UPPER       ') THEN
            BNDUPP = RVALUE
         ELSE IF (KEY .EQ. 'VARIABLES   ') THEN
            NN     = RVALUE
         ELSE IF (KEY .EQ. 'VERIFY      ') THEN
              IF (KEY2.EQ. 'OBJECTIVE   ') LVERFY =  1
              IF (KEY2.EQ. 'CONSTRAINTS ') LVERFY =  2
              IF (KEY2.EQ. 'NO          ') LVERFY = -1
              IF (KEY2.EQ. 'YES         ') LVERFY =  3
              IF (KEY2.EQ. 'GRADIENTS   ') LVERFY =  3
              IF (KEY2.EQ. 'LEVEL       ') LVERFY =  RVALUE
              IF (LOC2.EQ.  0            ) LVERFY =  3
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

*     End of NPKEY

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nploc ( n, nclin, ncnln, nctotl, litotl, lwtotl)

      implicit           double precision (a-h,o-z)

************************************************************************
*     nploc   allocates the addresses of the work arrays for npcore and
*     lscore.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version   14-February-1985.
*     This version of  NPLOC  dated 12-July-1986.
************************************************************************
      common    /sol1cm/ nout
      common    /sol3cm/ lennam, ldt   , ncolt , ldq

      parameter         (lenls = 20)
      common    /sol1ls/ locls(lenls)

      parameter         (lennp = 35)
      common    /sol1np/ locnp(lennp)

      logical            npdbg
      parameter         (ldbg = 5)
      common    /npdebg/ inpdbg(ldbg), npdbg

      miniw     = litotl + 1
      minw      = lwtotl + 1

*     Assign array lengths that depend upon the problem dimensions.

      if (nclin + ncnln .eq. 0) then
         lent      = 0
         lenzy     = 0
      else
         lent  = ldt *ncolt
         lenzy = ldq*ldq
      end if

      if (ncnln .eq. 0) then
         lenaqp = 0
      else
         lenaqp = (nclin + ncnln)*n
      end if

      lkactv    = miniw
      lkx       = lkactv + n
      lneedc    = lkx    + n
      liperm    = lneedc + ncnln
      miniw     = liperm + nctotl

      lhfrwd    = minw
      lhctrl    = lhfrwd + n
      lanorm    = lhctrl + n
      lqpgq     = lanorm + nclin + ncnln
      lgq       = lqpgq  + n
      lrlam     = lgq    + n
      lt        = lrlam  + n
      lq        = lt     + lent
      minw      = lq     + lenzy

      locls( 1) = lkactv
      locls( 2) = lanorm
      locls( 8) = lqpgq
      locls( 9) = lgq
      locls(10) = lrlam
      locls(11) = lt
      locls(12) = lq

*     Assign the addresses for the workspace arrays used by  NPIQP.

      lqpadx    = minw
      lqpdx     = lqpadx + nclin + ncnln
      lrpq      = lqpdx  + n
      lrpq0     = lrpq   + n
      lqphz     = lrpq0  + n
      lwtinf    = lqphz  + n
      lwrk1     = lwtinf + nctotl
      lqptol    = lwrk1  + nctotl
      minw      = lqptol + nctotl

      locls( 3) = lqpadx
      locls( 4) = lqpdx
      locls( 5) = lrpq
      locls( 6) = lrpq0
      locls( 7) = lqphz
      locls(13) = lwtinf
      locls(14) = lwrk1
      locls(15) = lqptol

*     Assign the addresses for arrays used in NPCORE.

      laqp      = minw
      ladx      = laqp   + lenaqp
      lbl       = ladx   + nclin  + ncnln
      lbu       = lbl    + nctotl
      ldx       = lbu    + nctotl
      lgq1      = ldx    + n
      lfeatl    = lgq1   + n
      lx1       = lfeatl + nctotl
      lwrk2     = lx1    + n
      minw      = lwrk2  + nctotl

      locnp( 1) = lkx
      locnp( 2) = liperm
      locnp( 3) = laqp
      locnp( 4) = ladx
      locnp( 5) = lbl
      locnp( 6) = lbu
      locnp( 7) = ldx
      locnp( 8) = lgq1
      locnp(10) = lfeatl
      locnp(11) = lx1
      locnp(12) = lwrk2

      lcs1      = minw
      lcs2      = lcs1   + ncnln
      lc1mul    = lcs2   + ncnln
      lcmul     = lc1mul + ncnln
      lcjdx     = lcmul  + ncnln
      ldlam     = lcjdx  + ncnln
      ldslk     = ldlam  + ncnln
      lrho      = ldslk  + ncnln
      lwrk3     = lrho   + ncnln
      lslk1     = lwrk3  + ncnln
      lslk      = lslk1  + ncnln
      minw      = lslk   + ncnln

      locnp(13) = lcs1
      locnp(14) = lcs2
      locnp(15) = lc1mul
      locnp(16) = lcmul
      locnp(17) = lcjdx
      locnp(18) = ldlam
      locnp(19) = ldslk
      locnp(20) = lrho
      locnp(21) = lwrk3
      locnp(22) = lslk1
      locnp(23) = lslk
      locnp(24) = lneedc

      lcjac     = minw
      lgrad     = lcjac  + ncnln*n
      minw      = lgrad  + n

      locnp(25) = lhfrwd
      locnp(26) = lhctrl
      locnp(27) = lcjac
      locnp(28) = lgrad

      litotl    = miniw - 1
      lwtotl    = minw  - 1

      return

*     end of  nploc
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npmrt ( feasqp, n, nclin, ncnln,
     $                   objalf, grdalf, qpcurv,
     $                   istate,
     $                   cjdx, cmul, cs,
     $                   dlam, rho, violn,
     $                   work1, work2 )

      implicit           double precision (a-h,o-z)

      logical            feasqp

      integer            istate(*)

      double precision   cjdx(*), cmul(*), cs(*),
     $                   dlam(*), rho(*), violn(*)
      double precision   work1(*), work2(*)

************************************************************************
*     npmrt   computes the value and directional derivative of the
*     augmented Lagrangian merit function.  The penalty parameters
*     rho(j) are boosted if the directional derivative of the resulting
*     augmented Lagrangian function is not sufficiently negative.  If
*     rho needs to be increased,  the perturbation with minimum two-norm
*     is found that gives a directional derivative equal to  - p'Hp.
*
*     Systems Optimization Laboratory, Stanford University, California.
*     Original version written  27-May-1985.
*     This version of  NPMRT  dated 14-November-1985.
************************************************************************
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      common    /sol1cm/ nout

      logical            incrun
      common    /sol6np/ rhomax, rhonrm, rhodmp, scale, incrun

      logical            npdbg
      parameter         (ldbg = 5)
      common    /npdebg/ inpdbg(ldbg), npdbg

      logical            boost , overfl
      external           ddiv  , ddot  , dnrm2
      intrinsic          abs   , max   , min   , sqrt
      parameter        ( zero   = 0.0d+0, half = 0.5d+0, one = 1.0d+0 )
      parameter        ( two    = 2.0d+0                              )

      if (ncnln .eq. 0) return

      rtmin  = wmach(6)

      objalf = objalf - ddot  ( ncnln, cmul, 1, cs, 1 )
      grdalf = grdalf - ddot  ( ncnln, dlam, 1, cs, 1 )

      call dcopy ( ncnln, cs, 1, work1, 1 )

      if (.not. feasqp) then
         nplin  = n + nclin

         do 100 i = 1, ncnln
            if (istate(nplin+i) .lt. 0  .or.  violn(i) .ne. zero)
     $         work1(i) = - cjdx(i)
  100    continue
      end if

      grdalf = grdalf + ddot  ( ncnln, work1, 1, cmul, 1 )

      if (npdbg  .and.  inpdbg(1) .gt. 0)
     $   write (nout, 1000) qpcurv, grdalf

      if (feasqp) then

*        Find the quantities that define  rhomin, the vector of minimum
*        two-norm such that the directional derivative is one half of
*        approximate curvature   - (dx)'H(dx).

         do 350, i = 1, ncnln
            if (abs( cs(i) ) .le. rtmin) then
               work2(i) = zero
            else
               work2(i) = cs(i)**2
            end if
  350    continue

         qnorm  = dnrm2 ( ncnln, work2, 1 )
         tscl   = ddiv  ( grdalf + half*qpcurv, qnorm, overfl )
         if (abs( tscl ) .le. rhomax  .and.  .not. overfl) then
*           ------------------------------------------------------------
*           Bounded  rhomin  found.  The final value of  rho(J)  will
*           never be less than  rhomin(j).  If the  QP  was feasible,  a
*           trial value  rhonew  is computed that is equal to the
*           geometric mean of the previous  rho  and a damped value of
*           rhomin.  The new  rho  is defined as  rhonew  if it is less
*           than half the previous  rho  and greater than  rhomin.
*           ------------------------------------------------------------
            scale  = one
            do 400, i = 1, ncnln
               rhomin = max(  (work2(i)/qnorm)*tscl, zero )
               rhoi   = rho(i)

               rhonew = sqrt( rhoi*(rhodmp + rhomin) )
               if (rhonew .lt. half*rhoi  ) rhoi = rhonew
               if (rhoi   .lt.      rhomin) rhoi = rhomin
               rho(i) = rhoi
  400       continue

            rho1   = rhonrm
            rhonrm = dnrm2 ( ncnln, rho, 1 )

*           ------------------------------------------------------------
*           If  incrun = true,  there has been a run of iterations in
*           which the norm of  rho  has not decreased.  Conversely,
*           incrun = false  implies that there has been a run of
*           iterations in which the norm of rho has not increased.  If
*           incrun changes during this iteration the damping parameter
*           rhodmp is increased by a factor of two.  This ensures that
*           rho(j) will oscillate only a finite number of times.
*           ------------------------------------------------------------
            boost  = .false.
            if (      incrun  .and.  rhonrm .lt. rho1) boost = .true.
            if (.not. incrun  .and.  rhonrm .gt. rho1) boost = .true.
            if (boost) then
               rhodmp = two*rhodmp
               incrun = .not. incrun
            end if
         end if

         if (npdbg  .and.  inpdbg(2) .gt. 0)
     $      write (nout, 1200) (rho(l), l=1,ncnln)

      else

*        The  QP  was infeasible.  Do not alter the penalty parameters,
*        but compute the scale factor so that the constraint violations
*        are reduced.

         call ddscl ( ncnln, rho, 1, work1, 1 )
         pterm2 = ddot  ( ncnln, work1, 1, cs, 1 )

         scale  = rhomax
         tscl   = ddiv  ( grdalf, pterm2, overfl )
         if (tscl .gt. scale  .and.  tscl .le. rhomax/(one+rhonrm)
     $                        .and.  .not. overfl)
     $      scale = tscl

         call dcopy ( ncnln, cs, 1, work1, 1 )
      end if

*     ------------------------------------------------------------------
*     Compute the new value and directional derivative of the
*     merit function.
*     ------------------------------------------------------------------
      call ddscl ( ncnln, rho, 1, work1, 1 )

      pterm  = ddot  ( ncnln, work1, 1, cs, 1 )
      objalf = objalf + half*scale*pterm

      if (feasqp)
     $  pterm2 = pterm

      grdalf = grdalf -      scale*pterm2

      if (npdbg  .and.  inpdbg(1) .gt. 0)
     $   write (nout, 1100) scale, rhonrm, grdalf

      return

 1000 format(/ ' //npmrt //        qpcurv        grdalf '
     $       / ' //npmrt //', 1p, 2e14.2 )
 1100 format(/ ' //npmrt //         scale        rhonrm        grdalf '
     $       / ' //npmrt //', 1p, 3e14.2 )
 1200 format(/ ' //npmrt //  Penalty parameters =       '/ (1p, 5e15.6))

*     end of  npmrt
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npoptn( string )
      character*(*)      string

************************************************************************
*     npoptn  loads the option supplied in string into the relevant
*     element of iprmls, rprmls, iprmnp or rprmnp.
************************************************************************

      logical             newopt
      common     /sol7np/ newopt
      save       /sol7np/

      double precision    wmach(15)
      common     /solmch/ wmach
      save       /solmch/

      external            mchpar
      character*16        key
      character*72        buffer
      logical             first , prnt
      save                first , nout  , prnt
      data                first /.true./

*     If first time in, set NOUT.
*     NEWOPT is true first time into NPFILE or NPOPTN
*     and just after a call to an optimization routine.
*     PRNT is set to true whenever NEWOPT is true.

      if (first) then
         first  = .false.
         newopt = .true.
         call mchpar()
         nout   =  wmach(11)
      end if
      buffer = string

*     Call NPKEY to decode the option and set the parameter value.
*     If NEWOPT is true, reset PRNT and test specially for NOLIST.

      if (newopt) then
         newopt = .false.
         prnt   = .true.
         call npkey ( nout, buffer, key )

         if (key .eq. 'NOLIST') then
            prnt   = .false.
         else
            write (nout, '(// a / a /)')
     $         ' Calls to NPOPTN',
     $         ' ---------------'
            write (nout, '( 6x, a )') buffer
         end if
      else
         if (prnt)
     $      write (nout, '( 6x, a )') buffer
         call npkey ( nout, buffer, key )

         if (key .eq.   'LIST') prnt = .true.
         if (key .eq. 'NOLIST') prnt = .false.
      end if

      return

*     end of  npoptn
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE NPPRT ( KTCOND, CONVRG, LSUMRY, MSGNP, MSGQP,
     $                   LDR, LDT, N, NCLIN, NCNLN,
     $                   NCTOTL, NACTIV, LINACT, NLNACT, NZ, NFREE,
     $                   MAJIT0, MAJITS, MINITS, ISTATE, ALFA, NFUN,
     $                   CONDHZ, CONDH, CONDT, OBJALF, OBJF,
     $                   GFNORM, GZNORM, CVNORM,
     $                   AX, C, R, T, VIOLN, X, WORK )

      IMPLICIT           DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*25       LSUMRY
      LOGICAL            KTCOND(2), CONVRG
      INTEGER            ISTATE(NCTOTL)
      DOUBLE PRECISION   AX(*), C(*), R(LDR,*), T(LDT,*), VIOLN(*)
      DOUBLE PRECISION   X(N)
      DOUBLE PRECISION   WORK(N)

************************************************************************
*     NPPRT  prints various levels of output for NPCORE.
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
*        ge  20        objective function,  x,  Ax  and  c.
*
*        ge  30        diagonals of  T  and  R.
*
*     Debug print is performed depending on the logical variable NPDBG.
*     NPDBG is set true when IDBG major iterations have been performed.
*     At this point,  printing is done according to a string of binary
*     digits of the form CLSVT (stored in the integer array INPDBG).
*
*     C  set 'on' gives detailed information from the checking routines.
*     L  set 'on' gives information from the linesearch.
*     S  set 'on' gives information from the maximum step routine NPALF.
*     V  set 'on' gives various vectors in  NPCORE  and its auxiliaries.
*     T  set 'on' gives a trace of which routine was called and an
*                 indication of the progress of the run.
*
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 66 version written November-1982.
*     This version of  NPPRT  dated  25-Nov-1989.
************************************************************************
      COMMON    /SOL1CM/ NOUT

      LOGICAL            INCRUN
      COMMON    /SOL6NP/ RHOMAX, RHONRM, RHODMP, SCALE, INCRUN

      LOGICAL            NPDBG
      PARAMETER         (LDBG = 5)
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG

      LOGICAL            PHEAD
      EXTERNAL           DNRM2

      IF (MSGNP .GE. 20) WRITE (NOUT, 1000) MAJITS

      IF (MSGNP  .GE. 5) THEN
*        ---------------------------------------------------------------
*        Print heading and terse line.
*        ---------------------------------------------------------------
         PHEAD = MSGQP .GT. 0  .OR.  MAJITS .EQ. MAJIT0

         IF (NCNLN .EQ. 0) THEN
            IF (PHEAD) WRITE (NOUT, 1100)
            WRITE (NOUT, 1300) MAJITS, MINITS, ALFA, NFUN, OBJALF,
     $                         N-NFREE, LINACT, NZ,
     $                         GFNORM, GZNORM, CONDH, CONDHZ, CONDT,
     $                         CONVRG, KTCOND(1), KTCOND(2), LSUMRY

         ELSE
            IF (PHEAD) WRITE (NOUT, 1110)
            WRITE (NOUT, 1310) MAJITS, MINITS, ALFA, NFUN, OBJALF,
     $                         N-NFREE, LINACT, NLNACT, NZ,
     $                         GFNORM, GZNORM, CONDH, CONDHZ, CONDT,
     $                         CVNORM, SCALE*RHONRM,
     $                         CONVRG, KTCOND(1), KTCOND(2), LSUMRY
         END IF

         IF (MSGNP .GE. 20) THEN
            IF (NCNLN .EQ. 0) THEN
               WRITE (NOUT, 1400) OBJF
            ELSE
               CVIOLS = DNRM2 ( NCNLN, VIOLN, 1 )
               WRITE (NOUT, 1410) OBJF, CVIOLS
            END IF

*           ------------------------------------------------------------
*           Print the constraint values.
*           ------------------------------------------------------------
            WRITE (NOUT, 2000)
            WRITE (NOUT, 2100) (X(J), ISTATE(J), J=1,N)
            IF (NCLIN .GT. 0)
     $         WRITE (NOUT, 2200) (AX(K), ISTATE(N+K),       K=1,NCLIN )
            IF (NCNLN .GT. 0)
     $         WRITE (NOUT, 2300) (C(K) , ISTATE(N+NCLIN+K), K=1,NCNLN )

            IF (MSGNP .GE. 30) THEN
*              ---------------------------------------------------------
*              Print the diagonals of  T  and  R.
*              ---------------------------------------------------------
               INCT   = LDT - 1
               IF (NACTIV .GT. 0) THEN
                  CALL DCOPY( NACTIV, T(NACTIV,NZ+1), INCT, WORK, 1 )
                  WRITE (NOUT, 3000) (WORK(J), J=1,NACTIV)
               END IF
               WRITE (NOUT, 3100) (R(J,J), J=1,N)
            END IF
         END IF
      END IF

      IF (MSGNP .GE. 20) WRITE (NOUT, 5000)

      LSUMRY(1:2) = '  '
      LSUMRY(4:5) = '  '

      RETURN

 1000 FORMAT(/// ' Major iteration', I5
     $       /   ' ====================' )
 1100 FORMAT(//  '  Itn', ' ItQP', '     Step',
     $           '  Nfun', '     Objective', ' Bnd', ' Lin', '  Nz',
     $           '  Norm Gf', '  Norm Gz', '  Cond H', ' Cond Hz',
     $           '  Cond T', ' Conv' )
 1110 FORMAT(//  '  Itn', ' ItQP', '     Step',
     $           '  Nfun', '         Merit', ' Bnd', ' Lin',
     $           ' Nln', '  Nz',
     $           '  Norm Gf', '  Norm Gz', '  Cond H', ' Cond Hz',
     $           '  Cond T' , '   Norm C', '  Penalty', ' Conv' )
 1300 FORMAT(2I5, 1PE9.1, I6, 1PE14.6, 3I4, 1P, 2E9.1, 1P, 3E8.0,
     $                        1X, L1, 1X, 2L1, A5 )
 1310 FORMAT(2I5, 1PE9.1, I6, 1PE14.6, 4I4, 1P, 2E9.1, 1P, 3E8.0,
     $            1P, 2E9.1,    1X, L1, 1X, 2L1, A5 )
 1400 FORMAT(/ ' Nonlinear objective value = ', 1PE15.6 )
 1410 FORMAT(/ ' Nonlinear objective value = ', 1PE15.6, '   Norm of',
     $         ' the nonlinear constraint violations = ', 1PE15.6 )
 2000 FORMAT(/ ' Values of the constraints and their predicted status'
     $       / ' ----------------------------------------------------')
 2100 FORMAT(/ ' Variables                  '/ (1X, 5(1PE15.6, I4)))
 2200 FORMAT(/ ' General linear constraints '/ (1X, 5(1PE15.6, I4)))
 2300 FORMAT(/ ' Nonlinear constraints      '/ (1X, 5(1PE15.6, I4)))
 3000 FORMAT(/ ' Diagonals of  T  =         '
     $       /   (1P, 5E15.6))
 3100 FORMAT(/ ' Diagonals of  R  =         '
     $       /   (1P, 5E15.6))
 5000 FORMAT(  ' ==================================================',
     $         '======================================='///)

*     End of  NPPRT .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE NPRSET( UNITQ,
     $                   N, NFREE, NZ, LDQ, LDR,
     $                   IPERM, KX,
     $                   GQ, R, ZY, WORK, QRWORK )

      IMPLICIT           DOUBLE PRECISION (A-H,O-Z)
      LOGICAL            UNITQ
      INTEGER            IPERM(N), KX(N)
      DOUBLE PRECISION   GQ(N), R(LDR,*), ZY(LDQ,*)
      DOUBLE PRECISION   WORK(N), QRWORK(2*N)

************************************************************************
*     NPRSET  bounds the condition estimator of the transformed Hessian.
*     On exit, R is of the form
*                  ( DRz   0     )
*                  (  0  sigma*I )
*     where D is a diagonal matrix such that DRz has a bounded condition
*     number,  I is the identity matrix and sigma  is the geometric mean
*     of the largest and smallest elements of DRz. The QR factorization
*     with interchanges is used to give diagonals of DRz that are
*     decreasing in modulus.
*
*     Systems Optimization Laboratory, Stanford University.
*
*     Original version of NPRSET dated  4-August-1986.
*     This version dated  26-Jun-1987.
************************************************************************

      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL6CM/ RCNDBD, RFROBN, DRMAX, DRMIN

      LOGICAL            NPDBG
      PARAMETER         (LDBG   = 5)
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG

      INTRINSIC          MAX   , MIN   , LOG   , REAL  , SQRT
      EXTERNAL           DDIV  , DDOT  , DNORM , DNRM2 , IDRANK
      PARAMETER        ( ZERO   =0.0D+0, HALF =0.5D+0, ONE    =1.0D+0 )

*     ==================================================================
*     Bound the condition estimator of Q'HQ.
*     The scheme used here reduces the modulus of the larger
*     diagonals while increasing the modulus of the smaller ones.
*     ==================================================================
      IF (NZ .GT. 1) THEN
*        ---------------------------------------------------------------
*        Refactorize Rz.  Interchanges are used to give diagonals
*        of decreasing magnitude.
*        ---------------------------------------------------------------
         DO 100, J = 1, NZ-1
            CALL DLOAD ( NZ-J, ZERO, R(J+1,J), 1 )
  100    CONTINUE

         CALL DGEQRP( 'Column iterchanges', NZ, NZ, R, LDR,
     $                WORK, IPERM, QRWORK, INFO )

         DO 110, J = 1, NZ
            JMAX = IPERM(J)
            IF (JMAX .GT. J) THEN
               IF (UNITQ) THEN
                  JSAVE    = KX(JMAX)
                  KX(JMAX) = KX(J)
                  KX(J)    = JSAVE
               ELSE
                  CALL DSWAP ( NFREE, ZY(1,JMAX), 1, ZY(1,J), 1 )
               END IF

               GJMAX    = GQ(JMAX)
               GQ(JMAX) = GQ(J)
               GQ(J)    = GJMAX
            END IF
  110    CONTINUE
      END IF

      DRGM  = ONE

      IF (NZ .GT. 0) THEN
         NRANK = IDRANK( NZ, R, LDR+1, ONE/RCNDBD )
         DRGM  = HALF*SQRT(ABS( R(1,1)*R(NRANK,NRANK) ))
         DRGS  = ABS( R(1,1) ) / RCNDBD

         IF (NZ .GT. NRANK) THEN
            DO 120, J = NRANK+1, NZ
               CALL DLOAD ( J-1, ZERO, R(1,J), 1 )
  120       CONTINUE
            CALL DLOAD ( NZ-NRANK, DRGS, R(NRANK+1,NRANK+1), LDR+1 )
         END IF
      END IF

*     ------------------------------------------------------------------
*     Reset the range-space partition of the Hessian.
*     ------------------------------------------------------------------
      IF (NZ .LT. N) THEN
         DO 130, J = NZ+1, N
            CALL DLOAD ( J, ZERO, R(1,J), 1 )
  130    CONTINUE
         CALL DLOAD ( N-NZ, DRGM, R(NZ+1,NZ+1), LDR+1 )
      END IF

*     Recompute the Frobenius norm of R.

      SCLE  = SQRT(REAL(N - NZ))*DRGM
      SUMSQ = ONE
      DO 140, J = 1, NZ
         CALL DSSQ  ( J, R(1,J), 1, SCLE, SUMSQ )
  140 CONTINUE
      RFROBN = DNORM( SCLE, SUMSQ )

      RETURN

*     End of  NPRSET.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE NPSAVR( UNITQ, N, NCLIN, NCNLN, LDR, LDQ,
     $                   NFREE, NSAVE, ITER, ISTATE, KX,
     $                   HFORWD, HCNTRL,
     $                   CMUL, R, RHO, X, ZY )

      IMPLICIT           DOUBLE PRECISION (A-H,O-Z)
      LOGICAL            UNITQ
      INTEGER            ISTATE(N+NCLIN+NCNLN), KX(N)
      DOUBLE PRECISION   R(LDR,*), X(N), ZY(LDQ,*)
      DOUBLE PRECISION   HFORWD(*), HCNTRL(*)
      DOUBLE PRECISION   CMUL(*), RHO(*)

C***********************************************************************
C     NPSAVR  saves the details of this run on unit NSAVE.
C
C     Original version   24-Nov-89.
C     This version of  NPSAVR  dated  1-Dec-89.
C***********************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4NP/ LVLDIF, NCDIFF, NFDIFF, LFDSET
      LOGICAL            INCRUN
      COMMON    /SOL6NP/ RHOMAX, RHONRM, RHODMP, SCALE, INCRUN

      IF (NSAVE .LE. 0) RETURN

*     We always save the computed difference intervals, if defined.

      IF (LVLDIF .GT. 0  .AND.  LFDSET .EQ. 0) THEN
         LFDNEW = 2
      ELSE 
         LFDNEW = LFDSET
      END IF
 
      WRITE( NSAVE, 1000 ) ITER, NFREE, LFDNEW, LVLDIF, UNITQ
      DO 110, J = 1, N
         WRITE( NSAVE, 1010 ) J, KX(J), ISTATE(J), X(J)
  110 CONTINUE

      DO 120, J = N+1, N+NCLIN
         WRITE( NSAVE, 1020 ) J, ISTATE(J)
  120 CONTINUE

      IF (NCNLN .GT. 0) THEN
         K = 1
         DO 130, J = N+NCLIN+1, N+NCLIN+NCNLN
            WRITE( NSAVE, 1030 ) J, ISTATE(J), CMUL(K), RHO(K)
            K = K + 1
  130    CONTINUE

         WRITE( NSAVE, 1040 ) RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
      END IF

*     ------------------------------------------------------------------
*     Write  Q(free)  and the factor of  Q'HQ.
*     ------------------------------------------------------------------
      IF (.NOT. UNITQ) THEN
         DO 210, J = 1, NFREE
            DO 200, I = 1, NFREE
              WRITE( NSAVE, 2000 ) I, J, ZY(I,J)
  200       CONTINUE
  210    CONTINUE
      END IF

      DO 240, J = 1, N
         DO 230, I = 1, J
           WRITE( NSAVE, 2000 ) I, J, R(I,J)
  230    CONTINUE
  240 CONTINUE

*     ------------------------------------------------------------------
*     Read the finite-difference intervals.  
*     ------------------------------------------------------------------
      IF (LVLDIF .GT. 0) THEN
         IF (LFDNEW .NE. 1) THEN
            DO 250, J = 1, N
               WRITE( NSAVE, 1050 ) J, HFORWD(J), HCNTRL(J)
  250       CONTINUE
         END IF
      END IF                                            
             
      WRITE( NOUT, 4000 ) NSAVE, ITER
      IF (NSAVE .NE. NOUT) REWIND NSAVE
      RETURN

 1000 FORMAT(2I8, 1X, 2I2, 1X, L1 )
 1010 FORMAT(2I8, 1X,  I2, 1P, 2E24.14 )
 1020 FORMAT( I8, 1X,  I2 )
 1030 FORMAT( I8, 1X,  I2, 1P, 2E24.14 )
 1040 FORMAT( 1P, 4E24.14, 1X, L1 )
 1050 FORMAT( I8, 1X, 1P, 2E24.14 )
 2000 FORMAT(2I8, 1X, 1P,  E24.14 )
 4000 FORMAT(  ' CURRENT RUN saved on file', I4, '    Itn = ', I7 )

*     End of  NPSAVR.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE NPSETX( UNITQ,
     $                   NCQP, NACTIV, NFREE, NZ,
     $                   N, NLNX, NCTOTL, LDQ, LDAQP, LDR, LDT,
     $                   ISTATE, KACTIV, KX,
     $                   DXNORM, GDX,
     $                   AQP, ADX, BL, BU, RPQ, RPQ0, DX, GQ,
     $                   R, T, ZY, WORK )

      IMPLICIT           DOUBLE PRECISION (A-H,O-Z)
      LOGICAL            UNITQ
      INTEGER            ISTATE(NCTOTL), KACTIV(N), KX(N)
      DOUBLE PRECISION   AQP(LDAQP,*), ADX(*), BL(NCTOTL), BU(NCTOTL),
     $                   RPQ(NLNX), RPQ0(NLNX), GQ(N), R(LDR,*),
     $                   T(LDT,*), ZY(LDQ,*), DX(N), WORK(N)

************************************************************************
*     NPSETX   defines a point which lies on the initial working set for
*     the QP subproblem.  This routine is a similar to LSSETX except
*     that advantage is taken of the fact that the initial estimate of
*     the solution of the least-squares subproblem is zero.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     Level 2 BLAS added 12-June-1986.
*     This version of NPSETX dated 11-June-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT

      LOGICAL            NPDBG
      PARAMETER         (LDBG = 5)
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG

      EXTERNAL           DDOT, DNRM2
      INTRINSIC          ABS , MIN
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      NFIXED = N - NFREE

      GDX    = ZERO
      CALL DLOAD ( N   , ZERO, DX  , 1 )
      CALL DLOAD ( NLNX, ZERO, RPQ , 1 )
      CALL DLOAD ( NLNX, ZERO, RPQ0, 1 )

      IF (NACTIV + NFIXED .GT. 0) THEN

*        Set  work = residuals for constraints in the working set.
*        Solve for  dx,  the smallest correction to  x  that gives a
*        point on the constraints in the working set.
*        Set the fixed variables on their bounds,  solve the triangular
*        system  T*(dxy) = residuals,  and define  dx = Y*(dxy).
*        Use  (dxy)  to update  d(=Pr)  as  d = d - R'(  0  ).
*                                                     ( dxy )

         DO 100 I = 1, NFIXED
            J   = KX(NFREE+I)
            IF (ISTATE(J) .LE. 3) THEN
               BND   = BL(J)
               IF (ISTATE(J) .EQ. 2) BND = BU(J)
               DX(J) = BND
               WORK(NFREE+I) = BND
            ELSE
               WORK(NFREE+I) = ZERO
            END IF
  100    CONTINUE

         DO 110 I = 1, NACTIV
            K   = KACTIV(I)
            J   = N + K
            BND = BL(J)
            IF (ISTATE(J) .EQ. 2) BND = BU(J)
            WORK(NZ+I) = BND - DDOT  ( N, AQP(K,1), LDAQP, DX, 1 )
  110    CONTINUE

         IF (NACTIV .GT. 0)
     $      CALL CMTSOL( 1, LDT, NACTIV, T(1,NZ+1), WORK(NZ+1) )
         CALL DCOPY ( NACTIV+NFIXED, WORK(NZ+1), 1, DX(NZ+1), 1 )
         IF (NZ .GT. 0)
     $      CALL DLOAD ( NZ, ZERO, DX, 1 )

         GDX  = DDOT  ( NACTIV+NFIXED, GQ(NZ+1), 1, DX(NZ+1), 1 )

         IF (NZ .LT. N) THEN
            CALL DGEMV ('N', NZ, N-NZ, -ONE, R(1,NZ+1), LDR,
     $                  DX(NZ+1), 1, ONE, RPQ, 1 )
            IF (NZ .LT. NLNX) THEN
               NR  = LDR
               IF (NZ+1 .EQ. N) NR = 1
               CALL DCOPY ( NLNX-NZ, DX(NZ+1), 1, RPQ(NZ+1), 1 )
               CALL DSCAL ( NLNX-NZ, (-ONE),      RPQ(NZ+1), 1 )
               CALL DTRMV ( 'U', 'N', 'N', NLNX-NZ, R(NZ+1,NZ+1), NR,
     $                      RPQ(NZ+1), 1 )
               IF (NLNX .LT. N) THEN
                  NR = LDR
                  IF (NLNX+1 .EQ. N) NR = N - NZ
                  CALL DGEMV( 'N', NLNX-NZ, N-NLNX, -ONE,R(NZ+1,NLNX+1),
     $                        NR, DX(NLNX+1), 1, ONE, RPQ(NZ+1), 1 )
               END IF
            END IF
         END IF

         CALL CMQMUL( 2, N, NZ, NFREE, LDQ, UNITQ, KX, DX, ZY, WORK )
      END IF

*     ------------------------------------------------------------------
*     Compute the 2-norm of  DX.
*     Initialize  A*DX.
*     ------------------------------------------------------------------
      DXNORM  = DNRM2 ( N, DX, 1 )
      IF (NCQP .GT. 0)
     $   CALL DGEMV ( 'N', NCQP, N, ONE, AQP, LDAQP, DX, 1, ZERO,ADX,1)

      IF (NPDBG  .AND.  INPDBG(2) .GT. 0)
     $   WRITE (NOUT, 1200) (DX(J), J = 1, N)

      RETURN

 1200 FORMAT(/ ' //NPSETX// Variables after NPSETX ... '/ (5G12.3))

*     End of  NPSETX.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npsrch( needfd, inform, n, ncnln,
     $                   ldcj, ldcju, nfun, ngrad,
     $                   needc, confun, objfun,
     $                   alfa, alfbnd, alfmax, alfsml, dxnorm,
     $                   epsrf, eta, gdx, grdalf, glf1, glf,
     $                   objf, objalf, qpcurv, xnorm,
     $                   c, c2, cjac, cjacu, cjdx, cjdx2,
     $                   cmul1, cmul, cs1, cs, dx, dlam, dslk,
     $                   grad, gradu, qpmul, rho, slk1, slk, x1, x,
     $                   work, w, lenw )

      implicit           double precision (a-h,o-z)
      logical            needfd
      integer            needc(*)
      double precision   dx(n), grad(n), gradu(n), x1(n), x(n)
      double precision   c(*), c2(*), cjac(ldcj,*), cjacu(ldcju,*),
     $                   cjdx(*), cjdx2(*), cmul1(*), cmul(*)
      double precision   cs1(*), cs(*), dlam(*), dslk(*), qpmul(*),
     $                   rho(*), slk1(*), slk(*),  work(*)
      double precision   w(lenw)
      external           objfun, confun

C***********************************************************************
C     npsrch finds the steplength alfa that gives sufficient decrease in
C     the augmented Lagrangian merit function.
C
C     On exit, if inform = 1, 2 or 3,  alfa will be a nonzero steplength
C     with an associated merit function value  objalf  which is lower
C     than that at the base point. If  inform = 4, 5, 6, 7 or 8,  alfa
C     is zero and  objalf will be the merit value at the base point.
C
C     Original version written  27-May-1985.
C     Level 2 BLAS added 12-June-1986.
C     This version of npsrch dated  22-Oct-91.
C***********************************************************************

      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      common    /sol1cm/ nout
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset
      logical            incrun
      common    /sol6np/ rhomax, rhonrm, rhodmp, scale, incrun

      logical            npdbg
      parameter        ( ldbg = 5 )
      common    /npdebg/ inpdbg(ldbg), npdbg

      logical            debug , done  , first , imprvd
      external           ddot  , dnrm2
      intrinsic          abs   , max   , min   , sqrt
      parameter        ( zero   =0.0d+0, half   =0.5d+0, one   =1.0d+0 )
      parameter        ( two    =2.0d+0                                )
      parameter        ( tolg   =1.0d-1                                )
      parameter        ( rmu    =1.0d-4                                )

      epsmch = wmach(3)

      if (.not. needfd  .and.  ncnln .gt. 0)
     $   cs1jdx = ddot( ncnln, cs1, 1, cjdx, 1 )

*     ------------------------------------------------------------------
*     Set the input parameters and tolerances for srchc and srchq.
*
*     tolrx   is the tolerance on relative changes in dx resulting from
*             changes in alfa.
*
*     tolax   is the tolerance on absolute changes in dx resulting from
*             changes in alfa.
*
*     tolabs  is the tolerance on absolute changes in alfa.
*
*     tolrel  is the tolerance on relative changes in alfa.
*
*     toltny  is the magnitude of the smallest allowable value of alfa.
*             if  M(tolabs) - M(0) .gt. epsaf,  the linesearch tries
*             steps in the range  toltny .le. alfa .le. tolabs.
*     ------------------------------------------------------------------
      nstate = 0
      debug  = npdbg  .and.  inpdbg(4) .gt. 0

      if (needfd) then
         maxf = 15
      else
         maxf = 10
      end if

      epsaf  = epsrf*(one + abs( objalf ))
      tolax  = epspt8
      tolrx  = epspt8

      if (tolrx*xnorm + tolax .lt. dxnorm*alfmax) then
         tolabs = (tolrx*xnorm + tolax) /  dxnorm
      else
         tolabs = alfmax
      end if
      tolrel = max( tolrx , epsmch )

      t      = zero
      do 10, j = 1, n
         s = abs( dx(j) )
         q = abs( x(j) )*tolrx + tolax
         if (s .gt. t*q) t = s / q
   10 continue

      if (t*tolabs .gt. one) then
         toltny = one / t
      else
         toltny = tolabs
      end if

      oldf   = objalf
      oldg   = grdalf
      alfbst = zero
      fbest  = zero
      gbest  = (one - rmu)*oldg
      targtg = (rmu - eta)*oldg
      g0     = gbest

      if (ncnln .gt. 0) call iload ( ncnln, (1), needc, 1 )

      if (needfd) then
         mode = 0
      else
         mode = 2
      end if

      first  = .true.

*     ------------------------------------------------------------------
*     Commence main loop, entering srchc or srchq two or more times.
*     first = true for the first entry,  false for subsequent entries.
*     done  = true indicates termination, in which case the value of
*     inform gives the result of the search.
*     inform = 1 if the search is successful and alfa < alfmax.
*            = 2 if the search is successful and alfa = alfmax.
*            = 3 if a better point was found but too many functions
*                were needed (not sufficient decrease).
*            = 4 if alfmax < tolabs (too small to do a search).
*            = 5 if alfa < alfsml (srchq only -- maybe want to switch
*                to central differences to get a better direction).
*            = 6 if the search found that there is no useful step.
*                The interval of uncertainty is less than 2*tolabs.
*                The minimizer is very close to alfa = zero
*                or the gradients are not sufficiently accurate.
*            = 7 if there were too many function calls.
*            = 8 if the input parameters were bad
*                (alfmax le toltny  or  oldg ge 0).
*     ------------------------------------------------------------------
*+    repeat
  100    if (needfd) then
            call srchq ( first , debug , done  , imprvd, inform,
     $                   maxf  , numf  , nout  , 
     $                   alfmax, alfsml, epsaf , 
     $                   g0    , targtg, ftry  ,         
     $                   tolabs, tolrel, toltny,
     $                   alfa  , alfbst, fbest  )
         else
            call srchc ( first , debug , done  , imprvd, inform,
     $                   maxf  , numf  , nout  ,
     $                   alfmax,         epsaf , 
     $                   g0    , targtg, ftry  , gtry  , 
     $                   tolabs, tolrel, toltny,
     $                   alfa  , alfbst, fbest , gbest )
         end if

         if (imprvd) then
            objf   = tobj
            objalf = tobjM

            if (ncnln .gt. 0)
     $         call dcopy ( ncnln, c2, 1, c, 1 )

            if (.not. needfd) then
               call dcopy ( n, gradu, 1, grad, 1 )
               gdx    = tgdx
               glf    = tglf

               if (ncnln .gt. 0) then
                  call dcopy ( ncnln, cjdx2, 1, cjdx, 1 )
                  call f06qff( 'General', ncnln, n, cjacu, ldcju,
     $                         cjac, ldcj )
               end if
            end if
         end if

*        ---------------------------------------------------------------
*        If done = false,  the problem functions must be computed for
*        the next entry to srchc or srchq.
*        If done = true,   this is the last time through.
*        ---------------------------------------------------------------
         if (.not. done) then

            call dcopy ( n,       x1, 1, x, 1 )
            call daxpy ( n, alfa, dx, 1, x, 1 )

            if (ncnln .gt. 0) then

*              Compute the new estimates of the multipliers and slacks.
*              If the step length is greater than one,  the multipliers
*              are fixed as the QP-multipliers.

               if (alfa .le. one) then
                  call dcopy ( ncnln,       cmul1, 1, cmul, 1 )
                  call daxpy ( ncnln, alfa, dlam , 1, cmul, 1 )
               end if
               call dcopy ( ncnln,       slk1, 1, slk, 1 )
               call daxpy ( ncnln, alfa, dslk, 1, slk, 1 )

*              ---------------------------------------------------------
*              Compute the new constraint vector and Jacobian.
*              ---------------------------------------------------------
               call confun( mode, ncnln, n, ldcju,
     $                      needc, x, c2, cjacu, nstate )
               if (mode .lt. 0) go to 999

               call dcopy ( ncnln,         c2 , 1, cs, 1 )
               call daxpy ( ncnln, (-one), slk, 1, cs, 1 )

               call dcopy ( ncnln, cs , 1,  work, 1 )
               call ddscl ( ncnln, rho, 1,  work, 1 )

               fterm  =            ddot( ncnln, cmul, 1, cs, 1 ) -
     $                  half*scale*ddot( ncnln, work, 1, cs, 1 )
            end if

*           ------------------------------------------------------------
*           Compute the value and gradient of the objective function.
*           ------------------------------------------------------------
            call objfun( mode, n, x, tobj, gradu, nstate )
            if (mode .lt. 0) go to 999

            if (ncnln .gt. 0) then
               tobjM = tobj  - fterm
            else
               tobjM = tobj
            end if

            ftry  = tobjM - oldf - rmu*oldg*alfa 

*           ------------------------------------------------------------
*           Compute auxiliary gradient information.
*           ------------------------------------------------------------
            if (.not. needfd) then
               gtry   = ddot( n, gradu, 1, dx, 1 )
               tgdx   = gtry
               tglf   = gtry
               if (ncnln .gt. 0) then

*                 Compute the Jacobian times the search direction.

                  call dgemv ( 'n', ncnln, n, one, cjacu, ldcju, dx, 1,
     $                         zero, cjdx2, 1 )

                  call dcopy ( ncnln,         cjdx2, 1, work, 1 )
                  call daxpy ( ncnln, (-one), dslk , 1, work, 1 )

                  gtry   = gtry - ddot( ncnln, cmul, 1, work, 1 )
                  if (alfa .le. one)
     $               gtry   = gtry - ddot( ncnln, dlam, 1, cs      , 1 )

                  call ddscl ( ncnln, rho , 1, work, 1 )
                  gtry = gtry  +
     $                     scale*ddot( ncnln, work, 1, cs   , 1 )

                  tglf = tgdx  - ddot( ncnln, cjdx2, 1, qpmul, 1 )

*                 ------------------------------------------------------
*                 If alfbnd .le. alfa .lt. alfmax and the norm of the
*                 quasi-Newton update is bounded, set alfmax to be alfa.
*                 This will cause the linesearch to stop if the merit
*                 function is decreasing at the boundary.
*                 ------------------------------------------------------
                  if (alfbnd .le. alfa  .and.  alfa .lt. alfmax) then

                     csjdx  = ddot   ( ncnln, cs, 1, cjdx2, 1 )

                     if (npdbg  .and.  inpdbg(1) .gt. 0)
     $                  write (nout, 1400) csjdx, cs1jdx, curvlf

                     curvlf = tglf  - glf1
                     curvc  = abs( csjdx - cs1jdx )
                     rhobfs = max( qpcurv*tolg - curvlf, zero )
                     if (rhobfs .le. curvc*rhomax) then
                        alfmax = alfa
                     else
                        alfbnd = min( two*alfa, alfmax )
                     end if
                     if (npdbg  .and.  inpdbg(1) .gt. 0)
     $                  write(nout,1300) alfbnd, alfa, alfmax
                  end if
               end if

               gtry = gtry  - rmu*oldg

            end if
         end if
*+    until (      done)
      if    (.not. done) go to 100

      nfun = nfun + numf
      if (.not. needfd) ngrad = ngrad + numf
      alfa = alfbst

      if (.not. imprvd) then
         call dcopy ( n,       x1, 1, x, 1 )
         call daxpy ( n, alfa, dx, 1, x, 1 )
         if (ncnln .gt. 0) then
            if (alfa .le. one) then
               call dcopy ( ncnln,       cmul1, 1, cmul, 1 )
               call daxpy ( ncnln, alfa, dlam , 1, cmul, 1 )
            end if
            call dcopy ( ncnln,         slk1 , 1, slk, 1 )
            call daxpy ( ncnln,   alfa, dslk , 1, slk, 1 )
            call dcopy ( ncnln,         c    , 1, cs , 1 )
            call daxpy ( ncnln, (-one), slk  , 1, cs , 1 )
         end if
      end if

      if (npdbg  .and.  inpdbg(1) .gt. 0)
     $   write (nout, 1200) inform

      return

*     The user wants to stop.  Who am I to object?

  999 inform = mode
      return

 1200 format(/ ' //npsrch// inform  = ', i4 )
 1300 format(/ ' //npsrch//        alfbnd          alfa        alfmax'
     $       / ' //npsrch//', 1p, 3e14.2 )
 1400 format(/ ' //npsrch//         csjdx        cs1jdx        curvlf'
     $       / ' //npsrch//', 1p, 3e14.2 )

*     end of  npsrch
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npupdt( lsumry, unitq,
     $                   n, ncnln, nfree, nz,
     $                   ldcj1, ldcj2, ldq, ldr, kx,
     $                   alfa, glf1, glf2, qpcurv,
     $                   cjac1, cjac2, cjdx1, cjdx2,
     $                   cs1, cs2, gq1, gq2, hpq, rpq,
     $                   qpmul, r, omega, zy, wrk1, wrk2 )

      implicit           double precision (a-h,o-z)
      character*25       lsumry
      logical            unitq
      integer            kx(n)
      double precision   cjac1(ldcj1,*), cjac2(ldcj2,*),
     $                   cjdx1(*), cjdx2(*), cs1(*), cs2(*),
     $                   gq1(n), gq2(n), hpq(n), rpq(n), qpmul(*),
     $                   r(ldr,*), omega(*), zy(ldq,*)
      double precision   wrk1(n+ncnln), wrk2(n)

************************************************************************
*     npupdt  computes the BFGS update for the approximate Hessian of
*     the Lagrangian.  If the approximate curvature of the Lagrangian
*     function is negative,  a nonnegative penalty vector OMEGA(i) of
*     minimum two norm is computed such that the approximate curvature
*     of the augmented Lagrangian will be positive. If no finite penalty
*     vector exists,  the BFGS update is performed with the approximate
*     curvature modified to be a small positive value.
*
*     On entry,  gq1 and gq2 contain the transformed objective gradients
*     at x1 and x2,  Hpq contains  R'R(pq), the transformed Hessian
*     times the transformed search direction.  The vectors gq1 and Hpq
*     are not saved.  If the regular BFGS quasi-Newton update could not
*     be performed, the first character of lsumry is loaded with 'M'.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 66 version written April 1984.
*     Level 2 BLAS added 12-June-1986.
*     Level 2 matrix routines added 22-Apr-1988.
*     This version of NPUPTD dated  22-Apr-1988.
************************************************************************
      common    /sol1cm/ nout
      common    /sol6cm/ rcndbd, rfrobn, drmax, drmin

      logical            incrun
      common    /sol6np/ rhomax, rhonrm, rhodmp, scale, incrun

      logical            npdbg
      parameter        ( ldbg   = 5 )
      common    /npdebg/ inpdbg(ldbg), npdbg

      logical            overfl, ssbfgs
      intrinsic          max   , min   , sqrt
      external           idamax, ddiv  , ddot  , dnrm2
      parameter        ( zero   = 0.0d+0, one    = 1.0d+0 )
      parameter        ( tolg   = 1.0d-1                  )

      if (ncnln .gt. 0) call dload ( ncnln, zero, omega, 1 )

*     ------------------------------------------------------------------
*     Set CURVL = (G2 - G1)'DX,  the approximate curvature along DX of
*     the (augmented) Lagrangian.  At first, the curvature is not scaled
*     by the steplength ALFA.
*     ------------------------------------------------------------------
      curvl  = glf2 -   glf1
      tinycl =        qpcurv * tolg
      ssbfgs = curvl .le. alfa*tinycl
      if (npdbg  .and.  inpdbg(1) .gt. 0)
     $   write (nout, 1000) ssbfgs, tinycl, curvl

*     ------------------------------------------------------------------
*     Test if CURVL is sufficiently positive.  If there are no nonlinear
*     constraints,  no update can be performed.
*     ------------------------------------------------------------------
      if (curvl  .lt. tinycl) then
         lsumry(1:13) = 'Modified BFGS'
         if (ncnln .gt. 0) then
            qmax = zero
            do 200 i = 1, ncnln
               qi  = cjdx2(i)*cs2(i) - cjdx1(i)*cs1(i)
               qmax = max( qmax, qi )
               if (qi .le. zero) wrk1(i) = zero
               if (qi .gt. zero) wrk1(i) = qi
  200       continue

            qnorm = dnrm2 ( ncnln, wrk1, 1 )

            test  = max( tinycl - curvl, zero )
            beta  = ddiv  ( qmax*test, qnorm*qnorm, overfl )
            if (beta .lt. rhomax  .and.  .not. overfl) then
               lsumry(1:1) = ' '
               beta  = test/(qnorm*qnorm)
               do 210 i = 1, ncnln
                  qi       = wrk1(i)
                  omega(i) =            beta*qi
                  curvl    = curvl    + beta*qi*qi
  210          continue

               if (npdbg) then
                  imax = idamax( ncnln, omega, 1 )
                  if (inpdbg(1) .gt. 0)
     $               write (nout, 1250) omega(imax)

                  if (inpdbg(2) .gt. 0)
     $               write (nout, 1300) (omega(i), i=1,ncnln)
               end if
            end if
         end if
      end if

*     ------------------------------------------------------------------
*     Compute the difference in the augmented Lagrangian gradient.
*     ------------------------------------------------------------------
*     Update gq1 to include the augmented Lagrangian terms.

      if (ncnln .gt. 0) then

         do 310 i = 1, ncnln
            wrk1(i) = - qpmul(i) + omega(i) * cs1(i)
  310    continue
         call dgemv ( 'T', ncnln, n, one, cjac1, ldcj1, wrk1, 1,
     $                zero, wrk2, 1 )

         do 320 i = 1, ncnln
            wrk1(i) =   qpmul(i) - omega(i) * cs2(i)
  320    continue
         call dgemv ( 'T', ncnln, n, one, cjac2, ldcj2, wrk1, 1,
     $                one, wrk2, 1 )

         call cmqmul( 6, n, nz, nfree, ldq, unitq, kx, wrk2, zy, wrk1 )
         call daxpy ( n, one, wrk2, 1, gq1, 1 )
      end if

      if (npdbg  .and.  inpdbg(1) .gt. 0)
     $   write (nout, 1100) alfa  , curvl

      if (curvl .lt. tinycl) curvl  = tinycl

      do 330 j = 1, n
         wrk2(j) = gq2(j) - gq1(j)
  330 continue

      rtgtp  = sqrt(qpcurv)
      rtyts  = sqrt(alfa*curvl)
      eta    = one
      if (ssbfgs)
     $   eta = rtyts / (rtgtp*alfa)

      trace1 = dnrm2 ( n,  hpq, 1 ) /  rtgtp
      trace2 = dnrm2 ( n, wrk2, 1 ) / (rtyts*eta)
      rfrobn = eta*sqrt( abs(  (rfrobn - trace1)*(rfrobn + trace1)
     $                                 + trace2**2) )

*     ==================================================================
*     Update the Cholesky factor of  Q'HQ.
*     ==================================================================
*     Normalize the vector  RPQ ( = R(pq) ).

      call dscal ( n, (one / rtgtp), rpq, 1 )

*     Do the self-scaled or regular BFGS update.
*     Form the vector WRK1 = gamma * (GQ2 - GQ1) - beta * R'R*PQ,
*     where  gamma = 1/SQRT( CURV ) = 1/SQRT( (GQ2 - GQ1)'SQ )

      call dscal ( n, (one / rtgtp), hpq, 1 )

      if (ssbfgs) then
         do 410 j   = 1, n
            call dscal ( j, eta, r(1,j), 1 )
            wrk1(j) = wrk2(j)/rtyts  -  eta * hpq(j)
  410    continue
      else
         do 420 j   = 1, n
            wrk1(j) = wrk2(j)/rtyts  -        hpq(j)
  420    continue
      end if

*     Perform the update to  R = R + RPQ*WRK1'.
*     RPQ is overwritten. Arrays  GQ1  and  HPQ  are used to store the
*     sines and cosines defined by the plane rotations.

      call cmr1md( n, 0, n, ldr, n, n, r, hpq, rpq, wrk1, gq1, hpq )

      return

 1000 format(/ ' //npupdt// ssbfgs    min. curvl         curvl '
     $       / ' //npupdt//   ', l4, 1p, 2e14.2 )
 1100 format(/ ' //npupdt//          alfa         curvl '
     $       / ' //npupdt//', 1p, 2e14.2 )
 1250 format(/ ' //npupdt//   omega(imax)'
     $       / ' //npupdt//', 1pe14.2 )
 1300 format(/ ' //npupdt//  Penalty parameters = '  / (1p, 5e15.6))

*     end of  npupdt
      end
