*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        
      subroutine chfgrd( inform, msglvl, n,
     $                   bigbnd, epsrf, oktol, fdchk, objf, xnorm,
     $                   objfun,
     $                   bl, bu, grad, gradu, dx, x, y, w, lenw )

      implicit           double precision(a-h,o-z)

      double precision   bl(n), bu(n), grad(n), gradu(n), dx(n)
      double precision   x(n), y(n), w(lenw)
      external           objfun

C***********************************************************************
C     CHFGRD  checks if the gradients of the objective function have
C     been coded correctly.
C
C     On input,  the value of the objective function at the point X is
C     stored in OBJF.  The corresponding gradient is stored in GRADU.
C     If any gradient component has not been specified,  it will have a
C     dummy value.  Missing values are not checked.
C
C     A cheap test is first undertaken by calculating the directional
C     derivative using two different methods.  If this proves 
C     satisfactory and no further information is desired, CHFGRD is 
C     terminated.  Otherwise, the routine CHCORE is called to give 
C     optimal step-sizes and a forward-difference approximation to each 
C     component of the gradient for which a test is deemed necessary,
C     either by the program or the user.
C
C     Other inputs:
C
C           X         The n-dimensional point at which the
C                     gradient is to be verified.
C           EPSRF     The positive bound on the relative error
C                     associated with computing the function at
C                     the point x.
C           OKTOL     The desired relative accuracy which the
C                     components of the gradient should satisfy.
C
C     LVRFYC has the following meaning...
C
C       -1        do not perform any check.
C        0        do the cheap test only.
C        1 or 3   do both cheap and full test.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written  19-May-1985.
C     This version of CHFGRD dated  19-Nov-91.  
C***********************************************************************
      common    /sol1cm/ nout
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      common    /sol5np/ lvrfyc, jverfy(4)

      logical            npdbg        
      parameter         (ldbg = 5)
      common    /npdebg/ inpdbg(ldbg), npdbg

      logical            const , debug , done  , first , headng
      logical            needed, ok
      character*4        key   , lbad  , lgood
      character*18       result(0:4)
      intrinsic          abs   , max   , min   , sqrt
      external           ddot
      parameter         (rdummy =-11111.0d+0              )
      parameter         (zero   =0.0d+0, half  = 0.5d+0, point9 =0.9d+0)
      parameter         (one    =1.0d+0, two   = 2.0d+0, ten    =1.0d+1)
      parameter         (lbad   ='BAD?', lgood = '  OK')
      data               result
     $                 / '                 ', 'Constant?      ',
     $                   'Linear or odd?   ', 'Too nonlinear?',
     $                   'Small derivative?'                   /

      inform = 0
      needed = lvrfyc .eq. 0  .or.  lvrfyc .eq. 1  .or.  lvrfyc .eq. 3
      if (.not. needed) return
 
      write (nout, 1000)
      debug  = npdbg  .and.  inpdbg(5) .gt. 0
      nstate = 0

      biglow = - bigbnd
      bigupp =   bigbnd

*     ==================================================================
*     Perform the cheap test.
*     ==================================================================
      h =     (one + xnorm)*fdchk

      if (     n .le. 100) then
         dxmult = 0.9
      else if (n .le. 250) then
         dxmult = 0.99
      else 
         dxmult = 0.999
      end if

      dxj  = one / n
      do 110, j = 1, n
         dx(j) =   dxj
         dxj   = - dxj*dxmult
  110 continue

*     ------------------------------------------------------------------
*     Do not perturb X(J) if the  J-th  element is missing.
*     Compute the directional derivative.
*     ------------------------------------------------------------------
      ncheck = 0
      do 120 j = 1, n
         if (grad(j) .eq. rdummy) then
            dx(j) = zero
         else
            ncheck = ncheck + 1

            xj     =   x(j)                             
            stepbl = - one       
            stepbu =   one
            if (bl(j) .gt. biglow)
     $         stepbl = max( stepbl, bl(j) - xj )
            if (bu(j) .lt. bigupp  .and.  bu(j) .gt. bl(j))
     $         stepbu = min( stepbu, bu(j) - xj )

            if (half*(stepbl + stepbu) .lt. zero) then
               dx(j) = dx(j)*stepbl
            else
               dx(j) = dx(j)*stepbu
            end if
         end if
  120 continue

      if (ncheck .eq. 0) then
         write (nout, 3500)
         return
      end if
      gdx    = ddot  ( n, gradu, 1, dx, 1 )

*     ------------------------------------------------------------------
*     Make forward-difference approximation along  p.
*     ------------------------------------------------------------------
      call dcopy ( n,     x, 1, y, 1 )
      call daxpy ( n, h, dx, 1, y, 1 )
      
      mode   = 0
      call objfun( mode, n, y, objf1, gradu, nstate )
      if (mode .lt. 0) go to 9999

      gdiff =    (objf1 - objf) / h
      error = abs(gdiff - gdx ) / (abs(gdx) + one)

      ok    = error .le. oktol
      if (ok) then
         if (msglvl .gt. 0) write (nout, 1100)
      else
         write (nout, 1200)
         if (error .ge. point9) inform = 1
      end if

      if (msglvl .gt. 0) write (nout, 1300) gdx, gdiff

*     ==================================================================
*     Component-wise check.
*     ==================================================================
      if (lvrfyc .eq. 1  .or.  lvrfyc .eq. 3) then
         headng = .true.
         itmax  = 3
         nwrong = 0
         ngood  = 0
         jmax   = 0
         emax   = zero
         ncheck = 0
         j1     = jverfy(1)
         j2     = jverfy(2)

*        ---------------------------------------------------------------
*        Loop over each of the components of  x.
*        ---------------------------------------------------------------
         do 500 j = j1, j2

            if (grad(j) .ne. rdummy) then
*              ---------------------------------------------------------
*              Check this gradient component.
*              ---------------------------------------------------------
               ncheck = ncheck + 1
               gj     = grad(j)
               gsize  = abs( gj )
               xj     = x(j)
*              ---------------------------------------------------------
*              Find a finite-difference interval by iteration.
*              ---------------------------------------------------------
               iter   = 0
               epsa   = epsrf*(one + abs( objf ))
               cdest  = zero
               sdest  = zero
               first  = .true.

               stepbl = biglow
               stepbu = bigupp
               if (bl(j) .gt. biglow) stepbl = bl(j) - xj
               if (bu(j) .lt. bigupp) stepbu = bu(j) - xj

               hopt   = two*(one + abs( xj ))*sqrt( epsrf )
               h      = ten*hopt
               if (half*(stepbl + stepbu) .lt. zero) h =  - h

*+             repeat
  400             x(j)  = xj + h
                  call objfun( mode, n, x, f1, gradu, nstate )
                  if (mode .lt. 0) go to 9999

                  x(j)  = xj + h + h
                  call objfun( mode, n, x, f2, gradu, nstate )
                  if (mode .lt. 0) go to 9999
                              
                  call chcore( debug, done, first, epsa, epsrf, objf,xj,
     $                         info, iter, itmax,
     $                         cdest, fdest, sdest, errbnd, f1,
     $                         f2, h, hopt, hphi )

*+             until     done
               if (.not. done) go to 400

*              ---------------------------------------------------------
*              Exit for this variable.
*              ---------------------------------------------------------
               gdiff = cdest
               x(j)  = xj

               error = abs(gdiff - gj) / (gsize + one)
               if (error .ge. emax) then
                  emax  = error
                  jmax  = j
               end if

               ok =  error .le. oktol
               if (ok) then
                  key    = lgood
                  ngood  = ngood  + 1
               else
                  key    = lbad
                  nwrong = nwrong + 1
               end if

*              Zero components are not printed.

               const = ok .and. info .eq. 1 .and. abs(gj) .lt. epspt8
               if (.not. const) then
                  if (headng) write (nout, 3000)
                  if (ok) then
                     write (nout, 3100) j, xj, hopt, gj, gdiff,
     $                                  key, iter
                  else
                     write (nout, 3110) j, xj, hopt, gj, gdiff,
     $                                  key, iter, result(info)
                  end if
                  headng = .false.
               end if
            end if
  500    continue

*        ===============================================================
*        Done.
*        ===============================================================
         inform = 0
         if (nwrong .eq. 0) then
            write (nout, 3200) ngood , ncheck, j1    , j2
         else
            write (nout, 3300) nwrong, ncheck, j1    , j2
            if (error .ge. point9) inform = 1
         end if
         write (nout, 3400) emax, jmax
      end if

      call dcopy ( n, grad, 1, gradu, 1 )

      return

 9999 inform = mode
      return

 1000 format(/// ' Verification of the objective gradients.'
     $       /   ' ----------------------------------------' )
 1100 format(/   ' The objective gradients seem to be ok.')
 1200 format(/   ' XXX  The objective gradients seem to be incorrect.')
 1300 format(/   ' Directional derivative of the objective', 1p, e18.8/
     $           ' Difference approximation               ', 1p, e18.8 )
 3000 format(// 4x, 'j', 4x, 'x(j)', 5x, 'dx(j)', 11x,
     $           'g(j)', 9x, '  Difference approxn  Itns' /)
 3100 format(  i5, 1p, 2e10.2,      2e18.8, 2x, a4, i6          )
 3110 format(  i5, 1p, 2e10.2,      2e18.8, 2x, a4, i6, 2x, a18 )
 3200 format(/ i7, '  Objective gradients out of the', i6,
     $             '  set in cols', i6, '  through', i6,
     $             '  seem to be ok.')
 3300 format(/   ' XXX  There seem to be', i6,
     $           '  incorrect objective gradients out of the', i6,
     $           '  set in cols', i6, '  through', i6 )
 3400 format(/   ' The largest relative error was', 1p, e12.2,
     $           '   in element', i6 /)
 3500 format(/   ' No gradient elements assigned.' )

*     end of  chfgrd.
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine chcjac( inform, lvlder, msglvl,
     $                   ncset, n, ncnln, ldcj, ldcju,
     $                   bigbnd, epsrf, oktol, fdchk, xnorm,
     $                   confun, needc,
     $                   bl, bu, c, c1, cjac, cjacu, cjdx,
     $                   dx, err, x, y, w, lenw )

      implicit           double precision(a-h,o-z)
      integer            needc(*)
      double precision   bl(n), bu(n), c(*), c1(*), cjdx(*),
     $                   cjac(ldcj,*), cjacu(ldcju,*), err(*)
      double precision   dx(n), x(n), y(n), w(lenw)
      external           confun

C***********************************************************************
C     CHCJAC  checks if the gradients of the constraints have been coded
C     correctly.
C
C     On input,  the values of the constraints at the point X are stored
C     in C.  Their corresponding gradients are stored in CJACU.  If any
C     Jacobian component has not been specified,  it will have a dummy
C     value.  Missing values are not checked.
C
C     A cheap test is first undertaken by calculating the directional
C     derivative using two different methods.  If this proves 
C     satisfactory and no further information is desired, CHCJAC is 
C     terminated.  Otherwise, CHCORE is called to give optimal stepsizes
C     and a central-difference approximation to each component of the 
C     Jacobian for which a test is deemed necessary, either by the 
C     program or the user.
C
C     LVRFYC has the following meaning...
C
C       -1        do not perform any check.
C        0        do the cheap test only.
C        2 or 3   do both cheap and full test.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written  19-May-1985.
C     This version of CHCJAC dated 19-Nov-91.  
C***********************************************************************
      common    /sol1cm/ nout
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      common    /sol5np/ lvrfyc, jverfy(4)

      logical            npdbg
      parameter         (ldbg = 5)
      common    /npdebg/ inpdbg(ldbg), npdbg

      logical            const , debug , done  , first , headng
      logical            needed, ok
      character*4        key   , lbad  , lgood
      character*18       result(0:4)
      intrinsic          abs   , max   , min   , sqrt
      external           ddot  , idamax
      parameter         (rdummy =-11111.0d+0              )
      parameter         (zero   =0.0d+0, half   =0.5d+0, point9 =0.9d+0)
      parameter         (one    =1.0d+0, two    =2.0d+0, ten    =1.0d+1)
      parameter         (lbad   ='BAD?', lgood  ='  OK')
      data               result
     $                 / '                 ', 'Constant?      ',
     $                   'Linear or odd?   ', 'Too nonlinear?',
     $                   'small derivative?'                   /

      inform = 0
      needed = ncnln  .gt. 0  .and.
     $         lvrfyc .eq. 0  .or.   lvrfyc .eq. 2  .or.  lvrfyc .eq. 3
      if (.not. needed) return

      write (nout, 1000)
      debug  = npdbg  .and.  inpdbg(5) .gt. 0
      nstate = 0

      biglow = - bigbnd
      bigupp =   bigbnd

*     ==================================================================
*     Perform the cheap test.
*     ==================================================================
      h = (one + xnorm)*fdchk

      if (     n .le. 100) then
         dxmult = 0.9
      else if (n .le. 250) then
         dxmult = 0.99
      else 
         dxmult = 0.999
      end if

      dxj  = one / n
      do 110, j = 1, n
         dx(j) =   dxj
         dxj   = - dxj*dxmult
  110 continue

*     ------------------------------------------------------------------
*     Do not perturb  X(J)  if the  J-th  column contains any
*     unknown elements.  Compute the directional derivative for each
*     constraint gradient.
*     ------------------------------------------------------------------
      ncheck = 0
      do 140 j = 1, n
         do 130 i = 1, ncnln
            if (cjac(i,j) .eq. rdummy) then
               dx(j) = zero      
               go to 140
            end if           
  130    continue
         ncheck = ncheck + 1

         xj     =   x(j)
         stepbl = - one
         stepbu =   one
         if (bl(j) .gt. biglow)
     $      stepbl = max( stepbl, bl(j) - xj )
         if (bu(j) .lt. bigupp  .and.  bu(j) .gt. bl(j))
     $      stepbu = min( stepbu, bu(j) - xj )

         if (half*(stepbl + stepbu) .lt. zero) then
            dx(j) = dx(j)*stepbl
         else
            dx(j) = dx(j)*stepbu
         end if
  140 continue

      if (ncheck .eq. 0) then
         write (nout, 2300)
      else

*        compute  (jacobian)*dx.

         call dgemv ( 'Normal', ncnln, n, one, cjacu, ldcju,
     $                dx, 1, zero, cjdx, 1 )

*        ---------------------------------------------------------------
*        Make forward-difference approximation along DX.
*        ---------------------------------------------------------------
         call dcopy ( n,     x, 1, y, 1 )
         call daxpy ( n, h, dx, 1, y, 1 )

         call iload ( ncnln, (1), needc, 1 )

         mode   = 0
         call confun( mode, ncnln, n, ldcju,
     $                needc, y, c1, cjacu, nstate )
         if (mode .lt. 0) go to 9999

*        Set  err = (c1 - c)/h  - Jacobian*dx.  This should be small.

         do 170 i = 1, ncnln
            err(i) = (c1(i) - c(i)) / h  -  cjdx(i)
  170    continue            
         imax  = idamax( ncnln, err, 1 )
         emax  = abs(err(imax)) / (abs(cjdx(imax)) + one)

         if (emax .le. oktol) then
            if (msglvl .gt. 0) write (nout, 2000)
         else
            write (nout, 2100)
            if (emax .ge. point9) inform = 1
         end if
         if (msglvl .ge. 0) write (nout, 2200) emax, imax
      end if

*     ==================================================================
*     Component-wise check.
*     ==================================================================
      if (lvrfyc .ge. 2) then
         if (lvlder .eq. 3) then

*           Recompute the Jacobian to find the non-constant elements.

            call f06qhf( 'General', ncnln, n, rdummy, rdummy, 
     $                   cjacu, ldcju )

            call iload ( ncnln, (1), needc, 1 )
            nstate = 0
            mode   = 2

            call confun( mode, ncnln, n, ldcju,
     $                   needc, x, c1, cjacu, nstate )
            if (mode .lt. 0) go to 9999
         end if

         call iload ( ncnln, (0), needc, 1 )

         itmax  =   3
         ncheck =   0
         nwrong =   0
         ngood  =   0
         colmax = - one
         jcol   =   0
         irow   =   0
         mode   =   0
         j3     =   jverfy(3)
         j4     =   jverfy(4)

*        ---------------------------------------------------------------
*        Loop over each column.
*        ---------------------------------------------------------------
         do 600 j = j3, j4

            call dload ( ncnln, zero, err, 1 )
            headng = .true.
            xj     = x(j)

            stepbl = biglow
            stepbu = bigupp
            if (bl(j) .gt. biglow) stepbl = bl(j) - xj
            if (bu(j) .lt. bigupp) stepbu = bu(j) - xj

            signh  = one
            if (half*(stepbl + stepbu) .lt. zero) signh =  - one

            do 500, i = 1, ncnln
               epsaci   = epsrf*(one + abs( c(i) ))

               if (cjacu(i,j) .ne. rdummy) then
*                 ------------------------------------------------------
*                 Check this Jacobian element.
*                 ------------------------------------------------------
                  ncheck   = ncheck + 1
                  needc(i) = 1

                  cij    = cjac(i,j)
                  cjsize = abs( cij )
*                 ------------------------------------------------------
*                 Find a finite-difference interval by iteration.
*                 ------------------------------------------------------
                  iter   = 0
                  hopt   = two*(one + abs( xj ))*sqrt( epsrf )
                  h      = ten*hopt*signh
                  cdest  = zero
                  sdest  = zero
                  first  = .true.

*+                repeat
  400                x(j)  = xj + h
                     call confun( mode, ncnln, n, ldcju,
     $                            needc, x, c1, cjacu, nstate )
                     if (mode .lt. 0) go to 9999
                     f1    = c1(i)

                     x(j)  = xj + h + h
                     call confun( mode, ncnln, n, ldcju,
     $                            needc, x, c1, cjacu, nstate )
                     if (mode .lt. 0) go to 9999
                     f2    = c1(i)
                                                     
                     call chcore( debug,done,first,epsaci,epsrf,c(i),xj,
     $                            info, iter, itmax,
     $                            cdest, fdest, sdest, errbnd, f1,
     $                            f2, h, hopt, hphi )

*+                until     done
                  if (.not. done) go to 400

*                 ------------------------------------------------------
*                 Exit for this element.
*                 ------------------------------------------------------
                  cjdiff   = cdest
                  err(i)   = abs(cjdiff - cij) / (cjsize + one)

                  ok       = err(i) .le. oktol
                  if (ok) then
                     key    = lgood
                     ngood  = ngood  + 1
                  else
                     key    = lbad
                     nwrong = nwrong + 1
                  end if

                  const = ok .and. info       .eq. 1
     $                       .and. abs( cij ) .lt. epspt8
                  if (.not. const) then
                     if (headng) then
                        write (nout, 4000)
                        if (ok)
     $                     write (nout, 4100)   j, xj    , hopt, i,
     $                                        cij, cjdiff, key , iter
                        if (.not. ok)
     $                     write (nout, 4110)   j, xj    , hopt, i,
     $                                        cij, cjdiff, key , iter,
     $                                        result(info)
                        headng = .false.
                     else
                        if (ok)
     $                     write (nout, 4200)              hopt, i,
     $                                        cij, cjdiff, key , iter
                        if (.not. ok)
     $                     write (nout, 4210)              hopt, i,
     $                                        cij, cjdiff, key , iter,
     $                                        result(info)
                     end if
                  end if
                  needc(i) = 0
               end if
  500       continue

*           ------------------------------------------------------------
*           Finished with this column.
*           ------------------------------------------------------------
            if (.not. headng) then
               imax = idamax( ncnln, err, 1 )
               emax = abs( err(imax) )

               if (emax .ge. colmax) then
                  irow   = imax
                  jcol   = j
                  colmax = emax
               end if
            end if
            x(j) = xj
  600    continue

         inform = 0
         if (ncheck .eq. 0) then
            write (nout, 4600) ncset
         else
            if (nwrong .eq. 0) then
               write (nout, 4300) ngood , ncheck, j3, j4
            else
               write (nout, 4400) nwrong, ncheck, j3, j4
               if (colmax .ge. point9) inform = 1
            end if
            write (nout, 4500) colmax, irow, jcol
         end if
      end if

*     Copy  ( constants + gradients + dummy values )  back into CJACU.
      
      call f06qff( 'General', ncnln, n, cjac, ldcj, cjacu, ldcju )

      return            

 9999 inform = mode
      return

 1000 format(/// ' Verification of the constraint gradients.'
     $       /   ' -----------------------------------------' )
 2000 format(/   ' The constraint Jacobian seems to be ok.')
 2100 format(/   ' XXX  The constraint Jacobian seems to be incorrect.')
 2200 format(/   ' The largest relative error was', 1p, e12.2,
     $           '  in row', i5 /)
 2300 format(/   ' Every column contains a constant or',
     $           ' missing element.')
 4000 format(// ' Column    x(j)     dx(j)    Row   ',
     $          ' Jacobian Value      Difference Approxn  Itns'    )
 4100 format(/ i7,      1p, 2e10.2, i5, 2e18.8, 2x, a4, i6         )
 4110 format(/ i7,      1p, 2e10.2, i5, 2e18.8, 2x, a4, i6, 2x, a18)
 4200 format(  7x, 10x, 1p,  e10.2, i5, 2e18.8, 2x, a4, i6         )
 4210 format(  7x, 10x, 1p,  e10.2, i5, 2e18.8, 2x, a4, i6, 2x, a18)
 4300 format(/ i7, '  constraint Jacobian elements out of the', i6,
     $             '  set in cols', i6, '  through', i6,
     $             '  seem to be ok.')
 4400 format(/   ' XXX  There seem to be', i6,
     $           '  incorrect Jacobian elements out of the', i6,
     $           '  set in cols', i6, '  through', i6 )
 4500 format(/ ' The largest relative error was', 1p, e12.2,
     $         '  in row', i5, ',  column', i5 /)
 4600 format(  ' All', i6, '   assigned Jacobian elements are',
     $         ' constant.' )

*     end of  chcjac
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine chfjac( inform, lvlder, msglvl,
     $                   nfset, m, n, ldfj, ldfju,
     $                   bigbnd, epsrf, oktol, fdchk, xnorm,
     $                   objfun,
     $                   bl, bu, f, f1, fjac, fjacu, fjdx,
     $                   dx, err, x, y, w, lenw )

      implicit           double precision(a-h,o-z)
      double precision   bl(n), bu(n), f(*), f1(*), fjdx(*),
     $                   fjac(ldfj,*), fjacu(ldfju,*), err(*)
      double precision   dx(n), x(n), y(n), w(lenw)
      external           objfun

C***********************************************************************
C     CHFJAC  checks if the objective Jacobian matrix has been coded
C     correctly.
C
C     On input,  the values of the objective vector at the point X are
C     stored in F.  Their corresponding gradients are stored in FJACU.
C     If any Jacobian component has not been specified,  it will have a
C     dummy value.  Missing values are not checked.
C
C     A cheap test is first undertaken by calculating the directional
C     derivative using two different methods.  If this proves
C     satisfactory and no further information is desired, CHFJAC is
C     terminated.  Otherwise, CHCORE is called to give optimal stepsizes
C     and a central-difference approximation to each component of the
C     Jacobian for which a test is deemed necessary, either by the
C     program or the user.
C
C     LVRFYC has the following meaning...
C
C       -1        do not perform any check.
C        0        do the cheap test only.
C        2 or 3   do both cheap and full test.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written  10-May-1988.
C     This version of CHFJAC dated  19-Nov-91.
C***********************************************************************
      common    /sol1cm/ nout
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      common    /sol5np/ lvrfyc, jverfy(4)

      logical            npdbg
      parameter         (ldbg = 5)
      common    /npdebg/ inpdbg(ldbg), npdbg

      logical            const , debug , done  , first , headng
      logical            needed, ok
      character*4        key   , lbad  , lgood
      character*18       result(0:4)
      intrinsic          abs   , max   , min   , sqrt
      external           ddot  , idamax
      parameter         (rdummy =-11111.0d+0              )
      parameter         (zero   =0.0d+0, half   =0.5d+0, point9 =0.9d+0)
      parameter         (one    =1.0d+0, two    =2.0d+0, ten    =1.0d+1)
      parameter         (lbad   ='BAD?', lgood  ='  OK')
      data               result
     $                 / '                 ', 'Constant?      ',
     $                   'Linear or odd?   ', 'Too nonlinear?',
     $                   'Small derivative?'                   /

      inform = 0
      needed = lvrfyc .eq. 0  .or.  lvrfyc .eq. 1  .or.  lvrfyc .eq. 3
      if (.not. needed) return

      write (nout, 1000)
      debug  = npdbg  .and.  inpdbg(5) .gt. 0
      nstate = 0

      biglow = - bigbnd
      bigupp =   bigbnd

*     ==================================================================
*     Perform the cheap test.
*     ==================================================================
      h = (one + xnorm)*fdchk

      if (     n .le. 100) then
         dxmult = 0.9
      else if (n .le. 250) then
         dxmult = 0.99
      else 
         dxmult = 0.999
      end if

      dxj  = one / n
      do 110, j = 1, n
         dx(j) =   dxj
         dxj   = - dxj*dxmult
  110 continue

*     ------------------------------------------------------------------
*     Do not perturb  X(J)  if the  J-th  column contains any
*     unknown elements.  Compute the directional derivative for each
*     objective gradient.
*     ------------------------------------------------------------------
      ncheck = 0
      do 140, j = 1, n
         do 130, i = 1, m
            if (fjac(i,j) .eq. rdummy) then
               dx(j) = zero
               go to 140
            end if
  130    continue
         ncheck = ncheck + 1

         xj     =   x(j)
         stepbl = - one
         stepbu =   one
         if (bl(j) .gt. biglow)
     $      stepbl = max( stepbl, bl(j) - xj )
         if (bu(j) .lt. bigupp  .and.  bu(j) .gt. bl(j))
     $      stepbu = min( stepbu, bu(j) - xj )

         if (half*(stepbl + stepbu) .lt. zero) then
            dx(j) = dx(j)*stepbl
         else
            dx(j) = dx(j)*stepbu
         end if
  140 continue

      if (ncheck .eq. 0) then
         write (nout, 2300)
      else

*        Compute  (Jacobian)*dx.

         call dgemv ( 'Normal', m, n, one, fjacu, ldfju,
     $                 dx, 1, zero, fjdx, 1 )

*        ---------------------------------------------------------------
*        Make forward-difference approximation along DX.
*        ---------------------------------------------------------------
         call dcopy ( n,     x, 1, y, 1 )
         call daxpy ( n, h, dx, 1, y, 1 )

         mode   = 0
         call objfun( mode, m, n, ldfju,
     $                y, f1, fjacu, nstate )
         if (mode .lt. 0) go to 9999

*        Set  err = (f1 - f)/h  - Jacobian*dx.  This should be small.

         do 170, i = 1, m
            err(i) = (f1(i) - f(i)) / h  -  fjdx(i)
  170    continue
         imax  = idamax( m, err, 1 )
         emax  = abs(err(imax)) / (abs(fjdx(imax)) + one)

         if (emax .le. oktol) then
            if (msglvl .gt. 0) write (nout, 2000)
         else
            write (nout, 2100)
            if (emax .ge. point9) inform = 1
         end if
         if (msglvl .ge. 0) write (nout, 2200) emax, imax
      end if

*     ==================================================================
*     Component-wise check.
*     ==================================================================
      if (lvrfyc .ge. 2) then
         if (lvlder .eq. 3) then

*           Recompute the Jacobian to find the non-constant elements.

            call f06qhf( 'General', m, n, rdummy, rdummy,
     $                   fjacu, ldfju )

            nstate = 0
            mode   = 2

            call objfun( mode, m, n, ldfju,
     $                   x, f1, fjacu, nstate )
            if (mode .lt. 0) go to 9999

         end if

         itmax  =   3
         ncheck =   0
         nwrong =   0
         ngood  =   0
         colmax = - one
         jcol   =   0
         irow   =   0
         mode   =   0
         j3     =   jverfy(3)
         j4     =   jverfy(4)

*        ---------------------------------------------------------------
*        Loop over each column.
*        ---------------------------------------------------------------
         do 600, j = j3, j4

            call dload ( m, zero, err, 1 )
            headng = .true.
            xj     = x(j)

            stepbl = biglow
            stepbu = bigupp
            if (bl(j) .gt. biglow) stepbl = bl(j) - xj
            if (bu(j) .lt. bigupp) stepbu = bu(j) - xj   
                                                 
            signh  = one
            if (half*(stepbl + stepbu) .lt. zero) signh =  - one

            do 500, i = 1, m
               epsafi   = epsrf*(one + abs( f(i) ))

               if (fjacu(i,j) .ne. rdummy) then
*                 ------------------------------------------------------
*                 Check this Jacobian element.
*                 ------------------------------------------------------
                  ncheck   = ncheck + 1

                  fij    = fjac(i,j)
                  fjsize = abs( fij )
*                 ------------------------------------------------------
*                 Find a finite-difference interval by iteration.
*                 ------------------------------------------------------
                  iter   = 0
                  hopt   = two*(one + abs( xj ))*sqrt( epsrf )
                  h      = ten*hopt*signh
                  cdest  = zero
                  sdest  = zero
                  first  = .true.

*+                repeat
  400                x(j)  = xj + h
                     call objfun( mode, m, n, ldfju,
     $                            x, f1, fjacu, nstate )
                     if (mode .lt. 0) go to 9999
                     fforw = f1(i)

                     x(j)  = xj + h + h
                     call objfun( mode, m, n, ldfju,
     $                            x, f1, fjacu, nstate )
                     if (mode .lt. 0) go to 9999
                     fback = f1(i)

                     call chcore( debug,done,first,epsafi,epsrf,f(i),xj,
     $                            info, iter, itmax,
     $                            cdest, fdest, sdest, errbnd, fforw,
     $                            fback, h, hopt, hphi )

*+                until     done
                  if (.not. done) go to 400

*                 ------------------------------------------------------
*                 Exit for this element.
*                 ------------------------------------------------------
                  fjdiff   = cdest
                  err(i)   = abs(fjdiff - fij) / (fjsize + one)

                  ok       = err(i) .le. oktol
                  if (ok) then
                     key    = lgood
                     ngood  = ngood  + 1
                  else
                     key    = lbad
                     nwrong = nwrong + 1
                  end if

                  const = ok .and. info       .eq. 1
     $                       .and. abs( fij ) .lt. epspt8
                  if (.not. const) then
                     if (headng) then
                        write (nout, 4000)
                        if (ok)
     $                     write (nout, 4100)   j, xj    , hopt, i,
     $                                        fij, fjdiff, key , iter
                        if (.not. ok)
     $                     write (nout, 4110)   j, xj    , hopt, i,
     $                                        fij, fjdiff, key , iter,
     $                                        result(info)
                        headng = .false.
                     else
                        if (ok)
     $                     write (nout, 4200)              hopt, i,
     $                                        fij, fjdiff, key , iter
                        if (.not. ok)
     $                     write (nout, 4210)              hopt, i,
     $                                        fij, fjdiff, key , iter,
     $                                        result(info)
                     end if
                  end if
               end if
  500       continue

*           ------------------------------------------------------------
*           Finished with this column.
*           ------------------------------------------------------------
            if (.not. headng) then
               imax = idamax( m, err, 1 )
               emax = abs( err(imax) )

               if (emax .ge. colmax) then
                  irow   = imax
                  jcol   = j
                  colmax = emax
               end if
            end if
            x(j) = xj
  600    continue

         inform = 0
         if (ncheck .eq. 0) then
            write (nout, 4600) nfset
         else
            if (nwrong .eq. 0) then
               write (nout, 4300) ngood , ncheck, j3, j4
            else
               write (nout, 4400) nwrong, ncheck, j3, j4
               if (colmax .ge. point9) inform = 1
            end if
            write (nout, 4500) colmax, irow, jcol
         end if
      end if

*     Copy  ( constants + gradients + dummy values )  back into fjacu.

      call f06qff( 'General', m, n, fjac, ldfj, fjacu, ldfju )

      return

 9999 inform = mode
      return

 1000 format(/// ' Verification of the objective gradients.'
     $       /   ' ----------------------------------------' )
 2000 format(/   ' The objective Jacobian seems to be ok.')
 2100 format(/   ' XXX  The objective Jacobian seems to be incorrect.')
 2200 format(/   ' The largest relative error was', 1p, e12.2,
     $           '  in row', i5 /)
 2300 format(/   ' Every column contains a constant or',
     $           ' missing element.')
 4000 format(// ' Column    x(j)     dx(j)    Row   ',
     $          ' Jacobian Value      Difference Approxn  Itns'    )
 4100 format(/ i7,      1p, 2e10.2, i5, 2e18.8, 2x, a4, i6         )
 4110 format(/ i7,      1p, 2e10.2, i5, 2e18.8, 2x, a4, i6, 2x, a18)
 4200 format(  7x, 10x, 1p,  e10.2, i5, 2e18.8, 2x, a4, i6         )
 4210 format(  7x, 10x, 1p,  e10.2, i5, 2e18.8, 2x, a4, i6, 2x, a18)
 4300 format(/ i7, '  Objective Jacobian elements out of the', i6,
     $             '  set in cols', i6, '  through', i6,
     $             '  seem to be ok.')
 4400 format(/   ' XXX  There seem to be', i6,
     $           '  incorrect objective Jacobian elements out of the',
     $             i6, '  set in cols', i6, '  through', i6 )
 4500 format(/ ' The largest relative error was', 1p, e12.2,
     $         '  in row', i5, ',  column', i5 /)
 4600 format(  ' All', i6, '   assigned Jacobian elements are',
     $         ' constant.' )

*     end of  chfjac
      end
