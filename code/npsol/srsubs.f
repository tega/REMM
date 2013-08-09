*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     srsubs.f
*
*     srchc    srchq
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine srchc ( first , debug , done  , imprvd, inform,
     $                   maxf  , numf  , nout  ,
     $                   alfmax,         epsaf , 
     $                   g0    , targtg, ftry  , gtry  , 
     $                   tolabs, tolrel, toltny,
     $                   alfa  , alfbst, fbest , gbest )

      implicit           double precision (a-h,o-z)
      logical            first , debug , done  , imprvd

************************************************************************
*     srchc  finds a sequence of improving estimates of a minimizer of
*     the univariate function f(alpha) in the interval (0,alfmax].
*     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
*     srchc  requires both  f(alpha)  and  f'(alpha)  to be evaluated at
*     points in the interval.  Estimates of the minimizer are computed
*     using safeguarded cubic interpolation.
*
*     Reverse communication is used to allow the calling program to
*     evaluate f and f'.  Some of the parameters must be set or tested
*     by the calling program.  The remainder would ordinarily be local
*     variables.
*
*     Input parameters (relevant to the calling program)
*     --------------------------------------------------
*
*     first         must be true on the first entry. It is subsequently
*                   altered by srchc.
*
*     debug         specifies whether detailed output is wanted.
*
*     maxf          is an upper limit on the number of times srchc is to
*                   be entered consecutively with done = false
*                   (following an initial entry with first = true).
*
*     alfa          is the first estimate of a minimizer.  alfa is
*                   subsequently altered by srchc (see below).
*
*     alfmax        is the upper limit of the interval to be searched.
*
*     epsaf         is an estimate of the absolute precision in the
*                   computed value of f(0).
*
*     ftry, gtry    are the values of f, f'  at the new point
*                   alfa = alfbst + xtry.
*
*     g0            is the value of f'(0).  g0 must be negative.
*
*     tolabs,tolrel define a function tol(alfa) = tolrel*alfa + tolabs
*                   such that if f has already been evaluated at alfa,
*                   it will not be evaluated closer than tol(alfa).
*                   These values may be reduced by srchc.
*
*     targtg        is the target value of abs(f'(alfa)). The search
*                   is terminated when 
*                    abs(f'(alfa)) le targtg and f(alfa) lt 0.
*
*     toltny        is the smallest value that tolabs is allowed to be
*                   reduced to.
*
*     Output parameters (relevant to the calling program)
*     ---------------------------------------------------
*
*     imprvd        is true if the previous alfa was the best point so
*                   far.  Any related quantities should be saved by the
*                   calling program (e.g., gradient arrays) before
*                   paying attention to the variable done.
*
*     done = false  means the calling program should evaluate
*                      ftry = f(alfa),  gtry = f'(alfa)
*                   for the new trial alfa, and re-enter srchc.
*
*     done = true   means that no new alfa was calculated.  The value
*                   of inform gives the result of the search as follows
*
*                   inform = 1 means the search has terminated
*                              successfully with alfbst < alfmax.
*
*                   inform = 2 means the search has terminated
*                              successfully with alfbst = alfmax.
*
*                   inform = 3 means that the search failed to find a 
*                              point of sufficient decrease in maxf
*                              functions, but a lower point was found.
*
*                   inform = 4 means alfmax is so small that a search
*                              should not have been attempted.
*
*                   inform = 5 is never set by srchc.
*
*                   inform = 6 means the search has failed to find a
*                              useful step.  The interval of uncertainty 
*                              is [0,b] with b < 2*tolabs. A minimizer
*                              lies very close to alfa = 0, or f'(0) is
*                              not sufficiently accurate.
*
*                   inform = 7 if no better point could be found after 
*                              maxf  function calls.
*
*                   inform = 8 means the input parameters were bad.
*                              alfmax le toltny  or g0 ge zero.
*                              No function evaluations were made.
*
*     numf          counts the number of times srchc has been entered
*                   consecutively with done = false (i.e., with a new
*                   function value ftry).
*
*     alfa          is the point at which the next function ftry and
*                   derivative gtry must be computed.
*
*     alfbst        should be accepted by the calling program as the
*                   approximate minimizer, whenever srchc returns
*                   inform = 1 or 2 (and possibly 3).
*
*     fbest, gbest  will be the corresponding values of f, f'.
*
*
*     The following parameters retain information between entries
*     -----------------------------------------------------------
*
*     braktd        is false if f and f' have not been evaluated at
*                   the far end of the interval of uncertainty.  In this
*                   case, the point b will be at alfmax + tol(alfmax).
*
*     crampd        is true if alfmax is very small (le tolabs).  If the
*                   search fails, this indicates that a zero step should
*                   be taken.
*
*     extrap        is true if xw lies outside the interval of
*                   uncertainty.  In this case, extra safeguards are
*                   applied to allow for instability in the polynomial
*                   fit.
*
*     moved         is true if a better point has been found, i.e.,
*                   alfbst gt 0.
*
*     wset          records whether a second-best point has been
*                   determined it will always be true when convergence
*                   is tested.
*
*     nsamea        is the number of consecutive times that the 
*                   left-hand end point of the interval of uncertainty
*                   has remained the same.
*
*     nsameb        similarly for the right-hand end.
*
*     a, b, alfbst  define the current interval of uncertainty. 
*                   A minimizer lies somewhere in the interval 
*                   [alfbst + a, alfbst + b].
*
*     alfbst        is the best point so far.  It is always at one end 
*                   of the interval of uncertainty.  hence we have
*                   either  a lt 0,  b = 0  or  a = 0,  b gt 0.
*
*     fbest, gbest  are the values of f, f' at the point alfbst.
*
*     factor        controls the rate at which extrapolated estimates 
*                   of alfa may expand into the interval of uncertainty.
*                   factor is not used if a minimizer has been bracketed
*                   (i.e., when the variable braktd is true).
*
*     fw, gw        are the values of f, f' at the point alfbst + xw.
*                   they are not defined until wset is true.
*
*     xtry          is the trial point in the shifted interval (a, b).
*
*     xw            is such that  alfbst + xw  is the second-best point.
*                   it is not defined until  wset  is true.
*                   in some cases,  xw  will replace a previous  xw  
*                   that has a lower function but has just been excluded
*                   from the interval of uncertainty.
*
*
*     Systems Optimization Laboratory, Stanford University, California.
*     Original version February 1982.  Rev. May 1983.
*     Original f77 version 22-August-1985.
*     This version of srchc dated  24-Oct-91.
************************************************************************

      logical            braktd, crampd, extrap, moved , wset
      save               braktd, crampd, extrap, moved , wset

      save               nsamea, nsameb
      save               a     , b     , factor
      save               xtry  , xw    , fw    , gw    , tolmax

      logical            badfun, closef, found 
      logical            quitF , quitI
      logical            fitok , setxw 
      intrinsic          abs   , sqrt

      parameter        ( zero  =0.0d+0, point1 =0.1d+0, half   =0.5d+0 )
      parameter        ( one   =1.0d+0, two    =2.0d+0, three  =3.0d+0 )
      parameter        ( five  =5.0d+0, ten    =1.0d+1, eleven =1.1d+1 )

*     ------------------------------------------------------------------
*     Local variables
*     ===============
*
*     closef     is true if the new function ftry is within epsaf of
*                fbest (up or down).
*
*     found      is true if the sufficient decrease conditions hold at
*                alfbst.
*
*     quitF      is true when  maxf  function calls have been made.
*
*     quitI      is true when the interval of uncertainty is less than
*                2*tol.
*  ---------------------------------------------------------------------

      badfun = .false.
      quitF  = .false.
      quitI  = .false.
      imprvd = .false.

      if (first) then
*        ---------------------------------------------------------------
*        First entry.  Initialize various quantities, check input data
*        and prepare to evaluate the function at the initial alfa.
*        ---------------------------------------------------------------
         first  = .false.
         numf   = 0
         alfbst = zero
         badfun = alfmax .le. toltny  .or.  g0 .ge. zero
         done   = badfun
         moved  = .false.

         if (.not. done) then
            braktd = .false.
            crampd = alfmax .le. tolabs
            extrap = .false.
            wset   = .false.
            nsamea = 0
            nsameb = 0

            tolmax = tolabs + tolrel*alfmax
            a      = zero
            b      = alfmax + tolmax
            factor = five
            tol    = tolabs
            xtry   = alfa
            if (debug) then
               write (nout, 1000) g0    , tolabs, alfmax, 
     $                            targtg, tolrel, epsaf , crampd
            end if
         end if
      else
*        ---------------------------------------------------------------
*        Subsequent entries. The function has just been evaluated at
*        alfa = alfbst + xtry,  giving ftry and gtry.
*        ---------------------------------------------------------------
         if (debug) write (nout, 1100) alfa, ftry, gtry

         numf   = numf   + 1
         nsamea = nsamea + 1
         nsameb = nsameb + 1

         if (.not. braktd) then
            tolmax = tolabs + tolrel*alfmax
            b      = alfmax - alfbst + tolmax
         end if

*        See if the new step is better.  If alfa is large enough that
*        ftry can be distinguished numerically from zero,  the function
*        is required to be sufficiently negative.

         closef = abs( ftry - fbest ) .le.  epsaf
         if (closef) then
            imprvd =  abs( gtry ) .le. abs( gbest )
         else
            imprvd = ftry .lt. fbest
         end if

         if (imprvd) then

*           We seem to have an improvement.  The new point becomes the
*           origin and other points are shifted accordingly.

            fw     = fbest
            fbest  = ftry
            gw     = gbest
            gbest  = gtry
            alfbst = alfa
            moved  = .true.

            a      = a    - xtry
            b      = b    - xtry
            xw     = zero - xtry
            wset   = .true.
            extrap =       xw .lt. zero  .and.  gbest .lt. zero
     $               .or.  xw .gt. zero  .and.  gbest .gt. zero

*           Decrease the length of the interval of uncertainty.

            if (gtry .le. zero) then
               a      = zero
               nsamea = 0
            else
               b      = zero
               nsameb = 0
               braktd = .true.
            end if
         else

*           The new function value is not better than the best point so
*           far.  The origin remains unchanged but the new point may
*           qualify as xw.  xtry must be a new bound on the best point.

            if (xtry .le. zero) then
               a      = xtry
               nsamea = 0
            else
               b      = xtry
               nsameb = 0
               braktd = .true.
            end if

*           If xw has not been set or ftry is better than fw, update the
*           points accordingly.

            if (wset) then
               setxw = ftry .lt. fw  .or.  .not. extrap
            else
               setxw = .true.
            end if

            if (setxw) then
               xw     = xtry
               fw     = ftry
               gw     = gtry
               wset   = .true.
               extrap = .false.
            end if
         end if

*        ---------------------------------------------------------------
*        Check the termination criteria.  wset will always be true.
*        ---------------------------------------------------------------
         tol    = tolabs + tolrel*alfbst
         truea  = alfbst + a
         trueb  = alfbst + b

         found  = abs(gbest) .le. targtg
         quitF  = numf       .ge. maxf
         quitI  = b - a      .le. tol + tol

         if (quitI  .and. .not. moved) then

*           The interval of uncertainty appears to be small enough,
*           but no better point has been found.  Check that changing 
*           alfa by b-a changes f by less than epsaf.

            tol    = tol/ten
            tolabs = tol
            quitI  = abs(fw) .le. epsaf  .or.  tol .le. toltny
         end if

         done  = quitF  .or.  quitI  .or.  found

         if (debug) then
            write (nout, 1200) truea    , trueb , b - a , tol   ,
     $                         nsamea   , nsameb, numf  , 
     $                         braktd   , extrap, closef, imprvd,
     $                         found    , quitI ,
     $                         alfbst   , fbest , gbest ,
     $                         alfbst+xw, fw    , gw
         end if

*        ---------------------------------------------------------------
*        Proceed with the computation of a trial steplength.
*        The choices are...
*        1. Parabolic fit using derivatives only, if the f values are
*           close.
*        2. Cubic fit for a minimizer, using both f and f'.
*        3. Damped cubic or parabolic fit if the regular fit appears to
*           be consistently overestimating the distance to a minimizer.
*        4. Bisection, geometric bisection, or a step of  tol  if
*           choices 2 or 3 are unsatisfactory.
*        ---------------------------------------------------------------
         if (.not. done) then
            xmidpt = half*(a + b)
            s      = zero
            q      = zero

            if (closef) then
*              ---------------------------------------------------------
*              Fit a parabola to the two best gradient values.
*              ---------------------------------------------------------
               s      = gbest
               q      = gbest - gw
               if (debug) write (nout, 2200)
            else
*              ---------------------------------------------------------
*              Fit cubic through  fbest  and  fw.
*              ---------------------------------------------------------
               if (debug) write (nout, 2100)
               fitok  = .true.
               r      = three*(fbest - fw)/xw + gbest + gw
               absr   = abs( r )
               s      = sqrt( abs( gbest ) ) * sqrt( abs( gw ) )

*              Compute  q =  the square root of  r*r - gbest*gw.
*              The method avoids unnecessary underflow and overflow.

               if ((gw .lt. zero  .and.  gbest .gt. zero) .or.
     $             (gw .gt. zero  .and.  gbest .lt. zero)) then
                  scale  = absr + s
                  if (scale .eq. zero) then
                     q  = zero
                  else
                     q  = scale*sqrt( (absr/scale)**2 + (s/scale)**2 )
                  end if
               else if (absr .ge. s) then
                  q     = sqrt(absr + s)*sqrt(absr - s)
               else
                  fitok = .false.
               end if

               if (fitok) then

*                 Compute a minimizer of the fitted cubic.

                  if (xw .lt. zero) q = - q
                  s  = gbest -  r - q
                  q  = gbest - gw - q - q
               end if
            end if
*           ------------------------------------------------------------
*           Construct an artificial interval  (artifa, artifb)  in which
*           the new estimate of a minimizer must lie.  Set a default
*           value of xtry that will be used if the polynomial fit fails.
*           ------------------------------------------------------------
            artifa = a
            artifb = b
            if (.not. braktd) then

*              A minimizer has not been bracketed.  Set an artificial
*              upper bound by expanding the interval  xw  by a suitable
*              factor.

               xtry   = - factor*xw
               artifb =   xtry
               if (alfbst + xtry .lt. alfmax) factor = five*factor

            else if (extrap) then

*              The points are configured for an extrapolation.
*              Set a default value of  xtry  in the interval  (a, b)
*              that will be used if the polynomial fit is rejected.  In
*              the following,  dtry  and  daux  denote the lengths of
*              the intervals  (a, b)  and  (0, xw)  (or  (xw, 0),  if
*              appropriate).  The value of  xtry is the point at which
*              the exponents of  dtry  and  daux  are approximately
*              bisected.

               daux = abs( xw )
               dtry = b - a
               if (daux .ge. dtry) then
                  xtry = five*dtry*(point1 + dtry/daux)/eleven
               else
                  xtry = half * sqrt( daux ) * sqrt( dtry )
               end if
               if (xw .gt. zero)   xtry = - xtry
               if (debug) write (nout, 2400) xtry, daux, dtry

*              Reset the artificial bounds.  If the point computed by
*              extrapolation is rejected,  xtry will remain at the
*              relevant artificial bound.

               if (xtry .le. zero) artifa = xtry
               if (xtry .gt. zero) artifb = xtry
            else

*              The points are configured for an interpolation.  The
*              default value xtry bisects the interval of uncertainty.
*              the artificial interval is just (a, b).

               xtry   = xmidpt
               if (debug) write (nout, 2300) xtry
               if (nsamea .ge. 3  .or.  nsameb .ge. 3) then

*                 If the interpolation appears to be overestimating the
*                 distance to a minimizer,  damp the interpolation.

                  factor = factor / five
                  s      = factor * s
               else
                  factor = one
               end if
            end if
*           ------------------------------------------------------------
*           The polynomial fits give  (s/q)*xw  as the new step.
*           Reject this step if it lies outside  (artifa, artifb).
*           ------------------------------------------------------------
            if (q .ne. zero) then
               if (q .lt. zero) s = - s
               if (q .lt. zero) q = - q
               if (s*xw .ge. q*artifa  .and.  s*xw .le. q*artifb) then

*                 Accept the polynomial fit.

                  if (abs( s*xw ) .ge. q*tol) then
                     xtry = (s/q)*xw
                  else
                     xtry = zero
                  end if
                  if (debug) write (nout, 2500) xtry
               end if
            end if
         end if
      end if

*     ==================================================================

      if (.not. done) then
         alfa  = alfbst + xtry
         if (braktd  .or.  alfa .lt. alfmax - tolmax) then

*           The function must not be evaluated too close to a or b.
*           (It has already been evaluated at both those points.)

            if (xtry .le. a + tol  .or.  xtry .ge. b - tol) then
               if (half*(a + b) .le. zero) then
                  xtry = - tol
               else
                  xtry =   tol
               end if
               alfa = alfbst + xtry
            end if
         else

*           The step is close to, or larger than alfmax, replace it by
*           alfmax to force evaluation of  f  at the boundary.

            braktd = .true.
            xtry   = alfmax - alfbst
            alfa   = alfmax
         end if
      end if

*     ------------------------------------------------------------------
*     Exit.
*     ------------------------------------------------------------------
      if (done) then
         if      (badfun) then
            inform = 8
         else if (found) then
            if (alfbst .lt. alfmax) then
               inform = 1
            else
               inform = 2
            end if
         else if (moved ) then
            inform = 3
         else if (quitF) then
            inform = 7
         else if (crampd) then
            inform = 4
         else
            inform = 6
         end if
      end if

      if (debug) write (nout, 3000)
      return

 1000 format(/'     g0  tolabs  alfmax        ', 1p, 2e22.14,   e16.8
     $       /' targtg  tolrel   epsaf        ', 1p, 2e22.14,   e16.8
     $       /' crampd                        ',  l3)
 1100 format(/' alfa    ftry    gtry          ', 1p, 2e22.14,   e16.8)
 1200 format(/' a       b       b - a   tol   ', 1p, 2e22.14,  2e16.8
     $       /' nsamea  nsameb  numf          ', 3i3
     $       /' braktd  extrap  closef  imprvd', 4l3
     $       /' found   quitI                 ', 2l3
     $       /' alfbst  fbest   gbest         ', 1p, 3e22.14          
     $       /' alfaw   fw      gw            ', 1p, 3e22.14)
 2100 format( ' Cubic.   ')
 2200 format( ' Parabola.')
 2300 format( ' Bisection.              xmidpt', 1p,  e22.14)
 2400 format( ' Geo. bisection. xtry,daux,dtry', 1p, 3e22.14)
 2500 format( ' Polynomial fit accepted.  xtry', 1p,  e22.14)
 3000 format( ' ----------------------------------------------------'/)

*     End of  srchc .
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine srchq ( first , debug , done  , imprvd, inform,
     $                   maxf  , numf  , nout  , 
     $                   alfmax, alfsml, epsaf , 
     $                   g0    , targtg, ftry  , 
     $                   tolabs, tolrel, toltny,
     $                   alfa  , alfbst, fbest  )

      implicit           double precision (a-h,o-z)
      logical            first , debug , done  , imprvd

************************************************************************
*     srchq  finds a sequence of improving estimates of a minimizer of
*     the univariate function f(alpha) in the interval (0,alfmax].
*     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
*     srchq  requires  f(alpha) (but not f'(alpha)) to be evaluated
*     in the interval.  New estimates of a minimizer are computed using
*     safeguarded quadratic interpolation.
*
*     Reverse communication is used to allow the calling program to
*     evaluate f.  Some of the parameters must be set or tested by the
*     calling program.  The remainder would ordinarily be local
*     variables.
*
*     Input parameters (relevant to the calling program)
*     --------------------------------------------------
*
*     first         must be true on the first entry.  It is subsequently
*                   altered by srchq.
*
*     debug         specifies whether detailed output is wanted.
*
*     maxf          is an upper limit on the number of times srchq is to
*                   be entered consecutively with done = false 
*                   (following an initial entry with first = true).
*
*     alfa          is the first estimate of a minimizer.  alfa is
*                   subsequently altered by srchq (see below).
*
*     alfmax        is the upper limit of the interval to be searched.
*
*     alfsml        is intended to prevent inefficiency when a minimizer
*                   is very small, for cases where the calling program
*                   would prefer to redefine f'(alfa).  alfsml is
*                   allowed to be zero.  Early termination will occur if
*                   srchq determines that a minimizer lies somewhere in
*                   the interval [0, alfsml) (but not if alfmax is 
*                   smaller that alfsml).
*
*     epsaf         is an estimate of the absolute precision in the
*                   computed value of f(0).
*
*     ftry          the value of f at the new point
*                   alfa = alfbst + xtry.
*
*     g0            is the value of f'(0).  g0 must be negative.
*
*     tolabs,tolrel define a function tol(alfa) = tolrel*alfa + tolabs
*                   such that if f has already been evaluated at alfa,
*                   it will not be evaluated closer than tol(alfa).
*                   These values may be reduced by srchc.
*
*     targtg        is the target value of abs(f'(alfa)). The search
*                   is terminated when 
*                    abs(f'(alfa)) le targtg and f(alfa) lt 0.
*
*     toltny        is the smallest value that tolabs is allowed to be
*                   reduced to.
*
*     Output parameters (relevant to the calling program)
*     ---------------------------------------------------
*
*     imprvd        is true if the previous alfa was the best point so
*                   far.  Any related quantities should be saved by the
*                   calling program (e.g., arrays) before paying
*                   attention to the variable done.
*
*     done = false  means the calling program should evaluate ftry
*                   for the new trial step alfa, and reenter srchq.
*
*     done = true   means that no new alfa was calculated.  The value
*                   of inform gives the result of the search as follows
*
*                   inform = 1 means the search has terminated 
*                              successfully with alfbst < alfmax.
*
*                   inform = 2 means the search has terminated
*                              successfully with alfbst = alfmax.
*
*                   inform = 3 means that the search failed to find a 
*                              point of sufficient decrease in maxf
*                              functions, but a lower point was found.
*
*                   inform = 4 means alfmax is so small that a search
*                              should not have been attempted.
*
*                   inform = 5 means that the search was terminated
*                              because of alfsml (see above).
*
*                   inform = 6 means the search has failed to find a
*                              useful step.  The interval of uncertainty 
*                              is [0,b] with b < 2*tolabs. A minimizer
*                              lies very close to alfa = 0, or f'(0) is
*                              not sufficiently accurate.
*
*                   inform = 7 if no better point could be found after 
*                              maxf  function calls.
*
*                   inform = 8 means the input parameters were bad.
*                              alfmax le toltny  or  g0 ge zero.
*                              No function evaluations were made.
*
*     numf          counts the number of times srchq has been entered
*                   consecutively with done = false (i.e., with a new
*                   function value ftry).
*
*     alfa          is the point at which the next function ftry must 
*                   be computed.
*
*     alfbst        should be accepted by the calling program as the
*                   approximate minimizer, whenever srchq returns
*                   inform = 1, 2 or 3.
*
*     fbest         will be the corresponding value of f.
*
*     The following parameters retain information between entries
*     -----------------------------------------------------------
*
*     braktd        is false if f has not been evaluated at the far end
*                   of the interval of uncertainty.  In this case, the
*                   point b will be at alfmax + tol(alfmax).
*
*     crampd        is true if alfmax is very small (le tolabs).  If the
*                   search fails, this indicates that a zero step should
*                   be taken.
*
*     extrap        is true if alfbst has moved at least once and xv 
*                   lies outside the interval of uncertainty.  In this
*                   case, extra safeguards are applied to allow for
*                   instability in the polynomial fit.
*
*     moved         is true if a better point has been found, i.e., 
*                   alfbst gt 0.
*
*     vset          records whether a third-best point has been defined.
*
*     wset          records whether a second-best point has been 
*                   defined.  It will always be true by the time the
*                   convergence test is applied.
*
*     nsamea        is the number of consecutive times that the
*                   left-hand end point of the interval of uncertainty
*                   has remained the same.
*
*     nsameb        similarly for the right-hand end.
*
*     a, b, alfbst  define the current interval of uncertainty.
*                   A minimizer lies somewhere in the  interval
*                   [alfbst + a, alfbst + b].
*
*     alfbst        is the best point so far.  It lies strictly within
*                   [atrue,btrue]  (except when alfbst has not been
*                   moved, in which case it lies at the left-hand end
*                   point).  Hence we have a .le. 0 and b .gt. 0.
*
*     fbest         is the value of f at the point alfbst.
*
*     fa            is the value of f at the point alfbst + a.
*
*     factor        controls the rate at which extrapolated estimates of
*                   alfa  may expand into the interval of uncertainty.
*                   Factor is not used if a minimizer has been bracketed
*                   (i.e., when the variable braktd is true).
*
*     fv, fw        are the values of f at the points alfbst + xv  and
*                   alfbst + xw.  They are not defined until  vset  or
*                   wset  are true.
*
*     xtry          is the trial point within the shifted interval
*                   (a, b).  The new trial function value must be
*                   computed at the point alfa = alfbst + xtry.
*
*     xv            is such that alfbst + xv is the third-best point. 
*                   It is not defined until vset is true.
*
*     xw            is such that alfbst + xw is the second-best point. 
*                   It is not defined until wset is true.  In some
*                   cases,  xw will replace a previous xw that has a
*                   lower function but has just been excluded from 
*                   (a,b).
*
*     Systems Optimization Laboratory, Stanford University, California.
*     Original version February 1982.  Rev. May 1983.
*     Original F77 version 22-August-1985.
*     This version of srchq dated  24-Oct-91.
************************************************************************

      logical            braktd, crampd, extrap, moved , vset  , wset
      save               braktd, crampd, extrap, moved , vset  , wset

      save               nsamea, nsameb
      save               a     , b     , fa    , factor
      save               xtry  , xw    , fw    , xv    , fv    , tolmax

      logical            badfun, closef, found 
      logical            quitF , quitFZ, quitI , quitS 
      logical            setxv , xinxw
      intrinsic          abs   , sqrt

      parameter        ( zero  =0.0d+0, point1 =0.1d+0, half   =0.5d+0 )
      parameter        ( one   =1.0d+0, two    =2.0d+0, five   =5.0d+0 )
      parameter        ( ten   =1.0d+1, eleven =1.1d+1                 )

*     ------------------------------------------------------------------
*     Local variables
*     ===============
*
*     closef     is true if the worst function fv is within epsaf of
*                fbest (up or down).
*
*     found      is true if the sufficient decrease conditions holds at
*                alfbst.
*
*     quitF      is true when  maxf  function calls have been made.
*
*     quitFZ     is true when the three best function values are within
*                epsaf of each other, and the new point satisfies
*                fbest le ftry le fbest+epsaf.
*
*     quitI      is true when the interval of uncertainty is less than
*                2*tol.
*
*     quitS      is true as soon as alfa is too small to be useful;
*                i.e., btrue le alfsml.
*
*     xinxw      is true if xtry is in (xw,0) or (0,xw).
*     ------------------------------------------------------------------

      imprvd = .false.
      badfun = .false.
      quitF  = .false.
      quitFZ = .false.
      quitS  = .false.
      quitI  = .false.

      if (first) then
*        ---------------------------------------------------------------
*        First entry.  Initialize various quantities, check input data
*        and prepare to evaluate the function at the initial step alfa.
*        ---------------------------------------------------------------
         first  = .false.
         numf   = 0
         alfbst = zero
         badfun = alfmax .le. toltny  .or.  g0 .ge. zero
         done   = badfun
         moved  = .false.

         if (.not. done) then
            braktd = .false.
            crampd = alfmax .le. tolabs
            extrap = .false.
            vset   = .false.
            wset   = .false.
            nsamea = 0
            nsameb = 0

            tolmax = tolrel*alfmax + tolabs
            a      = zero
            b      = alfmax + tolmax
            fa     = zero
            factor = five
            tol    = tolabs
            xtry   = alfa
            if (debug) then
               write (nout, 1000) g0    , tolabs, alfmax, 
     $                            targtg, tolrel, epsaf , crampd
            end if
         end if
      else
*        ---------------------------------------------------------------
*        Subsequent entries.  The function has just been evaluated at
*        alfa = alfbst + xtry,  giving ftry.
*        ---------------------------------------------------------------
         if (debug) write (nout, 1100) alfa, ftry

         numf   = numf   + 1
         nsamea = nsamea + 1
         nsameb = nsameb + 1

         if (.not. braktd) then
            tolmax = tolabs + tolrel*alfmax
            b      = alfmax - alfbst + tolmax
         end if

*        Check if xtry is in the interval (xw,0) or (0,xw).

         if (wset) then
            xinxw =        zero .lt. xtry  .and.  xtry .le. xw
     $               .or.    xw .le. xtry  .and.  xtry .lt. zero
         else
            xinxw = .false.
         end if

         imprvd = ftry .lt. fbest
         if (vset) then
            closef = abs( fbest - fv ) .le. epsaf
         else
            closef = .false.
         end if

         if (imprvd) then

*           We seem to have an improvement.  The new point becomes the
*           origin and other points are shifted accordingly.

            if (wset) then
               xv     = xw - xtry
               fv     = fw
               vset   = .true.
            end if

            xw     = zero - xtry
            fw     = fbest
            wset   = .true.
            fbest  = ftry
            alfbst = alfa
            moved  = .true.

            a      = a    - xtry
            b      = b    - xtry
            extrap = .not. xinxw

*           Decrease the length of (a,b).

            if (xtry .ge. zero) then
               a      = xw
               fa     = fw
               nsamea = 0
            else
               b      = xw
               nsameb = 0
               braktd = .true.
            end if
         else if (closef  .and.  ftry - fbest .lt. epsaf) then

*           Quit if there has been no progress and ftry, fbest, fw
*           and fv are all within epsaf of each other.

            quitFZ = .true.
         else

*           The new function value is no better than the current best
*           point.  xtry must an end point of the new (a,b).

            if (xtry .lt. zero) then
               a      = xtry
               fa     = ftry
               nsamea = 0
            else
               b      = xtry
               nsameb = 0
               braktd = .true.
            end if

*           The origin remains unchanged but xtry may qualify as xw.

            if (wset) then
               if (ftry .lt. fw) then
                  xv     = xw
                  fv     = fw
                  vset   = .true.

                  xw     = xtry
                  fw     = ftry
                  if (moved) extrap = xinxw
               else if (moved) then
                  if (vset) then
                     setxv = ftry .lt. fv  .or.  .not. extrap
                  else
                     setxv = .true.
                  end if

                  if (setxv) then
                     if (vset  .and.  xinxw) then
                        xw = xv
                        fw = fv
                     end if
                     xv   = xtry
                     fv   = ftry
                     vset = .true.
                  end if
               else
                  xw  = xtry
                  fw  = ftry
               end if
            else
               xw     = xtry
               fw     = ftry
               wset   = .true.
            end if
         end if

*        ---------------------------------------------------------------
*        Check the termination criteria.
*        ---------------------------------------------------------------
         tol    = tolabs + tolrel*alfbst
         truea  = alfbst + a
         trueb  = alfbst + b

         found  = moved  .and.  abs(fa - fbest) .le. -a*targtg
         quitF  = numf  .ge. maxf
         quitI  = b - a .le. tol + tol
         quitS  = trueb .le. alfsml

         if (quitI  .and.  .not. moved) then

*           The interval of uncertainty appears to be small enough,
*           but no better point has been found.  Check that changing 
*           alfa by b-a changes f by less than epsaf.

            tol    = tol/ten
            tolabs = tol
            quitI  = abs(fw) .le. epsaf  .or.  tol .le. toltny
         end if

         done  = quitF  .or.  quitFZ  .or.  quitS  .or.  quitI
     $                  .or.  found

         if (debug) then
            write (nout, 1200) truea    , trueb , b-a   , tol   ,
     $                         nsamea   , nsameb, numf  ,
     $                         braktd   , extrap, closef, imprvd,
     $                         found    , quitI , quitFZ, quitS ,
     $                         alfbst   , fbest ,
     $                         alfbst+xw, fw
            if (vset) then
               write (nout, 1300) alfbst + xv, fv
            end if
         end if

*        ---------------------------------------------------------------
*        Proceed with the computation of an estimate of a minimizer.
*        The choices are...
*        1. Parabolic fit using function values only.
*        2. Damped parabolic fit if the regular fit appears to be
*           consistently overestimating the distance to a minimizer.
*        3. Bisection, geometric bisection, or a step of tol if the
*           parabolic fit is unsatisfactory.
*        ---------------------------------------------------------------
         if (.not. done) then
            xmidpt = half*(a + b)
            s      = zero
            q      = zero

*           ============================================================
*           Fit a parabola.
*           ============================================================
*           See if there are two or three points for the parabolic fit.

            gw = (fw - fbest)/xw
            if (vset  .and.  moved) then

*              Three points available.  Use fbest, fw and fv.

               gv = (fv - fbest)/xv
               s  = gv - (xv/xw)*gw
               q  = two*(gv - gw)
               if (debug) write (nout, 2200)
            else

*              Only two points available.  Use fbest, fw and g0.

               if (moved) then
                  s  = g0 - two*gw
               else
                  s  = g0
               end if
               q = two*(g0 - gw)
               if (debug) write (nout, 2100)
            end if

*           ------------------------------------------------------------
*           Construct an artificial interval (artifa, artifb) in which 
*           the new estimate of the steplength must lie.  Set a default
*           value of  xtry  that will be used if the polynomial fit is
*           rejected. In the following, the interval (a,b) is considered
*           the sum of two intervals of lengths  dtry  and  daux, with
*           common end point the best point (zero).  dtry is the length
*           of the interval into which the default xtry will be placed
*           and endpnt denotes its non-zero end point.  The magnitude of
*           xtry is computed so that the exponents of dtry and daux are
*           approximately bisected.
*           ------------------------------------------------------------
            artifa = a
            artifb = b
            if (.not. braktd) then

*              A minimizer has not yet been bracketed.  
*              Set an artificial upper bound by expanding the interval
*              xw  by a suitable factor.

               xtry   = - factor*xw
               artifb =   xtry
               if (alfbst + xtry .lt. alfmax) factor = five*factor
            else if (vset .and. moved) then

*              Three points exist in the interval of uncertainty.
*              Check if the points are configured for an extrapolation
*              or an interpolation.

               if (extrap) then

*                 The points are configured for an extrapolation.

                  if (xw .lt. zero) endpnt = b
                  if (xw .gt. zero) endpnt = a
               else

*                 If the interpolation appears to be overestimating the
*                 distance to a minimizer,  damp the interpolation step.

                  if (nsamea .ge. 3  .or.   nsameb .ge. 3) then
                     factor = factor / five
                     s      = factor * s
                  else
                     factor = one
                  end if

*                 The points are configured for an interpolation.  The
*                 artificial interval will be just (a,b).  Set endpnt so
*                 that xtry lies in the larger of the intervals (a,b) 
*                 and  (0,b).

                  if (xmidpt .gt. zero) then
                     endpnt = b
                  else
                     endpnt = a
                  end if

*                 If a bound has remained the same for three iterations,
*                 set endpnt so that  xtry  is likely to replace the
*                 offending bound.

                  if (nsamea .ge. 3) endpnt = a
                  if (nsameb .ge. 3) endpnt = b
               end if

*              Compute the default value of  xtry.

               dtry = abs( endpnt )
               daux = b - a - dtry
               if (daux .ge. dtry) then
                  xtry = five*dtry*(point1 + dtry/daux)/eleven
               else
                  xtry = half*sqrt( daux )*sqrt( dtry )
               end if
               if (endpnt .lt. zero) xtry = - xtry
               if (debug) write (nout, 2500) xtry, daux, dtry

*              If the points are configured for an extrapolation set the
*              artificial bounds so that the artificial interval lies
*              within (a,b).  If the polynomial fit is rejected,  xtry 
*              will remain at the relevant artificial bound.

               if (extrap) then
                  if (xtry .le. zero) then
                     artifa = xtry
                  else
                     artifb = xtry
                  end if
               end if
            else

*              The gradient at the origin is being used for the
*              polynomial fit.  Set the default xtry to one tenth xw.

               if (extrap) then
                  xtry = - xw
               else
                  xtry   = xw/ten
               end if
               if (debug) write (nout, 2400) xtry
            end if

*           ------------------------------------------------------------
*           The polynomial fits give (s/q)*xw as the new step.  Reject
*           this step if it lies outside (artifa, artifb).
*           ------------------------------------------------------------
            if (q .ne. zero) then
               if (q .lt. zero) s = - s
               if (q .lt. zero) q = - q
               if (s*xw .ge. q*artifa   .and.   s*xw .le. q*artifb) then
 
*                 Accept the polynomial fit.

                  if (abs( s*xw ) .ge. q*tol) then
                     xtry = (s/q)*xw
                  else
                     xtry = zero
                  end if
                  if (debug) write (nout, 2600) xtry
               end if
            end if
         end if
      end if
*     ==================================================================

      if (.not. done) then
         alfa  = alfbst + xtry
         if (braktd  .or.  alfa .lt. alfmax - tolmax) then

*           The function must not be evaluated too close to a or b.
*           (It has already been evaluated at both those points.)

            xmidpt = half*(a + b)
            if (xtry .le. a + tol  .or.  xtry .ge. b - tol) then
               if (xmidpt .le. zero) then
                  xtry = - tol
               else
                  xtry =   tol
               end if
            end if

            if (abs( xtry ) .lt. tol) then
               if (xmidpt .le. zero) then
                  xtry = - tol
               else
                  xtry =   tol
               end if
            end if
            alfa  = alfbst + xtry
         else

*           The step is close to or larger than alfmax, replace it by
*           alfmax to force evaluation of the function at the boundary.

            braktd = .true.
            xtry   = alfmax - alfbst
            alfa   = alfmax
         end if
      end if
*     ------------------------------------------------------------------
*     Exit.
*     ------------------------------------------------------------------
      if (done) then
         if      (badfun) then
            inform = 8
         else if (quitS ) then
            inform = 5
         else if (found) then
            if (alfbst .lt. alfmax) then
               inform = 1
            else
               inform = 2
            end if
         else if (moved ) then
            inform = 3
         else if (quitF) then
            inform = 7
         else if (crampd) then
            inform = 4
         else
            inform = 6
         end if
      end if

      if (debug) write (nout, 3000)
      return

 1000 format(/'     g0  tolabs  alfmax        ', 1p, 2e22.14,   e16.8
     $       /' targtg  tolrel   epsaf        ', 1p, 2e22.14,   e16.8
     $       /' crampd                        ',  l3)
 1100 format(/' alfa    ftry                  ', 1p,2e22.14          )
 1200 format(/' a       b       b - a   tol   ', 1p,2e22.14,   2e16.8
     $       /' nsamea  nsameb  numf          ', 3i3
     $       /' braktd  extrap  closef  imprvd', 4l3
     $       /' found   quitI   quitFZ  quitS ', 4l3
     $       /' alfbst  fbest                 ', 1p,2e22.14          
     $       /' alfaw   fw                    ', 1p,2e22.14)
 1300 format( ' alfav   fv                    ', 1p,2e22.14 /)
 2100 format( ' Parabolic fit,    two points. ')
 2200 format( ' Parabolic fit,  three points. ')
 2400 format( ' Exponent reduced.  Trial point', 1p,  e22.14)
 2500 format( ' Geo. bisection. xtry,daux,dtry', 1p, 3e22.14)
 2600 format( ' Polynomial fit accepted.  xtry', 1p,  e22.14)
 3000 format( ' ----------------------------------------------------'/)

*     End of  srchq .
      end
