*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  mcsubs.f
*
*     mchpar   mcenvn   mcenv2   mcstor   mcmin    mcclos   mcopen
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mchpar()

C***********************************************************************
C     MCHPAR  must define certain machine parameters as follows:
C
C     wmach(1)  = NBASE  = base of floating-point arithmetic.
C     wmach(2)  = NDIGIT = no. of base wmach(1) digits of precision.
C     wmach(3)  = EPS    = floating-point precision.
C     wmach(4)  = RTEPS  = sqrt(EPS).
C     wmach(5)  = RMIN   = smallest positive normalized floating-point
C                          number.
C     wmach(6)  = RTRMIN = sqrt(RMIN).
C     wmach(7)  = BIG    = a very large number that can be represented
C                          without overflow. BIG is equal to 1.0/SMALL,
C                          where SMALL is a small number for which the
C                          machine can evaluate  1.0/SMALL  without
C                          causing overflow.
C     wmach(8)  = RTBIG  = sqrt(BIG).
C     wmach(9)  = UNDFLW = 0 if underflow is not fatal, +ve otherwise.
C                          (not implemented in post-1982 versions).
C     wmach(10) = NIN    = standard file number of the input stream.
C     wmach(11) = NOUT   = standard file number of the output stream.
C***********************************************************************
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      logical            first , hdwire
      integer            emin  , nbase , ndigit, nin   , nout

      double precision   base  , eps   , big   , rmin
      double precision   small , undflw
      external           mcenv2
      intrinsic          sqrt
      save               first
      data               first / .true. /


      if (first) then
         first = .false.

*        ---------------------------------------------------------------
*        Machine-dependent code.
*        1. Set UNDFLW, NIN, NOUT, HDWIRE as desired.
*        2. If  HDWIRE = .TRUE.  set the machine constants
*               NBASE, NDIGIT, EPS, RMIN, BIG
*           in-line.  Otherwise, they will be computed by MCENVN.
*        ---------------------------------------------------------------
         undflw = 0
         nin    = 5
*        nout   = 6
         nout   = 3
         hdwire = .true.

         if (hdwire) then

*           IEEE standard floating-point arithmetic.
*           (Rounded arithmetic is mandated).

            nbase  = 2
            ndigit = 52
            base   = nbase
            eps    = base**(- ndigit)
            rmin   = base**(- 126)
            big    = base**(+ 127)
         else
            call mcenvn( nbase, ndigit, eps, emin, rmin )
            small  = rmin*nbase**4
            big    = 1.0d0/small
         end if

         wmach( 1) = nbase
         wmach( 2) = ndigit
         wmach( 3) = eps
         wmach( 4) = sqrt( eps )
         wmach( 5) = rmin
         wmach( 6) = sqrt( rmin )
         wmach( 7) = big
         wmach( 8) = sqrt( big )
         wmach( 9) = undflw
         wmach(10) = nin
         wmach(11) = nout
      end if

      return

*     end of  mchpar.

      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mcenvn( beta, t, eps, emin, rmin )
C     ------------------------------------------------------------------
C     Based on NAG Mark 1.0 release ENVIRN.
C
C     MCENVN returns the machine parameters given by:
C
C        BETA - INTEGER.
C               The base of the machine.
C
C        T    - INTEGER.
C               The number of ( BETA ) digits in the mantissa.
C
C        EPS  - REAL.
C               The smallest positive number such that
C
C                  fl( 1.0 - EPS ) .lt. 1.0,
C
C               where fl denotes the computed value.
C
C        EMIN - INTEGER.
C               The minimum exponent before (gradual) underflow occurs.
C
C        RMIN - REAL.
C               The smallest normalized number for the machine given by
C               BASE**( EMIN - 1 ), where BASE is the floating point
C               value of BETA.
C
C
C     The computation of EPS, EMIN and RMIN is based on a routine,
C     PARANOIA by W. Kahan of the University of California at Berkeley.
C
C
C     Nag Fortran 77 O( 1 ) basic linear algebra routine (ENVIRN).
C
C     -- Written on 2-June-1987.
C     Sven Hammarling, Mick Pont and Janet Welding, Nag Central Office.
C     Modified by PEG, 7-Aug-1990.
C     ------------------------------------------------------------------
      double precision   eps, rmin
      integer            beta, emin, t

      double precision   a, b, leps, lrmin, one, rbase, small, two, zero
      integer            gnmin, gpmin, i, lbeta, lemin, lt, ngnmin,
     *                   ngpmin, nout
      logical            first, iwarn, lrnd

      double precision   mcstor
      external           mcstor

      external           mcenv2, mcmin

      intrinsic          abs, max, min

      common    /sol1cm/ nout

      save               first, iwarn, lbeta, lemin, leps, lrmin, lt

      data               first/.true./, iwarn/.false./

      if (first) then
         first = .false.
         zero  = 0
         one   = 1
         two   = 2

*        LBETA, LT, LEPS, LEMIN and LRMIN are the local values of BETA,
*        T, EPS, EMIN and RMIN.
*
*        Throughout this routine we use the function MCSTOR to ensure
*        that relevant values are stored and not held in registers, or
*        are not affected by optimizers.
*
*        MCENV2 returns the parameters LBETA and LT. ( LRND is not used
*        here. )

         call mcenv2( lbeta, lt, lrnd )

*        Start to find EPS.

         b = lbeta
         if (lrnd) then
            leps = (b**(1-lt))/two
         else
            leps = b**(1-lt)
         end if

*        Computation of EPS complete. Now find EMIN.
*        Let a = + or - 1, and + or - (1 + base**(-3)).
*        Keep dividing a by BETA until (gradual) underflow occurs.
*        This is detected when we cannot recover the previous a.

         rbase = one/lbeta
         small = one
         do 20, i = 1, 3
            small = mcstor( small*rbase, zero )
   20    continue
         a     = mcstor( one, small )
         call mcmin ( ngpmin,  one, lbeta )
         call mcmin ( ngnmin, -one, lbeta )
         call mcmin ( gpmin ,    a, lbeta )
         call mcmin ( gnmin ,   -a, lbeta )

         if (ngpmin .eq. ngnmin  .and.  gpmin .eq. gnmin) then

            if (ngpmin .eq. gpmin) then
               lemin = ngpmin
*              ( Non twos-complement machines, no gradual underflow;
*              eg VAX )
            else if (gpmin-ngpmin .eq. 3) then
               lemin = ngpmin - 1 + lt
*              ( Non twos-complement machines, with gradual underflow;
*              eg IEEE standard followers )
            else
               lemin = min( ngpmin, gpmin )
*              ( A guess; no known machine )
               iwarn = .true.
            end if

         else if (ngpmin .eq. gpmin .and. ngnmin .eq. gnmin) then
            if (abs( ngpmin-ngnmin ) .eq. 1) then
               lemin = max( ngpmin, ngnmin )
*              ( Twos-complement machines, no gradual underflow;
*              eg Cyber 205 )
            else
               lemin = min( ngpmin, ngnmin )
*              ( A guess; no known machine )
               iwarn = .true.
            end if

         else if (abs(ngpmin-ngnmin) .eq. 1 .and. gpmin .eq. gnmin) then
            if (gpmin-min( ngpmin, ngnmin ) .eq. 3) then
               lemin = max( ngpmin, ngnmin ) - 1 + lt
*              ( Twos-complement machines with gradual underflow;
*              no known machine )
            else
               lemin = max( ngpmin, ngnmin )
*              ( A guess; no known machine )
               iwarn = .true.
            end if
         else
            lemin = min( ngpmin, ngnmin, gpmin, gnmin )
*           ( A guess; no known machine )
            iwarn = .true.
         end if
*        **
*        Comment out this IF block if Emin is ok
         if (iwarn) then
            first = .true.
            write (nout,fmt=99999) lemin
         end if
*        **

*        Finally compute RMIN by successive division by BETA.
*        We could compute RMIN as base**( EMIN - 1 ), but some machines
*        underflow during this computation.

         lrmin = 1
         do 40, i = 1, 1 - lemin
            lrmin = lrmin/lbeta
   40    continue
      end if

      beta  = lbeta
      t     = lt
      eps   = leps
      emin  = lemin
      rmin  = lrmin
      return

*     End of mcenvn (envirn).

99999 format( // ' WARNING. The value Emin may be incorrect:-  Emin = ',
     $           I8 / ' If, after inspection, the value Emin looks',
     $           ' acceptable please comment out ' / ' the IF block',
     $           ' as marked within the code of routine MCENVN,' /
     $           ' otherwise contact UCSD. ' / )
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mcenv2( beta, t, rnd )
C     ------------------------------------------------------------------
C     Based on NAG Mark 1.0 release ENVRON.
C
C     MCENV2 returns the machine parameters given by:
C
C        BETA - INTEGER.
C               The base of the machine.
C
C        T    - INTEGER.
C               The number of ( BETA ) digits in the mantissa.
C
C        RND  - LOGICAL.
C               Whether proper rounding ( RND = .TRUE. ) or chopping
C               ( RND = .FALSE. ) occurs in addition. This may not be a
C               reliable guide to the way in which the machine perfoms
C               its arithmetic.
C
C     The routine is based on the routine of the same name by Malcolm
C     and incorporates suggestions by Gentleman and Marovich. See
C
C        Malcolm M. A. (1972) Algorithms to reveal properties of
C           floating-point arithmetic. Comms. of the ACM, 15, 949-951.
C
C        Gentleman W. M. and Marovich S. B. (1974) More on algorithms
C           that reveal properties of floating point arithmetic units.
C           Comms. of the ACM, 17, 276-277.
C
C
C     Nag Fortran 77 O( 1 ) basic linear algebra routine (envron).
C
C     -- Written on 26-November-1984.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C     Modified by PEG, 7-Aug-1990.
C     ------------------------------------------------------------------
      integer            beta, t
      logical            rnd

      double precision   a, b, c, c1, c2, f, one, qtr, sava
      integer            lbeta, lt
      logical            first, lrnd

      double precision   mcstor
      external           mcstor

      save               first, lbeta, lrnd, lt

      data               first/.true./


      if (first) then
         first = .false.
         one   = 1

*        LBETA, LT and LRND are the local values of BETA, T and RND.
*
*        Throughout this routine we use the function MCSTOR to ensure
*        that relevant values are stored and not held in registers, or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.

         a  = 1
         c  = 1

*       +       while( c .eq. one )loop
   20    if (c .eq. one) then
            a  = 2*a
            c  = mcstor( a, one )
            c  = mcstor( c, -a  )
            go to 20
         end if
*       +       end while

*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.

         b  = 1
         c  = mcstor( a, b )

*       +       while( c .eq. a )loop
   40    if (c .eq. a) then
            b  = 2*b
            c  = mcstor( a, b )
            go to 40
         end if
*       +       end while

*        Now compute the base. a and b are neighbouring floating point
*        numbers in the interval ( beta**t, beta**( t + 1 ) ) and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).

         qtr   = one/4
         c     = mcstor( c, -a )
         lbeta = c + qtr

*        Now determine whether rounding or chopping occurs, by adding
*        a bit less than beta/2 and a bit more than beta/2 to a.

         b    = lbeta
         f    = mcstor( b/2, -b/100 )
         c1   = mcstor( f  ,  a     )
         f    = mcstor( b/2,  b/100 )
         c2   = mcstor( f  ,  a     )
         sava = a

*        Now find the mantissa, t. It should be the integer part of
*        log to the base beta of a, however it is safer to determine t
*        by powering. So we find t as the smallest positive integer
*        for which
*
*           fl( beta**t + 1.0 ) = 1.0.

         lt = 0
         a  = 1
         c  = 1

*        +       while( c .eq. one )loop
   60    if (c .eq. one) then
            lt  = lt + 1
            a   = a*lbeta
            c   = mcstor( a, one )
            c   = mcstor( c,  -a )
            go to 60
         end if
*        +       end while

         lrnd = c1 .eq. sava  .and.  c2 .ne. sava
      end if

      beta = lbeta
      t    = lt
      rnd  = lrnd

      return

*     End of mcenv2 (envron).

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function mcstor( a, b )
C     ------------------------------------------------------------------
C     Based on NAG Mark 1.0 release.
C
C     MCSTOR is intended to force A and B to be stored prior to doing the
C     addition of A and B. For use in situations where optimizers might
C     hold one of these in a register.
C
C
C     Nag Fortran 77 O( 1 ) basic linear algebra routine (mcstor).
C
C     -- Written on 28-November-1984.
C     Sven Hammarling, Nag Central Office.
C     ------------------------------------------------------------------
      double precision                a, b

      mcstor = a + b

      return

C     end of mcstor.

      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mcmin ( emin, start, base )
C     ------------------------------------------------------------------
C     Based on NAG Mark 1.0 release.
C
C     Service routine for ENVIRN (mcenv2).
C
C
C     Nag Fortran 77 O( 1 ) basic linear algebra routine (getmin).
C
C     -- Written on 2-June-1987.
C     Mick Pont, Nag Central Office.
C     ------------------------------------------------------------------
      double precision   start
      integer            base, emin

      double precision   a, b1, b2, c1, c2, d1, d2, one, rbase, zero
      integer            i

      double precision   mcstor
      external           mcstor

      a     = start
      one   = 1
      rbase = one/base
      zero  = 0
      emin  = 1
      b1    = mcstor( a*rbase, zero )
      c1    = a
      c2    = a
      d1    = a
      d2    = a


   20 if ((c1 .eq. a)  .and.  (c2 .eq. a)  .and.
     $    (d1 .eq. a)  .and.  (d2 .eq. a)      ) then
         emin = emin - 1
         a    = b1
         b1   = mcstor(  a/base, zero )
         c1   = mcstor( b1*base, zero )
         d1   = zero
         do 40, i = 1, base
            d1 = d1 + b1
   40    continue
         b2   = mcstor(  a*rbase, zero )
         c2   = mcstor( b2/rbase, zero )
         d2   = zero
         do 60, i = 1, base
            d2 = d2 + b2
   60    continue
         go to 20
      end if

      return
      
*     End of mcmin (getmin).

      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mcclos( lun )

C***********************************************************************
C     MCCLOS  closes the file with logical unit number LUN.
C
C***********************************************************************

      close ( lun )

*     end of mcclos.
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mcopen( lun )

C***********************************************************************
C     MCOPEN  opens the file with logical unit number LUN.
C
C***********************************************************************

      open( lun, status='unknown' )

*     end of mcopen.
      end
