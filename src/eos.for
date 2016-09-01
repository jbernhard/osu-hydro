#include "defs.h"

      Subroutine InputRegulatedEOS
      Implicit none

!=======list of parameters===================================================
      Integer, Parameter :: RegEOSMudatasize = 501   !converted EOS Mu table data size
      Integer, Parameter :: IMax_Mu = 40        !maximum allowed stable particles for partically chemical equilibrium EOS
!=======list of parameters end===============================================

      Integer :: I
      double precision :: EOSe0=0.0d0   !lowest energy density
      double precision :: EOSde=0.0d0   !spacing of energy density
      double precision :: EOSEend   !the maximum energy density in the table
      Integer :: EOSne=EOSDATALENGTH
      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    TEOSdata(EOSDATALENGTH)

      double precision :: e0, p0, cs2
      double precision :: escale, epow, sscale, spow

      Integer :: Inumparticle = 0   !number of stable particles in PCE, 0 for chemical equilibrium EOS
      Integer :: EOS_Mu_ne = RegEOSMudatasize   !total rows in mu table
      double precision :: EOS_Mu_e0 = 0.0d0   !lowest energy density in mu table
      double precision :: EOS_Mu_de = 0.0d0   !spacing of energy density in mu table
      double precision :: MuEOSdata(RegEOSMudatasize, IMax_Mu)
      Integer :: IMuflag

      common /EOSdata/PEOSdata, SEOSdata, TEOSdata
      common /EOSdatastructure/ EOSe0, EOSde, EOSne, EOSEend

      common /EOSMudata/MuEOSdata, IMuflag
      common /EOSMudatastructure/ EOS_Mu_e0, EOS_Mu_de, EOS_Mu_ne,
     &                            Inumparticle
      common /EOScoeffs/ escale, epow, sscale, spow


      open(5,FILE='eos.dat',STATUS='OLD')

      ! read first row, save energy density as EOSe0
      I = 1
      read(5,*) EOSe0, PEOSdata(I), SEOSdata(I), TEOSdata(I)

      ! read rest of table, save last energy density as EOSEend
      do I=2,EOSne
        read(5,*) EOSEend, PEOSdata(I), SEOSdata(I), TEOSdata(I)
      enddo

      close(5)

      ! save energy density step size
      EOSde = (EOSEend - EOSe0)/(EOSne - 1)

      ! Prepare for extrapolating EOS beyond end of table.
      ! Assume a modified power-law ansatz for pressure in the high
      ! energy density limit:
      !
      !   p(e) = e/3 + a*e^b
      !
      ! For b < 1, this converges to the ideal gas limit as e -> inf and
      ! the speed of sound c_s^2 = dp/de -> 1/3 as e -> inf.
      ! The coefficients (a, b) may be determined from the pressure and
      ! its derivative (i.e. the speed of sound) at some reference
      ! energy density e0.  Use the second-to-last table entry as the
      ! reference point and estimate the derivative using a finite
      ! difference.
      e0 = EOSEend - EOSde
      p0 = PEOSdata(EOSne - 1)
      cs2 = (PEOSdata(EOSne) - PEOSdata(EOSne - 2)) / (2*EOSde)

      ! Use variable 'epow' for the exponent 'b' and 'escale' for the
      ! scale factor 'a'.
      epow = e0*(1 - 3*cs2)/(e0 - 3*p0)
      escale = (p0 - e0/3) * e0**(-epow)

      ! Having extrapolated p(e), the entropy density follows from the
      ! differential equation
      !
      !   (e + p(e))/s(e) = 1/(ds(e)/de)
      !
      ! (Note: both sides of this equation are the temperature.)
      ! For the chosen p(e) ansatz this diffeq can be analytically
      ! integrated for s(e).  These coefficients (spow, sscale) are used
      ! in SEOSL7, below.
      spow = 3/(4*(1 - epow))
      sscale = SEOSdata(EOSne - 1) * (e0**epow/(e0 + p0))**spow

      ! disable PCE
      !open(5,FILE='EOS/EOS_tables/EOS_particletable.dat', STATUS='OLD')
      !read(5,*) Inumparticle
      !close(5)
      !if(Inumparticle.ne.0) then
      !   open(5,FILE='EOS/EOS_tables/EOS_Mu.dat', STATUS='OLD')
      !   do I=1,EOS_Mu_ne
      !      read(5,*) ee, (MuEOSdata(I,J), J=1, Inumparticle)
      !      if(I.eq.1) EOS_Mu_e0 = ee
      !      if(I.eq.2) EOS_Mu_de = ee - EOS_Mu_e0
      !   enddo
      !endif

      end


C====EOS from table===================================================
      Double Precision Function PEOSL7(e)  ! for lattice P(e)
      Implicit none
      double precision :: e
      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    TEOSdata(EOSDATALENGTH)
      double precision :: EOSe0         !lowest energy density
      double precision :: EOSde         !spacing of energy density
      Integer :: EOSne                  !total rows of energy density
      double precision :: EOSEend   !the maximum energy density in the table
      double precision :: escale, epow, sscale, spow

      common /EOSdata/PEOSdata, SEOSdata, TEOSdata
      common /EOSdatastructure/ EOSe0, EOSde, EOSne, EOSEend
      common /EOScoeffs/ escale, epow, sscale, spow

      e = abs(e)

      if (e.lt.EOSe0) then
        PEOSL7 = PEOSdata(1)*e/EOSe0
      else if (e.lt.EOSEend) then
        call interpCubic(PEOSdata, EOSne, EOSe0, EOSde, e, PEOSL7)
      else
        ! extrapolate, see notes above
        PEOSL7 = e/3 + escale*e**epow
      endif

      return
      end

      Double Precision Function SEOSL7(e)  ! for lattice S(e)
      Implicit none
      double precision :: e
      double precision :: PEOSL7
      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    TEOSdata(EOSDATALENGTH)
      double precision :: EOSe0         !lowest energy density
      double precision :: EOSde         !spacing of energy density
      Integer :: EOSne                  !total rows of energy density
      double precision :: EOSEend   !the maximum energy density in the table
      double precision :: escale, epow, sscale, spow

      common /EOSdata/PEOSdata, SEOSdata, TEOSdata
      common /EOSdatastructure/ EOSe0, EOSde, EOSne, EOSEend
      common /EOScoeffs/ escale, epow, sscale, spow

      e = abs(e)

      if (e.lt.EOSe0) then
        SEOSL7 = SEOSdata(1)*e/EOSe0
      else if (e.lt.EOSEend) then
        call interpCubic(SEOSdata, EOSne, EOSe0, EOSde, e, SEOSL7)
      else
        ! extrapolate, see notes above
        SEOSL7 = sscale * ((e + PEOSL7(e))/e**epow)**spow
      endif

      return
      end

      Double Precision Function TEOSL7(e)  ! for lattice T(e)
      Implicit none
      double precision :: e
      double precision :: PEOSL7, SEOSL7
      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    TEOSdata(EOSDATALENGTH)
      double precision :: EOSe0         !lowest energy density
      double precision :: EOSde         !spacing of energy density
      Integer :: EOSne                  !total rows of energy density
      double precision :: EOSEend   !the maximum energy density in the table
      double precision :: escale, epow, sscale, spow

      common /EOSdata/PEOSdata, SEOSdata, TEOSdata
      common /EOSdatastructure/ EOSe0, EOSde, EOSne, EOSEend
      common /EOScoeffs/ escale, epow, sscale, spow

      e = abs(e)

      if (e.lt.EOSe0) then
        TEOSL7 = TEOSdata(1)*e/EOSe0
      else if (e.lt.EOSEend) then
        call interpCubic(TEOSdata, EOSne, EOSe0, EOSde, e, TEOSL7)
      else
        ! use extrapolated pressure and entropy density
        ! T = (e + p)/s
        TEOSL7 = (e + PEOSL7(e))/SEOSL7(e)
      endif

      return
      end
