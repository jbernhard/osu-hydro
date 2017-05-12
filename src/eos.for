#include "defs.h"

      Subroutine InputRegulatedEOS
      Implicit none

      double precision :: e0, p0, cs2

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    TEOSdata(EOSDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata

      double precision :: EOSe0         ! min (first) energy density
      double precision :: EOSde         ! energy density step
      double precision :: EOSEend       ! max (last) energy density
      Integer :: EOSne = EOSDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOSe0, EOSde, EOSEend, EOSne

      double precision :: escale, epow, sscale, spow
      common /EOScoeffs/ escale, epow, sscale, spow

      character(len=1000) :: find_data_file

      ! read EOS data from binary file
      open(5, file=find_data_file('eos.dat'), status='old',
     &     access='stream')
      read(5) EOSe0, EOSEend, PEOSdata, SEOSdata, TEOSdata
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

      end


C====EOS from table===================================================
      Double Precision Function PEOSL7(e)  ! for lattice P(e)
      Implicit none

      double precision :: e

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    TEOSdata(EOSDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata

      double precision :: EOSe0         ! min (first) energy density
      double precision :: EOSde         ! energy density step
      double precision :: EOSEend       ! max (last) energy density
      Integer :: EOSne                  ! number of table rows
      common /EOSdatastructure/ EOSe0, EOSde, EOSEend, EOSne

      double precision :: escale, epow, sscale, spow
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
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata

      double precision :: EOSe0         ! min (first) energy density
      double precision :: EOSde         ! energy density step
      double precision :: EOSEend       ! max (last) energy density
      Integer :: EOSne                  ! number of table rows
      common /EOSdatastructure/ EOSe0, EOSde, EOSEend, EOSne

      double precision :: escale, epow, sscale, spow
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
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata

      double precision :: EOSe0         ! min (first) energy density
      double precision :: EOSde         ! energy density step
      double precision :: EOSEend       ! max (last) energy density
      Integer :: EOSne                  ! number of table rows
      common /EOSdatastructure/ EOSe0, EOSde, EOSEend, EOSne

      double precision :: escale, epow, sscale, spow
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
