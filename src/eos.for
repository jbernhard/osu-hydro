#include "defs.h"

! Ver 1.2: EOS extended to infinite energy density by extrapolation
      Subroutine InputRegulatedEOS
      Implicit none

!=======list of parameters===================================================
      Integer, Parameter :: RegEOSMudatasize = 501   !converted EOS Mu table data size
      Integer, Parameter :: IMax_Mu = 40        !maximum allowed stable particles for partically chemical equilibrium EOS
!=======list of parameters end===============================================

!=======declear variables for reading equation of state======================
      Integer :: I, J, K, L, M, N    !loop variables
      double precision :: EOSe0=0.0d0   !lowest energy density
      double precision :: EOSde=0.0d0   !spacing of energy density
      Integer :: EOSne=EOSDATALENGTH   !total rows of energy density
      double precision :: EOSEend   !the maximum energy density in the table
      double precision :: ee
      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    TEOSdata(EOSDATALENGTH)

      double precision :: logx0, logy0, logx1, logy1
      double precision :: Pcoeff1, Pcoeff2, Scoeff1, Scoeff2,
     &                    Tcoeff1, Tcoeff2

      Integer :: Inumparticle = 0   !number of stable particles in PCE, 0 for chemical equilibrium EOS
      Integer :: EOS_Mu_ne = RegEOSMudatasize   !total rows in mu table
      double precision :: EOS_Mu_e0 = 0.0d0   !lowest energy density in mu table
      double precision :: EOS_Mu_de = 0.0d0   !spacing of energy density in mu table
      double precision :: MuEOSdata(RegEOSMudatasize, IMax_Mu)
      Integer :: IMuflag
!
!=======declear variables for reading equation of state end==================

!=======common blocks========================================================
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata
      common /EOSdatastructure/ EOSe0, EOSde, EOSne, EOSEend

      common /EOSMudata/MuEOSdata, IMuflag
      common /EOSMudatastructure/ EOS_Mu_e0, EOS_Mu_de, EOS_Mu_ne,
     &                            Inumparticle
      common /EOScoeffs/ Pcoeff1, Pcoeff2, Scoeff1, Scoeff2,
     &                   Tcoeff1, Tcoeff2
!=======common blocks end====================================================

      open(5,FILE='eos.dat',STATUS='OLD')
      do I=1,EOSne
        read(5,*) ee, PEOSdata(I), SEOSdata(I), TEOSdata(I)
        if(I.eq.1) EOSe0 = ee
        ! This is very dangerous!  Computing the step size like this,
        ! from two adjacent entries, could cause significant error.
        ! if(I.eq.2) EOSde = ee - EOSe0
        ! It's much safer to divide the full range by the number of steps:
        if(I.eq.EOSne) then
          EOSde = (ee - EOSe0)/(EOSne-1)
          EOSEend = ee
        end if
      enddo
      close(5)

      ! This is even crazier!  Just read it from the last table entry...
      ! EOSEend = EOSe0 + EOSde*(EOSne-1)

      ! Extrapolate the EOS table using a power law ansatz (a*x^b).
      ! This works exceptionally well for the thermodynamic quantities
      ! P, S, T at high energy density (well past the transition).
      !
      ! The coefficients (a, b) are easily derived by fitting the ansatz
      ! to two points, i.e.
      !
      !   y0 = a*x0^b,  y1 = a*x1^b
      !
      ! Compute the coefficients from the last table point and
      ! another somewhat earlier point
      logx0 = log(EOSEend - 10000*EOSde)
      logx1 = log(EOSEend)

      logy0 = log(PEOSdata(EOSne - 10000))
      logy1 = log(PEOSdata(EOSne))
      Pcoeff1 = exp((logx0*logy1 - logx1*logy0) / (logx0 - logx1))
      Pcoeff2 = (logy0 - logy1) / (logx0 - logx1)

      logy0 = log(SEOSdata(EOSne - 10000))
      logy1 = log(SEOSdata(EOSne))
      Scoeff1 = exp((logx0*logy1 - logx1*logy0) / (logx0 - logx1))
      Scoeff2 = (logy0 - logy1) / (logx0 - logx1)

      logy0 = log(TEOSdata(EOSne - 10000))
      logy1 = log(TEOSdata(EOSne))
      Tcoeff1 = exp((logx0*logy1 - logx1*logy0) / (logx0 - logx1))
      Tcoeff2 = (logy0 - logy1) / (logx0 - logx1)

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
      Double Precision Function PEOSL7(ee)  ! for lattice P(e)
      Implicit none
      double precision :: ee
      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    TEOSdata(EOSDATALENGTH)
      double precision :: EOSe0         !lowest energy density
      double precision :: EOSde         !spacing of energy density
      Integer :: EOSne                  !total rows of energy density
      double precision :: EOSEend   !the maximum energy density in the table
      double precision :: Pcoeff1, Pcoeff2, Scoeff1, Scoeff2,
     &                    Tcoeff1, Tcoeff2

      common /EOSdata/PEOSdata, SEOSdata, TEOSdata
      common /EOSdatastructure/ EOSe0, EOSde, EOSne, EOSEend
      common /EOScoeffs/ Pcoeff1, Pcoeff2, Scoeff1, Scoeff2,
     &                   Tcoeff1, Tcoeff2

      ee = abs(ee)

      if (ee.lt.EOSe0) then
        PEOSL7 = PEOSdata(1)*ee/EOSe0
      else if (ee.lt.EOSEend) then
        call interpCubic(PEOSdata, EOSne, EOSe0, EOSde, ee, PEOSL7)
      else
        PEOSL7 = Pcoeff1*(ee**Pcoeff2)
      endif

      return
      end

      Double Precision Function SEOSL7(ee)  ! for lattice S(e)
      Implicit none
      double precision :: ee
      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    TEOSdata(EOSDATALENGTH)
      double precision :: EOSe0         !lowest energy density
      double precision :: EOSde         !spacing of energy density
      Integer :: EOSne                  !total rows of energy density
      double precision :: EOSEend   !the maximum energy density in the table
      double precision :: Pcoeff1, Pcoeff2, Scoeff1, Scoeff2,
     &                    Tcoeff1, Tcoeff2

      common /EOSdata/PEOSdata, SEOSdata, TEOSdata
      common /EOSdatastructure/ EOSe0, EOSde, EOSne, EOSEend
      common /EOScoeffs/ Pcoeff1, Pcoeff2, Scoeff1, Scoeff2,
     &                   Tcoeff1, Tcoeff2

      ee = abs(ee)

      if (ee.lt.EOSe0) then
        SEOSL7 = SEOSdata(1)*ee/EOSe0
      else if (ee.lt.EOSEend) then
        call interpCubic(SEOSdata, EOSne, EOSe0, EOSde, ee, SEOSL7)
      else
        SEOSL7 = Scoeff1*(ee**Scoeff2)
      endif

      return
      end

      Double Precision Function TEOSL7(ee)  ! for lattice T(e)
      Implicit none
      double precision :: ee
      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    TEOSdata(EOSDATALENGTH)
      double precision :: EOSe0         !lowest energy density
      double precision :: EOSde         !spacing of energy density
      Integer :: EOSne                  !total rows of energy density
      double precision :: EOSEend   !the maximum energy density in the table
      double precision :: Pcoeff1, Pcoeff2, Scoeff1, Scoeff2,
     &                    Tcoeff1, Tcoeff2

      common /EOSdata/PEOSdata, SEOSdata, TEOSdata
      common /EOSdatastructure/ EOSe0, EOSde, EOSne, EOSEend
      common /EOScoeffs/ Pcoeff1, Pcoeff2, Scoeff1, Scoeff2,
     &                   Tcoeff1, Tcoeff2

      ee = abs(ee)

      if (ee.lt.EOSe0) then
        TEOSL7 = TEOSdata(1)*ee/EOSe0
      else if (ee.lt.EOSEend) then
        call interpCubic(TEOSdata, EOSne, EOSe0, EOSde, ee, TEOSL7)
      else
        TEOSL7 = Tcoeff1*(ee**Tcoeff2)
      endif

      return
      end
