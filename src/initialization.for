#include "defs.h"

      Subroutine InitializeAll(NX0,NY0,NZ0,NX,NY,NZ,
     &  NXPhy0,NYPhy0,NXPhy,NYPhy,T0,DX,DY,DZ,DT,
     &  TT00,TT01,TT02,ScT00,ScT01,ScT02,Vx,Vy,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22,
     &  PScT00,PScT01,PScT02,PScT33,
     &  PScT11,PScT12,PScT22,etaTtp0,etaTtp,PPI,PISc,XiTtP0,XiTtP,
     &  U0,U1,U2, PU0,PU1,PU2,SxyT,Stotal,StotalBv,StotalSv,
     &  Ed,PL,Sd,Time,Temp0,Temp,T00,T01,T02,IAA,CofAA,PNEW,
     &  TEM0,ATEM0,EPS0,V10,V20,AEPS0,AV10,AV20,TFREEZ)

      Implicit Double Precision (A-H, O-Z)

      Dimension TT00(NX0:NX, NY0:NY, NZ0:NZ)
      Dimension ScT00(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension TT01(NX0:NX, NY0:NY, NZ0:NZ)
      Dimension ScT01(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension TT02(NX0:NX, NY0:NY, NZ0:NZ)
      Dimension ScT02(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension T00(NX0:NX, NY0:NY, NZ0:NZ)! ideal T00  energy momentum tensor
      Dimension T01(NX0:NX, NY0:NY, NZ0:NZ)! ideal T01
      Dimension T02(NX0:NX, NY0:NY, NZ0:NZ)! ideal T02

      Dimension Vx(NX0:NX, NY0:NY, NZ0:NZ)
      Dimension Vy(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT00(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT01(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT02(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT33(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT11(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT12(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT22(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension PPI(NX0:NX, NY0:NY, NZ0:NZ)    !Bulk pressure
      Dimension PISc(NX0:NX, NY0:NY, NZ0:NZ)


      Dimension PU0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity from last time step
      Dimension PU1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Dimension PU2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

      Dimension U0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Dimension U1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Dimension U2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

      Dimension Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Dimension PL(NX0:NX, NY0:NY, NZ0:NZ) !pressure
      Dimension Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density

      Dimension Temp0(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature
      Dimension Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature
      Dimension IAA(NX0:NX, NY0:NY, NZ0:NZ)
      Dimension CofAA(0:2,NX0:NX, NY0:NY, NZ0:NZ)

      Dimension etaTtp(NX0:NX, NY0:NY, NZ0:NZ)  !extra (eta T)/tau_pi terms in I-S eqn 02/2008
      Dimension etaTtp0(NX0:NX, NY0:NY, NZ0:NZ)  !extra (eta T)/tau_pi terms in I-S eqn 02/2008

      Dimension XiTtP(NX0:NX, NY0:NY, NZ0:NZ)  !extra (Xi T)/tau_Pi terms in full I-S bulk eqn 08/2008
      Dimension XiTtP0(NX0:NX, NY0:NY, NZ0:NZ)  !extra (Xi T)/tau_Pi in last time step

CSHEN==========================================================================
C======output relaxation time for both shear and bulk viscosity================
      Dimension VRelaxT(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient relaxation time
      Dimension VRelaxT0(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient relaxation time
CSHEN==========================================================================

      ! energy density and temperature in previous step
      DIMENSION EPS0(NX0:NX,NY0:NY)
      DIMENSION TEM0(NX0:NX,NY0:NY)

      DIMENSION V10(NX0:NX,NY0:NY),V20(NX0:NX,NY0:NY)   !velocity in X Y in Previous step

      ! energy density, temperature, velocity in previous step
      DIMENSION AEPS0(NX0:NX, NY0:NY)
      DIMENSION ATEM0(NX0:NX, NY0:NY)
      DIMENSION AV10(NX0:NX, NY0:NY),AV20(NX0:NX, NY0:NY)

      ! stress tensor in previous step
      DIMENSION F0Pi00(NX0:NX,NY0:NY)
      DIMENSION F0Pi01(NX0:NX,NY0:NY)
      DIMENSION F0Pi02(NX0:NX,NY0:NY)
      DIMENSION F0Pi11(NX0:NX,NY0:NY)
      DIMENSION F0Pi12(NX0:NX,NY0:NY)
      DIMENSION F0Pi22(NX0:NX,NY0:NY)
      DIMENSION F0Pi33(NX0:NX,NY0:NY)

      double precision, parameter :: HbarC = M_HBARC

C--------------------------------------------------------------------------------------------------------
      parameter(NNEW=4)
      DIMENSION PNEW(NNEW)    !related to root finding
      Common /Tde/ Tde, Rdec1, Rdec2,TempIni !Decoupling Temperature !decoupling redious
      common/Edec/Edec
      common/Edec1/Edec1

      Common/R0Aeps/ R0,Aeps

      Integer InitialURead   ! specify if read in more profiles
      Common/LDInitial/ InitialURead

      Double Precision Time

      COMMON /IEin/ IEin     !  type of initialization  entropy/enrgy

      double precision :: VisHRG, VisMin, VisSlope, VisCurv, VisBeta
      common /VisShear/ VisHRG, VisMin, VisSlope, VisCurv, VisBeta

      Double Precision SEOSL7, PEOSL7, TEOSL7
      External SEOSL7

      Integer Initialpitensor
      Common/Initialpi/ Initialpitensor

      ! Freezeout energy density and temperature
      ee     = EDEC       ! GeV/fm^3
      TFREEZ = TEOSL7(ee) ! GeV

      Time = T0

      ! read initial energy/entropy density from file
      if (IEin == 0) then
        ! energy density
        open(10, file='ed.dat', status='old', access='stream')
        read(10) Ed(NXPhy0:NXPhy, NYPhy0:NYPhy, NZ0)
      else
        ! entropy density
        open(10, file='sd.dat', status='old', access='stream')
        read(10) Sd(NXPhy0:NXPhy, NYPhy0:NYPhy, NZ0)
        ! convert to energy density via equation of state
        do I=NXPhy0,NXPhy
          do J=NYPhy0,NYPhy
            call invertFunction_binary(
     &             SEOSL7, 0D0, 3D3, 1d-16, 1D-6,
     &             Sd(I, J, NZ0), Ed(I,J,NZ0))
          end do
        end do
      end if  ! IEin

      close(10)

      ! convert to fm^-4, add small constant (for numerical stability?)
      Ed = Ed/HbarC + 1d-10

      ! convert energy to pressure, entropy, temperature
      call EntropyTemp3 (Ed,PL, Temp,Sd,
     &         NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)

      Temp0 = Temp

      ! extra (eta T)/tau_pi terms in I-S eqn 02/2008
      etaTtp0 = (Ed + PL)*Temp / (6.0*VisBeta)

      ! initialize flow velocity
      if (InitialURead == 0) then
        ! no initial flow
        U0 = 1
        U1 = 0
        U2 = 0
      else
        ! read from files
        open(21, file='u1.dat', status='old', access='stream')
        read(21) U1(NXPhy0:NXPhy, NYPhy0:NYPhy, NZ0)
        close(21)

        open(22, file='u2.dat', status='old', access='stream')
        read(22) U2(NXPhy0:NXPhy, NYPhy0:NYPhy, NZ0)
        close(22)

        U0 = sqrt(1 + U1**2 + U2**2)
      end if  ! InitialURead

      PU0 = U0
      PU1 = U1
      PU2 = U2

      ! initialize viscous terms
      Pi00 = 0
      Pi01 = 0
      Pi02 = 0
      Pi11 = 0
      Pi12 = 0
      Pi22 = 0
      Pi33 = 0
      PPI = 0

      if (InitialURead /= 0) then
        ! read from files
        ! remember to convert to fm^-4
        open(111, file='pi11.dat', status='old', access='stream')
        read(111) Pi11(NXPhy0:NXPhy, NYPhy0:NYPhy, NZ0)
        close(111)
        Pi11 = Pi11/HbarC

        open(112, file='pi12.dat', status='old', access='stream')
        read(112) Pi12(NXPhy0:NXPhy, NYPhy0:NYPhy, NZ0)
        close(112)
        Pi12 = Pi12/HbarC

        open(122, file='pi22.dat', status='old', access='stream')
        read(122) Pi22(NXPhy0:NXPhy, NYPhy0:NYPhy, NZ0)
        close(122)
        Pi22 = Pi22/HbarC

        ! compute remaining components of the shear tensor using
        ! orgthogonality to the flow velocity and tracelessness
        Pi00 = (U1*U1*Pi11 + U2*U2*Pi22 + 2*U1*U2*Pi12)/(U0*U0)
        Pi01 = (U1*Pi11 + U2*Pi12)/U0
        Pi02 = (U1*Pi12 + U2*Pi22)/U0
        Pi33 = Pi00 - Pi11 - Pi22

        ! bulk pressure is the deviation from ideal pressure
        PPI = Ed/3 - PL

      else if (Initialpitensor == 1) then
        ! initialize to Navier-Stokes
        call TransportPi6(Pi00,Pi01,Pi02,Pi33, Pi11,Pi12,Pi22,
     &  PPI,Ed,Sd,PL,Temp,Temp0,U0,U1,U2,PU0,PU1,PU2, DX,DY,DZ,DT,
     &  NX0,NY0,NZ0, NX,NY,NZ, Time, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &  VRelaxT,VRelaxT0)
      end if

      ! initialize T^\mu\nu
      TT00 = ((Ed + PL + PPI)*U0*U0 + Pi00 - PL - PPI)*Time
      TT01 = ((Ed + PL + PPI)*U0*U1 + Pi01)*Time
      TT02 = ((Ed + PL + PPI)*U0*U2 + Pi02)*Time

       call dpSc8(TT00,TT01,TT02,ScT00,ScT01,ScT02,Vx,Vy,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22, PScT00,PScT01,PScT02,PScT33,
     &  PScT11,PScT12,PScT22,etaTtp0,etaTtp,  PPI,PISc, XiTtP0,XiTtP,
     &  U0,U1,U2, PU0,PU1,PU2,SxyT, Stotal,StotalBv,StotalSv,
     &  Ed,PL,Sd,Temp0,Temp, T00,T01,T02, IAA,CofAA,Time,DX,DY,
     &  DZ,DT,NXPhy0,NYPhy0,NXPhy,NYPhy,NX0,NX,NY0,NY,NZ0,NZ,PNEW,NNEW)

      DO 2600 J = NYPhy0-2,NYPhy+2
      DO 2600 I = NXPhy0-2,NXPhy+2
        EPS0(I,J) = Ed(I,J,NZ0)*HbarC
        V10(I,J)  = U1(I,J,NZ0)/U0(I,J,NZ0)
        V20(I,J)  = U2(I,J,NZ0)/U0(I,J,NZ0)
        TEM0(I,J) = Temp(I,J,NZ0)*HbarC

        F0Pi00(I,J) = Pi00(I,J,NZ0)
        F0Pi01(I,J) = Pi01(I,J,NZ0)
        F0Pi02(I,J) = Pi02(I,J,NZ0)
        F0Pi11(I,J) = Pi11(I,J,NZ0)
        F0Pi12(I,J) = Pi12(I,J,NZ0)
        F0Pi22(I,J) = Pi22(I,J,NZ0)
        F0Pi33(I,J) = Pi33(I,J,NZ0)

2600  CONTINUE

            DO 2610 J = NYPhy0-2,NYPhy+2
            DO 2610 I = NXPhy0-2,NXPhy+2
                AEPS0(I,J) = Ed(I,J,NZ0)
                AV10(I,J)  = U1(I,J,NZ0)/U0(I,J,NZ0)
                AV20(I,J)  = U2(I,J,NZ0)/U0(I,J,NZ0)
                ATEM0(I,J) = Temp(I,J,NZ0)
2610        CONTINUE

      End Subroutine
