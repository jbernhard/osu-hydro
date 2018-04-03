C*****************************************************************************
C                                                                            *
C         Viscous Israel Stewart Hydrodynamics in 2+1-dimention (VISH2+1)    *
C                                                                            *
C                           Huichao Song                                     *
C                                                                            *
C                    The Ohio State University                               *
C                                                                            *
C*****************************************************************************


C    REFERENCES

C   [1] H.Song and U.Heinz, Phys. Lett. B658, 279 (2008);
C   [2] H.Song and U.Heinz, Phys. Rev. C77, 064901 (2008);
C   [3] H.Song and U.Heinz, Phys. Rev. C78, 024902 (2008);
C   [4] H.Song and U.Heinz, arXiv:0909.1549 [nucl-th];
C   [5] H.Song, Ph.D thesis 2009, arXiv:0908.3656 [nucl-th].


#include "defs.h"

      Program Main

      Implicit Double Precision (A-H, O-Z)
      Integer LS
      Integer NXPhy0, NYPhy0
      Integer NXPhy, NYPhy
      Integer NX, NY, NZ
      Integer NX0,NY0, NZ0         ! dimension

      Integer InitialURead
      Common/LDInitial/ InitialURead  ! IintURead =1 read initial velocity profile

      Integer Initialpitensor
      Common/Initialpi/ Initialpitensor

      double precision :: VisT0, VisHRG, VisMin, VisSlope, VisCrv,
     &                    VisBeta
      common /VisShear/ VisT0, VisHRG, VisMin, VisSlope, VisCrv, VisBeta

      double precision :: VisBulkT0, VisBulkMax, VisBulkWidth, BulkTau
      integer :: IRelaxBulk
      common /VisBulk/ VisBulkT0, VisBulkMax, VisBulkWidth, BulkTau,
     &                 IRelaxBulk

      logical :: VisNonzero, VisBulkNonzero
      integer :: ViscousEqsType
      common /VisControl/ VisNonzero, VisBulkNonzero, ViscousEqsType

      Common /LS/ LS
      Common /R0Bdry/ R0Bdry

      COMMON /IEin/ IEin     !  type of initialization  entropy/enrgy

      double precision :: DT
      common /DT/ DT

      double precision :: DX, DY
      common /DXY/ DX, DY

      common/Edec/Edec    !decoupling energy density

      Integer MaxT

      Integer NDT, NDX, NDY
      Common /NDTXY/ NDT, NDX, NDY

      Double Precision T0 ! initial time tau_0
      Common /T0/ T0

      Double Precision R0, Aeps
      Common /R0Aeps/ R0,Aeps

      character(len=1000) :: find_data_file

      call prepareInputFun() ! this is the initialization function in InputFun.for

      ! read parameters from config file
      open(1, file=find_data_file('osu-hydro.conf'), status='old')

      ! initialization
      Read(1,*) T0               ! initial time [fm]
      Read(1,*) IEin             ! read initial condition as energy (0) or entropy (1) density
      Read(1,*) InitialURead     ! read initial flow and viscous terms if not 0
      Read(1,*) Initialpitensor  ! initialize shear tensor with zeros (0) or by Navier-Stokes (1)

      Read(1,*)

      ! grid
      Read(1,*) DT               ! timestep [fm]
      Read(1,*) DX               ! spatial step [fm]
      Read(1,*) LS               ! lattice size from origin (total size = 2*LS + 1)

      Read(1,*)

      ! freeze-out
      Read(1,*) Edec             ! decoupling energy density [GeV/fm^3]
      Read(1,*) NDT              ! freeze-out step in tau direction
      Read(1,*) NDX              ! freeze-out step in x, y directions

      Read(1,*)

      ! viscous equations
      Read(1,*) ViscousEqsType   ! Israel-Stewart (1) or 14-moment approximation (2)

      Read(1,*)

      ! shear viscosity
      Read(1,*) VisT0            ! temperature of minimum eta/s [GeV]
      Read(1,*) VisHRG           ! constant eta/s below T0
      Read(1,*) VisMin           ! eta/s at T0
      Read(1,*) VisSlope         ! slope of (eta/s)(T) above T0 [GeV^-1]
      Read(1,*) VisCrv           ! curvature of (eta/s)(T) above T0 (see readme)
      Read(1,*) VisBeta          ! shear relaxation time tau_pi = 6*VisBeta*eta/(sT)

      Read(1,*)

      ! bulk viscosity
      Read(1,*) VisBulkT0        ! peak location of Cauchy (zeta/s)(T) [GeV]
      Read(1,*) VisBulkMax       ! maximum value of zeta/s (at T0)
      Read(1,*) VisBulkWidth     ! width of (zeta/s)(T) [GeV]
      Read(1,*) IRelaxBulk       ! bulk relaxation time: critical slowing down (0), constant (1), 1.5/(2*pi*T) (2), ?? (3), ?? (4)
      Read(1,*) BulkTau          ! constant bulk relaxation time for IRelaxBulk == 1

      Close(1)

      ! read parameters from command line
      Call readInputFromCML2()

      ! copy X settings to Y
      DY = DX
      NDY = NDX

      ! radius and width of Fermi function for viscous regulation
      R0Bdry = LS*DX - 0.5
      Aeps = 0.05D0

      ! determine if shear and bulk viscosities are nonzero
      VisNonzero = (VisHRG > 1d-6) .or. (VisMin > 1d-6)
     &             .or.(VisSlope > 1d-6)
      VisBulkNonzero = (VisBulkMax > 1d-6) .and. (VisBulkWidth > 1d-6)

      DZ=0.01d0

      MaxT = int(40.0/DT)

      open(99, file='surface.dat', access='stream', status='replace')

      call InputRegulatedEOS

      NXPhy0=-LS
      NYPhy0=-LS
      NXPhy=LS
      NYPhy=LS
      NX=NXPhy+5
      NY=NYPhy+5
      NZ=1
      NX0=NXPhy0-5
      NY0=NYPhy0-5
      NZ0=1

      Call Mainpro(NX0,NY0,NZ0,NX,NY,NZ,NXPhy0,NYPhy0,
     &          NXPhy,NYPhy,T0,DX,DY,DZ,DT,MaxT,NDX,NDY,NDT)   ! main program

      close(99)

      End
!-----------------------------------------------------------------------



C######################################################################################
      Subroutine Mainpro(NX0,NY0,NZ0, NX,NY,NZ,NXPhy0,NYPhy0,
     &         NXPhy,NYPhy,T0,DX,DY,DZ,DT,MaxT,NDX,NDY,NDT)

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

      !Dimension Vperp(NX0:NX, NY0:NY, NZ0:NZ)

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

      ! backups

      Double Precision Pi00Regulated(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Pi01Regulated(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Pi02Regulated(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Pi11Regulated(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Pi12Regulated(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Pi22Regulated(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Pi33Regulated(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision PPIRegulated(NX0:NX,NY0:NY,NZ0:NZ)

C----------------------------------------------------------------------------------------------
      DIMENSION EPS0(NX0:NX,NY0:NY),EPS1(NX0:NX,NY0:NY) ! Energy density in previous and current step
      DIMENSION TEM0(NX0:NX,NY0:NY),TEM1(NX0:NX,NY0:NY) ! Temperature density in previous and current step

      DIMENSION V10(NX0:NX,NY0:NY),V20(NX0:NX,NY0:NY)   !velocity in X Y in Previous step
      DIMENSION V11(NX0:NX,NY0:NY),V21(NX0:NX,NY0:NY)   !velocity in X Y in current step

      DIMENSION AEPS0(NX0:NX, NY0:NY) ! Energy density in previous step
      DIMENSION ATEM0(NX0:NX, NY0:NY) ! Temperature in previous step
      DIMENSION AV10(NX0:NX, NY0:NY),AV20(NX0:NX, NY0:NY)   !velocity in X Y in Previous step

      DIMENSION F0Pi00(NX0:NX,NY0:NY),FPi00(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi01(NX0:NX,NY0:NY),FPi01(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi02(NX0:NX,NY0:NY),FPi02(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi11(NX0:NX,NY0:NY),FPi11(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi12(NX0:NX,NY0:NY),FPi12(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi22(NX0:NX,NY0:NY),FPi22(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi33(NX0:NX,NY0:NY),FPi33(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0PPI(NX0:NX,NY0:NY),FPPI(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step

      double precision, parameter :: HbarC = M_HBARC

      parameter(NNEW=4)
      DIMENSION PNEW(NNEW)    !related to root finding
      Common /Tde/ Tde, Rdec1, Rdec2,TempIni !Decoupling Temperature !decoupling radius
      common/Edec/Edec
      common/Edec1/Edec1

      Common /Nsm/ Nsm

      Common/R0Aeps/ R0,Aeps
      Common /R0Bdry/ R0Bdry

      COMMON /IEin/ IEin     !  type of initialization  entropy/energy

      double precision :: VisT0, VisHRG, VisMin, VisSlope, VisCrv,
     &                    VisBeta
      common /VisShear/ VisT0, VisHRG, VisMin, VisSlope, VisCrv, VisBeta

      double precision :: VisBulkT0, VisBulkMax, VisBulkWidth, BulkTau
      integer :: IRelaxBulk
      common /VisBulk/ VisBulkT0, VisBulkMax, VisBulkWidth, BulkTau,
     &                 IRelaxBulk

      logical :: VisNonzero, VisBulkNonzero
      integer :: ViscousEqsType
      common /VisControl/ VisNonzero, VisBulkNonzero, ViscousEqsType

      Double Precision SEOSL7, PEOSL7, TEOSL7
      External SEOSL7
      Integer iRegulateCounter, iRegulateCounterBulkPi


      Edec1 = Edec

!=======================================================================
!============ Initialization ===========================================
      Call InitializeAll(NX0,NY0,NZ0,NX,NY,NZ,
     &  NXPhy0,NYPhy0,NXPhy,NYPhy,T0,DX,DY,DZ,DT,
     &  TT00,TT01,TT02,ScT00,ScT01,ScT02,Vx,Vy,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22,
     &  PScT00,PScT01,PScT02,PScT33,
     &  PScT11,PScT12,PScT22,etaTtp0,etaTtp,PPI,PISc,XiTtP0,XiTtP,
     &  U0,U1,U2, PU0,PU1,PU2,SxyT,Stotal,StotalBv,StotalSv,
     &  Ed,PL,Sd,Time,Temp0,Temp,T00,T01,T02,IAA,CofAA,PNEW,
     &  TEM0,ATEM0,EPS0,V10,V20,AEPS0,AV10,AV20,TFREEZ)

       do 9999 ITime = 1,MaxT
!***********************  Begin  Time Loop ********************************

!   ---Zhi-Changes---
        Call determineR0(NX0,NY0,NZ0,NX,NY,NZ,Ed,PL,Sd,
     &  Pi00,Pi01,Pi02,Pi11,Pi12,Pi22,Pi33)  !fermi-dirac function to regulate boundary (by determine a suitable R0)

        DO 3006 K = NZ0,NZ !Store last step U0 U1 U2
        DO 3006 J = NY0,NY
        DO 3006 I = NX0,NX
            PU0(I,J,K)   = U0(I,J,K)
            PU1(I,J,K)   = U1(I,J,K)
            PU2(I,J,K)   = U2(I,J,K)
            Temp0(I,J,K)   = Temp(I,J,K)
            etaTtp0(I,J,K) = etaTtp(I,J,K)
            XiTtP0(I,J,K)  = XiTtP(I,J,K)
3006    Continue

      ! Initialize PiXXRegulated and backup oldTTXX
      if (VisNonzero .or. VisBulkNonzero) then

        If (VisNonzero) Then
          Pi00Regulated = Pi00
          Pi01Regulated = Pi01
          Pi02Regulated = Pi02
          Pi11Regulated = Pi11
          Pi22Regulated = Pi22
          Pi12Regulated = Pi12
          Pi33Regulated = Pi33
          iRegulateCounter = 0
        else
          Pi00 = 0D0
          Pi01 = 0D0
          Pi02 = 0D0
          Pi11 = 0D0
          Pi22 = 0D0
          Pi12 = 0D0
          Pi33 = 0D0
        endif

        if (VisBulkNonzero) then
          PPIRegulated = PPI
          iRegulateCounterBulkPi = 0
        else
          PPI = 0D0
        endif


      call dpSc8(TT00,TT01,TT02,ScT00,ScT01,ScT02,Vx,Vy,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22, PScT00,PScT01,PScT02,PScT33,
     &  PScT11,PScT12,PScT22,etaTtp0,etaTtp,  PPI,PISc, XiTtP0,XiTtP,
     &  U0,U1,U2, PU0,PU1,PU2,SxyT,Stotal,StotalBv,StotalSv,
     &  Ed,PL,Sd,Temp0,Temp, T00,T01,T02, IAA,CofAA,Time,DX,DY,
     &  DZ,DT,NXPhy0,NYPhy0,NXPhy,NYPhy,NX0,NX,NY0,NY,NZ0,NZ,PNEW,NNEW)  !PNEW NNEW  related to root finding

        DIFFC = 0.125D0
        !DIFFC = 0D0
        if (VisNonzero) then
         call UPShasta2 (Pi01,Vx,Vy,PScT01, NX0,NX,NY0,NY,NZ0,NZ, 0,0,  !Pi01
     &            DT,DX,DY, NXPhy0,NYPhy0,  NXPhy,NYPhy,-1,1, DIFFC)
         call UPShasta2 (Pi02,Vx,Vy,PScT02, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &            DT,DX,DY, NXPhy0,NYPhy0,  NXPhy,NYPhy,1,-1, DIFFC)
         call UPShasta2 (Pi00,Vx,Vy,PScT00, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &             DT,DX,DY, NXPhy0,NYPhy0,  NXPhy,NYPhy,1,1, DIFFC)
         call UPShasta2 (Pi33,Vx,Vy,PScT33, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &             DT,DX,DY, NXPhy0,NYPhy0,  NXPhy,NYPhy,1,1, DIFFC)
         call UPShasta2 (Pi11,Vx,Vy,PScT11, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &             DT,DX,DY, NXPhy0,NYPhy0,  NXPhy,NYPhy,1,1, DIFFC)
         call UPShasta2 (Pi12,Vx,Vy,PScT12, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &             DT,DX,DY, NXPhy0,NYPhy0,  NXPhy,NYPhy,1,1, DIFFC)
         call UPShasta2 (Pi22,Vx,Vy,PScT22, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &             DT,DX,DY, NXPhy0,NYPhy0,  NXPhy,NYPhy,1,1, DIFFC)
        endif

        if (VisBulkNonzero) then
          call UPShasta2 (PPI,Vx,Vy,PISc, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &             DT,DX,DY, NXPhy0,NYPhy0,  NXPhy,NYPhy,1,1, DIFFC)
        endif

!CHANGES
! Boundary regulation immediately
        DO K=NZ0,NZ
        DO J=NYPhy0,NYPhy
        DO I=NXPhy0,NXPhy
          xx=DX*I
          yy=DY*J
          rr=sqrt(xx**2+yy**2)
          ff=1.0/(Dexp((rr-R0Bdry)/Aeps)+1.0)
          Pi00(I,J,K) = Pi00(I,J,K)*ff
          Pi01(I,J,K) = Pi01(I,J,K)*ff
          Pi02(I,J,K) = Pi02(I,J,K)*ff
          Pi33(I,J,K) = Pi33(I,J,K)*ff
          Pi11(I,J,K) = Pi11(I,J,K)*ff
          Pi12(I,J,K) = Pi12(I,J,K)*ff
          Pi22(I,J,K) = Pi22(I,J,K)*ff
          PPI(I,J,K) = PPI(I,J,K)*ff
        End Do
        End Do
        End Do

        if (VisNonzero) then
         Do ! pi evolution
           call checkPiAll(iFailed, II, JJ, Time, Vx, Vy, Ed, PL,
     &     NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
     &     Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22,
     &     Pi00Regulated,Pi01Regulated,Pi02Regulated,Pi33Regulated,
     &     Pi11Regulated,Pi12Regulated,Pi22Regulated)
         If (iFailed==0) Exit
         call regulatePi(iRegulateCounter,Time,NX0,NY0,NZ0,NX,NY,NZ,
     &       NXPhy0,NXPhy,NYPhy0,NYPhy,
     &       Ed,PL,PPI,
     &       Pi00,Pi01,Pi02,Pi11,
     &       Pi12,Pi22,Pi33,Vx,Vy,II,JJ)
         iRegulateCounter = iRegulateCounter + 1
         End Do ! pi evolution
        endif

        if (VisBulkNonzero) then
          Do ! BulkPi evolution
            call checkBulkPi(iFailed, II, JJ, Time, PL,
     &      NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
     &      PPI)
          If (iFailed==0) Exit
          call regulateBulkPi(iRegulateCounterBulkPi,Time,NX0,NY0,NZ0,
     &        NX,NY,NZ,NXPhy0,NXPhy,NYPhy0,NYPhy,Ed,PL,PPI,II,JJ)
          iRegulateCounterBulkPi = iRegulateCounterBulkPi + 1
          End Do ! BulkPi evolution
        endif

      else ! ideal hydro

        Pi00 = 0D0
        Pi01 = 0D0
        Pi02 = 0D0
        Pi11 = 0D0
        Pi22 = 0D0
        Pi12 = 0D0
        Pi33 = 0D0
        PPI = 0D0

      call dpSc8(TT00,TT01,TT02,ScT00,ScT01,ScT02,Vx,Vy,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22, PScT00,PScT01,PScT02,PScT33,
     &  PScT11,PScT12,PScT22,etaTtp0,etaTtp,  PPI,PISc, XiTtP0,XiTtP,
     &  U0,U1,U2, PU0,PU1,PU2,SxyT,Stotal,StotalBv,StotalSv,
     &  Ed,PL,Sd,Temp0,Temp, T00,T01,T02, IAA,CofAA,Time,DX,DY,
     &  DZ,DT,NXPhy0,NYPhy0,NXPhy,NYPhy,NX0,NX,NY0,NY,NZ0,NZ,PNEW,NNEW)  !PNEW NNEW  related to root finding

      End If  ! VisNonzero .or. VisBulkNonzero


      ! T(mu,nu) evolution
      DIFFC = 0.125
      !DIFFC = 0D0
      call UPShasta2 (TT01,Vx,Vy,ScT01, NX0,NX,NY0,NY,NZ0,NZ, 0,0,   !T00
     &       DT,DX,DY,NXPhy0-3,NYPhy0-3,NXPhy+3,NYPhy+3,-1,1,DIFFC)
      call UPShasta2 (TT02,Vx,Vy,ScT02, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &       DT,DX,DY,NXPhy0-3,NYPhy0-3,NXPhy+3,NYPhy+3,1,-1,DIFFC)
      call UPShasta2 (TT00,Vx,Vy,ScT00, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &       DT,DX,DY,NXPhy0-3,NYPhy0-3,NXPhy+3,NYPhy+3,1,1,DIFFC)

      DO K=NZ0,NZ
      DO J=NYPhy0,NYPhy
      DO I=NXPhy0,NXPhy
        xx=DX*I
        yy=DY*J
        rr=sqrt(xx**2+yy**2)
        ff=1.0/(Dexp((rr-R0Bdry)/Aeps)+1.0)
        TT00(I,J,K) = TT00(I,J,K)*ff
        TT01(I,J,K) = TT01(I,J,K)*ff
        TT02(I,J,K) = TT02(I,J,K)*ff
      End Do
      End Do
      End Do

      Hc = HbarC


      Print '(I4, 2X, F7.4, 2X, F12.6, 2X, F9.6, 2X, I2, 2X, I2)',
     &      ITime, Time, maxval(Ed)*HBarC, maxval(Temp)*HBarC,
     &      iRegulateCounter, iRegulateCounterBulkPi


C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C~~~~     Freezeout Procedure (rewritten from Petor's code azhydro0p2)  A   ~~~~
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C      NDT = 5
      N   = ITime
      T   = Time

      IF (MOD(N+1,NDT).EQ.0) THEN
         X0 = 0.0
         Y0 = 0.0
         DO 5300 J = NY0,NY
         DO 5300 I = NX0,NX
            EPS1(I,J) = Ed(I,J,NZ0)*HbarC
            V11(I,J)  = Vx(I,J,NZ0)
            V21(I,J)  = Vy(I,J,NZ0)
            TEM1(I,J) = Temp(I,J,NZ0)*HbarC

            If (VisNonzero) Then
              FPi00(I,J) = Pi00(I,J,NZ0)
              FPi01(I,J) = Pi01(I,J,NZ0)
              FPi02(I,J) = Pi02(I,J,NZ0)
              FPi11(I,J) = Pi11(I,J,NZ0)
              FPi12(I,J) = Pi12(I,J,NZ0)
              FPi22(I,J) = Pi22(I,J,NZ0)
              FPi33(I,J) = Pi33(I,J,NZ0)
            End If
            if (VisBulkNonzero) then
              FPPI(I,J) = PPI(I,J,NZ0)
            endif

5300  CONTINUE

      Call FreezeoutPro10 (EDEC, NDX, NDY, NDT,
     &     EPS0,EPS1,V10,V20,V11,V21, NINT,
     &     F0Pi00,F0Pi01,F0Pi02,F0Pi33,F0Pi11,F0Pi12,F0Pi22,
     &     FPi00,FPi01,FPi02,FPi33,FPi11,FPi12,FPi22,
     &     F0PPI, FPPI,
     &     N,T,DX,DY,DT,NXPhy,NYPhy,NX0,NX,NY0,NY)

      DO 5400 J=NY0,NY
      DO 5400 I=NX0,NX
          TEM0(I,J) = TEM1(I,J)
          EPS0(I,J) = EPS1(I,J)
          V10(I,J)  = V11(I,J)
          V20(I,J)  = V21(I,J)
          TEM0(I,J) = TEM1(I,J)

          If (VisNonzero) Then
            F0Pi00(I,J) = Pi00(I,J,NZ0)
            F0Pi01(I,J) = Pi01(I,J,NZ0)
            F0Pi02(I,J) = Pi02(I,J,NZ0)
            F0Pi11(I,J) = Pi11(I,J,NZ0)
            F0Pi12(I,J) = Pi12(I,J,NZ0)
            F0Pi22(I,J) = Pi22(I,J,NZ0)
            F0Pi33(I,J) = Pi33(I,J,NZ0)
          End If
          if (VisBulkNonzero) then
            F0PPI(I,J) = PPI(I,J,NZ0)
          endif
5400  CONTINUE

      IF (NINT.EQ.0) THEN
        goto 10000
      END IF

      END IF    !  IF (MOD(N+1,NDT).EQ.0) THEN

      Time=Time+DT

9999  Continue  !***********************  End Time Loop  ******************************************
10000 Continue

      Return
      End


C##################################################################################
      Subroutine FreezeoutPro10(Edec, NDX, NDY, NDT,
     &     EPS0,EPS1,V10,V20,V11,V21, NINT,
     &     F0Pi00,F0Pi01,F0Pi02,F0Pi33,F0Pi11,F0Pi12,F0Pi22,
     &     FPi00,FPi01,FPi02,FPi33,FPi11,FPi12,FPi22,
     &     F0PPI, FPPI,
     &     N,T,DX,DY,DT,NXPhy,NYPhy,NX0,NX,NY0,NY)
*     a subroutine to calculate the freeze-out surface using cornelius from P. Huovinen
*     T=Time   N, Timestep for the largest Loop.
      Implicit none

      double precision :: Edec
      integer :: NDX, NDY, NDT, NINT
      integer :: N
      double precision :: T, X, Y
      double precision :: DX, DY, DT
      integer :: NXPHY, NYPHY, NX0, NX, NY0, NY

      double precision, Dimension(NX0:NX,NY0:NY) :: EPS0
      double precision, Dimension(NX0:NX,NY0:NY) :: EPS1 ! Energy density in previous and current step
      double precision, Dimension(NX0:NX,NY0:NY) :: V10
      double precision, Dimension(NX0:NX,NY0:NY) :: V20  !velocity in X Y in Previous step
      double precision, Dimension(NX0:NX,NY0:NY) :: V11
      double precision, Dimension(NX0:NX,NY0:NY) :: V21  !velocity in X Y in current step

      double precision, Dimension(NX0:NX,NY0:NY) :: F0Pi00
      double precision, Dimension(NX0:NX,NY0:NY) :: FPi00   !Stress Tensor in previous and current step
      double precision, Dimension(NX0:NX,NY0:NY) :: F0Pi01
      double precision, Dimension(NX0:NX,NY0:NY) :: FPi01   !Stress Tensor in previous and current step
      double precision, Dimension(NX0:NX,NY0:NY) :: F0Pi02
      double precision, Dimension(NX0:NX,NY0:NY) :: FPi02   !Stress Tensor in previous and current step
      double precision, Dimension(NX0:NX,NY0:NY) :: F0Pi11
      double precision, Dimension(NX0:NX,NY0:NY) :: FPi11   !Stress Tensor in previous and current step
      double precision, Dimension(NX0:NX,NY0:NY) :: F0Pi12
      double precision, Dimension(NX0:NX,NY0:NY) :: FPi12   !Stress Tensor in previous and current step
      double precision, Dimension(NX0:NX,NY0:NY) :: F0Pi22
      double precision, Dimension(NX0:NX,NY0:NY) :: FPi22   !Stress Tensor in previous and current step
      double precision, Dimension(NX0:NX,NY0:NY) :: F0Pi33
      double precision, Dimension(NX0:NX,NY0:NY) :: FPi33   !Stress Tensor in previous and current step
      double precision, Dimension(NX0:NX,NY0:NY) :: F0PPI
      double precision, Dimension(NX0:NX,NY0:NY) :: FPPI   !Bulk Stress Tensor in previous and current step

      double precision, Dimension(0:1,0:1,0:1) :: Cube
      double precision, Dimension(0:2,4) :: dSigma
      double precision, Dimension(0:2,4) :: Vmid
      double precision, Dimension(0:2) :: Vmidpoint
      integer :: Nsurf, Nambi, Ndisc
      integer :: iSurf

      double precision, parameter :: HbarC = M_HBARC

      Logical :: intersect
      double precision :: DTFreeze, DXFreeze, DYFreeze
      integer :: NXFUL, NYFUL
      double precision :: Tmid, Xmid, Ymid
      double precision :: v1mid, v2mid
      double precision :: CPi00, CPi01, CPi02
      double precision :: CPi11, CPi12, CPi22, CPi33
      double precision :: CPPI

!** Zhi ***
      Integer :: absI, absJ ! abs(I) and abs(J) used in the loop
      Integer :: I, J
      Integer :: tmpI

      DTFreeze = NDT*DT
      DXFreeze = NDX*DX
      DYFreeze = NDY*DY

      NYFUL = NYPhy + 2
      NXFUl = NXPhy + 2

      NINT = 0

      DO 500 absJ=0,NYFUL-NDY,NDY
      DO 510 absI=0,NXFUL-NDX,NDX
      Do 1911 tmpI=1,4 ! change this to tmpI=1,1 will output freeze-out data only in the 1st quadrant
       If (tmpI == 1) Then
         I = absI
         J = absJ
       ElseIf (tmpI == 2) Then
         I = -absI - NDX
         J = absJ
       ElseIf (tmpI == 3) Then
         I = -absI - NDX
         J = -absJ - NDY
       ElseIf (tmpI == 4) Then
         I = absI
         J = -absJ - NDY
       Else
         Print *, "You kidding. tmpI=", tmpI
         call exit(1)
       End If
       Y = J*DY
       X = I*DX

       CALL SECTIONCornelius(Edec, intersect, Cube, EPS0, EPS1,
     &                       I,J,NDX,NDY,NX0,NY0,NX,NY)

       IF (intersect) THEN
         NINT = NINT+1
         dSigma = 0.0d0
         Nsurf = 0
         Vmid = 0.0d0
         Nambi = 0
         Ndisc = 0
         CALL Cornelius2(Edec, Cube, dSigma, Nsurf, Vmid,
     &                   DTFreeze, DXFreeze, DYFreeze, Nambi, Ndisc)
         Do iSurf = 1, Nsurf  ! loop over the freeze out surface in the cube
           Tmid = T - DTFreeze + DT + Vmid(0,iSurf)
           Xmid = X + Vmid(1,iSurf)
           Ymid = Y + Vmid(2,iSurf)
           Vmidpoint = Vmid(:,iSurf)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,V10,V11,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,v1mid)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,V20,V21,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,v2mid)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,F0Pi00,FPi00,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,CPi00)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,F0Pi01,FPi01,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,CPi01)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,F0Pi02,FPi02,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,CPi02)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,F0Pi11,FPi11,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,CPi11)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,F0Pi12,FPi12,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,CPi12)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,F0Pi22,FPi22,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,CPi22)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,F0Pi33,FPi33,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,CPi33)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,F0PPI,FPPI,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,CPPI)

           write(99)
     &       Tmid, Xmid, Ymid,
     &       dSigma(:, iSurf),
     &       v1mid, v2mid,
     &       CPi00*HbarC, CPi01*HbarC, CPi02*HbarC,
     &       CPi11*HbarC, CPi12*HbarC, CPi22*HbarC,
     &       CPi33*HbarC,
     &       CPPI*HbarC

         Enddo  ! Nsurf
       ENDIF  ! intersect

1911    Continue
 510    CONTINUE
 500    CONTINUE

      Return
      End



C####################################################################
      Subroutine TransportPi6( Pi00,Pi01,Pi02,Pi33, Pi11,Pi12,Pi22,
     &  PPI,Ed, Sd, PL, Temp, Temp0, U0,U1,U2,PU0,PU1,PU2, DX,DY,DZ,DT,
     &  NX0,NY0,NZ0, NX,NY,NZ, Time, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &  VRelaxT, VRelaxT0)
C------- Transport Pi00,Pi01,Pi02,Pi03,Pi11, Pi12,Pi22 by first order theory
        Implicit Double Precision (A-H, O-Z)
        Dimension U0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U3(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

        Dimension Vx(NX0:NX, NY0:NY, NZ0:NZ) !used for checkPi
        Dimension Vy(NX0:NX, NY0:NY, NZ0:NZ) !used for checkPi

        Dimension PU0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity from last time step
        Dimension PU1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension PU2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension PU3(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

       Dimension Pi00(NX0:NX, NY0:NY, NZ0:NZ) ! Stress Tensor
       Dimension Pi01(NX0:NX, NY0:NY, NZ0:NZ)
       Dimension Pi02(NX0:NX, NY0:NY, NZ0:NZ)
       Dimension Pi33(NX0:NX, NY0:NY, NZ0:NZ)

       Dimension Pi11(NX0:NX, NY0:NY, NZ0:NZ) ! Stress Tensor
       Dimension Pi12(NX0:NX, NY0:NY, NZ0:NZ)
       Dimension Pi22(NX0:NX, NY0:NY, NZ0:NZ)

       Dimension PPI(NX0:NX, NY0:NY, NZ0:NZ)

        Dimension DPc00(NX0:NX, NY0:NY, NZ0:NZ) ! Differential part of Pi source term
        Dimension DPc01(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc02(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc33(NX0:NX, NY0:NY, NZ0:NZ) !

        Dimension DPc11(NX0:NX, NY0:NY, NZ0:NZ) ! Differential part of Pi source term
        Dimension DPc12(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc22(NX0:NX, NY0:NY, NZ0:NZ) !

        Dimension DDU0(NX0:NX, NY0:NY, NZ0:NZ) ! Differential part of Pi source term
        Dimension DDU1(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DDU2(NX0:NX, NY0:NY, NZ0:NZ) !

       Dimension Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
       Dimension Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density
       Dimension PL(NX0:NX, NY0:NY, NZ0:NZ) !pressure density
       Dimension Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature
       Dimension Temp0(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature

       Dimension VCoefi(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient
       Dimension VCBeta(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient Beta2
       Dimension VRelaxT(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient relaxation time
       Dimension SiLoc(NX0:NX, NY0:NY, NZ0:NZ) ! Local expansion rate \sita
       Dimension DLnT(NX0:NX, NY0:NY, NZ0:NZ) ! DlnT(x,y) terms

        Dimension etaTtp(NX0:NX, NY0:NY, NZ0:NZ)  !extra (eta T)/tau_pi terms in I-S eqn 02/2008

       Dimension VBulk(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient-bulk vicosity \xi
       Dimension VRelaxT0(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient relaxation time \tau_PI
       Dimension XiTtP(NX0:NX, NY0:NY, NZ0:NZ)  !extra (Xi T)/tau_Pi terms in full I-S bulk eqn 08/2008

       COMMON /IEin/ IEin  !parameter for initializtion by entropy 1 or energy 0


!   ---Changes-Zhi--- ***
      Double Precision R0,Aeps
      Common/R0Aeps/ R0,Aeps

      Vx = U1/U0
      Vy = U2/U0

       call getInitialR0(PU0,PU1,PU2,PU3,U0,U1,U2,U3,DX,DY,DZ,DT,
     & DPc00,DPc01,DPc02,DPc33, DPc11,DPc22,DPc12, DDU0,DDU1,DDU2,
     & Temp,Temp0,  SiLoc,DLnT, Time, NXPhy0,NYPhy0,NXPhy,NYPhy,
     & NX0,NX,NY0,NY,NZ0,NZ, Ed,Sd,PL,VCoefi)
!   ---Zhi-End---

        call ViscousCoefi8(Ed,Sd,PL,Temp,
     &  VCoefi,VCBeta,VRelaxT,etaTtp, VBulk, VRelaxT0,XiTtP,
     &  Time,DX,DY,DT,NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)

         do 10 k=NZ0,NZ
         do 10 j=NY0,NY
         do 10 i=NX0,NX
              Temp0(i,j,k)=Temp(i,j,k)
10       continue

        call PiS4U5(PU0,PU1,PU2,PU3,U0,U1,U2,U3, DX,DY,DZ, DT,
     & DPc00,DPc01,DPc02,DPc33, DPc11,DPc22,DPc12, DDU0,DDU1,DDU2,
     & Temp,Temp0,  SiLoc,DLnT,  Time, NXPhy0,NYPhy0,NXPhy,NYPhy,
     & NX0,NX,NY0,NY,NZ0,NZ)
C
      do 30 k=NZ0,NZ
      do 30 j=NYPhy0,NYPhy
      do 30 i=NXPhy0,NXPhy

        ConsP=VCoefi(i,j,k)*2


       Pi01(i,j,k)=ConsP*DPc01(i,j,k)
       Pi02(i,j,k)=ConsP*DPc02(i,j,k)
       Pi33(i,j,k)=ConsP*DPc33(i,j,k)

       !Pi00(i,j,k)=ConsP*DPc00(i,j,k)
       Pi11(i,j,k)=ConsP*DPc11(i,j,k)
       Pi12(i,j,k)=ConsP*DPc12(i,j,k)
       Pi22(i,j,k)=ConsP*DPc22(i,j,k)

       Pi00(i,j,k)=(Pi01(i,j,k)*U1(I,j,k)+Pi02(i,j,k)*U2(I,j,k))
     &      /DMax1(1d-34,U0(i,j,k))

        PPI(I,J,K)=(-1.0)*VBulk(i,j,k)*SiLoc(i,j,k)
30    continue

       call TriSembdary3(Pi00,Pi01, Pi02,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
       call TriSembdary3(Pi11,Pi22, Pi33,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
       call Sembdary3(Pi12, NX0,NY0,NZ0,NX,NY,NZ,
     &                         NXPhy0,NYPhy0,NXPhy,NYPhy)


       Return
       End





C####################555555555555555555555555555555555555#########################
C--------------------------------------------------------------------------------
       Subroutine PiS4U5(PU0,PU1,PU2,PU3,U0,U1,U2,U3, DX,DY,DZ, DT,
     & DPc00,DPc01,DPc02,DPc33, DPc11,DPc22,DPc12, DDU0,DDU1,DDU2,
     & Temp,Temp0,  SiLoc,DLnT,  Time, NXPhy0,NYPhy0,NXPhy,NYPhy,
     & NX0,NX,NY0,NY,NZ0,NZ)
        Implicit Double Precision (A-H, O-Z)

        Dimension PU0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity from last time step
        Dimension PU1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension PU2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension PU3(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

        Dimension U0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U3(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

        Dimension DPc00(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc01(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc02(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc33(NX0:NX, NY0:NY, NZ0:NZ) !

        Dimension DPc11(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc12(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc22(NX0:NX, NY0:NY, NZ0:NZ) !

        Dimension DDU0(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DDU1(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DDU2(NX0:NX, NY0:NY, NZ0:NZ) !

        Dimension Temp0(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature  in last time step
        Dimension Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature
        Dimension SiLoc(NX0:NX, NY0:NY, NZ0:NZ) ! Local expansion rate \sita
        Dimension DLnT(NX0:NX, NY0:NY, NZ0:NZ) ! DlnT(x,y) terms

        Common/R0Aeps/ R0,Aeps
        Common /R0Bdry/ R0Bdry

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C       Dimension A00(NX0:NX, NY0:NY, NZ0:NZ)
C       Dimension A01(NX0:NX, NY0:NY, NZ0:NZ)
C       Dimension A02(NX0:NX, NY0:NY, NZ0:NZ)

C       Dimension A10(NX0:NX, NY0:NY, NZ0:NZ)
C       Dimension A11(NX0:NX, NY0:NY, NZ0:NZ)
C       Dimension A12(NX0:NX, NY0:NY, NZ0:NZ)

C       Dimension A20(NX0:NX, NY0:NY, NZ0:NZ)
C       Dimension A21(NX0:NX, NY0:NY, NZ0:NZ)
C       Dimension A22(NX0:NX, NY0:NY, NZ0:NZ)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c       goto 199

        DO 100 K=NZ0,NZ
        DO 100 J=NYPhy0,NYPhy
        DO 100 I=NXPhy0,NXPhy

          xx=DX*I
          yy=DY*J
          rr=sqrt(xx**2+yy**2)
          ff=1.0/(Dexp((rr-R0Bdry)/Aeps)+1.0)


         D0U0=(U0(I,J,K)-PU0(I,J,K))/DT
         D0U1=(U1(I,J,K)-PU1(I,J,K))/DT
         D0U2=(U2(I,J,K)-PU2(I,J,K))/DT

#if DERIV_ORDER == 3
         D1U0=(U0(I+1,J,K)-U0(I-1,J,K))/(2.0*DX)
         D1U1=(U1(I+1,J,K)-U1(I-1,J,K))/(2.0*DX)
         D1U2=(U2(I+1,J,K)-U2(I-1,J,K))/(2.0*DX)
         D2U0=(U0(I,J+1,K)-U0(I,J-1,K))/(2.0*DY)
         D2U1=(U1(I,J+1,K)-U1(I,J-1,K))/(2.0*DY)
         D2U2=(U2(I,J+1,K)-U2(I,J-1,K))/(2.0*DY)
#elif DERIV_ORDER == 5
         D1U0=(U0(I+1,J,K)*2.0d0/3.0d0-U0(I-1,J,K)*2.0d0/3.0d0
     &            -U0(I+2,J,K)/12.0d0+U0(I-2,J,K)/12.0d0)/DX
         D1U1=(U1(I+1,J,K)*2.0d0/3.0d0-U1(I-1,J,K)*2.0d0/3.0d0
     &            -U1(I+2,J,K)/12.0d0+U1(I-2,J,K)/12.0d0)/DX
         D1U2=(U2(I+1,J,K)*2.0d0/3.0d0-U2(I-1,J,K)*2.0d0/3.0d0
     &            -U2(I+2,J,K)/12.0d0+U2(I-2,J,K)/12.0d0)/DX

         D2U0=(U0(I,J+1,K)*2.0d0/3.0d0-U0(I,J-1,K)*2.0d0/3.0d0
     &           -U0(I,J+2,K)/12.0d0+U0(I,J-2,K)/12.0d0)/DY
         D2U1=(U1(I,J+1,K)*2.0d0/3.0d0-U1(I,J-1,K)*2.0d0/3.0d0
     &           -U1(I,J+2,K)/12.0d0+U1(I,J-2,K)/12.0d0)/DY
         D2U2=(U2(I,J+1,K)*2.0d0/3.0d0-U2(I,J-1,K)*2.0d0/3.0d0
     &           -U2(I,J+2,K)/12.0d0+U2(I,J-2,K)/12.0d0)/DY
#else
#error "DERIV_ORDER must be 3 or 5"
#endif

         CS=(D0U0+D1U1+D2U2+U0(I,J,K)/Time)/3.0

      If (Time > 0.8) Then
      If (isnan(CS)) Then
        Print *, "Invalid CS!"
        Print *, "i,j,k=", i,j,k
        Print *, "D0U0=", D0U0, "D1U1=", D1U1, "D2U2=", D2U2
        Print *, "D2U1=", D2U1, "D2U0=", D2U0
        Print *, "U0=", U0(I,J,K), "Time=",Time
        Print *, "U0(I,J+1,K), U0(I,J-1,K)=", U0(I,J+1,K), U0(I,J-1,K)
        Print *, "U1(I,J+1,K), U1(I,J-1,K)=", U1(I,J+1,K), U1(I,J-1,K)
        Print *, "U2(I,J+1,K), U2(I,J-1,K)=", U2(I,J+1,K), U2(I,J-1,K)
      EndIf
      EndIf

      DU0=U0(I,J,K)*D0U0+U1(I,J,K)*D1U0+U2(I,J,K)*D2U0
      DU1=U0(I,J,K)*D0U1+U1(I,J,K)*D1U1+U2(I,J,K)*D2U1
      DU2=U0(I,J,K)*D0U2+U1(I,J,K)*D1U2+U2(I,J,K)*D2U2

      DDU0(i,j,k)=DU0
      DDU1(i,j,k)=DU1
      DDU2(i,j,k)=DU2

       DPc00(I,J,K)=D0U0-U0(I,J,K)*DU0+CS*(U0(I,J,K)**2-1.0)
       DPc01(I,J,K)=0.5*(D0U1-D1U0)-0.5*(U1(I,J,K)*DU0+U0(I,J,K)*DU1)
     &                      +CS*(U1(I,J,K)*U0(I,J,K))
       DPc02(I,J,K)=0.5*(D0U2-D2U0)-0.5*(U2(I,J,K)*DU0+U0(I,J,K)*DU2)
     &                      +CS*(U2(I,J,K)*U0(I,J,K))
       DPc33(I,J,K)=CS-U0(I,J,K)/Time

       DPc11(I,J,K)=(-1.0)*D1U1-U1(I,J,K)*DU1+CS*(U1(I,J,K)**2+1.0)
       DPc22(I,J,K)=(-1.0)*D2U2-U2(I,J,K)*DU2+CS*(U2(I,J,K)**2+1.0)
       DPc12(I,J,K)=(-0.5)*(D2U1+D1U2)-0.5*(U1(I,J,K)*DU2+U2(I,J,K)*DU1)
     &                      +CS*(U1(I,J,K)*U2(I,J,K))

         SiLoc(I,J,K)=CS*3.0  !   CS=(D0U0+D1U1+D2U2+U0(I,J,K)/Time)/3.0

         D0LnT=(DLog(Temp(I,J,K))-Dlog(Temp0(I,J,K)))/DT

C       If(abs(Accu-3.0).le.0.00001) then  !3pt formula
          D1LnT=(DLog(Temp(I+1,J,K))-Dlog(Temp(I-1,J,K)))/(2.0*DX)
          D2LnT=(DLog(Temp(I,J+1,K))-Dlog(Temp(I,J-1,K)))/(2.0*DY)
C       else if (abs(Accu-5.0).le.0.00001) then  !5pt formula
C         D1LnT=(DLog(Temp(I+1,J,K))*2.0d0/3.0d0
C     &      -Dlog(Temp(I-1,J,K))*2.0d0/3.0d0
C     &      -DLog(Temp(I+2,J,K))/12.0d0+DLog(Temp(I-2,J,K))/12.0d0 )/DX
C         D2LnT=(DLog(Temp(I,J+1,K))*2.0d0/3.0d0
C     &      -Dlog(Temp(I,J-1,K))*2.0d0/3.0d0
C     &      -DLog(Temp(I,J+2,K))/12.0d0+DLog(Temp(I,J-2,K))/12.0d0 )/DY
C        else
C           Print*, "wrong input for Accu,
C     &         Accu=3or5 for 3pt or 5pt cal of deriv. "
C        end if
       DlnT(I,J,K)=U0(I,J,K)*D0LnT +U1(I,J,K)*D1LnT +U2(I,J,K)*D2LnT         !DLnT term, added 02/2008
       DlnT(I,J,K)= DlnT(I,J,K)*ff  !make boundary works
 100   Continue

      Return
      End




       Subroutine StressTensorZero(Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22,
     &                         PPI, NX0,NY0,NZ0,  NX,NY,NZ) !Set Zero ST for Ideal Hydro
C-----------Used for Ideal Hydro related stress tensors are set zero-------------------
        Implicit Double Precision (A-H, O-Z)
        Dimension Pi00(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension Pi01(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension Pi02(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension Pi33(NX0:NX, NY0:NY, NZ0:NZ)

        Dimension Pi11(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension Pi12(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension Pi22(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension PPI(NX0:NX, NY0:NY, NZ0:NZ)
c        Dimension Pi03(NX0:NX, NY0:NY, NZ0:NZ)
c        Dimension Pi13(NX0:NX, NY0:NY, NZ0:NZ)

          do 10 I=NX0,NX
          do 10 J=NY0,NY
          do 10 K=NZ0,NZ
             Pi00(I,J,K)=0.0
             Pi01(I,J,K)=0.0
             Pi02(I,J,K)=0.0
             Pi33(I,J,K)=0.0
             Pi11(I,J,K)=0.0
             Pi12(I,J,K)=0.0
             Pi22(I,J,K)=0.0
             PPI(I,J,K)=0.0
c             Pi03(I,J,K)=0.0
c             Pi13(I,J,K)=0.0
 10      Continue

       Return
       End

C#####################################################
       Subroutine ViscousCoefi8(Ed,Sd,PL,Temp,
     &  VCoefi,VCBeta,VRelaxT,etaTtp, VBulk, VRelaxT0,XiTtP,
     &  Time,DX,DY,DT,NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
!--------call visous coeficient from from entropy
!--------bulk viscosity are included
      Implicit Double Precision (A-H, O-Z)

      Dimension Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Dimension Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density
      Dimension PL(NX0:NX, NY0:NY, NZ0:NZ) !pressure density
      Dimension Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature

      Dimension VCoefi(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient shear viscosity eta
      Dimension VCBeta(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient Beta2
      Dimension VRelaxT(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient relaxation time
      Dimension etaTtp(NX0:NX, NY0:NY, NZ0:NZ)  !extra (eta T)/tau_pi terms in I-S eqn 02/2008

      Dimension VBulk(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient-bulk vicosity \xi
      Dimension VRelaxT0(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient relaxation time \tau_PI
      Dimension XiTtP(NX0:NX, NY0:NY, NZ0:NZ)  !extra (Xi T)/tau_Pi terms in full I-S bulk eqn 08/2008

      double precision :: VisT0, VisHRG, VisMin, VisSlope, VisCrv,
     &                    VisBeta
      common /VisShear/ VisT0, VisHRG, VisMin, VisSlope, VisCrv, VisBeta

      double precision :: VisBulkT0, VisBulkMax, VisBulkWidth, BulkTau
      integer :: IRelaxBulk
      common /VisBulk/ VisBulkT0, VisBulkMax, VisBulkWidth, BulkTau,
     &                 IRelaxBulk

      logical :: VisNonzero, VisBulkNonzero
      integer :: ViscousEqsType
      common /VisControl/ VisNonzero, VisBulkNonzero, ViscousEqsType

      Common /Tde/ Tde, Rdec1, Rdec2,TempIni !Decoupling Temperature !decoupling radius
      Common/R0Aeps/ R0,Aeps
      Common/R0Bdry/ R0Bdry
      double precision, parameter :: HbarC = M_HBARC

      do 10 k=1,1
      do 10 j=NYPhy0-2,NYPhy+2 ! -2,NYPhy+2
      do 10 i=NXPhy0-2,NXPhy+2
      xx=DX*I
      yy=DY*J
      rr=sqrt(xx**2+yy**2)
      ff=1.0/(Dexp((rr-R0Bdry)/Aeps)+1.0)

      If (VisNonzero) then
        VCoefi(i,j,k) = ViscousCTemp(Temp(i,j,k)*HbarC) * Sd(i,j,k)
        VCBeta(i,j,k) = 1.D0/dmax1(Ed(i,j,k)+PL(i,j,k),1e-30)

        if(ViscousEqsType .eq. 1) then
          VRelaxT(i,j,k)=1.0
     &      /dmax1(5.0*VCoefi(i,j,k)*VCBeta(i,j,k),1e-30)
        else if(ViscousEqsType .eq. 2) then   ! 14-moment result
          VRelaxT(i,j,k)=(Ed(i,j,k)+PL(i,j,k))
     &      /dMax1(VCoefi(i,j,k)*5.D0,1e-30)
        else
          write(*, *) "No such viscous equation type:",ViscousEqsType
          call exit(1)
        end if

        etaTtp(i,j,k)=(VCoefi(i,j,k)*Temp(i,j,k))*VRelaxT(i,j,k)  ! A(e+p)T/6 !(eta T/tau_Pi) for extra term in full I-S
        etaTtp(i,j,k)=dmax1(etaTtp(i,j,k),1e-30)

        VCoefi(i,j,k)=VCoefi(i,j,k)*ff
        VRelaxT(i,j,k)=VRelaxT(i,j,k)*ff
        If (Time > 0.8) Then
        If (isnan(etaTtp(i,j,k))) Then
          Print *, "Invalid etaTtp!"
          Print *, "i,j,k=",i,j,k
          Print *, "VRelaxT=",VRelaxT(i,j,k)
          Print *, "VCoefi=", VCoefi(i,j,k), "Temp=", Temp(i,j,k)
          Print *, "VCBeta=", VCBeta(i,j,k)
          Print *, "VisBeta=", VisBeta
          Print *, "Sd=", Sd(i,j,k), "Ed=", Ed(i,j,k)
        EndIf
        EndIf
      else
        VCoefi(i,j,k)=0.0
        VCBeta(i,j,k)=0.0
        VRelaxT(i,j,k)=0.0
        etaTtp(i,j,k)=0.0
      end if

!------------ for bulk pressure------------
      If (VisBulkNonzero) then
        VBulk(i,j,k) = ViscousBulkTemp(Temp(i,j,k)*HbarC) * Sd(i,j,k)

        If (IRelaxBulk.eq.0) then
          TTpi=DMax1(0.1d0, 120* VBulk(i,j,k)/DMax1(Sd(i,j,k),0.1d0))
          VRelaxT0(i,j,k)=1.0/TTpi
        else if (IRelaxBulk.eq.1) then
          VRelaxT0(i,j,k)=1.0/BulkTau
        else if (IRelaxBulk.eq.2) then
          VRelaxT0(i,j,k)=2*M_PI*Temp(i,j,k)/1.5
        else if (IRelaxBulk .eq. 3) then
          TauPi = 9.0*VBulk(i,j,k)/(Ed(i,j,k) - 3.*PL(i,j,k))
          VRelaxT0(i,j,k) = 1.0/DMax1(0.1d0, TauPi)
        else if (IRelaxBulk .eq. 4) then
          ! lower bound of bulk relaxation time
          delta_cs2 = (1.0/3.0 - getCS2(Ed(i,j,k)*HbarC))**2
          Taupi = VBulk(i,j,k)/(14.55*delta_cs2*(Ed(i,j,k)+PL(i,j,k)))
          VRelaxT0(i,j,k) = 1.D0/DMax1(0.1D0, TauPi)
        else
          Print*,'This option is not supported by this version'
          Print*,'IRelaxBulk'
          call exit(1)
        end if

        XiTtP(i,j,k)=(VBulk(i,j,k)*Temp(i,j,k))*VRelaxT0(i,j,k)  !(Xi T/tau_Pi)  for extra term in full I-S

        VBulk(i,j,k)=VBulk(i,j,k)*ff
        VRelaxT0(i,j,k)=VRelaxT0(i,j,k)*ff
      else
        VBulk(i,j,k)=0.0
        VRelaxT0(i,j,k)=0.0
        XiTtP(i,j,k)=0.0
      end if

10    continue

      Return
      End

CSHEN======================================================================
C====eta/s dependent on local temperature==================================

      double precision function ViscousCTemp(T)
      double precision :: T

      double precision :: VisT0, VisHRG, VisMin, VisSlope, VisCrv,
     &                    VisBeta
      common /VisShear/ VisT0, VisHRG, VisMin, VisSlope, VisCrv, VisBeta

      if (T > VisT0) then
        ViscousCTemp = VisMin + VisSlope*(T - VisT0) * (T/VisT0)**VisCrv
      else
        ViscousCTemp = VisHRG
      endif

      return
      end

CSHEN=====end==============================================================

      double precision function ViscousBulkTemp(T)
      double precision :: T

      double precision :: VisBulkT0, VisBulkMax, VisBulkWidth, BulkTau
      integer :: IRelaxBulk
      common /VisBulk/ VisBulkT0, VisBulkMax, VisBulkWidth, BulkTau,
     &                 IRelaxBulk

      ! Cauchy distribution for zeta/s
      ViscousBulkTemp = VisBulkMax /
     &                  (1 + ((T - VisBulkT0)/VisBulkWidth)**2)

      return
      end


C---J.Liu------------------------------------------------------------------------

      double precision function getCS2(Ed)
      Implicit double precision (A-H, O-Z)
      ! unit of Ed should be GeV/fm^3
      de = 0.05*Ed
      p1 = PEOSL7(Ed - de/2.)
      p2 = PEOSL7(Ed + de/2.)
      cs2 = (p2 - p1)/de   !cs^2 = dP/de

      getCS2 = cs2
      Return
      End

      Subroutine ViscousShearTransCoefs(Ed, PL, taupi_inverse,
     &        deltaSpiSpi, lambdaSpiBPi, phi7, taupipi)
      ! calculate the shear transport coefficients according to
      ! PL input should be in unit of fm^-4
      Implicit Double Precision (A-H, O-Z)
      Double precision :: lambdaSpiBPi

      taupi = 1.D0/dMax1(taupi_inverse, 1e-30)

      deltaSpiSpi = 4.0*taupi/3.0
      lambdaSpiBPi= 6.0*taupi/5.0
      phi7 = 4.0*9.0/70.0/DMax1(Ed+PL, 1e-30)
      taupipi = 10.0*taupi/7.0

      Return
      End


      Subroutine ViscousBulkTransCoefs(Ed, tauBPi_inverse,
     &        deltaBPiBPi, lambdaBPiSpi)
      ! calculate the bulk transport coefficients according to
      ! Ed input should be in unit of GeV/fm^3
      Implicit Double Precision (A-H, O-Z)
      double precision lambdaBPiSpi
      tauBPi = 1.D0/dMax1(tauBPi_inverse, 1e-18)
      cs2 = getCS2(Ed)

      deltaBPiBPi = 2.D0/3.D0*tauBPi
      lambdaBPiSpi = 8.D0/5.D0*(1.D0/3.D0 - cs2)*tauBPi

      Return
      End

C---J.Liu----end-----------------------------------------------------------------

C----------------------------------------------------
       Double Precision FUNCTION BulkAdSH0(eta,tt)  !Ads-CFT Bulk-Vis=0 for HRG
       Implicit Double Precision (A-H, O-Z)
          !eta shear viscosity  Cs speed of Sound:Katz05 Lattice EOS4 here
          !tt temperature !GeV

           Cper=1.0d0/(exp((0.17d0-tt)/0.001d0)+1.0d0)
           Bper=1.0d0/(exp((0.180d0-tt)/0.01d0)+1.0d0)

       Cs1=(0.3d0*(1.0d0/(exp((0.180d0-TT)/0.025d0)+1.0d0))*Bper
     &  +(0.150d0*(1.0-exp(-((0.185d0-TT)/0.037d0)**2))+0.035d0)
     &  *(1-Bper))*Cper    !QGP

CCcc       Width=0.02d0
C       Width=0.008d0      !sharp head
C       CC0=0.07458569d0
C       Cs2=((1.0d0/3.0d0-CC0)*(1.0-exp(-((0.172d0-TT)/Width)**2))+CC0)
C     &          *(1-Cper)  !HRG

       P2=0.07458569d0    !round head
       Cs2=((1.0d0/3.0d0-P2)/(exp((tt-0.157d0)/0.002d0)+1.0d0)+P2)
     &        *(1-Cper)

       Cs=Cs1+Cs2

           Zeta=eta*2.0*(1.0/3.0-Cs)
           BulkAdSH0=Zeta
       return
       end


C###########################################################################
       Subroutine EntropyTemp3 (Ed,PL, Temp,Sd,
     &       NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
C------call Temperatuen mu  Entropy from ee pp-----------
        Implicit Double Precision (A-H, O-Z)
       Dimension Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
       Dimension PL(NX0:NX, NY0:NY, NZ0:NZ) !Pressure

       Dimension Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature
       Dimension Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density

      double precision, parameter :: HbarC = M_HBARC

      Do 1001 k=1,1                           !do 10 is in the unit of fm-1
      Do 1001 j=NYPhy0-2,NYPhy+2 ! -2,NYPhy
      Do 1001 i=NXPhy0-2,NXPhy+2
        Ed(i,j,k) = max(1e-30, Ed(i,j,k))
 1001 Continue

      do k = 1,1
      do j = NYPhy0-2,NYPhy+2
      do i = NXPhy0-2,NXPhy+2
       ee = Ed(i,j,k)*HbarC     !ee in Gev/fm^3
       PL(i,j,k)   = PEOSL7(ee)/HbarC    !PEOSL7 in GeV/fm^3, PL in fm^-4
       Temp(i,j,k) = TEOSL7(ee)/HbarC    !TEOSL7 in GeV, Temp in fm^-1
       Sd(i,j,k)   = SEOSL7(ee)          !SEOSL7 in fm^-3, Sd in fm^-3
      enddo   ! i loop
      enddo   ! j loop
      enddo   ! k loop

      Return
      end


C###################################################################
      SUBROUTINE FYSCOR(T00,T01,T02,C)
*************************************************************
** THIS ROUTINE ASSERTS THAT THE VELOCITY OF THE ENERGY     *
** FLOW CAN'T EXCEED 1 AND THAT ENERGY DENSITY IS POSITIVE  *
** OR ZERO.                                                 *
*************************************************************
      Implicit Double Precision (A-H, O-Z)
      C = 1.0
      T00=dMAX1(0.0d0,T00)
      W=SQRT(T01*T01+T02*T02)
      IF (W.GT.T00) THEN
        C=0.999*T00/W
        T01=C*T01
        T02=C*T02
      END IF

      RETURN
      END

C###################################################################
      SUBROUTINE FYSCOR_velocity(T00,T01,T02,Ed,PL,BulkPi)
*************************************************************
** THIS ROUTINE ASSERTS THAT THE VELOCITY OF THE ENERGY     *
** FLOW CAN'T EXCEED 1 AND THAT ENERGY DENSITY IS POSITIVE  *
** OR ZERO.                                                 *
*************************************************************
      Implicit none
      double precision T00, T01, T02, Ed, PL, BulkPi
      double precision cstilde2
      double precision M0, M
      double precision scaleFactor
      double precision temp1, temp2

      cstilde2 = PL/Ed
      T00 = dMax1(0.0d0, T00)
      M0 = T00
      M = sqrt(T01**2 + T02**2)
      if (M0 .lt. M) then
        scaleFactor = 0.999*M0/M
        T01 = scaleFactor*T01
        T02 = scaleFactor*T02
      endif

      temp1 = 2.*sqrt(cstilde2)*M - M0*(1.+cstilde2)
      if (BulkPi .lt. temp1) then
        BulkPi = temp1 + 1D-3*abs(temp1)
      endif

      temp2 = (M - M0)
      if (BulkPi .lt. temp2) then
        BulkPi = 0.999*temp2
      endif

      RETURN
      END

C###################################################################
      SUBROUTINE FYSCOR2(T00,T01,T02,PL,C)
*************************************************************
** THIS ROUTINE ASSERTS THAT THE VELOCITY OF THE ENERGY     *
** FLOW CAN'T EXCEED 1 AND THAT ENERGY DENSITY IS POSITIVE  *
** OR ZERO.                                                 *
*************************************************************
        Implicit Double Precision (A-H, O-Z)
      C = 1.0
      T00=dMAX1(0.0d0,T00)
      W=SQRT(T01*T01+T02*T02)
      IF (W.GT.(T00+PL)) THEN
        C=0.999*(T00+PL)/W
        T01=C*T01
        T02=C*T02
      END IF

      RETURN
      END

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C###################################################################
      SUBROUTINE FYSCOR_withBulkvis(T00,T01,T02,BulkPi)
*************************************************************
** THIS ROUTINE ASSERTS THAT THE VELOCITY OF THE ENERGY     *
** FLOW CAN'T EXCEED 1 AND THAT ENERGY DENSITY IS POSITIVE  *
** OR ZERO.                                                 *
*************************************************************
      Implicit Double Precision (A-H, O-Z)
      temp = dmax1(0.0d0, T00)
      temp = temp*(temp + BulkPi)
      temp = dmax1(0.0d0, temp)
      W = sqrt(T01*T01+T02*T02)
      IF (W .GT. sqrt(temp)) THEN ! regulate the violation
        C=0.999*sqrt(temp)/W
        T01=C*T01
        T02=C*T02
      END IF

      RETURN
      END

C###################################################################
      Subroutine dpSc8(TT00,TT01,TT02,ScT00,ScT01,ScT02,Vx,Vy,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22, PScT00,PScT01,PScT02,PScT33,
     &  PScT11,PScT12,PScT22,etaTtp0,etaTtp,  PPI,PISc, XiTtP0,XiTtP,
     &  U0,U1,U2, PU0,PU1,PU2,SxyT,Stotal,StotalBv,StotalSv,
     &  Ed,PL,Sd,Temp0,Temp, T00,T01,T02, IAA,CofAA,Time,DX,DY,
     &  DZ,DT,NXPhy0,NYPhy0,NXPhy,NYPhy,NX0,NX,NY0,NY,NZ0,NZ, PNEW,NNEW)
      !PNEW NNEW related to root finding

        Implicit Double Precision (A-H, O-Z)
        Dimension TT00(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension ScT00(NX0:NX, NY0:NY, NZ0:NZ)

        Dimension TT01(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension ScT01(NX0:NX, NY0:NY, NZ0:NZ)

        Dimension TT02(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension ScT02(NX0:NX, NY0:NY, NZ0:NZ)

        Dimension Vx(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension Vy(NX0:NX, NY0:NY, NZ0:NZ)

        Dimension U0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U3(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

        Dimension PU0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity from last time step
        Dimension PU1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension PU2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension PU3(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

        Dimension Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
        Dimension PScT00(NX0:NX, NY0:NY, NZ0:NZ)  !
        Dimension Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
        Dimension PScT01(NX0:NX, NY0:NY, NZ0:NZ)  !
        Dimension Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
        Dimension PScT02(NX0:NX, NY0:NY, NZ0:NZ)  !
        Dimension Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
        Dimension PScT33(NX0:NX, NY0:NY, NZ0:NZ)  !

        Dimension Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
        Dimension PScT11(NX0:NX, NY0:NY, NZ0:NZ)  !
        Dimension Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
        Dimension PScT12(NX0:NX, NY0:NY, NZ0:NZ)  !
        Dimension Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
        Dimension PScT22(NX0:NX, NY0:NY, NZ0:NZ)  !

       Dimension PPI(NX0:NX, NY0:NY, NZ0:NZ)    !Bulk pressure
       Dimension PISc(NX0:NX, NY0:NY, NZ0:NZ)  !

C--------------------------------------------
        Dimension T00(NX0:NX, NY0:NY, NZ0:NZ)! ideal T00  energy momentum tensor
        Dimension T01(NX0:NX, NY0:NY, NZ0:NZ)! ideal T01
        Dimension T02(NX0:NX, NY0:NY, NZ0:NZ)! ideal T02
C-------------------------------------------
        Dimension Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
        Dimension PL(NX0:NX, NY0:NY, NZ0:NZ) !local pressure
        Dimension Bd(NX0:NX, NY0:NY, NZ0:NZ) !net baryon density
        Dimension Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density

        Dimension Temp0(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature  in last time step
        Dimension Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature

        Dimension VCoefi(NX0:NX, NY0:NY, NZ0:NZ) !viscous coefficient
        Dimension VCBeta(NX0:NX, NY0:NY, NZ0:NZ) !viscous coefficient  Beta 2
        Dimension VRelaxT(NX0:NX, NY0:NY, NZ0:NZ) !viscous coefficient relaxation Time
        Dimension SiLoc(NX0:NX, NY0:NY, NZ0:NZ) ! Local expansion rate \theta
        Dimension DLnT(NX0:NX, NY0:NY, NZ0:NZ) ! DlnT(x,y) terms  !added 02/2008

       Dimension VBulk(NX0:NX, NY0:NY, NZ0:NZ) !viscous coefficient-bulk viscosity \xi
       Dimension VRelaxT0(NX0:NX, NY0:NY, NZ0:NZ) !viscous coefficient relaxation time \tau_PI

        Dimension etaTtp(NX0:NX, NY0:NY, NZ0:NZ)  !extra (eta T)/tau_pi terms in I-S eqn 02/2008
        Dimension etaTtp0(NX0:NX, NY0:NY, NZ0:NZ)  !extra (eta T)/tau_pi terms in I-S eqn 02/2008

        Dimension XiTtP(NX0:NX, NY0:NY, NZ0:NZ)  !extra (Xi T)/tau_Pi terms in full I-S bulk eqn 08/2008
        Dimension XiTtP0(NX0:NX, NY0:NY, NZ0:NZ)  !extra (Xi T)/tau_Pi in last time step

        Dimension DDU0(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DDU1(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DDU2(NX0:NX, NY0:NY, NZ0:NZ) !


        Dimension DPc00(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc01(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc02(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc33(NX0:NX, NY0:NY, NZ0:NZ) !

        Dimension DPc11(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc12(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc22(NX0:NX, NY0:NY, NZ0:NZ) !

        Dimension IAA(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension CofAA(0:2,NX0:NX, NY0:NY, NZ0:NZ)

      double precision :: VisT0, VisHRG, VisMin, VisSlope, VisCrv,
     &                    VisBeta
      common /VisShear/ VisT0, VisHRG, VisMin, VisSlope, VisCrv, VisBeta

      double precision :: VisBulkT0, VisBulkMax, VisBulkWidth, BulkTau
      integer :: IRelaxBulk
      common /VisBulk/ VisBulkT0, VisBulkMax, VisBulkWidth, BulkTau,
     &                 IRelaxBulk

      logical :: VisNonzero, VisBulkNonzero
      integer :: ViscousEqsType
      common /VisControl/ VisNonzero, VisBulkNonzero, ViscousEqsType

        DIMENSION PNEW(NNEW)!something related to root finding
       Parameter(XC0=1.0d-15, AAC=1.0d-16) !root finding accuracy

        Common/R0Aeps/ R0,Aeps
        Common /R0Bdry/ R0Bdry

      Double Precision T0
      Common /T0/ T0

       double precision, parameter :: HbarC = M_HBARC

       Common /Nsm/ Nsm
       common/Edec1/Edec1

      ! ----- Use in root search -----
      Double Precision :: RSDM0, RSDM, RSPPI !, RSee
      Common /findEdHookData/ RSDM0, RSDM, RSPPI ! M0, M, Pi (see 0510014)

      Double precision :: deltaBPiBPi, lambdaBPiSpi ! bulk transport coefficients
      Double precision :: deltaSpiSpi, lambdaSpiBPi, phi7, taupipi ! shear transport coefficients
      Double precision :: piSigma
         Nsm=5

        if(NZ0.ne.NZ)then
          Print *,' VSc2d , a 2-dimensinal subroutine '
          call exit(1)
        end if


      PNEW(1)=XC0   !dimension related to newton root finding
      PNEW(2)=33335.5
      PNEW(3)=0.3
      PNEW(4)=0.3
      AAC0=AAC
      AAC10=AAC*0.1


!   ---Zhi-Changes---
!-------Regulate Pi(mu,nu) before adding it to T tensor
!      If (VisNonzero) Then
!        call regulatePi(Time,NX0,NY0,NZ0,NX,NY,NZ,
!     &  NXPhy0,NXPhy,NYPhy0,NYPhy,
!     &  Ed,PL,PPI,
!     &  Pi00,Pi01,Pi02,Pi11,Pi12,Pi22,Pi33)
!      End If
!-------End of regulation---------
!       ---Zhi-End---

! Ver 1.6.29RC5: NXPhy0->NX0, etc

       DO 100 K=NZ0,NZ
       DO 100 J=NY0,NY
       DO 100 I=NX0,NX

        T00M=TT00(I,J,K)/Time-Pi00(I,J,K)   ! ---Zhi-Changes---
        T01M=TT01(I,J,K)/Time-Pi01(I,J,K)   ! ---Zhi-Changes---
        T02M=TT02(I,J,K)/Time-Pi02(I,J,K)   ! ---Zhi-Changes---
        PLM=PL(I,J,K)
        EdM = Ed(I,J,K)
        BulkPi = PPI(I,J,K)

        !*** ABORT SMALL NUMBERS TO AVOID UNDERFLOW **
        !CALL FYSCOR(T00M,T01M,T02M,CC)  !T00=AMIN1(T00,0.0)and  other treatment ???????
        !call FYSCOR_withBulkvis(T00M, T01M, T02M, BulkPi)
        call FYSCOR_velocity(T00M, T01M, T02M, EdM, PLM, BulkPi)

        if(T00M .gt. AAC) then
          T00(I,J,K) = T00M
          PPI(I,J,K) = BulkPi
          IF (ABS(T01M).GT.AAC) THEN
            T01(I,J,K)=T01M
          ELSE
            T01(I,J,K)=0.0
            Pi01(I,J,k)=0.0   !song viscous
          END IF
          IF (ABS(T02M).GT.AAC) THEN
            T02(I,J,K)=T02M
          ELSE
            T02(I,J,K)=0.0
            Pi02(I,J,K)=0.0
          END IF
        else
          T00(I,J,K) = 1D-30
          T01(I,J,K) = 0.0
          T02(I,J,K) = 0.0
          Pi01(I,J,k)=0.0   !song viscous
          Pi02(I,J,K)=0.0
          PPI(I,J,K) = 0.0
          PL(I,J,K) = 0.0
        endif
 100  CONTINUE

C---------------------------------------------------------------
      DO 300 K=NZ0,NZ
      DO 300 J=NYPhy0,NYPhy
      DO 310 I=NXPhy0,NXPhy !changed by song

        Bd(I,J,K)=0
        T00IJ=T00(I,J,K)
        T01IJ=T01(I,J,K)
        T02IJ=T02(I,J,K)
        DM0 = T00IJ
        DM = sqrt(T01IJ*T01IJ + T02IJ*T02IJ)

        RSDM0 = DM0
        RSDM = DM
        RSPPI = PPI(I,J,K)

        !eeH = findEdHook(1D0)
        !Call invertFunctionD(findEdHook,0D0,5D3,1D-3,ED(I,J,K),0D0,eeH)
        VP_local = findvHook(0.0D0)
        Call invertFunctionH(findvHook, 0.0D0, 1.0D0, 0.0, 1D-6,
     &                       VP_local)
        U0_guess = 1./sqrt(1. - VP_local*VP_local)
        U0_low_boundary = dmax1(1.0, 0.5*U0_guess)
        U0_upper_boundary = 1.5*U0_guess
        U0_local = findU0Hook(1D0)
        Call invertFunctionH(findU0Hook, U0_low_boundary,
     &                       U0_upper_boundary, 0.0, 1D-6,
     &                       U0_local)
        U0_critial = 1.21061
        if(U0_local .gt. U0_critial) then
          VP_local = sqrt(1. - 1./U0_local/U0_local)
        else
          U0_local = 1./sqrt(1. - VP_local*VP_local)
        endif

        eeH = DM0 - VP_local*DM

        !eeH = quadESolver(ED(I,J,K),DM0,DM,PPI(I,J,K),I,J) ! use Ed from last time step as a starting point
        ED(I,J,K) = dmax1(eeH, 1D-10)
        DM = dmax1(DM, 1D-10)
        !VP = (DM0-eeH)/DM
        VP = VP_local

        ee=Ed(I,J,NZ0)*HbarC        !fm-4 to GeV/fm3  peter unit
        pp=DMAX1(0.0D0,PEOSL7(ee))
        PL(I,J,K)=pp/HbarC           !GeV/fm3 to fm-4

        Vx(I,J,K) = VP*T01IJ/DM
        Vy(I,J,K) = VP*T02IJ/DM
        if (isnan(Vx(I,J,K))) then
          print*, 'vx', VP_local, U0_local, VP, T01IJ, DM
          call exit(1)
        endif
        if (isnan(Vy(I,J,K))) then
          print*, 'vy', VP_local, U0_local, VP, T02IJ, DM
          call exit(1)
        endif
!----------------------------------------------------------------------
 310   CONTINUE
 300   CONTINUE


!---------- End of root finding ------------------------------------------------------
! CAREFUL

      TT00=(T00+Pi00)*Time ! ---Zhi-Changes---
      TT01=(T01+Pi01)*Time ! ---Zhi-Changes---
      TT02=(T02+Pi02)*Time ! ---Zhi-Changes---


!=============== Extrapolation Process ======================================

      DO K=NZ0,NZ
      DO J=NYPhy0,NYPhy
      DO I=NXPhy0,NXPhy
        U0(I,J,K)=1D0/sqrt(1D0-Vx(I,J,K)**2-Vy(I,J,k)**2+1D-10) ! also calculate U0 to save another separated loop
        U1(I,J,K)=U0(I,J,K)*Vx(I,J,K)
        U2(I,J,K)=U0(I,J,K)*Vy(I,J,K)
      End Do
      End Do
      End Do

      call VSBdary3(U1,U2,NX0,NY0,NZ0,NX,NY,NZ, ! smearing and extrapolation
     & NXPhy0,NYPhy0, NXPhy,NYPhy,AAC0)

      DO K=NZ0,NZ
      DO J=NYPhy0-3,NYPhy+3
      DO I=NXPhy0-3,NXPhy+3
        U0(I,J,K)=sqrt(1D0+U1(I,J,K)**2+U2(I,J,k)**2) ! also calculate U0 to save another separated loop
        Vx(I,J,K)=U1(I,J,K)/U0(I,J,K)
        Vy(I,J,K)=U2(I,J,K)/U0(I,J,K)
        If (Vx(I,J,K)**2+Vy(I,J,K)**2>1D0) Then
          Print*, "Error:"
          Print*, "I,J=",I,J
          Print*, "Vx,Vy=",Vx(I,J,K),Vy(I,J,K)
          Print*, "norm=",Vx(I,J,K)**2+Vy(I,J,K)**2
          Print*, "U0,U1,U2=",U0(I,J,K),U1(I,J,K),U2(I,J,K)
          call exit(1)
        End If
      End Do
      End Do
      End Do

!      call VSBdary3(Vx,Vy,NX0,NY0,NZ0,NX,NY,NZ, ! smearing and extrapolation
!     & NXPhy0,NYPhy0, NXPhy,NYPhy,AAC0)
      !Do K=NZ0,NZ
      !Do I=NXPhy0-3,NXPhy+3
      !Do J=NYPhy0-3,NYPhy+3
      !  dlength = sqrt(Vx(I,J,K)*Vx(I,J,K) + Vy(I,J,K)*Vy(I,J,K))
      !  If (dlength>1D0) Then
      !    Vx(I,J,K) = Vx(I,J,K)/dlength*0.99999999999D0
      !    Vy(I,J,K) = Vy(I,J,K)/dlength*0.99999999999D0
      !  End If
      !  U0(I,J,K)=1D0/max(sqrt(1D0-Vx(I,J,K)**2-Vy(I,J,k)**2),1D-30) ! also calculate U0 to save another separated loop
      !  U1(I,J,K)=U0(I,J,K)*Vx(I,J,K)
      !  U2(I,J,K)=U0(I,J,K)*Vy(I,J,K)
      !End Do
      !End Do
      !End Do

!=======================================================================

       call TriSembdary3(Bd,PL, Ed,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)

      If (VisNonzero) Then

       call TriSembdary3(Pi00,Pi01, Pi02,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
       call TriSembdary3(Pi11,Pi22, Pi33,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
       call Sembdary3(Pi12, NX0,NY0,NZ0,NX,NY,NZ,
     &                         NXPhy0,NYPhy0,NXPhy,NYPhy)

      End If

       call Sembdary3(PPI, NX0,NY0,NZ0,NX,NY,NZ,
     &                         NXPhy0,NYPhy0,NXPhy,NYPhy)


C-----------------------------------------------------------------
      Call CoefASM2d(CofAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &               NX0,NY0,NZ0, NX,NY,NZ)
c------------------------------------------------------------------

      If (VisNonzero) Then

       DO 670 K=NZ0,NZ
       DO 670 J=NYPhy0,NYPhy
       DO 670 I=NXPhy0,NXPhy

        SDX1=(0.0 + Pi11(I+1,J,K) - Pi11(I-1,J,K)
     &   -Vx(I+1,J,K)*Pi01(I+1,J,K)
     &   +Vx(I-1,J,K)*Pi01(I-1,J,K))/(2.0*DX)
        SDY1=(Pi12(I,J+1,K)-Pi12(I,J-1,K)
     &   -Vy(I,J+1,K)*Pi01(I,J+1,K)
     &   +Vy(I,J-1,K)*Pi01(I,J-1,K))/(2.0*DY)
        ScT01(i,j,K)=(SDX1+SDY1)*Time  !ScT01

        SDX2=(Pi12(I+1,J,K)-Pi12(I-1,J,K)
     &   -Vx(I+1,J,K)*Pi02(I+1,J,K)
     &   +Vx(I-1,J,K)*Pi02(I-1,J,K))/(2.0*DX)
        SDY2=(0.0 + Pi22(I,J+1,K)-Pi22(I,J-1,K)
     &   -Vy(I,J+1,K)*Pi02(I,J+1,K)
     &   +Vy(I,J-1,K)*Pi02(I,J-1,K))/(2.0*DY)
        ScT02(i,j,K)=(SDX2+SDY2)*Time  !ScT02

        SDX0=(Pi01(I+1,J,K)-Pi01(I-1,J,K)+Vx(I+1,J,K)*(0.0
     & -Pi00(I+1,J,K))-Vx(I-1,J,K)*(0.0-Pi00(I-1,J,K)))/(2.0*DX)
        SDY0=(Pi02(I,J+1,K)-Pi02(I,J-1,K)+Vy(I,J+1,K)*(0.0
     & -Pi00(I,J+1,K))-Vy(I,J-1,K)*(0.0-Pi00(I,J-1,K)))/(2.0*DY)
        ScT00(i,j,K)=0.0+Time*(SDX0+SDY0)+0.0 !ScT00

 670   continue
      Else
        ScT01 = 0D0
        ScT02 = 0D0
        ScT00 = 0D0
      End If


      If (VisNonzero) Then

       call TriSembdary3(ScT00,ScT01,ScT02,
     &        NX0,NY0,NZ0, NX,NY,NZ,NXPhy0,NYPhy0, NXPhy,NYPhy)
      End If


       DO 690 K=NZ0,NZ
       DO 690 J=NYPhy0,NYPhy
       DO 690 I=NXPhy0,NXPhy

        ScT01(i,j,K)=ScT01(i,j,k)
     &    +Time*(PL(I+1,J,K)-PL(I-1,J,K)
     &    +PPi(I+1,J,K)-PPi(I-1,J,K))/(2.0*DX)!ScT01
        ScT02(i,j,K)=ScT02(i,j,K)
     &    +Time*(PL(I,J+1,K)-PL(I,J-1,K)
     &    +PPi(I,J+1,K)-PPi(I,J-1,K))/(2.0*DY)

        ScT00(i,j,K)=ScT00(i,j,K)+PL(i,j,K)+PPi(i,j,K)+Pi33(I,J,K)
     & +Time*(Vx(I+1,J,K)*PL(I+1,J,K)-Vx(I-1,J,K)*PL(I-1,J,K))/(2.0*DX)
     &+Time*(Vx(I+1,J,K)*PPi(I+1,J,K)-Vx(I-1,J,K)*PPi(I-1,J,K))/(2.0*DX)
     & +Time*(Vy(I,J+1,K)*PL(I,J+1,K)-Vy(I,J-1,K)*PL(I,J-1,K))/(2.0*DY)
     &+Time*(Vy(I,J+1,K)*PPi(I,J+1,K)-Vy(I,J-1,K)*PPi(I,J-1,K))/(2.0*DY)
 690   continue

        call TriSembdary3(ScT00,ScT01, ScT02,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)

       If(Nsm.gt.0) then
       do MM=1,NSm
       call ASM2D6(ScT00,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                         NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(ScT01,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(ScT02,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                         NX0,NY0,NZ0, NX,NY,NZ)
       end do
       end if

C---------------------------------------------------------------------
        call EntropyTemp3 (Ed,PL, Temp,Sd,
     &         NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)

        call ViscousCoefi8(Ed,Sd,PL,Temp,
     &  VCoefi,VCBeta,VRelaxT,etaTtp, VBulk, VRelaxT0,XiTtP,
     &  Time,DX,DY,DT,NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
C--------------------

      If (VisNonzero .or. VisBulkNonzero) Then
        call PiS4U5(PU0,PU1,PU2,PU3,U0,U1,U2,U3, DX,DY,DZ, DT,
     & DPc00,DPc01,DPc02,DPc33, DPc11,DPc22,DPc12, DDU0,DDU1,DDU2,
     & Temp,Temp0,  SiLoc,DLnT,  Time, NXPhy0,NYPhy0,NXPhy,NYPhy,
     & NX0,NX,NY0,NY,NZ0,NZ)
      Else ! Ver 1.6.19RC5: set to 0
        DPc00=0D0
        DPc01=0D0
        DPc02=0D0
        DPc11=0D0
        DPc12=0D0
        DPc22=0D0
        DPc33=0D0
      End If


      If (VisNonzero) Then

       DO 700 K=NZ0,NZ
       DO 700 J=NYPhy0,NYPhy
       DO 700 I=NXPhy0,NXPhy

        Tr0=Pi00(i,j,k)*DDU0(i,j,k)-Pi01(i,j,k)*DDU1(i,j,k)
     &                           -Pi02(i,j,k)*DDU2(i,j,k)   !additonal term to keep transversality of pi
        Tr1=Pi01(i,j,k)*DDU0(i,j,k)-Pi11(i,j,k)*DDU1(i,j,k)
     &                           -Pi12(i,j,k)*DDU2(i,j,k)   !additonal term to keep transversality of pi
        Tr2=Pi02(i,j,k)*DDU0(i,j,k)-Pi12(i,j,k)*DDU1(i,j,k)
     &                           -Pi22(i,j,k)*DDU2(i,j,k)   !additonal term to keep transversality of pi

        TPi00=2.0*Tr0
        TPi01=Tr1+Vx(i,j,k)*Tr0
        TPi02=Tr2+Vy(i,j,k)*Tr0
        TPi11=2.0*Vx(i,j,k)*Tr1
        TPi12=Vx(i,j,k)*Tr2+Vy(i,j,k)*Tr1
        TPi22=2.0*Vy(i,j,k)*Tr2
        TPi33=0.0

        PScT00(i,j,K)=(0.0-TPi00)*(-1.0)
        PScT01(i,j,K)=(0.0-TPi01)*(-1.0)
        PScT02(i,j,K)=(0.0-TPi02)*(-1.0)
        PScT33(i,j,K)=(0.0-TPi33)*(-1.0)

        PScT11(i,j,K)=(0.0-TPi11)*(-1.0)
        PScT12(i,j,K)=(0.0-TPi12)*(-1.0)
        PScT22(i,j,K)=(0.0-TPi22)*(-1.0)
 700   continue
       Else ! Ver 1.6.19RC5: set to 0
        PScT00=0D0
        PScT01=0D0
        PScT02=0D0
        PScT11=0D0
        PScT12=0D0
        PScT22=0D0
        PScT33=0D0
      End If
C--------------------------------------------------------------------------------


      If (VisNonzero) Then

       call TriSembdary3(PScT00,PScT01, PScT02,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
       call TriSembdary3(PScT11,PScT22, PScT33,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
       call Sembdary3(PScT12, NX0,NY0,NZ0,NX,NY,NZ,
     &                         NXPhy0,NYPhy0,NXPhy,NYPhy)
      End If

      DO 710 K=NZ0,NZ
      DO 710 J=NYPhy0,NYPhy
      DO 710 I=NXPhy0,NXPhy
        xx=DX*I
        yy=DY*J
        rr=sqrt(xx**2+yy**2)
        ff=1.0/(Dexp((rr-R0Bdry)/Aeps)+1.0)

#if DERIV_ORDER == 3
          PA=(Vx(I+1,J,K)-Vx(I-1,J,K))/(2*DX)
     &      +(Vy(I,J+1,K)-Vy(I,J-1,K))/(2*DY)
#elif DERIV_ORDER == 5
          PA=( Vx(I+1,J,K)*2.0d0/3.0d0-Vx(I-1,J,K)*2.0d0/3.0d0
     &       -Vx(I+2,J,K)/12.0d0+Vx(I-2,J,K)/12.0d0 )/DX
     &      +( Vy(I,J+1,K)*2.0d0/3.0d0-Vy(I,J-1,K)*2.0d0/3.0d0
     &       -Vy(I,J+2,K)/12.0d0+Vy(I,J-2,K)/12.0d0 )/DY
#else
#error "DERIV_ORDER must be 3 or 5"
#endif

        PT=VRelaxT(I,J,K)/U0(I,J,K)
        PS=2.0*VCoefi(I,J,K)

        PT0=VRelaxT0(I,J,K)/U0(I,J,K)
        PS0=VBulk(I,J,K)

        if (VisNonzero) then

          if(ViscousEqsType .eq. 1) then
          ! old version: shear tensor pressure terms from entropy generation
            D0Ln=(DLog(etaTtp0(I,J,K))-DLog(etaTtp(I,J,K)))/DT
            D1Ln=(DLog(etaTtp(I-1,J,K))-DLog(etaTtp(I+1,J,K)))/(2.0*DX)
            D2Ln=(DLog(etaTtp(I,J-1,K))-DLog(etaTtp(I,J+1,K)))/(2.0*DY)

            ADLnT=-0.5*(SiLoc(I,J,K)/U0(I,J,K)         !high precision version to reduce round-off error
     &           +(D0Ln+ Vx(I,J,K)*D1Ln +Vy(I,J,K)*D2Ln)*ff)

            PScT00(i,j,K)=PScT00(i,j,K) +( Pi00(I,J,K)*ADLnT+
     &         PA*Pi00(I,J,K)-PT*(Pi00(I,J,K)-PS*DPc00(i,j,K)))*(-1.0)
            PScT01(i,j,K)=PScT01(i,j,K) +( Pi01(I,J,K)*ADLnT+
     &         PA*Pi01(I,J,K)-PT*(Pi01(I,J,K)-PS*DPc01(i,j,K)))*(-1.0)
            PScT02(i,j,K)=PScT02(i,j,K) +( Pi02(I,J,K)*ADLnT+
     &         PA*Pi02(I,J,K)-PT*(Pi02(I,J,K)-PS*DPc02(i,j,K)))*(-1.0)
            PScT33(i,j,K)=PScT33(i,j,K) +( Pi33(I,J,K)*ADLnT+
     &         PA*Pi33(I,J,K)-PT*(Pi33(I,J,K)-PS*DPc33(i,j,K)))*(-1.0)
            PScT11(i,j,K)=PScT11(i,j,K) +( Pi11(I,J,K)*ADLnT+
     &         PA*Pi11(I,J,K)-PT*(Pi11(I,J,K)-PS*DPc11(i,j,K)))*(-1.0)
            PScT22(i,j,K)=PScT22(i,j,K) +( Pi22(I,J,K)*ADLnT+
     &         PA*Pi22(I,J,K)-PT*(Pi22(I,J,K)-PS*DPc22(i,j,K)))*(-1.0)
            PScT12(i,j,K)=PScT12(i,j,K) +( Pi12(I,J,K)*ADLnT+
     &         PA*Pi12(I,J,K)-PT*(Pi12(I,J,K)-PS*DPc12(i,j,K)))*(-1.0)

          else if(ViscousEqsType .eq. 2) then
          ! shear terms from 14-moments expansion
            p00 = Pi00(I,J,K)    ! some short hand
            p01 = Pi01(I,J,K)
            p02 = Pi02(I,J,K)
            p11 = Pi11(I,J,K)
            p12 = Pi12(I,J,K)
            p22 = Pi22(I,J,K)
            p33 = Pi33(I,J,K)

            s00 = DPc00(I,J,K)
            s01 = DPc01(I,J,K)
            s02 = DPc02(I,J,K)
            s11 = DPc11(I,J,K)
            s12 = DPc12(I,J,K)
            s22 = DPc22(I,J,K)
            s33 = DPc33(I,J,K)

            u0_temp = U0(I,J,K)
            vx_temp = Vx(I,J,K)
            vy_temp = Vy(I,J,K)

            theta_temp = SiLoc(I,J,K)
            ppi_temp = PPI(I,J,K)
            taupi = 1.D0/VRelaxT(I,J,K)

            s0d = s00 - vx_temp*s01 - vy_temp*s02
            s1d = s01 - vx_temp*s11 - vy_temp*s12
            s2d = s02 - vx_temp*s12 - vy_temp*s22

            ! calculate Tr(pi^2)
            TrPi2 = p00*p00 + p11*p11+  p22*p22 + p33*p33
     &        -2*p01*p01-2*p02*p02+2*p12*p12
            ! calculate Tr(pi*sigma)
            TrPiSigma = p00*s00 + p11*s11 + p22*s22 + p33*s33
     &           - 2.D0*p01*s01 - 2.D0*p02*s02 + 2*p12*s12

            call ViscousShearTransCoefs(Ed(I,J,K), PL(I,J,K),
     &         VRelaxT(I,J,K), deltaSpiSpi, lambdaSpiBPi,
     &         phi7, taupipi)

            ! pi source 0,0
            phi7Add = p00**2.0 - p01**2.0 - p02**2.0
     &        -(1.0-u0_temp**2.0)*TrPi2/3.D0
            taupipiAdd = p00*s00 - p01*s01 - p02*s02
     &        -u0_temp**2.0*(p00*s0d-p01*s1d-p02*s2d)
     &        -(1.D0-u0_temp**2.0)*TrPiSigma/3.D0
            PScT00(i,j,K)=PScT00(i,j,K) + (-PT)
     &        *(taupi*u0_temp*PA*p00-(p00-PS*s00)
     &         -deltaSpiSpi*p00*theta_temp+lambdaSpiBPi*ppi_temp*s00
     &         +phi7*phi7Add - taupipi*taupipiAdd)

            ! pi source 0,1
            phi7Add = p00*p01 - p01*p11 - p12*p02
     &        +u0_temp**2.0*vx_temp*TrPi2/3.D0
            taupipiAdd = 0.5*(p00*s01+p01*s00-p01*s11-p02*s12-p11*s01
     &         -p12*s02)
     &        - 0.5*u0_temp**2.0*vx_temp*(p00*s0d-p01*s1d-p02*s2d)
     &        - 0.5*u0_temp**2.0*(p01*s0d-p11*s1d-p12*s2d)
     &        + u0_temp**2.0*vx_temp*TrPiSigma/3.D0
            PScT01(i,j,k)=PScT01(i,j,k) + (-PT)
     &        *(taupi*u0_temp*PA*p01-(p01-PS*s01)
     &         -deltaSpiSpi*p01*theta_temp+lambdaSpiBPi*ppi_temp*s01
     &         +phi7*phi7Add - taupipi*taupipiAdd)

            ! pi source 0,2
            phi7Add = p00*p02 - p01*p12 - p02*p22
     &        +u0_temp**2.0*vy_temp*TrPi2/3.D0
            taupipiAdd = 0.5*(p00*s02+p02*s00-p01*s12-p02*s22-p12*s01
     &         -p22*s02)
     &        -0.5*u0_temp**2.0*vy_temp*(p00*s0d-p01*s1d-p02*s2d)
     &        -0.5*u0_temp**2.0*(p02*s0d-p12*s1d-p22*s2d)
     &        +u0_temp**2.0*vy_temp*TrPiSigma/3.D0
            PScT02(i,j,k)=PScT02(i,j,k) + (-PT)
     &        *(taupi*u0_temp*PA*p02-(p02-PS*s02)
     &         -deltaSpiSpi*p02*theta_temp+lambdaSpiBPi*ppi_temp*s02
     &         +phi7*phi7Add - taupipi*taupipiAdd)

            ! pi source 1,1
            phi7Add = p01**2.0-p11**2.0-p12**2.0
     &        +(1.0+u0_temp**2.0*vx_temp**2.0)*TrPi2/3.D0
            taupipiAdd = p01*s01-p11*s11-p12*s12
     &        -u0_temp**2.0*vx_temp*(p01*s0d-p11*s1d-p12*s2d)
     &        +(1.0+u0_temp**2*vx_temp**2)*TrPiSigma/3.D0
            PScT11(i,j,k)=PScT11(i,j,k) + (-PT)
     &        *(taupi*u0_temp*PA*p11-(p11-PS*s11)
     &         -deltaSpiSpi*p11*theta_temp+lambdaSpiBPi*ppi_temp*s11
     &         +phi7*phi7Add - taupipi*taupipiAdd)

            ! pi source 1,2
            phi7Add = p01*p02-p11*p12-p12*p22
     &        +u0_temp**2.0*vx_temp*vy_temp*TrPi2/3.D0
            taupipiAdd = 0.5*(p01*s02-p11*s12-p12*s22+p02*s01-p12*s11
     &         -p22*s12)
     &        -0.5*u0_temp**2.0*vy_temp*(p01*s0d-p11*s1d-p12*s2d)
     &        -0.5*u0_temp**2.0*vx_temp*(p02*s0d-p12*s1d-p22*s2d)
     &        +u0_temp**2.0*vx_temp*vy_temp*TrPiSigma/3.D0
            PScT12(i,j,k)=PScT12(i,j,k) + (-PT)
     &        *(taupi*u0_temp*PA*p12-(p12-PS*s12)
     &         -deltaSpiSpi*p12*theta_temp+lambdaSpiBPi*ppi_temp*s12
     &         +phi7*phi7Add - taupipi*taupipiAdd)

            ! pi source 2,2
            phi7Add = p02**2.0-p12**2.0-p22**2.0
     &        +(1.0+u0_temp**2.0*vy_temp**2.0)*TrPi2/3.D0
            taupipiAdd = p02*s02-p12*s12-p22*s22
     &        -u0_temp**2.0*vy_temp*(p02*s0d-p12*s1d-p22*s2d)
     &        +(1.0+u0_temp**2.0*vy_temp**2.0)*TrPiSigma/3.0
            PScT22(i,j,k)=PScT22(i,j,k) + (-PT)
     &        *(taupi*u0_temp*PA*p22-(p22-PS*s22)
     &         -deltaSpiSpi*p22*theta_temp+lambdaSpiBPi*ppi_temp*s22
     &         +phi7*phi7Add - taupipi*taupipiAdd)

            ! pi source 3,3
            phi7Add = -p33**2.0+TrPi2/3.D0
            taupipiAdd = -p33*s33+TrPiSigma/3.D0
            PScT33(i,j,k)=PScT33(i,j,k) + (-PT)
     &        *(taupi*u0_temp*PA*p33-(p33-PS*s33)
     &         -deltaSpiSpi*p33*theta_temp+lambdaSpiBPi*ppi_temp*s33
     &         +phi7*phi7Add - taupipi*taupipiAdd)

          else
            write(*, *) "No such viscous equation type:",ViscousEqsType
            call exit(1)
          end if
        else
          ADLnT = 0.D0
        end if

        if (VisBulkNonzero) then

          if(ViscousEqsType .eq. 1) then
          ! old version: bulk pressure terms from entropy generation
            DB0Ln=(DLog(XiTtp0(I,J,K))-DLog(XiTtp(I,J,K)))/DT
            DB1Ln=(DLog(XiTtp(I-1,J,K))-DLog(XiTtp(I+1,J,K)))/(2.0*DX)
            DB2Ln=(DLog(XiTtp(I,J-1,K))-DLog(XiTtp(I,J+1,K)))/(2.0*DX)

            Badd=-0.5*(SiLoc(I,J,K)/U0(I,J,K)    !high precision version to reduce round-off error
     &          +(DB0Ln +Vx(I,J,K)*DB1Ln +Vy(I,J,K)*DB2Ln)*ff)

          elseif(ViscousEqsType .eq. 2) then
          ! Bulk pressure terms from 14-moments expansion
            call ViscousBulkTransCoefs(ED(i,j,k)*Hbarc, VRelaxT0(i,j,k),
     &          deltaBPiBPi, lambdaBPiSpi) ! calculate transport coefficients

            piSigma = Pi00(i,j,k)*DPc00(i,j,k)+
     &           Pi11(i,j,k)*DPc11(i,j,k) + Pi22(i,j,k)*DPc22(i,j,k)
     &           + Pi33(i,j,k)*DPc33(i,j,k)
     &           - 2.D0*Pi01(i,j,k)*DPc01(i,j,k)
     &           - 2.D0*Pi02(i,j,k)*DPc02(i,j,k)
     &           + 2.D0*Pi12(i,j,k)*DPc12(i,j,k) !Tr(pi*sigma)

            Badd = (-deltaBPiBPi*SiLoc(i,j,k)+sign(1.D0,PPI(i,j,k))*
     &       lambdaBPiSpi/DMax1(abs(PPI(i,j,k)), 1e-30)*piSigma)
     &       /U0(i,j,k)*VRelaxT0(i,j,k)
          else
            write(*, *) "No such viscous equation type:", ViscousEqsType
            call exit(1)
          endif

        else
          Badd=0.0
        end if

        PScTSum = PScT00(i,j,k) + PScT01(i,j,k) + PScT02(i,j,k)
     &        + PScT11(i,j,k) + PScT12(i,j,k) + PScT22(i,j,k)
     &        + PScT33(i,j,k)

        If (TIME > T0) Then
        If (isnan(PScTSum)) Then
          Print *, "Invalid PScT terms!"
          Print *, "i,j,k=", i,j,k
          Print *, "PScT00=", PScT00(i,j,K)
          Print *, "PScT01=", PScT01(i,j,K)
          Print *, "PScT02=", PScT02(i,j,K)
          Print *, "PScT11=", PScT11(i,j,K)
          Print *, "PScT12=", PScT12(i,j,K)
          Print *, "PScT22=", PScT22(i,j,K)
          Print *, "PScT33=", PScT33(i,j,K)
          Print *, "Pi00,ADLnT=", Pi00(I,J,K), ADLnT
          Print *, "PA", PA
          Print *, "Vx(I+1,J,K)", Vx(I+1,J,K),"Vx(I-1,J,K)",Vx(I-1,J,K)
          Print *, "Vy(I,J+1,K)", Vy(I,J+1,K),"Vy(I,J-1,K)",Vy(I,J-1,K)
          Print *, "PT=", PT, "Pi00=", Pi00(I,J,K)
          Print *, "PS=", PS, "DPc00=", DPc00(i,j,K)
          Print *, "SiLoc=", SiLoc(I,J,K), "U0=", U0(I,J,K)   !high precision version to reduce round-off error
          Print *, "D0Ln=", D0Ln, "Vx=", Vx(I,J,K), "Vy=", Vy(I,J,K)
          Print *, "D1Ln=", D1Ln, "D2Ln=", D2Ln, "ff=", ff
          Print *, "etaTtp(I,J-1,K)=", etaTtp(I,J-1,K)
          Print *, "etaTtp(I,J+1,K)=", etaTtp(I,J+1,K)

          call exit(1)
        EndIf
        EndIf

        PISc(i,j,K)=0.0 +( PPI(I,J,K)*BAdd+
     &      PA*PPI(I,J,K)-PT0*(PPI(I,J,K)+PS0*SiLoc(i,j,K)))*(-1.0)

 710   continue

       call TriSembdary3(PScT00,PScT01, PScT02,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
       call TriSembdary3(PScT11,PScT22, PScT33,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
       call Sembdary3(PScT12, NX0,NY0,NZ0,NX,NY,NZ,
     &                         NXPhy0,NYPhy0,NXPhy,NYPhy)
       call Sembdary3(PISc, NX0,NY0,NZ0,NX,NY,NZ,
     &                         NXPhy0,NYPhy0,NXPhy,NYPhy)

       If(Nsm.gt.0) then
       do MM=1,NSm
       call ASM2D6(PScT00,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                         NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(PScT01,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(PScT02,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                         NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(PScT11,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(PScT12,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(PScT22,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                         NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(PScT33,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(PISc,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)

       end do
       end if

C*********************************************************************
       call TriSembdary3(T00, T01, T02,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)

       call TriSembdary3(TT00, TT01, TT02,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)

C-------------------------------------------------------------------------------------


       If(Nsm.gt.0) then
        do MM=1,Nsm
        Call ASM2D6bsm(TT00,CofAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
        Call ASM2D6bsm(TT01,CofAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
        Call ASM2D6bsm(TT02,CofAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
        end do
       end if

         If(Time.lt.T0+dT) then
            StotalSv=0.0
            StotalBv=0.0
            Stotal=0.0
               K=1
            do 111 J=NYPhy0,NYPhy
            do 111 I=NXPhy0,NXPhy
             if(Ed(I-1,J-1,K).ge. Edec1/HbarC) Stotal=Stotal+Sd(i,j,K)
 111        continue
            Stotal=Stotal*time*dx*dy
            StotalSv=Stotal
            StotalBv=Stotal
         else
              DsP1=0.0
              DsP2=0.0
              K=1
            do 222 J=NYPhy0,NYPhy
            do 222 I=NXPhy0,NXPhy
            if(Ed(I-1,J-1,K).ge. Edec1/HbarC)  then
             if (VisNonzero) then
              Dsp1=Dsp1+(Pi00(I,J,K)**2 +Pi11(I,J,K)**2 +Pi22(I,J,K)**2
     &             +Pi33(I,J,K)**2 +2.0*Pi01(I,J,K)**2
     &             +2.0*Pi02(I,J,K)**2 +2.0*Pi12(I,J,K)**2)/
     &                  (2.0*VCoefi(I,J,K)*Temp(I,J,K))
             end if
             if (VisBulkNonzero .and. VBulk(I,J,K).ge.0.000001)then
             Dsp2=Dsp2+(PPI(I,J,K)*PPI(I,J,K))/
     &                    (VBulk(I,J,K)*Temp(I,J,K))
             end if
            end if
 222       continue
             StotalSv=StotalSv+ Dsp1*time*dx*dy*dt
             StotalBv=StotalBv+ Dsp2*time*dx*dy*dt

             Stotal=Stotal+(Dsp1+Dsp2)*time*dx*dy*dt
         end if

C            Print *, 'time',time,'Stotal', Stotal,StotalSv,StotalBv


         If(Time.le.T0+dT) SxyT=0.0d0

        call anisotrop10 (Vx,Vy,U0,U1,U2, Ed,PL,Sd, PPI,
     &  Pi00,Pi01,Pi02, Pi33, Pi11,Pi12,Pi22, ScT00,ScT01,ScT02,
     &  DX,DY,DZ,DT,Time,EpsX,EpsP,TEpsP, Vaver,VavX,VavY, STotal22,
     &  AE,ATemp, APL, APLdx,APLdy, AScT00, AScT01, AScT02,
     &  AP11dx,AP22dy,  AP33,AP00, AP01,AP02, AP11, Ap12, AP22,SxyT,
     &  PaScT00,PaScT01,PaScT02, difScT0102,difPLxy, difPixy, difVxVy,
     &  NXPhy0,NYPhy0, NXPhy,NYPhy, NX0,NX,NY0,NY,NZ0,NZ)  !PNEW NNEW something related to root finding



      AMV=Hbarc*1000.0 !fm-1 change to MEV
      GV=Hbarc !fm-1 change to GEV

      Do J=0,NXPhy,NXPhy+1
      Do I=NYPhy0,NYPhy,10
        PPI_NS = (-1.0)*VBulk(I,J,NZ0)*SiLoc(I,J,NZ0) ! Navier-Stokes limit
      enddo
      enddo

!     Checking codes deleted
!     See previous versions
!     Search "goto 888"

      Return
      End








C###########################################################################################

      subroutine anisotrop10 (Vx,Vy,U0,U1,U2, Ed,PL,Sd, PPI,
     &  Pi00,Pi01,Pi02, Pi33, Pi11,Pi12,Pi22, ScT00,ScT01,ScT02,
     &  DX,DY,DZ,DT,Time,EpsX,EpsP,TEpsP, Vaver,VavX,VavY, STotal,
     &  AE,ATemp, APL, APLdx,APLdy, AScT00, AScT01, AScT02,
     &  AP11dx,AP22dy,  AP33,AP00, AP01,AP02, AP11, Ap12, AP22,SxyT,
     &  PaScT00,PaScT01,PaScT02, difScT0102,difPLxy, difPixy, difVxVy,
     &  NXPhy0,NYPhy0, NXPhy,NYPhy, NX0,NX,NY0,NY,NZ0,NZ)  !PNEW NNEW  related to root finding

      Implicit Double Precision (A-H, O-Z)
      Dimension Vx(NX0:NX, NY0:NY, NZ0:NZ) !fluid velocity
      Dimension Vy(NX0:NX, NY0:NY, NZ0:NZ) !fluid velocity

      Dimension U0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Dimension U1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Dimension U2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

      Dimension Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Dimension PL(NX0:NX, NY0:NY, NZ0:NZ) !local pressure
      Dimension Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density
      Dimension Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature

      Dimension Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
       Dimension PPI(NX0:NX, NY0:NY, NZ0:NZ)    !Bulk pressure


      Dimension ScT00(NX0:NX, NY0:NY, NZ0:NZ) !Source Term  ScT=ScX+ScY+ScZ
      Dimension ScT01(NX0:NX, NY0:NY, NZ0:NZ) !Source Term  ScT=ScX+ScY+ScZ
      Dimension ScT02(NX0:NX, NY0:NY, NZ0:NZ) !Source Term  ScT=ScX+ScY+ScZ

      common/Edec/Edec    !decoupling temperature

       EpsX1=0.0  !relate spacial ellipticity
       EpsX2=0.0  !relate spacial ellipticity
       EpsP1=0.0  !relate momentum ellipticity
       EpsP2=0.0  !relate momentum ellipticity
       TEpsP1=0.0  !relate momentum ellipticity
       TEpsP2=0.0  !relate momentum ellipticity
       STotal=0.0
       Vaver1=0.0    !average velocity
       Vaver2=0.0
       VavX1=0.0
       VavY1=0.0

      DO 100 K=NZ0,NZ
      DO 100 J=-NYPhy,NYPhy
      DO 100 I=-NXPhy,NXPhy
        !If (Ed(I,J,K) < EDec) Cycle
        gamma=1.0D0/sqrt(1D0-Vx(I,J,NZ0)**2-Vy(I,J,NZ0)**2)
        !gamma=1.0D0/max(sqrt(abs(1D0-Vx(I,J,K)**2-Vy(I,J,K)**2)),1D-15)

        xx=I*DX
        yy=J*DY

        x2=xx*xx
        y2=yy*yy

        EpsX1=EpsX1+(y2-x2)*Ed(I,J,K)*gamma   !*
        EpsX2=EpsX2+(y2+x2)*Ed(I,J,K)*gamma   !*

        !angle = 2*atan2(YY,XX)
        !RR = sqrt(XX*XX+YY*YY)
        !
        !gammaE=1.0/sqrt(1-Vx(I,J,K)**2-Vy(I,J,K)**2)*Ed(I,J,K)

        ! eccentricity, defined using r^2 as weight function:
        !Weight=Weight+RR*RR*gammaE ! Note that this weight is repeatedly calculated. But since it is much cleaner written this way and it is very fast...
        !EpsX1=EpsX1+RR*RR*cos(angle)*gammaE
        !EpsX2=EpsX2+RR*RR*gammaE

        ep=Ed(I,J,K)+PL(I,J,K)
        TXX=ep*U1(I,J,K)*U1(I,J,K)+PL(I,J,K)
        TYY=ep*U2(I,J,K)*U2(I,J,K)+PL(I,J,K)

        EpsP1=EpsP1+(TXX-TYY)
        EpsP2=EpsP2+(TXX+TYY)

        TTXX=TXX+Pi11(I,J,K)+PPI(I,J,K)
     &   +PPI(I,J,K)*U1(I,J,K)*U1(I,J,K)
        TTYY=TYY+Pi22(I,J,K)+PPI(I,J,K)
     &   +PPI(I,J,K)*U2(I,J,K)*U2(I,J,K)

        TEpsP1=TEpsP1+(TTXX-TTYY)
        TEpsP2=TEpsP2+(TTXX+TTYY)

        Vr=sqrt(Vx(I,J,K)**2+Vy(I,J,K)**2)
        Vaver1=Vaver1+Ed(I,J,K)*Vr*gamma   !*
        VavX1=VavX1+Ed(I,J,K)*Vx(I,J,K)*gamma   !*
        VavY1=VavY1+Ed(I,J,K)*Vy(I,J,K)*gamma   !*
        Vaver2=Vaver2+Ed(I,J,K)*gamma   !*

        STotal=STotal+Sd(I,J,K)*U0(I,J,K)*Time*dx*dy
 100  continue

        EpsX=EpsX1/EpsX2  !spacial ellipticity
        EpsP=EpsP1/EpsP2  !Momentum  ellipticity
        TEpsP=TEpsP1/TEpsP2 ! total Momentum  ellipticity

        Vaver=Vaver1/Vaver2
        VavX=VaVX1/Vaver2
        VavY=VaVY1/Vaver2

      K=NZ0
      Sx=0.0d0
      I=NXPhy
      DO 110 J=NYPhy0,NYPhy
       Sx=Sx+1.0*Time*dy*dT*U1(I,J,K)*Sd(I,J,K)
 110  continue


      K=NZ0
      Sy=0.0d0
      J=NYPhy
      DO 120 I=NXPhy0,NXPhy
       Sy=Sy+1.0*Time*dx*dT*U2(I,J,K)*Sd(I,J,K)
 120  continue


         SxyT=SxyT+Sx+Sy
         STotal=(STotal+SxyT)*4.0
C        STotal=STotal*Time*dx*dy*DT


       AE=0.0
       ATemp=0.0
       AP00=0.0
       AP01=0.0
       AP02=0.0
       AP33=0.0
       AP11=0.0
       AP12=0.0
       AP22=0.0

       APLdx=0.0
       APLdy=0.0
       APL=0.0

       PaScT00=0.0
       PaSCT01=0.0
       PaSCT02=0.0

       AScT00=0.0
       ASCT01=0.0
       ASCT02=0.0


       difScT0102=0.0
       difPLxy=0.0
       difPixy=0.0
       difVxVy=0.0
C----------------------------------------------------------------
      DO 200 K=NZ0,NZ
      DO 200 J=0,NYPhy
      DO 200 I=0,NXPhy
      AE=AE+Ed(I,J,K)*gamma  !*
      ppe=PL(I,J,K)+Ed(I,J,K)
      ATemp=ATemp+Ed(I,J,K)*Temp(I,J,K)*gamma  !*

      AP00=AP00+Ed(I,J,K)*Pi00(I,J,K)*gamma/ppe  !*
      AP01=AP01+Ed(I,J,K)*Pi01(I,J,K)*gamma/ppe  !*
      AP02=AP02+Ed(I,J,K)*Pi02(I,J,K)*gamma/ppe  !*
      AP11=AP11+Ed(I,J,K)*Pi11(I,J,K)*gamma/ppe  !*
      AP12=AP12+Ed(I,J,K)*Pi12(I,J,K)*gamma/ppe  !*
      AP22=AP22+Ed(I,J,K)*Pi22(I,J,K)*gamma/ppe  !*
      AP33=AP33+Ed(I,J,K)*Pi33(I,J,K)*gamma/ppe  !*

      APL=APL+PL(I,J,K)*Ed(I,J,K)*gamma  !*
     & +Time*Ed(I,J,K)*gamma  !*
     & *((PL(I+1,J,K)*Vx(I+1,J,K) -PL(I-1,J,K)*Vx(I-1,J,K))/(2.0*DX)
     &    +(PL(I,J+1,K)*Vy(I,J+1,K)-PL(I,J-1,K)*Vy(I,J-1,K))/(2.0*DY))

      APLdx=APLdx+(PL(I+1,J,K)-PL(I-1,J,K))*Ed(I,J,K)/(2.0*DX)*gamma  !*
      APLdy=APLdy+(PL(I,J+1,K)-PL(I,J-1,K))*Ed(I,J,K)/(2.0*DY)*gamma  !*

      AP11dx=AP11dx+(Pi11(I+1,J,K)-Pi11(I-1,J,K))
     &                                  *Ed(I,J,K)/(2.0*DX) *gamma  !*
      AP22dy=AP22dy+(Pi22(I,J+1,K)-Pi22(I,J-1,K))
     &                                  *Ed(I,J,K)/(2.0*DY) *gamma  !*

      AScT00=ASCT00+ScT00(I,J,K)*Ed(I,J,K)*gamma  !*
      AScT01=ASCT01+ScT01(I,J,K)*Ed(I,J,K)*gamma  !*
      AScT02=ASCT02+ScT02(I,J,K)*Ed(I,J,K)*gamma  !*

      difScT0102=difScT0102+(ScT01(I,J,K)-ScT02(I,J,K))*Ed(I,J,K)*gamma  !*
      difPLxy=difPLxy
     &     +Time*(PL(I+1,J,K)-PL(I-1,J,K))*Ed(I,J,K)/(2.0*DX)*gamma  !*
     &     -Time*(PL(I,J+1,K)-PL(I,J-1,K))*Ed(I,J,K)/(2.0*DY)*gamma  !*
      difPixy=difPixy
     &      +Time*(PL(I+1,J,K)-PL(I-1,J,K))*Ed(I,J,K)/(2.0*DX)*gamma  !*
     &      +Time*(Pi11(I+1,J,K)-Pi11(I-1,J,K))*Ed(I,J,K)/(2.0*DX)*gamma  !*
     &      -Time*(PL(I,J+1,K)-PL(I,J-1,K))*Ed(I,J,K)/(2.0*DY)*gamma  !*
     &      -Time*(Pi22(I,J+1,K)-Pi22(I,J-1,K))*Ed(I,J,K)/(2.0*DY)*gamma  !*

      difVxVy=difVxVy+(Vx(I,J,K)-Vy(I,J,K))*Ed(I,J,K) *gamma  !*

      PaScT00=PaScT00+PL(I,J,K)*Ed(I,J,K)*gamma  !*
     &  +Time*Ed(I,J,K)*gamma  !*
     &  *((PL(I+1,J,K)*Vx(I+1,J,K) -PL(I-1,J,K)*Vx(I-1,J,K))/(2.0*DX)
     &    +(PL(I,J+1,K)*Vy(I,J+1,K)-PL(I,J-1,K)*Vy(I,J-1,K))/(2.0*DY))
     & +Ed(I,J,K)*Pi33(I,J,K)*gamma  !*

      PaSCT01=PaScT01
     &   +(PL(I+1,J,K)-PL(I-1,J,K))*Ed(I,J,K)/(2.0*DX)*gamma  !*
     &   +(Pi11(I+1,J,K)-Pi11(I-1,J,K))*Ed(I,J,K)/(2.0*DX)*gamma  !*
      PaSCT02=PaScT02
     &    +(PL(I,J+1,K)-PL(I,J-1,K))*Ed(I,J,K)/(2.0*DY)*gamma  !*
     &    +(Pi22(I,J+1,K)-Pi22(I,J-1,K))*Ed(I,J,K)/(2.0*DY)*gamma  !*
 200  continue

      ATemp=ATemp/AE

      AP00=AP00/AE
      AP01=AP01/AE
      AP02=AP02/AE
      AP11=AP11/AE
      AP12=AP12/AE
      AP22=AP22/AE
      AP33=AP33/AE

      APL=APL/AE
      APLdx=Time*APLdx/AE
      APLdy=Time*APLdy/AE
      AP11dx=Time*AP11dx/AE
      AP22dy=Time*AP22dy/AE

      AScT00=AScT00/AE
      AScT01=AScT01/AE
      AScT02=AScT02/AE

      PaScT00=PaScT00/AE
      PaScT01=Time*PaScT01/AE
      PaScT02=Time*PaScT02/AE

      difScT0102=difScT0102/AE
      difPLxy=difPLxy/AE
      difPixy=difPixy/AE
      difVxVy=difVxVy/AE

      return
      end

!======================================================================
      Double Precision Function findEdHook(ee)
      ! Root finding by iterating quadratic equation formula for the equation
      ! (M0+cs^2*e+Pi)(M0-e)=M^2

      Implicit None

      Double Precision PEOSL7

      Double Precision ee ! energy density

      Double Precision :: RSDM0, RSDM, RSPPI
      Common /findEdHookData/ RSDM0, RSDM, RSPPI ! M0, M, Pi (see 0510014)

      Double Precision cs2 ! energy density from previous iteration, p/e, pressure
      ! Note that cs2 is not dp/de, but rather p/e!!!

      Double Precision, Parameter :: zero=1e-30 ! a small number
      double precision, parameter :: HbarC = M_HBARC

      Double Precision A,B ! intermedia step variables

      cs2=PEOSL7(ee*Hbarc)/dmax1(abs(ee),zero)/Hbarc ! check dimension?
      A=RSDM0*(1-cs2)+RSPPI
      B=RSDM0*(RSDM0+RSPPI)-RSDM*RSDM
      findEdHook = ee-2*B/dmax1((sqrt(A*A+4*cs2*B)+A), zero)

      End Function
!----------------------------------------------------------------------

!======================================================================
      Double Precision Function findvHook(v)
      ! Root finding by iterating quadratic equation formula for the equation
      ! (M0+cstilde^2*e+Pi)v=M

      Implicit None

      Double Precision PEOSL7
      Double Precision v

      Double Precision :: RSDM0, RSDM, RSPPI, RSee
      Common /findEdHookData/ RSDM0, RSDM, RSPPI ! M0, M, Pi (see 0510014)

      Double Precision cstilde2 ! energy density from previous iteration, p/e, pressure
      ! Note that cstilde2 is NOT dp/de, but rather p/e!!!

      Double Precision, Parameter :: zero=1e-30 ! a small number
      double precision, parameter :: HbarC = M_HBARC

      Double Precision A ! temporary variables

      RSee = RSDM0 - v*RSDM
      cstilde2=PEOSL7(RSee*Hbarc)/dmax1(abs(RSee),zero)/Hbarc
      A=RSDM0*(1+cstilde2)+RSPPI

      findvHook = v - (2*RSDM)/(A + sqrt(dmax1(A*A
     &                  - 4*cstilde2*RSDM*RSDM, 1D-30)) + 1D-30)
      if(isnan(findvHook)) then
        print*, cstilde2, A, v
        print*, A*A - 4*cstilde2*RSDM*RSDM
        print*, A + sqrt(A*A - 4*cstilde2*RSDM*RSDM)
        print*, RSDM0, RSDM, RSPPI, RSee
        call exit(1)
      endif
      End Function
!----------------------------------------------------------------------

!======================================================================
      Double Precision Function findU0Hook(U0)
      ! Root finding by iterating quadratic equation formula for the equation
      Implicit None

      Double Precision PEOSL7
      Double Precision v, U0

      Double Precision :: RSDM0, RSDM, RSPPI, RSee
      Common /findEdHookData/ RSDM0, RSDM, RSPPI ! M0, M, Pi (see 0510014)

      Double Precision cstilde2 ! energy density from previous iteration, p/e, pressure
      ! Note that cstilde2 is NOT dp/de, but rather p/e!!!

      Double Precision, Parameter :: zero=1e-30 ! a small number
      double precision, parameter :: HbarC = M_HBARC

      Double Precision A ! temporary variables

      v = sqrt(1 - 1/max(U0, 1d0)**2)
      RSee = RSDM0 - v*RSDM
      cstilde2=PEOSL7(RSee*Hbarc)/dmax1(abs(RSee),zero)/Hbarc
      A=RSDM0*(1+cstilde2)+RSPPI
      v = (2*RSDM)/(A+sqrt(dmax1(A*A-4*cstilde2*RSDM*RSDM, 1D-30))
     &              + 1D-30)

      findU0Hook = U0 - 1/sqrt(1-v*v)
      End Function
!----------------------------------------------------------------------


      Subroutine checkPiAll(failed, II, JJ, Time, Vx, Vy, Ed, PL,
     &  NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22,
     &  Pi00Regulated,Pi01Regulated,Pi02Regulated,Pi33Regulated,
     &  Pi11Regulated,Pi12Regulated,Pi22Regulated)
      ! Check tracelessness + transversality + positiveness of Tr(pi^2)

      Implicit None

      Integer NXPhy0,NXPhy,NYPhy0,NYPhy,NX0,NX,NY0,NY,NZ0,NZ
      Integer failed, II, JJ ! (II,JJ): where it fails
      Double Precision Time

      Double Precision Ed(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision PL(NX0:NX,NY0:NY,NZ0:NZ)

      Double Precision Vx(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Vy(NX0:NX,NY0:NY,NZ0:NZ)

      Double Precision Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      Double Precision Pi00Regulated(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi01Regulated(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi02Regulated(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi33Regulated(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi11Regulated(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi12Regulated(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi22Regulated(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      Integer I,J,K
      Double Precision trace_pi, trans, TrPi2

      Double Precision :: Tideal_scale, pi_scale
      Double Precision :: absNumericalzero = 1D-2
      Double Precision :: relNumericalzero = 1D1  !Xsi_0 in Zhi's thesis

      Double Precision maxPiRatio
      Double Precision gamma_perp
      Common /maxPiRatio/ maxPiRatio

      failed = 0

      Do K=NZ0,NZ
      Do J=NYPhy0-3,NYPhy+3
      Do I=NXPhy0-3,NXPhy+3

        Tideal_scale = sqrt(Ed(I,J,K)**2 + 3*PL(I,J,K)**2)

        TrPi2 = Pi00(I,J,K)**2+Pi11(I,J,K)**2+Pi22(I,J,K)**2
     &          +Pi33(I,J,K)**2
     &          -2*Pi01(I,J,K)**2-2*Pi02(I,J,K)**2+2*Pi12(I,J,K)**2

        pi_scale = sqrt(abs(TrPi2))
        !Positivity of Tr(pi^2)
        if(TrPi2 < -relNumericalzero*pi_scale) then
#ifdef OUTPUT_VIOLATION
            Print*, "Time=", Time
            Print*, "Trace Pi^2 is negative!"
            Print*, "I,J=", I,J
            Print*, "Trace pi^2=", TrPi2
            Print*, "Ed,PL=", Ed(I,J,K),PL(I,J,K)
            Print*, "Before: Pi00,Pi11,Pi22,Pi33*T^2,Pi01,Pi02,Pi12=",
     &        Pi00Regulated(I,J,K),Pi11Regulated(I,J,K),
     &        Pi22Regulated(I,J,K),Pi33Regulated(I,J,K),
     &        Pi01Regulated(I,J,K),Pi02Regulated(I,J,K),
     &        Pi12Regulated(I,J,K)
            Print*, "Before: Trace pi^2=",
     &        Pi00Regulated(I,J,K)**2+Pi11Regulated(I,J,K)**2
     &        +Pi22Regulated(I,J,K)**2+Pi33Regulated(I,J,K)**2
     &        -2*Pi01Regulated(I,J,K)**2-2*Pi02Regulated(I,J,K)**2
     &        +2*Pi12Regulated(I,J,K)**2
#endif
          II = I
          JJ = J
          failed = 1
          return
        endif

        !Largeness of pi tensor Tr(pi^2)
        If (pi_scale > max(maxPiRatio*Tideal_scale,
     &      absNumericalzero)) Then
#ifdef OUTPUT_VIOLATION
          Print*, "Time=", Time
          Print*, "pi is too large!"
          Print*, "I,J=", I,J
          Print*, "Trace pi^2=", TrPi2
          Print*, "Ed,PL=", Ed(I,J,K),PL(I,J,K)
          Print*, "Before: Pi00,Pi11,Pi22,Pi33*T^2,Pi01,Pi02,Pi12=",
     &      Pi00Regulated(I,J,K),Pi11Regulated(I,J,K),
     &      Pi22Regulated(I,J,K),Pi33Regulated(I,J,K),
     &      Pi01Regulated(I,J,K),Pi02Regulated(I,J,K),
     &      Pi12Regulated(I,J,K)
          Print*, "Before: Trace pi^2=",
     &      Pi00Regulated(I,J,K)**2+Pi11Regulated(I,J,K)**2
     &      +Pi22Regulated(I,J,K)**2+Pi33Regulated(I,J,K)**2
     &      -2*Pi01Regulated(I,J,K)**2-2*Pi02Regulated(I,J,K)**2
     &      +2*Pi12Regulated(I,J,K)**2
#endif
          II = I
          JJ = J
          failed = 1
          return
        End If

        !trace of pi tensor
        trace_pi = Pi00(I,J,K)-Pi11(I,J,K)-Pi22(I,J,K)-Pi33(I,J,K)
        If(abs(trace_pi) > max(relNumericalzero*pi_scale,
     &     absNumericalzero) ) Then
#ifdef OUTPUT_VIOLATION
          Print*, "Time=", Time
          Print*, "Pi should be traceless!"
          Print*, "I,J=", I,J
          Print*, "Pi00-Pi11-Pi22-Time**2*Pi33=", trace_pi
          Print*, "Pi00,Pi11,Pi22,Time**2*Pi33=",
     &      Pi00(I,J,K),Pi11(I,J,K),Pi22(I,J,K),Pi33(I,J,K)
          Print*, "Ed,PL=", Ed(I,J,K),PL(I,J,K)
          Print*, "Before: Pi00,Pi11,Pi22,Pi33*T^2,Pi01,Pi02,Pi12=",
     &      Pi00Regulated(I,J,K),Pi11Regulated(I,J,K),
     &      Pi22Regulated(I,J,K),Pi33Regulated(I,J,K),
     &      Pi01Regulated(I,J,K),Pi02Regulated(I,J,K),
     &      Pi12Regulated(I,J,K)
          Print*, "Before: Pi00-Pi11-Pi22-Time**2*Pi33=",
     &      Pi00Regulated(I,J,K)-Pi11Regulated(I,J,K)
     &      -Pi22Regulated(I,J,K)-Pi33Regulated(I,J,K)
#endif
          II = I
          JJ = J
          failed = 1
          return
        End If

        !transversality of pi tensor  u_mu pi^{mu,x} = 0
        gamma_perp = sqrt(1./(1. - Vx(I,J,K)**2 - Vy(I,J,K)**2 + 1D-30))
        trans = gamma_perp*(Pi01(I,J,K) - Vx(I,J,K)*Pi11(I,J,K)
     &                      - Vy(I,J,K)*Pi12(I,J,K))
        If (abs(trans) > max(relNumericalzero*pi_scale,
     &      absNumericalzero)) Then
#ifdef OUTPUT_VIOLATION
          Print*, "Time=", Time
          Print*, "Tranversality violated!"
          Print*, "I,J=", I,J
          Print*, "Pi01-Vx*Pi11-Vy*Pi12=", trans
          Print*, "Pi01,Pi11,Pi12=",
     &      Pi01(I,J,K),Pi11(I,J,K),Pi12(I,J,K)
          Print*, "Vx,Vy=",Vx(I,J,K),Vy(I,J,K)
          Print*, "Ed,PL=", Ed(I,J,K),PL(I,J,K)
          Print*, "Before: Pi00,Pi11,Pi22,Pi33*T^2,Pi01,Pi02,Pi12=",
     &      Pi00Regulated(I,J,K),Pi11Regulated(I,J,K),
     &      Pi22Regulated(I,J,K),Pi33Regulated(I,J,K),
     &      Pi01Regulated(I,J,K),Pi02Regulated(I,J,K),
     &      Pi12Regulated(I,J,K)
          Print*, "Before: Pi01-Vx*Pi11-Vy*Pi12=",
     &      Pi01Regulated(I,J,K)-Vx(I,J,K)*Pi11Regulated(I,J,K)
     &      -Vy(I,J,K)*Pi12Regulated(I,J,K)
#endif
          II = I
          JJ = J
          failed = 1
          return
        End If

        !transversality of pi tensor  u_mu pi^{mu,y} = 0
        trans = gamma_perp*(Pi02(I,J,K) - Vx(I,J,K)*Pi12(I,J,K)
     &                      - Vy(I,J,K)*Pi22(I,J,K))
        If (abs(trans) > max(relNumericalzero*pi_scale,
     &      absNumericalzero)) Then
#ifdef OUTPUT_VIOLATION
          Print*, "Time=", Time
          Print*, "Tranversality violated!"
          Print*, "I,J=", I,J
          Print*, "Pi02-Vx*Pi12-Vy*Pi22=", trans
          Print*, "Pi02,Pi12,Pi22=",
     &      Pi02(I,J,K),Pi12(I,J,K),Pi22(I,J,K)
          Print*, "Vx,Vy=",Vx(I,J,K),Vy(I,J,K)
          Print*, "Ed,PL=", Ed(I,J,K),PL(I,J,K)
          Print*, "Before: Pi00,Pi11,Pi22,Pi33*T^2,Pi01,Pi02,Pi12=",
     &      Pi00Regulated(I,J,K),Pi11Regulated(I,J,K),
     &      Pi22Regulated(I,J,K),Pi33Regulated(I,J,K),
     &      Pi01Regulated(I,J,K),Pi02Regulated(I,J,K),
     &      Pi12Regulated(I,J,K)
          Print*, "Before: Pi02-Vx*Pi12-Vy*Pi22=",
     &      Pi02Regulated(I,J,K)-Vx(I,J,K)*Pi12Regulated(I,J,K)
     &      -Vy(I,J,K)*Pi22Regulated(I,J,K)
#endif
          II = I
          JJ = J
          failed = 1
          return
        End If

        !transversality of pi tensor  u_mu pi^{mu,tau} = 0
        trans = gamma_perp*(Pi00(I,J,K) - Vx(I,J,K)*Pi01(I,J,K)
     &                      - Vy(I,J,K)*Pi02(I,J,K))
        If (abs(trans) > max(relNumericalzero*pi_scale,
     &      absNumericalzero)) Then
#ifdef OUTPUT_VIOLATION
          Print*, "Time=", Time
          Print*, "Tranversality violated!"
          Print*, "I,J=", I,J
          Print*, "Pi00,Pi01,Pi02=",
     &      Pi00(I,J,K),Pi01(I,J,K),Pi02(I,J,K)
          Print*, "Vx,Vy=",Vx(I,J,K),Vy(I,J,K)
          Print*, "Ed,PL=", Ed(I,J,K),PL(I,J,K)
          Print*, "Before: Pi00,Pi11,Pi22,Pi33*T^2,Pi01,Pi02,Pi12=",
     &      Pi00Regulated(I,J,K),Pi11Regulated(I,J,K),
     &      Pi22Regulated(I,J,K),Pi33Regulated(I,J,K),
     &      Pi01Regulated(I,J,K),Pi02Regulated(I,J,K),
     &      Pi12Regulated(I,J,K)
          Print*, "Before: Pi00-Vx*Pi01-Vy*Pi02=",
     &      Pi00Regulated(I,J,K)-Vx(I,J,K)*Pi01Regulated(I,J,K)
     &      -Vy(I,J,K)*Pi02Regulated(I,J,K)
#endif
          II = I
          JJ = J
          failed = 1
          return
        End If


      End Do
      End Do
      End Do

      End Subroutine
!-----------------------------------------------------------------------


      Subroutine checkBulkPi(failed, II, JJ, Time, PL,
     &  NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
     &  PPI)
      ! Check the size of Bulk pressure

      Implicit None

      Integer NXPhy0,NXPhy,NYPhy0,NYPhy,NX0,NX,NY0,NY,NZ0,NZ
      Integer failed, II, JJ ! (II,JJ): where it fails
      Double Precision Time

      Double Precision PL(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision PPI(NX0:NX, NY0:NY, NZ0:NZ)     !Bulk pressure

      Integer I,J,K
      Double Precision :: pressure_scale, bulkPi_scale
      Double Precision :: absNumericalzero = 1D-2
      !Double Precision :: relNumericalzero = 1D-2  !Xsi_0 in Zhi's thesis

      Double Precision maxBulkPiRatio
      Common /maxBulkPiRatio/ maxBulkPiRatio

      failed = 0

      Do K=NZ0,NZ
      Do J=NYPhy0-3,NYPhy+3
      Do I=NXPhy0-3,NXPhy+3

        !Largeness of bulk pressure
        pressure_scale = abs(PL(I,J,K))
        bulkPi_scale = abs(PPI(I,J,K))

        If (bulkPi_scale > max(maxBulkPiRatio*pressure_scale,
     &      absNumericalzero)) Then
#ifdef OUTPUT_VIOLATION
            Print*, "Time=", Time
            Print*, "Bulk Pi is larger than pressure!"
            Print*, "I,J=", I,J
            Print*, "Bulk Pi=", PPI(I,J,K)
            Print*, "Pressure=", PL(I,J,K)
#endif
          II = I
          JJ = J
          failed = 1
          return
        End If

      End Do
      End Do
      End Do

      End Subroutine
