#include "defs.h"

      Subroutine prepareInputFun()
!     Purpose:
!       Do some common preparation for other functions in this file

      Implicit None

      Integer regMethod ! used to determine which method of regulation should be applied.
!     1: find R0 (see PiRatio) + Fermi-Dirac; 2: use tanh (see maxPiRatio)
      Common /regMethod/ regMethod

      Double Precision PiRatio ! used to determine R0; within r<R0, Pi/(e+p) < PiRatio
      Common /PiRatio/ PiRatio

      Double Precision maxPiRatio ! used in tanh regulation method: Pi=Pi_max*tanh(Pi/Pi_max), Pi_max=maxPiRatio*(e+p)
      Common /maxPiRatio/ maxPiRatio

      Double Precision maxBulkPiRatio ! used in tanh regulation method: Pi=Pi_max*tanh(Pi/Pi_max), Pi_max=maxPiRatio*(e+p)
      Common /maxBulkPiRatio/ maxBulkPiRatio

      ! hard-code rather than read from extra file
      regMethod = 2
      PiRatio = 0.141421
      maxPiRatio = 10.0
      maxBulkPiRatio = 10.0

      End Subroutine


      Subroutine readInputFromCML2()
!     Purpose:
!     Read inputs from command line.

      Implicit None ! enforce explicit variable declaration

      Integer IEin
      Common /IEin/ IEin     !  type of initialization  entropy/enrgy

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

      double precision :: DT
      common /DT/ DT

      double precision :: DX, DY
      common /DXY/ DX, DY

      Double Precision Edec
      Common/Edec/Edec    !decoupling temperature

      Integer InitialURead
      Common/LDInitial/ InitialURead  ! IintURead =1 read initial velocity profile

      Integer NDT, NDX, NDY
      Common /NDTXY/ NDT, NDX, NDY

      Double Precision T0
      Common /T0/ T0

      Integer LS
      Common /LS/ LS

      Integer QNum, ArgIndex ! QNum is the total number of arguments, ArgIndex gives the index to the one currently reading

      Character*60 :: buffer
      Character*20 :: varName
      Integer IResult
      Double Precision DResult

      QNum = iargc ()

      Do ArgIndex = 1, QNum
        Call getarg(ArgIndex, buffer)
        Call processAssignment(buffer, "=", varName, IResult, DResult)

        If (varName=="iein") IEin=IResult ! 0: initialize by energy density; 1: initialize by entropy density
        If (varName=="iin") IEin=IResult

        If (varName=="dt") DT=DResult
        If (varName=="dxy") DX=DResult  ! will be copied to DY

        If (varName=="edec") EDec=DResult ! decouple energy density, in GeV/fm^3
        If (varName=="e_dec") EDec=DResult
        If (varName=="e_d") EDec=DResult

        If (varName=="t0") T0=DResult ! initial proper time tau_0, in fm/c

        If (varName=="vist0") VisT0=DResult
        If (varName=="etast0") VisT0=DResult
        If (varName=="etas_t0") VisT0=DResult
        If (varName=="eta_s_t0") VisT0=DResult

        If (varName=="vishrg") VisHRG=DResult
        If (varName=="etashrg") VisHRG=DResult
        If (varName=="etas_hrg") VisHRG=DResult
        If (varName=="eta_s_hrg") VisHRG=DResult

        If (varName=="vismin") VisMin=DResult
        If (varName=="etasmin") VisMin=DResult
        If (varName=="etas_min") VisMin=DResult
        If (varName=="eta_s_min") VisMin=DResult

        If (varName=="visslope") VisSlope=DResult
        If (varName=="etasslope") VisSlope=DResult
        If (varName=="etas_slope") VisSlope=DResult
        If (varName=="eta_s_slope") VisSlope=DResult

        If (varName=="viscurv") VisCrv=DResult
        If (varName=="etascurv") VisCrv=DResult
        If (varName=="etas_curv") VisCrv=DResult
        If (varName=="eta_s_curv") VisCrv=DResult
        If (varName=="viscrv") VisCrv=DResult
        If (varName=="etascrv") VisCrv=DResult
        If (varName=="etas_crv") VisCrv=DResult
        If (varName=="eta_s_crv") VisCrv=DResult

        If (varName=="visbulkt0") VisBulkT0=DResult
        If (varName=="zetast0") VisBulkT0=DResult
        If (varName=="zetas_t0") VisBulkT0=DResult
        If (varName=="zeta_s_t0") VisBulkT0=DResult

        If (varName=="visbulkmax") VisBulkMax=DResult
        If (varName=="zetasmax") VisBulkMax=DResult
        If (varName=="zetas_max") VisBulkMax=DResult
        If (varName=="zeta_s_max") VisBulkMax=DResult

        If (varName=="visbulkwidth") VisBulkWidth=DResult
        If (varName=="zetaswidth") VisBulkWidth=DResult
        If (varName=="zetas_width") VisBulkWidth=DResult
        If (varName=="zeta_s_width") VisBulkWidth=DResult

        If (varName=="ils") LS=IResult
        If (varName=="nls") LS=IResult

        If (varName=="ndt") NDT=IResult
        If (varName=="ndxy") NDX=IResult  ! will be copied to NDY

        If (varName=="visbeta") VisBeta=DResult ! VisBeta, used for proper time tau_pi

        If (varName=="initialuread") InitialURead=IResult ! read in initial flow velocity profiles

        If (varName=="initialpitensor") Initialpitensor=IResult ! initialization of pi tensor

      End Do ! ArgIndex

      End Subroutine
!-----------------------------------------------------------------------




************************************************************************
      Subroutine processAssignment(string, separator,
     &                            varName, IResult, DResult)
!     This subroutine process a string assignment.
!     First it seprate string into LHS and RHS according to separator.
!     Then the LHS is converted into variable using only lower case
!     letters, and the RHS is converted into numerical values.
!     The variable IResult holds result for integer and DResult holds
!     one for double.
!     Convention: integer-valued variable should start with "I" or "N"
!     for its name.

      Implicit None

      Character (*) :: string, varName
      Character*60 :: LHS, RHS
      Character separator
      Integer IResult
      Double Precision DResult

      Integer break_here, I, cha

      varName = ""

      break_here = index(string, separator)
      LHS = adjustl(string(:break_here-1))
      RHS = adjustl(string(break_here+1:))

      ! convert LHS to lower case:
      Do I = 1,len_trim(LHS)
        cha = ichar(LHS(I:I))
        If (cha>=65 .and. cha<90) Then
          varName(I:I) = char(cha+32)
        Else
          varName(I:I) = LHS(I:I)
        EndIf
      EndDo

      ! convert RHS to numerics:
      If (varName(1:1)=="i" .or. varName(1:1)=="n") Then
        Read(RHS, *) IResult
      Else
        Read(RHS, *) DResult
      EndIf

      End Subroutine
!-----------------------------------------------------------------------




************************************************************************
      Subroutine getInitialR0(PU0,PU1,PU2,PU3,U0,U1,U2,U3,DX,DY,DZ,DT,
     &  DPc00,DPc01,DPc02,DPc33,DPc11,DPc22,DPc12,DDU0,DDU1,DDU2,
     &  Temp,Temp0,SiLoc,DLnT,Time, NXPhy0,NYPhy0,NXPhy,NYPhy,
     &  NX0,NX,NY0,NY,NZ0,NZ,Ed,Sd,PL,VCoefi)
!     Purpose:
!       Return a suitable initial (before the initialization of
!       Pi(mu,nu)) R0 (thru common) to regulate Pi(mu,nu)

      Implicit None

      Integer I,J,K
      Integer NX0,NY0,NZ0,NX,NY,NZ
      Integer NXPhy0,NYPhy0,NXPhy,NYPhy

      Double Precision PU0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity from last time step
      Double Precision PU1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Double Precision PU2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Double Precision PU3(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

      Double Precision U0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Double Precision U1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Double Precision U2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Double Precision U3(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

      Double Precision DPc00(NX0:NX, NY0:NY, NZ0:NZ) ! DIfferential part of Pi source term
      Double Precision DPc01(NX0:NX, NY0:NY, NZ0:NZ) !
      Double Precision DPc02(NX0:NX, NY0:NY, NZ0:NZ) !
      Double Precision DPc33(NX0:NX, NY0:NY, NZ0:NZ) !

      Double Precision DPc11(NX0:NX, NY0:NY, NZ0:NZ) ! DIfferential part of Pi source term
      Double Precision DPc12(NX0:NX, NY0:NY, NZ0:NZ) !
      Double Precision DPc22(NX0:NX, NY0:NY, NZ0:NZ) !

      Double Precision DDU0(NX0:NX, NY0:NY, NZ0:NZ) ! DIfferential part of Pi source term
      Double Precision DDU1(NX0:NX, NY0:NY, NZ0:NZ) !
      Double Precision DDU2(NX0:NX, NY0:NY, NZ0:NZ) !

      Double Precision Temp0(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature  in last time step
      Double Precision Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature
      Double Precision SiLoc(NX0:NX, NY0:NY, NZ0:NZ) ! Local expansion rate \sita
      Double Precision DLnT(NX0:NX, NY0:NY, NZ0:NZ) ! DlnT(x,y) terms

      Double Precision R0, Aeps
      Common /R0Aeps/ R0,Aeps

      Double Precision Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Double Precision Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density
      Double Precision PL(NX0:NX, NY0:NY, NZ0:NZ) !pressure density
      Double Precision VCoefi(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient shear viscosity eta
      Double Precision RMin, PiEPRatio, SigmaLargeness, EAndP

      double precision :: VisT0, VisHRG, VisMin, VisSlope, VisCrv,
     &                    VisBeta
      common /VisShear/ VisT0, VisHRG, VisMin, VisSlope, VisCrv, VisBeta

      Double Precision PiRatio ! used to determine R0; within r<R0, Pi/(e+p) < PiRatio
      Common /PiRatio/ PiRatio ! should already be setuped in prepareInputFun function

      Double Precision D0U0,D0U1,D0U2,D1U0,D1U1,D1U2,D2U0,D2U1,D2U2
      Double Precision CS,DT,DX,DY,DZ,Time,DU0,DU1,DU2

      double precision ViscousCTemp

      Integer regMethod
      Common /regMethod/ regMethod

      If (regMethod .eq. 1) Then
        DO 791 K=NZ0,NZ
        DO 791 J=NYPhy0,NYPhy
        DO 791 I=NXPhy0,NXPhy

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
     &        -U0(I+2,J,K)/12.0d0+U0(I-2,J,K)/12.0d0)/DX
          D1U1=(U1(I+1,J,K)*2.0d0/3.0d0-U1(I-1,J,K)*2.0d0/3.0d0
     &        -U1(I+2,J,K)/12.0d0+U1(I-2,J,K)/12.0d0)/DX
          D1U2=(U2(I+1,J,K)*2.0d0/3.0d0-U2(I-1,J,K)*2.0d0/3.0d0
     &        -U2(I+2,J,K)/12.0d0+U2(I-2,J,K)/12.0d0)/DX
          D2U0=(U0(I,J+1,K)*2.0d0/3.0d0-U0(I,J-1,K)*2.0d0/3.0d0
     &        -U0(I,J+2,K)/12.0d0+U0(I,J-2,K)/12.0d0)/DY
          D2U1=(U1(I,J+1,K)*2.0d0/3.0d0-U1(I,J-1,K)*2.0d0/3.0d0
     &        -U1(I,J+2,K)/12.0d0+U1(I,J-2,K)/12.0d0)/DY
          D2U2=(U2(I,J+1,K)*2.0d0/3.0d0-U2(I,J-1,K)*2.0d0/3.0d0
     &        -U2(I,J+2,K)/12.0d0+U2(I,J-2,K)/12.0d0)/DY
#else
#error "DERIV_ORDER must be 3 or 5"
#endif

        CS=(D0U0+D1U1+D2U2+U0(I,J,K)/Time)/3.0

        DU0=U0(I,J,K)*D0U0+U1(I,J,K)*D1U0+U2(I,J,K)*D2U0
        DU1=U0(I,J,K)*D0U1+U1(I,J,K)*D1U1+U2(I,J,K)*D2U1
        DU2=U0(I,J,K)*D0U2+U1(I,J,K)*D1U2+U2(I,J,K)*D2U2

        DPc00(I,J,K)=D0U0-U0(I,J,K)*DU0+CS*(U0(I,J,K)**2-1.0)
        DPc01(I,J,K)=0.5*(D0U1-D1U0)-0.5*(U1(I,J,K)*DU0+U0(I,J,K)*DU1)
     &              +CS*(U1(I,J,K)*U0(I,J,K))
        DPc02(I,J,K)=0.5*(D0U2-D2U0)-0.5*(U2(I,J,K)*DU0+U0(I,J,K)*DU2)
     &              +CS*(U2(I,J,K)*U0(I,J,K))
        DPc33(I,J,K)=CS-U0(I,J,K)/Time
        DPc11(I,J,K)=(-1.0)*D1U1-U1(I,J,K)*DU1+CS*(U1(I,J,K)**2+1.0)
        DPc22(I,J,K)=(-1.0)*D2U2-U2(I,J,K)*DU2+CS*(U2(I,J,K)**2+1.0)
        DPc12(I,J,K)=(-0.5)*(D2U1+D1U2)-0.5*(U1(I,J,K)*DU2
     &              +U2(I,J,K)*DU1)+CS*(U1(I,J,K)*U2(I,J,K))
 791    Continue

        RMin = NX*DX+NY*DY !---Upper-Bound-R0---

        DO 3007 K=NZ0,NZ !Check for Pi tensor
        DO 3007 J=NYPhy0,NYPhy
        DO 3007 I=NXPhy0,NXPhy
          EAndP = Abs(Ed(I,J,K)+PL(I,J,K))
          SigmaLargeness = 1/7.0*(Abs(DPc00(I,J,K))+
     &      Abs(DPc01(I,J,K))+Abs(DPc02(I,J,K))+Abs(DPc33(I,J,K))+
     &      Abs(DPc11(I,J,K))+Abs(DPc12(I,J,K))+Abs(DPc22(I,J,K)))
          PiEPRatio = 2*ViscousCTemp(Temp(I,J,K))*Sd(I,J,K)
     &                 *SigmaLargeness/EAndP

        If (PiEPRatio > PiRatio) Then
          If (sqrt(DX*DX*I*I+DY*DY*J*J) < RMin) Then
            RMin = sqrt(DX*DX*I*I+DY*DY*J*J)
          EndIf
        EndIf

 3007   Continue
        R0 = RMin
      ElseIf (regMethod .eq. 2) Then ! use maximun possible R0
        R0 = NX*DX+NY*DY
      Else  ! use R0=12
        R0 = 12.0
      EndIf ! corresponding to the one on variable "regMethod"

      End Subroutine
!-----------------------------------------------------------------------



************************************************************************
      Subroutine determineR0(NX0,NY0,NZ0,NX,NY,NZ,Ed,PL,Sd,
     &  Pi00,Pi01,Pi02,Pi11,Pi12,Pi22,Pi33)
!     Purpose:
!       Return a suitable R0 (thru common) to regulate Pi(mu,nu)

      Implicit None

      Integer NX0,NY0,NZ0,NX,NY,NZ
      Integer I,J,K

      double precision :: DX, DY
      common /DXY/ DX, DY

      Double Precision Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Double Precision PL(NX0:NX, NY0:NY, NZ0:NZ) !local pressure
      Double Precision Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density

      Double Precision Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      Double Precision R0, Aeps
      Common /R0Aeps/ R0, Aeps

      Double Precision PiRatio ! used to determine R0; within r<R0, Pi/(e+p) < PiRatio
      Common /PiRatio/ PiRatio ! should already be setuped in prepareInputFun function


      Double Precision PiLargeness ! as a measurement of how large Pi is
      Double Precision PiEPRatio, EAndP
      Double Precision RMin ! an intermedia variable that trace the radius of the largest possible region

      Integer regMethod
      Common /regMethod/ regMethod

      If (regMethod .eq. 1) Then
        RMin = NX*DX+NY*DY !---Upper-Bound-R0---
        DO 3007 K=NZ0,NZ !Check for Pi tensor
        DO 3007 J=NY0,NY
        DO 3007 I=NX0,NX
          EAndP = Abs(Ed(I,J,K)+PL(I,J,K))
          PiLargeness = 1/7.0*(Abs(Pi00(I,J,K))+
     &      Abs(Pi01(I,J,K))+Abs(Pi02(I,J,K))+Abs(Pi33(I,J,K))+
     &      Abs(Pi11(I,J,K))+Abs(Pi12(I,J,K))+Abs(Pi22(I,J,K)))
          PiEPRatio=PiLargeness/EAndP

          If (PiEPRatio > PiRatio) Then
              If (sqrt(DX*DX*I*I+DY*DY*J*J) < RMin) Then
                RMin = sqrt(DX*DX*I*I+DY*DY*J*J)
              EndIf
          EndIf
 3007   Continue
        R0 = RMin
      ElseIf (regMethod .eq. 2) Then ! use maximun possible R0
        R0 = (NX*DX+NY*DY)*2.0
      Else ! use R0=12.0
        R0 = 12.0
      EndIf

      End Subroutine
!-----------------------------------------------------------------------


************************************************************************
      Subroutine regulateBulkPi(regStr,Time,NX0,NY0,NZ0,NX,NY,NZ,
     &  NXPhy0,NXPhy,NYPhy0,NYPhy,
     &  Ed,PL,PPI,II,JJ)
!     Purpose:
!       Regulate Bulk pressure tensor by restrain it under a maximum
!       value using tanh function

      Implicit None

      Integer NX0,NY0,NZ0,NX,NY,NZ,NXPhy0,NXPhy,NYPhy0,NYPhy
      Integer I,J,K,II,JJ,regStr

      Double Precision Time

      Double Precision Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Double Precision PL(NX0:NX, NY0:NY, NZ0:NZ) !local pressure

      Double Precision PPI(NX0:NX, NY0:NY, NZ0:NZ) ! Bulk Pressure Tensor

      Double Precision BulkPi

      Integer regMethod
      Common /regMethod/ regMethod

      Double Precision :: Xsi0 = 1D0  !adaptive zero
      Double Precision :: pressure_scale, bulkPi_scale
      Double Precision regStrength

      Double Precision maxBulkPiRatio
      Common /maxBulkPiRatio/ maxBulkPiRatio

      Xsi0 = 1D-2/(regStr+1D0) ! VER-1.29RC: adaptive zero chooser VER-1.29RC4: bug fix: regStr -> regStr+1D0

      If (regMethod == 2) Then ! do tanh regulation

        DO 3019 K=NZ0,NZ
        DO 3019 J=NY0,NY
        DO 3018 I=NX0,NX

        regStrength = 1D-30

        pressure_scale = abs(PL(I,J,K))

        BulkPi = PPI(I,J,K)

        ! get Bulk pi scale
        bulkPi_scale = abs(BulkPi) + 1d-30
        if (isnan(bulkPi_scale)) then
           print*, "Bulk Pi is NaN, I,J =", I, J
           call exit(1)
        endif

        ! find regulation strength using largeness comparison
        regStrength = max(bulkPi_scale/(maxBulkPiRatio*pressure_scale),
     &                    regStrength)

        PPI(I,J,K)=PPI(I,J,K)*(tanh(regStrength)/regStrength) ! Bulk pressure PPI is regulated here

3018    Continue
3019    Continue

      EndIf ! on regMethod

      End Subroutine
!-----------------------------------------------------------------------------

************************************************************************
      Subroutine regulatePi(regStr,Time,NX0,NY0,NZ0,NX,NY,NZ,
     &  NXPhy0,NXPhy,NYPhy0,NYPhy,
     &  Ed,PL,PPI,
     &  Pi00,Pi01,Pi02,Pi11,Pi12,Pi22,Pi33,Vx,Vy,II,JJ)
!     Purpose:
!       Regulate Pi(mu,nu) tensor by restrain it under a maximum
!       value using tanh function

      Implicit None

      Integer NX0,NY0,NZ0,NX,NY,NZ,NXPhy0,NXPhy,NYPhy0,NYPhy
      Integer I,J,K,II,JJ,regStr

      Double Precision Time

      Double Precision Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Double Precision PL(NX0:NX, NY0:NY, NZ0:NZ) !local pressure

      Double Precision Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      Double Precision Vx(NX0:NX, NY0:NY, NZ0:NZ)
      Double Precision Vy(NX0:NX, NY0:NY, NZ0:NZ)

      Double Precision PPI(NX0:NX, NY0:NY, NZ0:NZ) ! Bulk Pressure Tensor

      Double Precision p00,p01,p02,p11,p12,p22,p33,vvx,vvy! These are just Pi at I,J,K (in a loop)

      Double Precision TrPi2 ! Tr(pi^2)

      Double Precision rTrPi2EAndP ! ratio between sqrt(Tr(pi^2)) and e+p

      Integer regMethod
      Common /regMethod/ regMethod

      Double Precision :: Xsi0 = 1D0  !adaptive zero
      Double Precision :: Tideal_scale, pi_scale
      Double Precision regStrength

      Double Precision maxPiRatio
      Common /maxPiRatio/ maxPiRatio
      Double Precision maxPi ! maxPi = maxPiRatio*(e+p)

      Double Precision PiAvg
      Double Precision gamma_perp

      Double Precision PiTr, trans

      Xsi0 = 1D-2/(regStr+1D0) ! VER-1.29RC: adaptive zero chooser VER-1.29RC4: bug fix: regStr -> regStr+1D0

      If (regMethod == 2) Then ! do tanh regulation
        DO 3009 K=NZ0,NZ
        DO 3009 J=NY0,NY
        DO 3008 I=NX0,NX

        regStrength = 1D-30

        vvx = Vx(I,J,K)
        vvy = Vy(I,J,K)
        gamma_perp = 1./sqrt(1. - vvx**2 - vvy**2 + 1D-30)
        Tideal_scale = sqrt(Ed(I,J,K)**2 + 3*PL(I,J,K)**2)

        p00 = Pi00(I,J,K)
        p01 = Pi01(I,J,K)
        p02 = Pi02(I,J,K)
        p11 = Pi11(I,J,K)
        p12 = Pi12(I,J,K)
        p22 = Pi22(I,J,K)
        p33 = Pi33(I,J,K)   ! pi are in t-xyz coordinate

        ! calculate Tr(pi^2)
        TrPi2 = p00*p00+p11*p11+p22*p22+p33*p33
     &    -2*p01*p01-2*p02*p02+2*p12*p12
        pi_scale = sqrt(abs(TrPi2)) + 1D-30

        ! find regulation strength

        ! first, tracelessness
        PiTr = p00-p11-p22-p33
        regStrength = max(abs(PiTr)/(Xsi0*MaxPiRatio*pi_scale),
     &                    regStrength)

        ! next transversality
        trans = gamma_perp*(p01-vvx*p11-vvy*p12)
        regStrength = max(abs(trans)/(Xsi0*MaxPiRatio*pi_scale),
     &                    regStrength)
        trans = gamma_perp*(p02-vvx*p12-vvy*p22)
        regStrength = max(abs(trans)/(Xsi0*MaxPiRatio*pi_scale),
     &                    regStrength)
        trans = gamma_perp*(p00-vvx*p01-vvy*p02)
        regStrength = max(abs(trans)/(Xsi0*MaxPiRatio*pi_scale),
     &                    regStrength)

        ! largeness comparision
        rTrPi2EAndP = pi_scale/(MaxPiRatio*Tideal_scale) + 1e-30
        regStrength = max(rTrPi2EAndP, regStrength)

        Pi00(I,J,K)=Pi00(I,J,K)*(tanh(regStrength)/regStrength) ! Pi## is regulated here
        Pi01(I,J,K)=Pi01(I,J,K)*(tanh(regStrength)/regStrength) ! Pi## is regulated here
        Pi02(I,J,K)=Pi02(I,J,K)*(tanh(regStrength)/regStrength) ! Pi## is regulated here
        Pi11(I,J,K)=Pi11(I,J,K)*(tanh(regStrength)/regStrength) ! Pi## is regulated here
        Pi12(I,J,K)=Pi12(I,J,K)*(tanh(regStrength)/regStrength) ! Pi## is regulated here
        Pi22(I,J,K)=Pi22(I,J,K)*(tanh(regStrength)/regStrength) ! Pi## is regulated here
        Pi33(I,J,K)=Pi33(I,J,K)*(tanh(regStrength)/regStrength) ! Pi## is regulated here

        PiAvg = 1.0d0/7.0d0*
     &      (abs(Pi00(I,J,K))+abs(Pi01(I,J,K))
     &      +abs(Pi02(I,J,K))
     &      +abs(Pi11(I,J,K))+abs(Pi12(I,J,K))
     &      +abs(Pi22(I,J,K))
     &      +abs(Pi33(I,J,K)))

        If (isnan(PiAvg)) Then
          Print *, "Invalid PiAvg"
          Print *, "(I,J,K)=",I,J,K
          Print *, "e=", Ed(I,J,K)
          Print *, "p=", PL(I,J,K)
          Print *, "maxPi=", maxPi
          Print *, "TrPi2=",TrPi2
          Print *, "rTrPi2EAndP=",rTrPi2EAndP
          Print *, "Pi00=", Pi00(I,J,K)
          Print *, "Pi01=", Pi01(I,J,K)
          Print *, "Pi02=", Pi02(I,J,K)
          Print *, "Pi11=", Pi11(I,J,K)
          Print *, "Pi12=", Pi12(I,J,K)
          Print *, "Pi22=", Pi22(I,J,K)
          Print *, "Pi33=", Pi33(I,J,K)
          call exit(1)
        EndIf

3008    Continue
3009    Continue

      EndIf ! on regMethod

      End Subroutine
!-----------------------------------------------------------------------------




!*****************************************************************************
      Subroutine regulateAllPi(NX0,NY0,NZ0,NX,NY,NZ,Ed,PL,U0,U1,U2,Time,
     &  Pi00,Pi01,Pi02,Pi11,Pi12,Pi22,Pi33,
     &  NXPhy0,NXPhy,NYPhy0,NYPhy,ratio)

      Implicit None

      Integer NX0,NY0,NZ0,NX,NY,NZ,NXPhy0,NXPhy,NYPhy0,NYPhy
      Double Precision ratio

      Double Precision Time

      Double Precision Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Double Precision PL(NX0:NX, NY0:NY, NZ0:NZ) !local pressure

      Double Precision U0(NX0:NX, NY0:NY, NZ0:NZ) !local pressure
      Double Precision U1(NX0:NX, NY0:NY, NZ0:NZ) !local pressure
      Double Precision U2(NX0:NX, NY0:NY, NZ0:NZ) !local pressure

      Double Precision Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      Double Precision CC(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      CC = 1D0
      CC = min(CC, ratio*abs(Pi00)/max(abs(Ed*U0*U0), 1D-30))
      CC = min(CC, ratio*abs(Pi11)/max(abs((Ed+PL)*U1*U1+PL), 1D-30))
      CC = min(CC, ratio*abs(Pi22)/max(abs((Ed+PL)*U2*U2+PL), 1D-30))
      CC = min(CC, ratio*abs(Pi33)/max(abs(PL), 1D-30))
      CC = min(CC, ratio*abs(Pi01)/max(abs((Ed+PL)*U0*U1), 1D-30))
      CC = min(CC, ratio*abs(Pi02)/max(abs((Ed+PL)*U0*U2), 1D-30))
      CC = min(CC, ratio*abs(Pi12)/max(abs((Ed+PL)*U1*U2), 1D-30))


      Pi00 = CC*Pi00
      Pi11 = CC*Pi11
      Pi22 = CC*Pi22
      Pi33 = CC*Pi33
      Pi01 = CC*Pi01
      Pi02 = CC*Pi02
      Pi12 = CC*Pi12

      End Subroutine
