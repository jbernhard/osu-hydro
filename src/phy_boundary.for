!###########################################################################
! Version 1.02
! -- For changes, see ChangeLogsPhyBDary.txt
!---------------------------------------------------------------------------

C###########################################################################
      Subroutine UPShasta2(TT,Vx,Vy,ScT, NX0,NX,NY0,NY,NZ0,NZ,  IPX,IPY,
     &  DT,DX,DY,  NXPhy0,NYPhy0, NXPhy,NYPhy, ISYM1,ISYM2, DIFFC )
C------ using P.Colb 2+1 Shastala to cal the tansport part
C------ Transport diffusion part:
C-------limiter (antidiffusion part)
C-----  IPX,IPY ! differential x(or y) for source term
C------ NXPhy, NYphy :Real Physics demension < NX <NY
C------ ISYM1,ISYM2 symetric boundary condition for 2+1
        Implicit Double Precision (A-H, O-Z)

        Dimension TT(NX0:NX, NY0:NY, NZ0:NZ)  !Tansport density
        Dimension Vx(NX0:NX, NY0:NY, NZ0:NZ)  !Velocity in X
        Dimension Vy(NX0:NX, NY0:NY, NZ0:NZ)  !Velocity in Y
        Dimension ScT(NX0:NX, NY0:NY, NZ0:NZ) !Source Term  ScT=ScX+ScY+ScZ

        Dimension PTT(NX0:NX, NY0:NY)
        Dimension PVx(NX0:NX, NY0:NY)
        Dimension PVy(NX0:NX, NY0:NY)
        Dimension PSc(NX0:NX, NY0:NY)
        DIMENSION WORK1(NX0:NX,NY0:NY),WORK2(NX0:NX,NY0:NY)

C        if(NX0.ne.-2.or.NY0.ne.-2.or.NX.ne.MAX .or. NY.ne.MAY ) then
C          print 'in UPShasta wrong dimension'
C          stop
C        end if

         do 100 I=NX0,NX   !(NX0,NX)
         do 100 J=NY0,NY   !(NY0,NY)
           PTT(I,J)=TT(I,J,NZ0)     !TT(I+3, J+3)
           PVx(I,J)=Vx(I,J,NZ0)
           PVy(I,J)=Vy(I,J,NZ0)
           PSc(I,J)=ScT(I,J,NZ0)
 100      continue

        CALL SHATRA(NXPhy0,NYPhy0,NXPhy,NYPhy,PTT,PVx,PVy,PSc, IPX, IPY,
     &         WORK1,  ISYM1,ISYM2, DT/DX,DT/DY,  NX0, NY0,  NX,NY)
        CALL SHACOR(NXPhy0,NYPhy0,NXPhy,NYPhy,PTT,PVx,PVy, WORK1,
     &    DT/DX,DT/DY, WORK2,WORK1,DIFFC, ISYM1,ISYM2, NX0, NY0,  NX,NY)

          If(IPX.eq.0.and.IPY.eq.0 ) then
C           DO 500 J=NYPhy0-2,NYPhy+2
C           DO 500 I=NXPhy0-2,NXPhy+2
           DO 500 J=NYPhy0,NYPhy
           DO 500 I=NXPhy0,NXPhy
           PTT(I,J)=PTT(I,J)-PSc(I,J)*DT
 500       Continue
          end if

         do 800 I=NX0,NX   !(NX0,NX)
         do 800 J=NY0,NY   !(NY0,NY)
           TT(I,J,NZ0)=PTT(I,J)     !TT(I+3, J+3)
 800      continue

      Return
      End

C##################################################################
*************************************************************

      SUBROUTINE SHATRA(NXPhy0,NYPhy0,NXPhy,NYPhy,A,U,V,B,ISL1,ISL2,AB,
     &        ISYM1,ISYM2,  TDX,TDY, NX0, NY0, NX,NY)
*     In order to get the symmetry relations right we need again
*     a routine to introduce the right conditions, especially for
*     the y-Component of the velocity, that is T02
        Implicit Double Precision (A-H, O-Z)
      DIMENSION A(NX0:NX,NY0:NY),U(NX0:NX,NY0:NY),
     &          V(NX0:NX,NY0:NY),B(NX0:NX,NY0:NY),
     &          AB(NX0:NX,NY0:NY)


      CALL SHAS1(NXPhy0,NYPhy0,NXPhy,NYPhy,A,U,V,B,ISL1,ISL2,AB,
     &                        TDX,TDY, NX0, NY0, NX,NY)

*      Ab  work1
*     Now AB is the transported and diffused field that will
*     replace A
*     Care must still be taken on the boundary conditions

      DO I=NXPhy0-1,NXPhy+1
cCC        AB(I,-1)=ISYM2*AB(I,1)
cCC        AB(I,-2)=ISYM2*AB(I,2)
!        AB(I,NYPhy+2) = 0.0
!        AB(I,NYPhy0-2) = 0.0
        AB(I,NYPhy+2) =2.0*AB(I,NYPhy+1)-AB(I,NYPhy)
        AB(I,NYPhy0-2)=2.0*AB(I,NYPhy0-1)-AB(I,NYPhy0)
      END DO

      DO J=NYPhy0-2,NYPhy+2
CCC        AB(-1,J)=ISYM1*AB(1,J)
CCC        AB(-2,J)=ISYM1*AB(2,J)
!        AB(NXPhy+2,J) = 0.0
!        AB(NXPhy0-2,J) = 0.0
        AB(NXPhy+2,J) =2.0*AB(NXPhy+1,J)-AB(NXPhy,J)
        AB(NXPhy0-2,J)=2.0*AB(NXPhy0-1,J)-AB(NXPhy0,J)
      END DO

      DO I=NXPhy0-2,NXPhy+2
        AB(I,NYPhy+2)=(2.0*AB(I,NYPhy+1)-AB(I,NYPhy))*0.5
     &                                      +0.5*AB(I,NYPhy+2)
        AB(I,NYPhy0-2)=(2.0*AB(I,NYPhy0-1)-AB(I,NYPhy0))*0.5
     &                                      +0.5*AB(I,NYPhy0-2)
      END DO

cCC        AB(-1,-1)=ISYM1*ISYM2*AB(1,1)
cCC        AB(-2,-1)=ISYM1*ISYM2*AB(2,1)
cCC        AB(-1,-2)=ISYM1*ISYM2*AB(1,2)
cCC        AB(-2,-2)=ISYM1*ISYM2*AB(2,2)

      RETURN
      END

*************************************************************

      SUBROUTINE SHACOR(NXPhy0,NYPhy0,NXPhy,NYPhy,A,U,V,AB,TDX,TDY,
     &                  AC1,GAB,DIFFC,ISYM1,ISYM2, NX0,NY0,NX,NY)
        Implicit Double Precision (A-H, O-Z)
*************************************************************
**  THIS SUBROUTINE CALLS THE ANTIDIFFUSION ALGORITHM TO    *
** DO THE CORRECTION ON THE TRANSPORTED AND DIFFUSED        *
** SOLUTION OF THE DIFFERENTIAL EQUATION                    *
**        DA/DT=-D(AU)/DX-D(AV)/DY-DB/DX                    *
**      USING 'SHASTA-MULTIDIMENSIONAL FCT' ALGORITHM       *
** ONE CALL OF THE ROUTINE ADVANCES THE DEPENDENT VARIABLE  *
** A BY ONE TIME STEP. THIS IS DONE IN TWO STEPS: THE FIRST *
** (SHAS1) IS THE 'TRANSPORT AND DIFFUSION', AND THE        *
** SECOND (SHAS2) IS THE 'ANTIDIFFUSION' STEP.              *
** AFTER EACH STEP THE BOUNDARIES MUST BE HANDLED           *
** SEPARATELY ACCORDING TO BOUNDARY CONDITIONS OR SYMMETRY. *
*************************************************************
** BOUNDARY IN THIS VERSION:                                *
** -SOLUTION IS ASSUMED TO BE SYMMETRIC OR ANTI-            *
**  SYMMETRIC ACCORDING TO PARAMETERS 'ISYM%' (=-1,0 OR 1)  *
**  WITH RESPECT TO X=X(MIN) (I=0).                         *
** -SOLUTION IS ASSUMED TO BE CONST. NEAR X=X(MAX).         *
*************************************************************

      DIMENSION A(NX0:NX,NY0:NY),U(NX0:NX,NY0:NY),
     &          V(NX0:NX,NY0:NY),
     &          AB(NX0:NX,NY0:NY),AC1(NX0:NX,NY0:NY),
     &          GAB(NX0:NX,NY0:NY)

      CALL SHAS2(NXPhy0,NYPhy0,NXPhy,NYPhy,A,U,V,AB,AC1,GAB,DIFFC,
     &                       TDX,TDY,NX0,NY0,NX,NY)

      DO 180 I=NXPhy0-2,NXPhy+2
        A(I,NYPhy+1)=0.0
        A(I,NYPhy+2)=0.0
        A(I,NYPhy0-1)=0.0
        A(I,NYPhy0-2)=0.0
CCC        A(I,-1)= ISYM2*A(I,1)
CCC        A(I,-2)= ISYM2*A(I,2)
180   CONTINUE

      DO 190 J=NYPhy0-2,NYPhy+2
        A(NXPhy+1,J)=0.0
        A(NXPhy+2,J)=0.0
        A(NXPhy0-1,J)=0.0
        A(NXPhy0-2,J)=0.0
CCC        A(-1,J)=ISYM1*A(1,J)
CCC        A(-2,J)=ISYM2*A(2,J)
190   CONTINUE

cCC      A(-1,-1)=ISYM1*ISYM2*A(1,1)
CCC      A(-2,-1)=ISYM1*ISYM2*A(2,1)
CCC      A(-1,-2)=ISYM1*ISYM2*A(1,2)
CCC      A(-2,-2)=ISYM1*ISYM2*A(2,2)

      RETURN
      END

*********************************************************

      SUBROUTINE SHAS1(NXPhy0,NYPhy0,NXPhy,NYPhy,A,U,V,B,ISL1,ISL2,AB,
     &                 TDX,TDY,NX0,NY0,NX,NY)
        Implicit Double Precision (A-H, O-Z)
*********************************************************
** THIS ROUTINE PRODUCES THE 'TRANSPORTED AND DIFFUSED' *
** SOLUTION AB FOR THE SHASTA ALGORITHM.                *
** NOTE: FOR INPUT: A(I,J),U(I,J),V(I,J),B(I,J)         *
**       FOR OUTPUT: AB(I,J)                            *
**            I=-1-N,...,N+1; J=-1-N,...,N+1            *
** NOTE: 'ISL' MUST BE EITHER 1 OR 0. (THIS IS NOT      *
**       CHECKED HERE!)                                 *
** NOTE: IF U(I,J)*DT/DX >0.5                           *
**  OR   IF V(I,J)*DT/DY >0.5 FOR SOME I, EXECUTION IS  *
**       TERMINATED AND AN ERROR MESSAGE PRINTED        *
*********************************************************

      DIMENSION A(NX0:NX,NY0:NY),U(NX0:NX,NY0:NY),
     &          V(NX0:NX,NY0:NY),B(NX0:NX,NY0:NY),
     &          AB(NX0:NX,NY0:NY)

          DO 100 I=NXPhy0-1,NXPhy+1
          DO 120 J=NYPhy0-1,NYPhy+1

          EIM1=U(I-1,J)*TDX
          EII=U(I,J)*TDX
          EIP1=U(I+1,J)*TDX
          EJM1=V(I,J-1)*TDY
          EJI=V(I,J)*TDY
          EJP1=V(I,J+1)*TDY

          IF (EJM1.GE.0.5) THEN
          WRITE(*,*) 'SHAS1: DT/DY TOO LARGE ! J=',J-1,'  I=',I
          WRITE(*,*) 'TDX=',TDX,'  V(I,J-1)=',V(I,J-1)
          call exit(1)
          END IF
          IF (EIM1.GE.0.5) THEN
          WRITE(*,*) 'SHAS1: DT/DX TOO LARGE ! I=',I-1,'   J=',J
          WRITE(*,*) 'TDX=',TDX,'  U(I-1,J)=',U(I-1,J)
          call exit(1)
          END IF


          QUP=(0.5-EII)/(1.0+(EIP1-EII))
          QUM=(0.5+EII)/(1.0-(EIM1-EII))
          QVP=(0.5-EJI)/(1.0+(EJP1-EJI))
          QVM=(0.5+EJI)/(1.0-(EJM1-EJI))
          AIM1=A(I,J)-A(I-1,J)
          AJM1=A(I,J)-A(I,J-1)
          AI=A(I+1,J)-A(I,J)
          AJ=A(I,J+1)-A(I,J)
          BIM1=ISL1*TDX*(B(I,J)-B(I-1,J)) !ISL1
          BJM1=ISL2*TDY*(B(I,J)-B(I,J-1)) !ISL2
          BI=ISL1*TDX*(B(I+1,J)-B(I,J))
          BJ=ISL2*TDY*(B(I,J+1)-B(I,J))


* This is the transportation step given by Boris and Book
* The twice appearing -0.5*A(I,J) is taking the 2 spacial directions
* into account

*          AB(I,J)=0.5*(QUP*QUP*AI-QUM*QUM*AIM1)
*     &             +QUP*(A(I,J)-BI)+QUM*(A(I,J)-BIM1)
*     &             -0.5*A(I,J)
*     &        +     0.5*(QVP*QVP*AJ-QVM*QVM*AJM1)
*     &             +QVP*(A(I,J)-BJ)+QVM*(A(I,J)-BJM1)
*     &             -0.5*A(I,J)

* The idea for the following transportation step is given
* in Rischkes Paper in Nucl Ph A595(1995)346


           AB(I,J)=0.5*(QUP*QUP*AI-QUM*QUM*AIM1)
     &             +QUP*A(I,J)+QUM*A(I,J)
     &             -0.5*ISL1*TDX*(B(I+1,J)-B(I-1,J))
     &             -0.5*A(I,J)
     &        +     0.5*(QVP*QVP*AJ-QVM*QVM*AJM1)
     &             +QVP*A(I,J)+QVM*A(I,J)
     &             -0.5*ISL2*TDY*(B(I,J+1)-B(I,J-1) )
     &             -0.5*A(I,J)


120    CONTINUE
100   CONTINUE


       IF (EJP1.GE.0.5) THEN
       WRITE(*,*) 'SHAS1: DT/DY TOO LARGE ! J=',J+1
       call exit(1)
       END IF
       IF (EIP1.GE.0.5) THEN
       WRITE(*,*) 'SHAS1: DT/DX TOO LARGE ! I=',I+1
       call exit(1)
       END IF

      RETURN
      END

******************************************************

      SUBROUTINE SHAS2(NXPhy0,NYPhy0,NXPhy,NYPhy,A,U,V,AB,AC1,
     &                    GAB,DIFFC,TDX,TDY,NX0,NY0,NX,NY)
        Implicit Double Precision (A-H, O-Z)
******************************************************
** THIS ROUTINE PRODUCES THE 'ANTIDIFFUSED' SOLUTION *
** A FROM THE 'TRANSPORTED AND DIFFUSED' SOLUTION AB *
** FOR THE SHASTA-FCT ALGORITHM USING ZALESAK:S      *
** METHOD                                            *
** NOTE: FOR INPUT: AB(I,J),U(I,J),V(I,J)            *
**       FOR OUTPUT: A(I,J)   I= 0,..,N; J= 0,..,N   *
** NOTE: IN THIS VERSION THE ANTIDIFFUSION CONSTANT  *
**       'DIFFC' IS FIXED NUMBER GIVEN TO ROUTINE    *
**       AS A PARAMETER.(SO U AND V ARE NOT NEEDED.) *
******************************************************
      DIMENSION AB(NX0:NX,NY0:NY),A(NX0:NX,NY0:NY),
     &           U(NX0:NX,NY0:NY),V(NX0:NX,NY0:NY),
     &         AC1(NX0:NX,NY0:NY),GAB(NX0:NX,NY0:NY)

      COMMON /PARAM/ DX, DY, DT, X0, Y0, ALPHA

      T1 = 0.0
      ACON = 0.0

      DO I=NXPhy0-2,NXPhy+2
       DO J=NYPhy0-2,NYPhy+2
       AC1(I,J)=A(I,J)
       END DO
      END DO



      DO I=NXPhy0,NXPhy
      IP=I-1
      IPP=I+1
      DO J=NYPhy0,NYPhy
        YY=0
        YPY1=0
        YMY1=0
      JP=J-1
      JPP=J+1

*     Here are some remnants to check the influence of the antidiffusion

      DIFF = DIFFC

*      EPSX = ABS(0.5*(U(IP+1,J)+U(IP,J)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AIPIP=DIFF*(AB(IP+1,J)-AB(IP,J))

*      EPSX = ABS(0.5*(U(IP,J)+U(IP-1,J)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AIPIM=DIFF*(AB(IP,J)-AB(IP-1,J))

*      EPSY = ABS(0.5*(V(IP,J+1)+V(IP,J)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AIPJP=DIFF*(AB(IP,J+1)-AB(IP,J))

*      EPSY = ABS(0.5*(V(IP,J)+V(IP,J-1)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AIPJM=DIFF*(AB(IP,J)-AB(IP,J-1))

*      EPSX = ABS(0.5*(U(I+1,J)+U(I,J)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AIJIP=DIFF*(AB(I+1,J)-AB(I,J))

*      EPSX = ABS(0.5*(U(I,J)+U(I-1,J)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AIJIM=DIFF*(AB(I,J)-AB(I-1,J))

*      EPSY = ABS(0.5*(V(I,J+1)+V(I,J)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AIJJP=DIFF*(AB(I,J+1)-AB(I,J))

*      EPSY = ABS(0.5*(V(I,J)+V(I,J-1)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AIJJM=DIFF*(AB(I,J)-AB(I,J-1))

*      EPSX = ABS(0.5*(U(IPP+1,J)+U(IPP,J)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AIPPIP=DIFF*(AB(IPP+1,J)-AB(IPP,J))

*      EPSX =ABS( 0.5*(U(IPP,J)+U(IPP-1,J)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AIPPIM=DIFF*(AB(IPP,J)-AB(IPP-1,J))

*      EPSY = ABS(0.5*(V(IPP,J+1)+V(IPP,J)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AIPPJP=DIFF*(AB(IPP,J+1)-AB(IPP,J))

*      EPSY = ABS(0.5*(V(IPP,J)+V(IPP,J-1)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AIPPJM=DIFF*(AB(IPP,J)-AB(IPP,J-1))

*      EPSX = ABS(0.5*(U(I+1,JP)+U(I,JP)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AJPIP=DIFF*(AB(I+1,JP)-AB(I,JP))

*      EPSX = ABS(0.5*(U(I,JP)+U(I-1,JP)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AJPIM=DIFF*(AB(I,JP)-AB(I-1,JP))

*      EPSY = ABS(0.5*(V(I,JP+1)+V(I,JP)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AJPJP=DIFF*(AB(I,JP+1)-AB(I,JP))

*     EPSY = ABS(0.5*(V(I,JP)+V(I,JP-1)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AJPJM=DIFF*(AB(I,JP)-AB(I,JP-1))

*      EPSX = ABS(0.5*(U(I+1,JPP)+U(I,JPP)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AJPPIP=DIFF*(AB(I+1,JPP)-AB(I,JPP))

*      EPSX = ABS(0.5*(U(I,JPP)+U(I-1,JPP)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AJPPIM=DIFF*(AB(I,JPP)-AB(I-1,JPP))

*      EPSY = ABS(0.5*(V(I,JPP+1)+V(I,JPP)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AJPPJP=DIFF*(AB(I,JPP+1)-AB(I,JPP))

*      EPSY = ABS(0.5*(V(I,JPP)+V(I,JPP-1)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AJPPJM=DIFF*(AB(I,JPP)-AB(I,JPP-1))


C **** SEE PG.349 EQ-14 IN ZALESAK FOR THE BELOW CONDITIONS ********

      PRPI=AIJIP*(AB(I+1,J)-AB(I,J))
      PRPIPP=AIJIP*(AB(I+2,J)-AB(I+1,J))
      PRPIP=AIJIP*(AB(I,J)-AB(I-1,J))
      PRMI=AIJIM*(AB(I+1,J)-AB(I,J))
      PRMIP=AIJIM*(AB(I,J)-AB(I-1,J))
      PRMIP1=AIJIM*(AB(I-1,J)-AB(I-2,J))

      PRPJ=AIJJP*(AB(I,J+1)-AB(I,J))
      PRPJPP=AIJJP*(AB(I,J+2)-AB(I,J+1))
      PRPJP=AIJJP*(AB(I,J)-AB(I,J-1))
      PRMJ=AIJJM*(AB(I,J+1)-AB(I,J))
      PRMJP=AIJJM*(AB(I,J)-AB(I,J-1))
      PRMJP1=AIJJM*(AB(I,J-1)-AB(I,J-2))



      UIPJ=dMAX1(AC1(IP,J),GAB(IP,J))
      UIP1J=dMAX1(AC1(IP-1,J),GAB(IP-1,J))
      UIJ=dMAX1(AC1(I,J),GAB(I,J))
      UIPPJ=dMAX1(AC1(IPP,J),GAB(IPP,J))
      UIPP1J=dMAX1(AC1(IPP+1,J),GAB(IPP+1,J))
      UIJP=dMAX1(AC1(I,JP),GAB(I,JP))
      UIJP1=dMAX1(AC1(I,JP-1),GAB(I,JP-1))
      UIJPP=dMAX1(AC1(I,JPP),GAB(I,JPP))
      UIJPP1=dMAX1(AC1(I,JPP+1),GAB(I,JPP+1))
      UIPJP=dMAX1(AC1(IP,JP),GAB(IP,JP))
      UIPJPP=dMAX1(AC1(IP,JPP),GAB(IP,JPP))
      UIPPJP=dMAX1(AC1(IPP,JP),GAB(IPP,JP))
      UIPPJPP=dMAX1(AC1(IPP,JPP),GAB(IPP,JPP))

      VIPJ=dMIN1(AC1(IP,J),GAB(IP,J))
      VIP1J=dMIN1(AC1(IP-1,J),GAB(IP-1,J))
      VIJ=dMIN1(AC1(I,J),GAB(I,J))
      VIPPJ=dMIN1(AC1(IPP,J),GAB(IPP,J))
      VIPP1J=dMIN1(AC1(IPP+1,J),GAB(IPP+1,J))
      VIJP=dMIN1(AC1(I,JP),GAB(I,JP))
      VIJP1=dMIN1(AC1(I,JP-1),GAB(I,JP-1))
      VIJPP=dMIN1(AC1(I,JPP),GAB(I,JPP))
      VIJPP1=dMIN1(AC1(I,JPP+1),GAB(I,JPP+1))
      VIPJP=dMIN1(AC1(IP,JP),GAB(IP,JP))
      VIPJPP=dMIN1(AC1(IP,JPP),GAB(IP,JPP))
      VIPPJP=dMIN1(AC1(IPP,JP),GAB(IPP,JP))
      VIPPJPP=dMIN1(AC1(IPP,JPP),GAB(IPP,JPP))

      WMXIPJ=dMAX1(UIP1J,UIPJ,UIJ,UIPJP,UIPJPP)
      WMXIJ=dMAX1(UIPJ,UIJ,UIPPJ,UIJP,UIJPP)
      WMXIPPJ=dMAX1(UIJ,UIPPJ,UIPP1J,UIPPJP,UIPPJPP)
      WMXIJP=dMAX1(UIPJP,UIJP,UIPPJP,UIJP1,UIJ)
      WMXIJPP=dMAX1(UIPJPP,UIJPP,UIPPJPP,UIJ,UIJPP1)



      WMNIPJ=dMIN1(VIP1J,VIPJ,VIJ,VIPJP,VIPJPP)
      WMNIJ=dMIN1(VIPJ,VIJ,VIPPJ,VIJP,VIJPP)
      WMNIPPJ=dMIN1(VIJ,VIPPJ,VIPP1J,VIPPJP,VIPPJPP)
      WMNIJP=dMIN1(VIPJP,VIJP,VIPPJP,VIJP1,VIJ)
      WMNIJPP=dMIN1(VIPJPP,VIJPP,VIPPJPP,VIJ,VIJPP1)

      PIPJP=dMAX1(0.0d0,AIPIM)-dMIN1(0.0d0,AIPIP)
     &     +dMAX1(0.0d0,AIPJM)-dMIN1(0.0d0,AIPJP)
      QIPJP=(WMXIPJ-GAB(IP,J))
      IF (PIPJP.GT.0.0) THEN
        RIPJP=dMIN1(1.0d0,QIPJP/PIPJP)
      ELSE
        RIPJP=0.0
      END IF

      PIJP=dMAX1(0.0d0,AIJIM)-dMIN1(0.0d0,AIJIP)
     &    +dMAX1(0.0d0,AIJJM)-dMIN1(0.0d0,AIJJP)
      QIJP=(WMXIJ-GAB(I,J))
      IF (PIJP.GT.0.0) THEN
        RIJP=dMIN1(1.0d0,QIJP/PIJP)
      ELSE
        RIJP=0.0
      END IF

      PIPPJP=dMAX1(0.0d0,AIPPIM)-dMIN1(0.0d0,AIPPIP)
     &      +dMAX1(0.0d0,AIPPJM)-dMIN1(0.0d0,AIPPJP)
      QIPPJP=(WMXIPPJ-GAB(IPP,J))
      IF (PIPPJP.GT.0.0) THEN
        RIPPJP=dMIN1(1.0d0,QIPPJP/PIPPJP)
      ELSE
        RIPPJP=0.0
      END IF

      PIJPP=dMAX1(0.0d0,AJPIM)-dMIN1(0.0d0,AJPIP)
     &     +dMAX1(0.0d0,AJPJM)-dMIN1(0.0d0,AJPJP)
      QIJPP=(WMXIJP-GAB(I,JP))
      IF (PIJPP.GT.0.0) THEN
        RIJPP=dMIN1(1.0d0,QIJPP/PIJPP)
      ELSE
        RIJPP=0.0
      END IF

      PIJPPP=dMAX1(0.0d0,AJPPIM)-dMIN1(0.0d0,AJPPIP)
     &      +dMAX1(0.0d0,AJPPJM)-dMIN1(0.0d0,AJPPJP)
      QIJPPP=(WMXIJPP-GAB(I,JPP))
      IF (PIJPPP.GT.0.0) THEN
        RIJPPP=dMIN1(1.0d0,QIJPPP/PIJPPP)
      ELSE
        RIJPPP=0.0
      END IF

      PIPJM=dMAX1(0.0d0,AIPIP)-dMIN1(0.0d0,AIPIM)
     &     +dMAX1(0.0d0,AIPJP)-dMIN1(0.0d0,AIPJM)
      QIPJM=(GAB(IP,J)-WMNIPJ)
      IF (PIPJM.GT.0.0) THEN
        RIPJM=dMIN1(1.0d0,QIPJM/PIPJM)
      ELSE
        RIPJM=0.0
      END IF

      PIJM=DMAX1(0.0d0,AIJIP)-dMIN1(0.0d0,AIJIM)
     &    +DMAX1(0.0d0,AIJJP)-DMIN1(0.0d0,AIJJM)
      QIJM=(GAB(I,J)-WMNIJ)
      IF (PIJM.GT.0.0) THEN
        RIJM=DMIN1(1.0D0,QIJM/PIJM)
      ELSE
        RIJM=0.0
      END IF

      PIPPJM=DMAX1(0.0D0,AIPPIP)-DMIN1(0.0D0,AIPPIM)
     &      +DMAX1(0.0D0,AIPPJP)-dMIN1(0.0D0,AIPPJM)
      QIPPJM=(GAB(IPP,J)-WMNIPPJ)
      IF (PIPPJM.GT.0.0) THEN
        RIPPJM=DMIN1(1.0D0,QIPPJM/PIPPJM)
      ELSE
        RIPPJM=0.0
      END IF

      PIJPM=DMAX1(0.0D0,AJPIP)-DMIN1(0.0D0,AJPIM)
     &     +DMAX1(0.0D0,AJPJP)-DMIN1(0.0D0,AJPJM)
      QIJPM=(GAB(I,JP)-WMNIJP)
      IF (PIJPM.GT.0.0) THEN
        RIJPM=DMIN1(1.0D0,QIJPM/PIJPM)
      ELSE
        RIJPM=0.0
      END IF

      PIJPPM=DMAX1(0.0D0,AJPPIP)-DMIN1(0.0D0,AJPPIM)
     &      +DMAX1(0.0D0,AJPPJP)-DMIN1(0.0D0,AJPPJM)
      QIJPPM=(GAB(I,JPP)-WMNIJPP)
      IF (PIJPPM.GT.0.0) THEN
        RIJPPM=DMIN1(1.0D0,QIJPPM/PIJPPM)
      ELSE
        RIJPPM=0.0
      END IF

      IF (AIJIP.GE.0.0) THEN
       CIJIP=DMIN1(RIPPJP,RIJM)
      ELSE
       CIJIP=DMIN1(RIJP,RIPPJM)
      END IF

      IF (AIJIM.GE.0.0) THEN
       CIJIM=DMIN1(RIJP,RIPJM)
      ELSE
       CIJIM=DMIN1(RIPJP,RIJM)
      END IF

      IF (AIJJP.GE.0.0) THEN
       CIJJP=DMIN1(RIJPPP,RIJM)
      ELSE
       CIJJP=DMIN1(RIJP,RIJPPM)
      END IF

      IF (AIJJM.GE.0.0) THEN
       CIJJM=DMIN1(RIJP,RIJPM)
      ELSE
       CIJJM=DMIN1(RIJPP,RIJM)
      END IF


      ACIJIP=CIJIP*AIJIP
      ACIJJP=CIJJP*AIJJP
      test1= 0.5*(0.125-CIJJP*DIFFC)*YPY1*SIJJP
*      ACIJJP=CIJJP*AIJJP-0.5*(0.125-CIJJP*DIFFC)*YPY1*SIJJP
      ACIJIM=CIJIM*AIJIM
      ACIJJM=CIJJM*AIJJM
      test2= 0.5*(0.125-CIJJM*DIFFC)*YMY1*SIJJM
*      ACIJJM=CIJJM*AIJJM+0.5*(0.125-CIJJM*DIFFC)*YMY1*SIJJM

*      IF (test1.NE.0.) WRITE(*,*) test1
*      IF (test2.NE.0.) WRITE(*,*) test2

      A(I,J)=AB(I,J)-(ACIJIP-ACIJIM
     &         +ACIJJP-ACIJJM)


      END DO
      END DO

      RETURN
      END


*******************************************************************
      SUBROUTINE P4(I,J,NDX,NDY,NDT,VMID,F0,F1,
     &    NX0,NY0,NX,NY,DT,DX,DY,F)

************************************************************
** THIS ROUTINE PERFORMS A THREE DIMENSIONAL INTERPOLATION
** USING DATA IN MATRICES F0 AND F1.
** Assumes NDX,NDY,NDT are all positive
************************************************************
      Implicit Double Precision (A-H, O-Z)
      Double Precision F0(NX0:NX,NY0:NY),F1(NX0:NX,NY0:NY)
      Double Precision V000,V100,V010,V001,V101,V011,V110,V111
      Double Precision VMID(0:2)
      Double Precision F

      Double Precision x,y,z

      x=VMID(1)/DX
      y=VMID(2)/DY
      z=VMID(0)/DT

      V000=F0(I,J)
      V100=F0(I+NDX,J)
      V010=F0(I,J+NDY)
      V110=F0(I+NDX,J+NDY)

      V001=F1(I,J)
      V101=F1(I+NDX,J)
      V011=F1(I,J+NDY)
      V111=F1(I+NDX,J+NDY)

      F = V000*(1-x)*(1-y)*(1-z) + V100*x*(1-y)*(1-z) +
     &    V010*(1-x)*y*(1-z) + V001*(1-x)*(1-y)*z +
     &    V101*x*(1-y)*z + V011*(1-x)*y*z +
     &    V110*x*y*(1-z) + V111*x*y*z

      RETURN
      END


*******************************************************************
      SUBROUTINE SECTIONCornelius(E0,INTERSECT,ECUBE,EPS0,EPS1,
     &                            I,J,NDX,NDY,NX0,NY0,NX,NY)

*****************************************************************
** THIS ROUTINE CHECKS WHETHER THE eps=E0 -SURFACE INTERSECTS   *
** WITH A 'CUBE' (I,J). IF THAT IS THE CASE, THE LOGICAL        *
** VARIABLE 'INTERSECT' IS RETURNED WITH A VALUE .TRUE. AND     *
** THE MATRIX ECUBE(2,2,2) WITH THE eps VALUES OF THE CUBE      *
** CORNERS for Cornelius routine.                               *
*****************************************************************
      Implicit Double Precision (A-H, O-Z)
      DIMENSION EPS0(NX0:NX, NY0:NY),EPS1(NX0:NX, NY0:NY)
      DIMENSION ECUBE(0:1,0:1,0:1)
      LOGICAL INTERSECT

      IF ((E0-EPS0(I,J))*(EPS1(I+NDX,J+NDY)-E0).LT.0.0) THEN
        IF ((E0-EPS0(I+NDX,J))*(EPS1(I,J+NDY)-E0).LT.0.0) THEN
          IF ((E0-EPS0(I+NDX,J+NDY))*(EPS1(I,J)-E0).LT.0.0) THEN
            IF ((E0-EPS0(I,J+NDY))*(EPS1(I+NDX,J)-E0).LT.0.0) THEN
              INTERSECT=.FALSE.
              GO TO 200
            END IF
          END IF
        END IF
      END IF

      INTERSECT=.TRUE.
      ECUBE(0,0,0)=EPS0(I,J)
      ECUBE(0,1,0)=EPS0(I+NDX,J)
      ECUBE(0,1,1)=EPS0(I+NDX,J+NDY)
      ECUBE(0,0,1)=EPS0(I,J+NDY)
      ECUBE(1,0,0)=EPS1(I,J)
      ECUBE(1,1,0)=EPS1(I+NDX,J)
      ECUBE(1,1,1)=EPS1(I+NDX,J+NDY)
      ECUBE(1,0,1)=EPS1(I,J+NDY)
200   continue

      RETURN
      END


C####################################################################
      Subroutine ASM2D6(AA,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
        Implicit Double Precision (A-H, O-Z)
       Dimension AA(NX0:NX, NY0:NY, NZ0:NZ)
       Dimension IAA(NX0:NX, NY0:NY, NZ0:NZ)
       Dimension FAA(NX0-3:NX+3, NY0-3:NY+3, NZ0:NZ)
        C0=0.486/1.514
        C1=0.343/1.514
        C2=-0.086/1.514

         do 100 k=NZ0,NZ
         do 100 i=NXPhy0-3,NXPhy+3
         do 100 j=NYPhy0-3,NYPhy+3
        FAA(i,j,k)=AA(i,j,k)
 100   continue
       do 300 i=NXPhy0-2,NXPhy+2
       do 300 j=NYPhy0-2,NYPhy+2
c       If(IAA(i,j,NZ0).eq.0) then
         AA(i,j,NZ0)=C0*FAA(i,j,NZ0) +C1*(FAA(i-1,j,NZ0)+FAA(i+1,j,NZ0)
     &     +FAA(i,j-1,NZ0)+FAA(i,j+1,NZ0)) +C2*(FAA(i-1,j-1,NZ0)
     &     +FAA(i+1,j+1,NZ0)+FAA(i+1,j-1,NZ0)+FAA(i-1,j+1,NZ0))
c       else
        AA(i,j,NZ0)=AA(i,j,NZ0)
c       end if
 300   continue
      Return
      End
C####################################################################
      Subroutine CoefASM2d(CofAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
      Implicit Double Precision (A-H, O-Z)
      Dimension CofAA(0:2,NX0:NX, NY0:NY, NZ0:NZ)

      double precision :: DX, DY
      common /DXY/ DX, DY

      Common /R0Bdry/ R0Bdry
      Double Precision R0Bdry, AepsBdry

      AepsBdry = 0.5
      !Print *, "R0Bdry=", R0Bdry

        C0=0.486d0/1.514d0
        C1=0.343d0/1.514d0
        C2=-0.086d0/1.514d0

       do 300 i=NXPhy0-2,NXPhy+2
       do 300 j=NYPhy0-2,NYPhy+2
         xx=DX*I
         yy=DY*J
         rr=sqrt(xx**2+yy**2)
         ff=1.0/(Dexp((rr-R0Bdry)/AepsBdry)+1.0)

        C0a=(1.0-C0)*ff+C0
        C1a=C1*(1-ff)
        C2a=C2*(1-ff)

        CofAA(0,I,J,NZ0)=C0a
        CofAA(1,I,J,NZ0)=C1a
        CofAA(2,I,J,NZ0)=C2a

300   continue
      Return
      End



C####################################################################
      Subroutine ASM2D6bsm(AA,CofAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
        Implicit Double Precision (A-H, O-Z)
       Dimension AA(NX0:NX, NY0:NY, NZ0:NZ)
       Dimension FAA(NX0-3:NX+3, NY0-3:NY+3, NZ0:NZ)
       Dimension CofAA(0:2,NX0:NX, NY0:NY, NZ0:NZ)

         do 100 k=NZ0,NZ
         do 100 i=NXPhy0-3,NXPhy+3
         do 100 j=NYPhy0-3,NYPhy+3
        FAA(i,j,k)=AA(i,j,k)
 100   continue

       do 300 i=NXPhy0-2,NXPhy+2
       do 300 j=NYPhy0-2,NYPhy+2

            C0=CofAA(0,I,J,NZ0)
            C1=CofAA(1,I,J,NZ0)
            C2=CofAA(2,I,J,NZ0)
          AA(i,j,NZ0)=C0*FAA(i,j,NZ0) +C1*(FAA(i-1,j,NZ0)+FAA(i+1,j,NZ0)
     &    +FAA(i,j-1,NZ0)+FAA(i,j+1,NZ0)) +C2*(FAA(i-1,j-1,NZ0)
     &    +FAA(i+1,j+1,NZ0)+FAA(i+1,j-1,NZ0)+FAA(i-1,j+1,NZ0))

300   continue
      Return
      End


C###############################################################################################
C###############################################################################################
      Subroutine VSBdary3(Vx,Vy, NX0,NY0,NZ0,NX,NY,NZ,
     &  NXPhy0,NYPhy0, NXPhy,NYPhy,AAC0)
        Implicit Double Precision (A-H, O-Z)
        Dimension Vx(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension Vy(NX0:NX, NY0:NY, NZ0:NZ)
        LOGICAL V1FOUND,V2FOUND,V3FOUND,V4FOUND
C-----------Boundary Treatment-for Velocity ------------------

      DO 1800 K=NZ0,NZ
      DO 1800 I=NXPhy0,NXPhy
        do kk=1,3
        Vx(I,NYPhy+kk,K)=2.0*Vx(I,NYPhy,K)-Vx(I,NYPhy-kk,K)
        Vx(I,NYPhy0-kk,K)=2.0*Vx(I,NYPhy0,K)-Vx(I,NYPhy0+kk,K)

        Vy(I,NYPhy+kk,K)=2.0*Vy(I,NYPhy,K)-Vy(I,NYPhy-kk,K)
        Vy(I,NYPhy0-kk,K)=2.0*Vy(I,NYPhy0,K)-Vy(I,NYPhy0+kk,K)
        end do
1800  CONTINUE

      DO 1810 K=NZ0,NZ
      DO 1810 J=NYPhy0-3,NYPhy+3
        do kk=1,3
        Vx(NXPhy+kk,J,K)=2.0*Vx(NXPhy,J,K)-Vx(NXPhy-kk,J,K)
        Vx(NXPhy0-kk,J,K)=2.0*Vx(NXPhy0,J,K)-Vx(NXPhy0+kk,J,K)

        Vy(NXPhy+kk,J,K)=2.0*Vy(NXPhy,J,K)-Vy(NXPhy-kk,J,K)
        Vy(NXPhy0-kk,J,K)=2.0*Vy(NXPhy0,J,K)-Vy(NXPhy0+kk,J,K)
        end do
1810  CONTINUE

      DO 1820 K=NZ0,NZ
      DO 1820 I=NXPhy0-3,NXPhy-3
        do kk=1,3
        Vx(I,NYPhy+kk,K)=(2.0*Vx(I,NYPhy,K)-Vx(I,NYPhy-kk,K))*0.5
     &                                      +0.5*Vx(I,NYPhy+kk,K)
        Vx(I,NYPhy0-kk,K)=(2.0*Vx(I,NYPhy0,K)-Vx(I,NYPhy0+kk,K))*0.5
     &                                    +0.5*Vx(I,NYPhy0-kk,K)

        Vy(I,NYPhy+kk,K)=(2.0*Vy(I,NYPhy,K)-Vy(I,NYPhy-kk,K))*0.5
     &                                    +0.5*Vy(I,NYPhy+kk,K)
        Vy(I,NYPhy0-kk,K)=(2.0*Vy(I,NYPhy0,K)-Vy(I,NYPhy0+kk,K))*0.5
     &                                    +0.5*Vy(I,NYPhy0-kk,K)
        end do
1820  CONTINUE

          DO 3310 K=NZ0,NZ
          DO 3310 I=NXPhy0-3,NXPhy+3
            V1FOUND=.TRUE.
            V2FOUND=.TRUE.
            V3FOUND=.TRUE.
            V4FOUND=.TRUE.

            DO 3311 J=NYPhy-1,0,-1
             AVy=Vy(I,J,K)
      IF (V1FOUND.AND.ABS(AVy).GT.AAC0) THEN
        V1FOUND=.FALSE.
          Vy(I,J+1,K)=2.0*Vy(I,J,K)-Vy(I,J-1,K)
        Vy(I,J+2,K)=2.0*Vy(I,J,K)-Vy(I,J-2,K)
        Vy(I,J+3,K)=2.0*Vy(I,J,K)-Vy(I,J-3,K)
        Vy(I,J+4,K)=2.0*Vy(I,J,K)-Vy(I,J-4,K)
      ENDIF
 3311     CONTINUE

            DO 3312 J=NYPhy0+1,0,+1
             AVy=Vy(I,J,K)
      IF (V2FOUND.AND.ABS(AVy).GT.AAC0) THEN
        V2FOUND=.FALSE.
          Vy(I,J-1,K)=2.0*Vy(I,J,K)-Vy(I,J+1,K)
        Vy(I,J-2,K)=2.0*Vy(I,J,K)-Vy(I,J+2,K)
        Vy(I,J-3,K)=2.0*Vy(I,J,K)-Vy(I,J+3,K)
        Vy(I,J-4,K)=2.0*Vy(I,J,K)-Vy(I,J+4,K)
      ENDIF
 3312     CONTINUE

 3310     CONTINUE

        DO 3320 K=NZ0,NZ
        DO 3320 J=NYPhy0-3,NYPhy+3

          V1FOUND=.TRUE.
          V2FOUND=.TRUE.
          V3FOUND=.TRUE.
          V4FOUND=.TRUE.

          DO 3321 I=NXPhy-1,0,-1
             AVx=Vx(I,J,K)
          IF (V3FOUND.AND.ABS(AVx).GT.AAC0) THEN
      V3FOUND=.FALSE.
C            Print *,I,J, AAC0,AVx
      Vx(I+1,J,K)=2.0*Vx(I,J,K)-Vx(I-1,J,K)
      Vx(I+2,J,K)=2.0*Vx(I,J,K)-Vx(I-2,J,K)
      Vx(I+3,J,K)=2.0*Vx(I,J,K)-Vx(I-3,J,K)
      Vx(I+4,J,K)=2.0*Vx(I,J,K)-Vx(I-4,J,K)
          END IF
 3321  CONTINUE

          DO 3322 I=NXPhy0+1,0,+1
             AVx=Vx(I,J,K)
          IF (V4FOUND.AND.ABS(AVx).GT.AAC0) THEN
      V4FOUND=.FALSE.
      Vx(I-1,J,K)=2.0*Vx(I,J,K)-Vx(I+1,J,K)
      Vx(I-2,J,K)=2.0*Vx(I,J,K)-Vx(I+2,J,K)
      Vx(I-3,J,K)=2.0*Vx(I,J,K)-Vx(I+3,J,K)
      Vx(I-4,J,K)=2.0*Vx(I,J,K)-Vx(I+4,J,K)
          END IF
 3322  CONTINUE

 3320  CONTINUE

      Return
      End





C###############################################################################################
C###############################################################################################
      Subroutine VSBdary3Single(VV,NX0,NY0,NZ0,NX,NY,NZ,
     &  NXPhy0,NYPhy0, NXPhy,NYPhy,AAC0)
        Implicit Double Precision (A-H, O-Z)
        Dimension VV(NX0:NX, NY0:NY, NZ0:NZ)
C-----------Boundary Treatment-for Velocity ------------------

      DO 1800 K=NZ0,NZ
      DO 1800 I=NXPhy0,NXPhy
        do kk=1,3
        VV(I,NYPhy+kk,K)=2.0*VV(I,NYPhy,K)-VV(I,NYPhy-kk,K)
        VV(I,NYPhy0-kk,K)=2.0*VV(I,NYPhy0,K)-VV(I,NYPhy0+kk,K)
        end do
1800  CONTINUE

      DO 1810 K=NZ0,NZ
      DO 1810 J=NYPhy0-3,NYPhy+3
        do kk=1,3
        VV(NXPhy+kk,J,K)=2.0*VV(NXPhy,J,K)-VV(NXPhy-kk,J,K)
        VV(NXPhy0-kk,J,K)=2.0*VV(NXPhy0,J,K)-VV(NXPhy0+kk,J,K)
        end do
1810  CONTINUE

      DO 1820 K=NZ0,NZ
      DO 1820 I=NXPhy0-3,NXPhy-3
        do kk=1,3
        VV(I,NYPhy+kk,K)=(2.0*VV(I,NYPhy,K)-VV(I,NYPhy-kk,K))*0.5
     &                                      +0.5*VV(I,NYPhy+kk,K)
        VV(I,NYPhy0-kk,K)=(2.0*VV(I,NYPhy0,K)-VV(I,NYPhy0+kk,K))*0.5
     &                                    +0.5*VV(I,NYPhy0-kk,K)
        end do
1820  CONTINUE

      Return
      End




C########################################################################
      Subroutine Sembdary3(PL,NX0,NY0,NZ0,NX,NY,NZ,
     &                         NXPhy0,NYPhy0,NXPhy,NYPhy)
        Implicit Double Precision (A-H, O-Z)
        Dimension PL(NX0:NX, NY0:NY, NZ0:NZ)
C-----------Boundary Treatment-------------------
      DO 1800 K=NZ0,NZ
      DO 1800 I=NXPhy0,NXPhy
        do kk=1,3
        PL(I,NYPhy+kk,K)=2.0*PL(I,NYPhy,K)-PL(I,NYPhy-kk,K)
        PL(I,NYPhy0-kk,K)=2.0*PL(I,NYPhy0,K)-PL(I,NYPhy0+kk,K)
        end do
1800  CONTINUE

      DO 1810 K=NZ0,NZ
      DO 1810 J=NYPhy0-3,NYPhy+3
        do kk=1,3
        PL(NXPhy+kk,J,K)=2.0*PL(NXPhy,J,K)-PL(NXPhy-kk,J,K)
        PL(NXPhy0-kk,J,K)=2.0*PL(NXPhy0,J,K)-PL(NXPhy0+kk,J,K)
        end do
1810  CONTINUE

      DO 1820 K=NZ0,NZ
      DO 1820 I=NXPhy0-3,NXPhy-3
        do kk=1,3
        PL(I,NYPhy+kk,K)=(2.0*PL(I,NYPhy,K)-PL(I,NYPhy-kk,K))*0.5
     &                                    +0.5*PL(I,NYPhy+kk,K)
        PL(I,NYPhy0-kk,K)=(2.0*PL(I,NYPhy0,K)-PL(I,NYPhy0+kk,K))*0.5
     &                                    +0.5*PL(I,NYPhy0-kk,K)
        end do
1820  CONTINUE

      Return
      End


C########################################################################
      Subroutine TriSembdary3(PL1, PL2, PL3,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
        Implicit Double Precision (A-H, O-Z)
        Dimension PL1(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension PL2(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension PL3(NX0:NX, NY0:NY, NZ0:NZ)
C-----------Boundary Treatment-------------------

      DO 1800 K=NZ0,NZ
      DO 1800 I=NXPhy0,NXPhy
        do kk=1,3
        PL1(I,NYPhy+kk,K)=2.0*PL1(I,NYPhy,K)-PL1(I,NYPhy-kk,K)
        PL1(I,NYPhy0-kk,K)=2.0*PL1(I,NYPhy0,K)-PL1(I,NYPhy0+kk,K)

        PL2(I,NYPhy+kk,K)=2.0*PL2(I,NYPhy,K)-PL2(I,NYPhy-kk,K)
        PL2(I,NYPhy0-kk,K)=2.0*PL2(I,NYPhy0,K)-PL2(I,NYPhy0+kk,K)

        PL3(I,NYPhy+kk,K)=2.0*PL3(I,NYPhy,K)-PL3(I,NYPhy-kk,K)
        PL3(I,NYPhy0-kk,K)=2.0*PL3(I,NYPhy0,K)-PL3(I,NYPhy0+kk,K)
        end do
1800  CONTINUE

      DO 1810 K=NZ0,NZ
      DO 1810 J=NYPhy0-3,NYPhy+3
        do kk=1,3
        PL1(NXPhy+kk,J,K)=2.0*PL1(NXPhy,J,K)-PL1(NXPhy-kk,J,K)
        PL1(NXPhy0-kk,J,K)=2.0*PL1(NXPhy0,J,K)-PL1(NXPhy0+kk,J,K)

        PL2(NXPhy+kk,J,K)=2.0*PL2(NXPhy,J,K)-PL2(NXPhy-kk,J,K)
        PL2(NXPhy0-kk,J,K)=2.0*PL2(NXPhy0,J,K)-PL2(NXPhy0+kk,J,K)

        PL3(NXPhy+kk,J,K)=2.0*PL3(NXPhy,J,K)-PL3(NXPhy-kk,J,K)
        PL3(NXPhy0-kk,J,K)=2.0*PL3(NXPhy0,J,K)-PL3(NXPhy0+kk,J,K)
        end do
1810  CONTINUE

      DO 1820 K=NZ0,NZ
      DO 1820 I=NXPhy0-3,NXPhy-3
        do kk=1,3
        PL1(I,NYPhy+kk,K)=(2.0*PL1(I,NYPhy,K)-PL1(I,NYPhy-kk,K))*0.5
     &                                    +0.5*PL1(I,NYPhy+kk,K)
        PL1(I,NYPhy0-kk,K)=(2.0*PL1(I,NYPhy0,K)-PL1(I,NYPhy0+kk,K))*0.5
     &                                      +0.5*PL1(I,NYPhy0-kk,K)

        PL2(I,NYPhy+kk,K)=(2.0*PL2(I,NYPhy,K)-PL2(I,NYPhy-kk,K))*0.5
     &                                    +0.5*PL2(I,NYPhy+kk,K)
        PL2(I,NYPhy0-kk,K)=(2.0*PL2(I,NYPhy0,K)-PL2(I,NYPhy0+kk,K))*0.5
     &                                      +0.5*PL2(I,NYPhy0-kk,K)

        PL3(I,NYPhy+kk,K)=(2.0*PL3(I,NYPhy,K)-PL3(I,NYPhy-kk,K))*0.5
     &                                    +0.5*PL3(I,NYPhy+kk,K)
        PL3(I,NYPhy0-kk,K)=(2.0*PL3(I,NYPhy0,K)-PL3(I,NYPhy0+kk,K))*0.5
     &                                      +0.5*PL3(I,NYPhy0-kk,K)
        end do
1820  CONTINUE

      Return
      End
