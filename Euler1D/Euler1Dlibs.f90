!######################################!
!#           LIBS Library             #!
!# Contains subroutines and functions #!
!# used in the Euler1D program.       #!
!######################################!

!############################!
!# Loads Initial Conditions #!
!############################!
SUBROUTINE INITCOND(ICS,NX,P,P0,NT)

  INTEGER, INTENT(IN) :: ICS
  INTEGER, INTENT(IN) :: NX
  REAL*8, INTENT(OUT) :: P(NX,3)
  REAL*8, INTENT(OUT) :: P0(8)
  REAL*8, INTENT(OUT) :: NT

  REAL*8 :: X0=0, UL0=0, UR0=0, DL0=0, DR0=0, PL0=0, PR0=0
  INTEGER :: MID

  IF (ICS==0) THEN
    ! Custom Initial Conditions
    ! The vectors of PRIMITIVES must be initialized here
  ELSE
    SELECT CASE (ICS)
    ! Test 1 of Toro
    CASE (1)
      X0 = 0.3
      UL0 = 0.75
      UR0 = 0.0
      DL0 = 1.0
      DR0 = 0.125
      PL0 = 1.0
      PR0 = 0.1
      NT = 0.2
    ! Test 2 of Toro
    CASE (2)
      X0 = 0.5
      UL0 = -2.0
      UR0 = 2.0
      DL0 = 1.0
      DR0 = 1.0
      PL0 = 0.4
      PR0 = 0.4
      NT = 0.15
    ! Test 3 of Toro
    CASE (3)
      X0 = 0.5
      UL0 = 0.0
      UR0 = 0.0
      DL0 = 1.0
      DR0 = 1.0
      PL0 = 1000.0
      PR0 = 0.01
      NT = 0.012
    ! Test 4 of Toro
    CASE (4)
      X0 = 0.4
      UL0 = 19.5975
      UR0 = -6.19633
      DL0 = 5.99924
      DR0 = 5.99242
      PL0 = 460.894
      PR0 = 46.0950
      NT = 0.035
    ! Test 5 of Toro
    CASE (5)
      X0 = 0.8
      UL0 = -19.59745
      UR0 = -19.59745
      DL0 = 1.0
      DR0 = 1.0
      PL0 = 1000.0
      PR0 = 0.01
      NT = 0.012
    ! Custom Riemann test
    CASE (6)
      X0 = 0.8
      UL0 = -19.59745
      UR0 = -19.59745
      DL0 = 1.0
      DR0 = 1.0
      PL0 = 1000.0
      PR0 = 0.01
      NT = 1.0
    END SELECT

    ! Set Riemann Problem
    ! Stored in primitives
    MID = X0*NX
    P(1:MID,1) = UL0
    P(MID+1:NX,1) = UR0
    P(1:MID,2) = DL0
    P(MID+1:NX,2) = DR0
    P(1:MID,3) = PL0
    P(MID+1:NX,3) = PR0

    ! Bundle initial conditions
    P0(1) = X0
    P0(2) = UL0
    P0(3) = UR0
    P0(4) = DL0
    P0(5) = DR0
    P0(6) = PL0
    P0(7) = PR0
    P0(8) = NT

  END IF

END SUBROUTINE


!#############################!
!# Calculates integration    #!
!# variables from primitives #!
!#############################!
SUBROUTINE GETVARS(NX,P,GAM,U)

  INTEGER, INTENT(IN) :: NX
  REAL*8, INTENT(IN) :: P(NX,3)
  REAL*8, INTENT(IN) :: GAM
  REAL*8, INTENT(OUT) :: U(NX,3)

  U(:,1) = P(:,2)
  U(:,2) = P(:,2)*P(:,1)
  U(:,3) = 0.5*P(:,2)*P(:,1)**2 + P(:,3)/(GAM-1)

END SUBROUTINE

!##############################!
!# Calculates primitives from #!
!# integration variables      #!
!##############################!
SUBROUTINE GETPRIMS(NX,U,GAM,P)

  INTEGER, INTENT(IN) :: NX
  REAL*8, INTENT(IN) :: U(NX,3)
  REAL*8, INTENT(IN) :: GAM
  REAL*8, INTENT(OUT) :: P(NX,3)

  P(:,1) = U(:,2)/U(:,1)
  P(:,2) = U(:,1)
  P(:,3) = (GAM-1)*(U(:,3)-0.5*U(:,2)**2/U(:,1))

END SUBROUTINE

!#######################!
!# Calculates Timestep #!
!#######################!
SUBROUTINE GETSTEP(P,NX,GAM,DX,CP,DT)

  INTEGER, INTENT(IN) :: NX
  REAL*8, INTENT(IN) :: P(NX,3)
  REAL*8, INTENT(IN) :: GAM
  REAL*8, INTENT(IN) :: DX
  REAL*8, INTENT(IN) :: CP
  REAL*8, INTENT(OUT) :: DT

  REAL*8 :: UMAX, ULOC
  INTEGER :: I

  ! The timestep is calculated following the
  ! Courant-Friedrichs-Lewis criterion:
  ! DT = DX / UMAX * CP
  ! where DX is the grid spacing, CP is the
  ! Courant parameter and UMAX is the maximum
  ! value of local flow speed + local sound speed

  UMAX = 0
  DO I=1,NX
    ULOC = ABS(P(I,1)) + SQRT(GAM*P(I,3)/P(I,2))
    IF (ULOC > UMAX) THEN
      UMAX = ULOC
    END IF
  END DO

  DT = DX / UMAX * CP

END SUBROUTINE

!##############!
!# Lax Method #!
!##############!
SUBROUTINE LAX(NX,U,P,GAM,DT,DX,UP,F)

  INTEGER, INTENT(IN) :: NX
  REAL*8, INTENT(IN) :: U(NX,3)
  REAL*8, INTENT(IN) :: P(NX,3)
  REAL*8, INTENT(IN) :: GAM
  REAL*8, INTENT(IN) :: DT
  REAL*8, INTENT(IN) :: DX
  REAL*8, INTENT(OUT) :: UP(NX,3)
  REAL*8, INTENT(OUT) :: F(NX,3)

  INTEGER :: I

  ! Calculate Fluxes
  F(:,1) = P(:,2)*P(:,1)
  F(:,2) = P(:,2)*P(:,1)**2 + P(:,3)
  F(:,3) = P(:,1)*(0.5*P(:,2)*P(:,1)**2 + GAM/(GAM-1)*P(:,3))

  ! Lax algorithm is applied to non-boundary grid points
  DO I=2,NX-1
    UP(I,:) =  0.5*(U(I+1,:) + U(I-1,:)) - 0.5 * DT/DX * (F(I+1,:) - F(I-1,:))
  END DO

END SUBROUTINE

!#####################!
!# MacCormack Method #!
!#####################!
SUBROUTINE MAC(NX,U,P,GAM,ETA,DT,DX,SOLVER,UP,F)

  INTEGER, INTENT(IN) :: NX
  REAL*8, INTENT(IN) :: U(NX,3)
  REAL*8, INTENT(IN) :: P(NX,3)
  REAL*8, INTENT(IN) :: GAM
  REAL*8, INTENT(IN) :: ETA
  REAL*8, INTENT(IN) :: DT
  REAL*8, INTENT(IN) :: DX
  INTEGER, INTENT(IN) :: SOLVER
  REAL*8, INTENT(OUT) :: UP(NX,3)
  REAL*8, INTENT(OUT) :: F(NX,3)

  INTEGER :: I
  REAL*8 :: UB(NX,3)
  REAL*8 :: FB(NX,3)

  ! Fluxes
  F(:,1) = P(:,2)*P(:,1)
  F(:,2) = P(:,2)*P(:,1)**2 + P(:,3)
  F(:,3) = P(:,1)*(0.5*P(:,2)*P(:,1)**2 + GAM/(GAM-1)*P(:,3))

  ! Predictor
  DO I=2,NX-1
    UB(I,:) = U(I,:) - DT/DX*(F(I+1,:) - F(I,:))
  END DO
  UB(1,:) = UB(2,:)
  UB(NX,:) = UB(NX-1,:)

  ! Update predicted fluxes
  FB(:,1) = UB(:,2)
  FB(:,2) = 0.5*(3-GAM)*UB(:,2)**2/UB(:,1) + (GAM-1)*UB(:,3)
  FB(:,3) = UB(:,2)/UB(:,1)*(GAM*UB(:,3) - 0.5*(GAM-1)*UB(:,2)**2/UB(:,1))

  ! Corrector
  DO I=2,NX-1
    UP(I,:) = 0.5*(UB(I,:) + U(I,:) - DT/DX*(FB(I,:) - FB(I-1,:)))
  END DO

  ! Apply artificial viscosity
  CALL VISCOSITY(NX,SOLVER,ETA,UP)

END SUBROUTINE

!########################!
!# Applies viscosity to #!
!# MacCormack algorithm #!
!########################!
SUBROUTINE VISCOSITY(NX,SOLVER,ETA,UP)

  INTEGER, INTENT(IN) :: NX
  INTEGER, INTENT(IN) :: SOLVER
  REAL*8, INTENT(IN) :: ETA
  REAL*8, INTENT(OUT) :: UP(NX,3)

  INTEGER :: I, J
  REAL*8, DIMENSION(2,3) :: Q

  DO I=2,NX-1
    SELECT CASE (SOLVER)
      CASE (2) ! Uniform vicosity
        UP(I,:) = UP(I,:) + ETA*(UP(I+1,:) + UP(I-1,:) - 2*UP(I,:))
      CASE (3) ! Flux-corrector criterion
        DO J=1,3
          Q(1,J) = UP(I,J) - UP(I-1,J)
          Q(2,J) = UP(I+1,J) - UP(I,J)
          IF (Q(1,J)*Q(2,J) < 0) THEN
            UP(I,J) = UP(I,J) + ETA*(UP(I+1,J) + UP(I-1,J) - 2*UP(I,J))
          END IF
        END DO
      END SELECT
  END DO

END SUBROUTINE

!################################!
!# Simple Godunov upwind scheme #!
!################################!
SUBROUTINE GOD(NX,U,P,GAM,DT,DX,UP,F)

  INTEGER, INTENT(IN) :: NX
  REAL*8, INTENT(IN) :: U(NX,3)
  REAL*8, INTENT(IN) :: P(NX,3)
  REAL*8, INTENT(IN) :: GAM
  REAL*8, INTENT(IN) :: DT
  REAL*8, INTENT(IN) :: DX
  REAL*8, INTENT(OUT) :: UP(NX,3)
  REAL*8, INTENT(OUT) :: F(NX,3)

  INTEGER :: I

  ! Calculate raw fluxes
  F(:,1) = P(:,2)*P(:,1)
  F(:,2) = P(:,2)*P(:,1)**2 + P(:,3)
  F(:,3) = P(:,1)*(0.5*P(:,2)*P(:,1)**2 + GAM/(GAM-1)*P(:,3))

  ! Obtain new integration variables
  DO I=2,NX-1
    UP(I,:) = U(I,:) + DT/DX*(F(I-1,:)-F(I,:))
  END DO

END SUBROUTINE

!############################!
!# Wavespeed estimators     #!
!# for HLL and HLLC solvers #!
!############################!
SUBROUTINE WAVESPEED(NX,P,GAM,WSPD,I,SL,SR,SP,SS)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NX
  REAL*8, INTENT(IN) :: P(NX,3)
  REAL*8, INTENT(IN) :: GAM
  INTEGER, INTENT(IN) :: WSPD
  INTEGER, INTENT(IN) :: I
  REAL*8, INTENT(OUT) :: SL, SR, SP, SS

  REAL*8 :: UL, UR, DL, DR, PL, PR
  REAL*8 :: CL, CR
  REAL*8 :: UA, DA, PA, CA
  REAL*8 :: US, PS, QL, QR
  REAL*8 :: Z, PLR
  REAL*8 :: PPV, GL, GR, P0
  REAL*8 :: AL, AR, BL, BR
  INTEGER :: SUBWSPD
  REAL*8 :: PMIN, PMAX, Q

  ! Flow variables at cells
  ! on left/right of interface
  UL = P(I,1)
  UR = P(I+1,1)
  DL = P(I,2)
  DR = P(I+1,2)
  PL = P(I,3)
  PR = P(I+1,3)

  ! Sound speeds
  CL = SQRT(GAM*PL/DL)
  CR = SQRT(GAM*PR/DR)

  SELECT CASE(WSPD)

  ! Davis direct, 10.37 of Toro
  CASE (1)
    SL = UL-CL
    SR = UR+CR

  ! Davis direct bounded, 10.38 of Toro
  CASE (2)
    SL = MIN( UL-CL, UR-CR )
    SR = MAX( UL+CL, UR+CR )

  ! Bounded single 'S+', 10.43 of Toro
  CASE (3)
    SP = MAX( ABS(UL)+CL, ABS(UR)+CR )

  ! Pressure-based estimates: PVRS, TRRS, TSRS
  CASE (4,5,6,7)

    ! Use Adaptive solver (7) to select the scheme,
    ! or just continue with selected scheme (4-6)
    IF (WSPD==7) THEN
      ! Averages
      CA = 0.5*(CL+CR)
      UA = 0.5*(UL+UR)
      DA = 0.5*(DL+DR)
      PA = 0.5*(PL+PR)
      ! PVRS star pressure
      PS = PA-0.5*(UR-UL)*DA*CA
      ! PMIN and PMAX
      IF (PL<=PR) THEN
        PMIN = PL
        PMAX = PR
      ELSE
        PMIN = PR
        PMAX = PL
      END IF
      ! Calculate Q
      Q = PMAX/PMIN
      ! Select wavespeed estimation method; see fig 9.4 of Toro
      ! This uses Q_user = 2.
      IF ((Q<2.0).AND.(PS>PMIN).AND.(PS<PMAX)) THEN
        ! Choose PVRS
        SUBWSPD = 4
      ELSEIF (PS<PMIN) THEN
        ! Choose TRRS
        SUBWSPD = 5
      ELSE
        ! Choose TSRS
        SUBWSPD = 6
      END IF
    ELSE
      ! Just continue with selected scheme
      SUBWSPD = WSPD
    END IF

    ! Specific part of each algorithm
    SELECT CASE (SUBWSPD)

    ! PVRS pressure-based estimate, 10.51 of Toro
    CASE (4)
      ! Averages
      CA = 0.5*(CL+CR)
      UA = 0.5*(UL+UR)
      DA = 0.5*(DL+DR)
      PA = 0.5*(PL+PR)
      ! Star pressure
      PS = PA-0.5*(UR-UL)*DA*CA
      ! Star velocity
      US = UA-0.5*(PR-PL)/(DA*CA)

    ! TRRS pressure-based estimate, 10.53 of Toro
    CASE (5)
      Z = (GAM-1)/(2*GAM)
      PLR = (PL/PR)**Z
      ! Star pressure
      PS = ((CL+CR-0.5*(GAM-1)*(UR-UL))/(CL/(PL**Z)+CR/(PR**Z)))**(1.0/Z)
      ! Star velocity
      US = (PLR*UL/CL+UR/CR+2.0*(PLR-1)/(GAM-1))/(PLR/CL+1.0/CR)

    ! TSRS pressure-based estimate, 10.55 of Toro
    CASE (6)
      ! Averages
      CA = 0.5*(CL+CR)
      UA = 0.5*(UL+UR)
      DA = 0.5*(DL+DR)
      PA = 0.5*(PL+PR)
      ! As and Bs
      AL = 2.0/(GAM+1)/DL
      AR = 2.0/(GAM+1)/DR
      BL = (GAM-1)/(GAM+1)*PL
      BR = (GAM-1)/(GAM+1)*PR
      ! P0
      PPV = PA-0.5*(UR-UL)*DA*CA
      IF (PPV>0) THEN
        P0 = PPV
      ELSE
        P0 = 0
      END IF
      ! Gs at P0
      GL = (AL/(P0+BL))**(0.5)
      GR = (AR/(P0+BR))**(0.5)
      ! Star pressure
      PS = (GL*PL+GR*PR-(UR-UL))/(GL+GR)
      ! Star velocity
      US = UA+0.5*((PS-PR)*GR-(PS-PL)*GL)

    CASE DEFAULT
      PS = 0
      US = 0

    END SELECT

    ! Calculate Qs
    IF (PS<=PL) THEN
      QL = 1
    ELSE
      QL = (1+(GAM+1)/(2*GAM)*(PS/PL-1))**0.5
    END IF
    IF (PS<=PR) THEN
      QR = 1
    ELSE
      QR = (1+(GAM+1)/(2*GAM)*(PS/PR-1))**0.5
    END IF

    ! Finally, wavespeed estimates
    SL = UL-CL*QL
    SR = UR+CR*QR
    SS = US

  END SELECT

END SUBROUTINE


!###########################!
!# Godunov upwind scheme   #!
!# with HLL Riemann solver #!
!###########################!
SUBROUTINE HLL(NX,U,P,GAM,DT,DX,WSPD,UP,F)

  INTEGER, INTENT(IN) :: NX
  REAL*8, INTENT(IN) :: U(NX,3)
  REAL*8, INTENT(IN) :: P(NX,3)
  REAL*8, INTENT(IN) :: GAM
  REAL*8, INTENT(IN) :: DT
  REAL*8, INTENT(IN) :: DX
  INTEGER, INTENT(IN) :: WSPD
  REAL*8, INTENT(OUT) :: UP(NX,3)
  REAL*8, INTENT(OUT) :: F(NX,3)

  INTEGER :: I
  REAL*8, DIMENSION(NX-1,3) :: FC
  REAL*8, DIMENSION(3) :: UL, UR, FL, FR
  REAL*8 :: SL, SR, SP, SS

  ! Calculate raw fluxes for all cells
  F(:,1) = P(:,2)*P(:,1)
  F(:,2) = P(:,2)*P(:,1)**2 + P(:,3)
  F(:,3) = P(:,1)*(0.5*P(:,2)*P(:,1)**2 + GAM/(GAM-1)*P(:,3))

  ! Obtain intercell fluxes FC(I)
  ! FC(I) is the flux at the interface
  ! between cells I and I+1
  DO I=1,NX-1

    ! Primitive variables and fluxes
    ! on left/right of interface
    UL(:) = U(I,:)
    UR(:) = U(I+1,:)
    FL(:) = F(I,:)
    FR(:) = F(I+1,:)

    ! Obtain wavespeed estimates
    CALL WAVESPEED(NX,P,GAM,WSPD,I,SL,SR,SP,SS)

    ! Obtain FC according to WSPD
    IF (WSPD==3) THEN
      ! Reduced intercell fluxes, 10.42 of Toro
      FC(I,:) = 0.5*(FL(:)+FR(:))-0.5*SP*(UR(:)-UL(:))
    ELSE
      ! Intercell fluxes given by 10.21 of Toro
      ! 'S star' is not needed in HLL, so it is ignored
      IF (SL>=0) THEN
        FC(I,:) = FL(:)
      ELSEIF ((SL<=0).AND.(SR>=0)) THEN
        FC(I,:) = (SR*FL(:)-SL*FR(:)+SL*SR*(UR(:)-UL(:)))/(SR-SL)
      ELSEIF (SR<=0) THEN
        FC(I,:) = FR(:)
      END IF
    END IF

  END DO

  ! Obtain new integration variables
  ! from intercell fluxes; exclude boundaries
  ! Standard Godunov upwind scheme, 10.2 of Toro
  DO I=2,NX-1
    UP(I,:) = U(I,:) + DT/DX*(FC(I-1,:)-FC(I,:))
  END DO

END SUBROUTINE


!##############################!
!# Returns average of A and B #!
!# depending on method AVG    #!
!##############################!
REAL*8 FUNCTION AVERAGE(A,B,AVG)

  REAL*8, INTENT(IN) :: A,B
  INTEGER, INTENT(IN) :: AVG

  REAL*8 :: E

  SELECT CASE (AVG)

  CASE (1)  ! Direct arithmetic average
    AVERAGE = (A+B)/2.0

  CASE (2)  ! Minmod
    IF (ABS(A)<ABS(B)) THEN
      AVERAGE = A
    ELSE
      AVERAGE = B
    END IF

  CASE (3)  ! van Albada
    E = 0.001
    AVERAGE = ((B**2+E**2)*A+(A**2+E**2)*B)/(A**2+B**2+E**2)

  CASE DEFAULT  ! Safety
    AVERAGE = (A+B)/2.0

  END SELECT

END FUNCTION

!#########################!
!# Domain Reconstruction #!
!#########################!
SUBROUTINE RECONSTRUCT()

END SUBROUTINE


!############################!
!# Godunov upwind scheme    #!
!# with HLLC Riemann solver #!
!############################!
SUBROUTINE HLLC(NX,U,P,GAM,DT,DX,SOLVER,WSPD,AVG,UP,F)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NX
  REAL*8, INTENT(IN) :: U(NX,3)
  REAL*8, INTENT(IN) :: P(NX,3)
  REAL*8, INTENT(IN) :: GAM
  REAL*8, INTENT(IN) :: DT
  REAL*8, INTENT(IN) :: DX
  INTEGER, INTENT(IN) :: SOLVER
  INTEGER, INTENT(IN) :: WSPD
  INTEGER, INTENT(IN) :: AVG
  REAL*8, INTENT(OUT) :: UP(NX,3)
  REAL*8, INTENT(OUT) :: F(NX,3)

  REAL*8 :: AVERAGE

  INTEGER :: I, J
  REAL*8, DIMENSION(3) :: FL, FR, UUL, UUR, USL, USR
  REAL*8, DIMENSION(NX-1,3) :: FC
  REAL*8, DIMENSION(NX,3) :: DU, DP
  REAL*8, DIMENSION(3) :: DUL, DUR, DPL, DPR
  REAL*8 :: UL, UR, DL, DR, PL, PR, EL, ER
  REAL*8 :: SL, SR, SP, SS

  ! Calculate raw fluxes for all cells
  F(:,1) = P(:,2)*P(:,1)
  F(:,2) = P(:,2)*P(:,1)**2 + P(:,3)
  F(:,3) = P(:,1)*(0.5*P(:,2)*P(:,1)**2 + GAM/(GAM-1)*P(:,3))

  ! Calculate averaged cell slopes
  ! (Linear Data Reconstruction)
  IF (SOLVER==7) THEN
    DO I=2,NX-1
      ! Slopes at left and right
      DUL(:) = U(I,:)-U(I-1,:)
      DUR(:) = U(I+1,:)-U(I,:)
      DPL(:) = P(I,:)-P(I-1,:)
      DPR(:) = P(I+1,:)-P(I,:)
      ! Cell slope is average of DUL and DUR
      DO J=1,3
        DU(I,J) = AVERAGE(DUL(J),DUR(J),AVG)
        DP(I,J) = AVERAGE(DPL(J),DPR(J),AVG)
      END DO
    END DO
    DU(1,:) = U(2,:)-U(1,:)
    DU(NX,:) = U(NX,:)-U(NX-1,:)
    DP(1,:) = P(2,:)-P(1,:)
    DP(NX,:) = P(NX,:)-P(NX-1,:)
  END IF

  ! Obtain intercell fluxes FC(I)
  ! FC(I) is the flux at the interface
  ! between cells I and I+1
  DO I=1,NX-1

    ! Variables on left/right
    ! of interface
    IF (SOLVER==6) THEN
      ! Piecewise Constant Data
      UL = P(I,1)
      UR = P(I+1,1)
      DL = P(I,2)
      DR = P(I+1,2)
      PL = P(I,3)
      PR = P(I+1,3)
      EL = U(I,3)
      ER = U(I+1,3)
      UUL(:) = U(I,:)
      UUR(:) = U(I+1,:)
      FL(:) = F(I,:)
      FR(:) = F(I+1,:)
    ELSE
      ! Linear Data Reconstruction
      UL = P(I,1)+DP(I,1)/2.0
      UR = P(I+1,1)-DP(I+1,1)/2.0
      DL = P(I,2)+DP(I,2)/2.0
      DR = P(I+1,2)-DP(I+1,2)/2.0
      PL = P(I,3)+DP(I,3)/2.0
      PR = P(I+1,3)-DP(I+1,3)/2.0
      EL = U(I,3)+DU(I,3)/2.0
      ER = U(I+1,3)-DU(I+1,3)/2.0
      UUL(:) = U(I,:)+DU(I,:)/2.0
      UUR(:) = U(I+1,:)-DU(I+1,:)/2.0
      ! Fluxes untouched
      FL(:) = F(I,:)
      FR(:) = F(I+1,:)
    END IF

    ! Obtain wavespeed estimates SL, SR and SS
    SELECT CASE (WSPD)
    CASE (4,5,6,7)
      CALL WAVESPEED(NX,P,GAM,WSPD,I,SL,SR,SP,SS)
    CASE DEFAULT
      WRITE(*,'(A,I1,A)') 'HLLC not compatible with wavespeed estimation method ',WSPD,'; use one of 4,5,6,7.'
      STOP
    END SELECT

    ! Compute star region states, 10.33 of Toro
    USL(1) = DL*(SL-UL)/(SL-SS)
    USL(2) = DL*(SL-UL)/(SL-SS)*SS
    USL(3) = DL*(SL-UL)/(SL-SS)*(EL/DL+(SS-UL)*(SS+PL/(DL*(SL-UL))))
    USR(1) = DR*(SR-UR)/(SR-SS)
    USR(2) = DR*(SR-UR)/(SR-SS)*SS
    USR(3) = DR*(SR-UR)/(SR-SS)*(ER/DR+(SS-UR)*(SS+PR/(DR*(SR-UR))))

    ! Computer intercell fluxes; 10.34 of Toro
    IF (SL>=0) THEN
      FC(I,:) = FL(:)
    ELSEIF ((SL<=0).AND.(SS>=0)) THEN
      FC(I,:) = FL(:)+SL*(USL(:)-UUL(:))
    ELSEIF ((SS<=0).AND.(SR>=0)) THEN
      FC(I,:) = FR(:)+SR*(USR(:)-UUR(:))
    ELSEIF (SR<=0) THEN
      FC(I,:) = FR(:)
    END IF

  END DO

  ! Obtain new integration variables
  ! from intercell fluxes; exclude boundaries
  ! Standard Godunov upwind scheme, 10.2 of Toro
  DO I=2,NX-1
    UP(I,:) = U(I,:) + DT/DX*(FC(I-1,:)-FC(I,:))
  END DO

END SUBROUTINE


!############################!
!# Sets Boundary Conditions #!
!############################!
SUBROUTINE BOUNDARY(NX, BCS, UP)

  INTEGER, INTENT(IN) :: NX
  INTEGER, INTENT(IN) :: BCS
  REAL*8, INTENT(INOUT) :: UP(NX,3)

  SELECT CASE (BCS)

  CASE (1)
    ! Transmission boundary condition
    UP(1,:) = UP(2,:)
    UP(NX,:) = UP(NX-1,:)

  CASE (2)
    ! Reflection boundary condition
    UP(1,:) = UP(2,:)
    UP(NX,:) = UP(NX-1,:)
    UP(1,2) = -1*UP(2,2)
    UP(NX,2) = -1*UP(NX-1,2)

  END SELECT

END SUBROUTINE

!#######################!
!# Advances Simulation #!
!#######################!
SUBROUTINE ADVANCE(NX,UP,DT,U,T,IT)

  INTEGER, INTENT(IN) :: NX
  REAL*8, INTENT(IN) :: UP(NX,3)
  REAL*8, INTENT(IN) :: DT
  REAL*8, INTENT(OUT) :: U(NX,3)
  REAL*8, INTENT(INOUT) :: T, IT

  U(:,:) = UP(:,:)
  T = T + DT
  IT = IT + 1

END SUBROUTINE

!################!
!# Update CVARS #!
!################!
SUBROUTINE GETCVARS(NDUMP,IT,T,DT,CVARS)

  REAL*8, INTENT(IN) :: NDUMP,IT,T,DT
  REAL*8, INTENT(OUT) :: CVARS(4)

  CVARS(1) = NDUMP
  CVARS(2) = IT
  CVARS(3) = T
  CVARS(4) = DT

END SUBROUTINE

!########################!
!# Writes primitives to #!
!# output data file     #!
!########################!
SUBROUTINE OUTPUT(DIR,ROOT,EXT,NX,P,METHOD,PARAMS,CVARS)

  CHARACTER(*), INTENT(IN) :: DIR
  CHARACTER(*), INTENT(IN) :: ROOT, EXT
  INTEGER, INTENT(IN) :: NX
  REAL*8, INTENT(IN) :: P(NX,3)
  INTEGER, INTENT(IN) :: METHOD
  REAL*8, INTENT(IN) :: PARAMS(5)
  REAL*8, INTENT(IN) :: CVARS(4)

  REAL*8 :: CP, ETA, GAM, L, NT
  REAL*8 :: IT, T, DT
  INTEGER :: I
  CHARACTER(30) :: FNAME
  INTEGER :: TIME(8)

  ! De-bundle simulation parameters and variables
  CP = PARAMS(1)
  ETA = PARAMS(2)
  GAM = PARAMS(3)
  L = PARAMS(4)
  NT = PARAMS(5)
  NDUMP = CVARS(1)
  IT = CVARS(2)
  T = CVARS(3)
  DT = CVARS(4)

  ! Build output filename
  100 FORMAT(A,A,A,I4.4,A,A)
  WRITE(FNAME,100) DIR, ROOT, '.', NDUMP, '.', EXT

  OPEN(UNIT=10,FILE=FNAME,STATUS='REPLACE')

  101 FORMAT('Dumping output ',I4.4,' to ', A)
  WRITE(*,101) INT(NDUMP), FNAME

  ! Header
  CALL DATE_AND_TIME(VALUES=TIME)
  201 FORMAT('# Dumped on ',I2.2,'/',I2.2,'/',I4,1X,I2.2,':',I2.2,':',I2.2)
  WRITE(10,201) TIME(3),TIME(2),TIME(1),TIME(5),TIME(6),TIME(7)
  202 FORMAT('# Using METHOD=',I1)
  WRITE(10,202) METHOD
  203 FORMAT('# CP=',F5.3,', ETA=',F5.3,', GAM=',F4.2)
  WRITE(10,203) CP,ETA,GAM
  204 FORMAT('# NX=',I4,', L=',F5.1,' NT=',F7.4)
  WRITE(10,204) NX,L,NT
  205 FORMAT('# IT= ',I5.5,',   T=',ES12.4,' (',F7.3,'%),   DT=',ES12.4)
  WRITE(10,205) INT(IT),T,T/NT*100,DT

  ! Data
  200 FORMAT(F5.3,3(1X,ES12.4))
  DX = L / NX
  DO I=1,NX
    WRITE(10,200) I*DX,P(I,:)
  END DO
  CLOSE(UNIT=10)

END SUBROUTINE

!#######################!
!# Plotting subroutine #!
!# Calls PGPLOT        #!
!#######################!
SUBROUTINE PLOT(NX,P,PARAMS,CVARS,P0)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NX
  REAL*8, INTENT(IN) :: P(NX,3)
  REAL*8, INTENT(IN) :: PARAMS(5)
  REAL*8, INTENT(IN) :: CVARS(4)
  REAL*8, INTENT(IN) :: P0(8)

  REAL*8 :: CP, ETA, GAM, L, NT, IT, T, DT
  REAL*8 :: X0, UL0, UR0, DL0, DR0, PL0, PR0

  INTEGER :: I
  CHARACTER(50) :: buff
  REAL :: XP(NX), DX
  REAL :: XMIN, XMAX, YMIN, YMAX

  ! De-bundle simulation parameters and variables
  CP = PARAMS(1)
  ETA = PARAMS(2)
  GAM = PARAMS(3)
  L = PARAMS(4)
  NT = PARAMS(5)
  IT = CVARS(2)
  T = CVARS(3)
  DT = CVARS(4)

  ! De-bundle initial conditions
  X0 = P0(1)
  UL0 = P0(2)
  UR0 = P0(3)
  DL0 = P0(4)
  DR0 = P0(5)
  PL0 = P0(6)
  PR0 = P0(7)
  NT = P0(8)

  DX = L/NX
  ! Generate x-axis
  DO I=1,NX
    XP(I)=I*DX
  END DO
  XMIN = XP(1)
  XMAX = XP(NX)

  ! Start plot buffer
  CALL PGBBUF()

  ! Velocity
  CALL PGPANL(1,1)
  CALL PGERAS()
  CALL GETRANGE(NX,P(:,1),YMIN,YMAX)
  CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
  CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
  CALL PGPOINT(NX,XP,REAL(P(:,1)),4)
  CALL PGLINE(NX,XP,REAL(P(:,1)))
  CALL PGLABEL('\Fix\Fn','\Fiu\Fn','Velocity')

  ! Density
  CALL GETRANGE(NX,P(:,2),YMIN,YMAX)
  CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
  CALL PGPANL(2,1)
  CALL PGERAS()
  CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
  CALL PGPOINT(NX,XP,REAL(P(:,2)),4)
  CALL PGLINE(NX,XP,REAL(P(:,2)))
  CALL PGLABEL('\Fix\Fn','\gr','Density')

  ! Pressure
  CALL GETRANGE(NX,P(:,3),YMIN,YMAX)
  CALL PGPANL(1,2)
  CALL PGERAS()
  CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
  CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
  CALL PGPOINT(NX,XP,REAL(P(:,3)),4)
  CALL PGLINE(NX,XP,REAL(P(:,3)))
  CALL PGLABEL('\Fix\Fn','\Fip\Fn','Pressure')

  ! Info panel
  CALL PGSWIN(0.0,100.0,0.0,100.0)
  CALL PGPANL(2,2)
  CALL PGERAS()
  CALL PGSCH(2.5)
  100 FORMAT(A,F8.5,A,F6.2,A)
  WRITE(buff,100) 'T=',T,' (',T/NT*100,'%)'
  CALL PGTEXT(0.,95.0,buff)
  CALL PGSCH(2.0)
  101 FORMAT(A,I5,A,F8.5,A)
  WRITE(buff,101) 'IT=',INT(IT),'  (DT=',DT,')'
  CALL PGTEXT(0.,87.0,buff)

  ! Parameters
  CALL PGTEXT(0.,72.0,'Simulation Parameters')
  200 FORMAT(A,I5.4,A,F5.2,A,F5.2)
  WRITE(buff,200) 'NX=',NX,', CP=',CP,', ETA=',ETA
  CALL PGTEXT(0.,63.0,buff)
  CALL PGTEXT(0.,48.0,'Initial Conditions')
  300 FORMAT(A,F9.5,A,F9.5)
  WRITE(buff,300) '\Fiu\dL\u\Fn=',UL0,', \Fiu\dR\u\Fn=',UR0
  CALL PGTEXT(0.,39.0,buff)
  400 FORMAT(A,F6.3,A,F6.3)
  WRITE(buff,400) '\Fi\gr\dL\u\Fn=',DL0,', \Fi\gr\dR\u\Fn=',DR0
  CALL PGTEXT(0.,30.0,buff)
  500 FORMAT(A,F7.1,A,F6.3)
  WRITE(buff,500) '\Fip\dL\u\Fn=',PL0,', \Fip\dR\u\Fn=',PR0
  CALL PGTEXT(0.,21.0,buff)
  CALL PGSCH(1.75)

  ! Plot buffer
  CALL PGEBUF()

END SUBROUTINE

!#####################!

!#####################!
!# UTILITY FUNCTIONS #!
!#####################!

!#######################!
!# Determines plotting #!
!# range for PLOT()    #!
!#######################!
SUBROUTINE GETRANGE(NX,P,YMIN,YMAX)

  INTEGER, INTENT(IN) :: NX
  REAL*8, INTENT(IN) :: P(NX)
  REAL, INTENT(OUT) :: YMIN, YMAX

  REAL :: MINV, MAXV
  REAL :: AMIN, AMAX, DIF

  AMIN = MINV(NX,P)
  AMAX = MAXV(NX,P)

  IF (AMIN==AMAX) THEN
    ! Treat empty range
    IF (AMIN<0) THEN
      YMIN = AMIN*1.1
      YMAX = 0
    ELSEIF (AMIN>0) THEN
      YMIN = 0
      YMAX = AMAX*1.1
    ELSE
      YMIN = -1
      YMAX = 1
    END IF
  ELSE
    ! Otherwise, expand range +/- 10%
    DIF =  AMAX - AMIN
    YMIN = AMIN - DIF/10
    YMAX = AMAX + DIF/10
  END IF

END SUBROUTINE

!########################!
!# Find maximum/minimum #!
!# value of an array    #!
!#########################!
REAL FUNCTION MAXV(N,A)

  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: A(N)
  REAL :: BUF
  INTEGER :: I

  BUF = REAL(A(1))

  DO I=2,N
    IF (REAL(A(I)) > BUF) BUF = REAL(A(I))
  END DO

  MAXV = BUF

  RETURN 

END FUNCTION

REAL FUNCTION MINV(N,A)

  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: A(N)
  REAL :: BUF
  INTEGER :: I

  BUF = REAL(A(1))
  DO I=2,N
    IF (REAL(A(I)) < BUF) BUF = REAL(A(I))
  END DO

  MINV = BUF

  RETURN

END FUNCTION