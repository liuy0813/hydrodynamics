!###########################################!
!#                EULER1D                  #!
!# ======================================= #!
!#              10 Mar 2009                #!
!#           by J.C. Toledo-Roy            #!
!###########################################!
!# This code implements various solution   #!
!# methods for the 1D hidrodynamic Euler   #!
!# equations:                              #!
!#  1) Lax method (first-order)            #!
!#  2) Macormack method (second-order)     #!
!#     with uniform artificial viscosity   #!
!#  3) Macormack as above, but with flux-  #!
!#     corrected adaptive viscosity        #!
!#  4) Godunov upwind scheme (first-order) #!
!#  5) Godunov as above, but paired with   #!
!#     the HLL Riemann solver,             #!
!#  6) Godunov as above, but paired with   #!
!#     the HLLC Riemann solver             #!
!#  7) Godunov + HLLC solver with a linear #!
!#     data reconstruction method          #!
!# Several wavespeed estimation methods    #!
!# can be used for the HLL/HLLC schemes.   #!
!#                                         #!
!# This program uses no global variables!  #!
!###########################################!

PROGRAM Euler1D

  ! Load Parameters
  USE PARAMETERS

!###############################################################!

  IMPLICIT NONE

  !################!
  !# DECLARATIONS #!
  !################!

  ! Simulation Variables
  !  U(:,1)= rho      P(:,1)= u       F(:,1)=rho*u
  !  U(:,2)= rho*u    P(:,2)= rho     F(:,2)=rho*u^2+P
  !  U(:,3)= E        P(:,3)= P       F(:,3)=u(E+P)
  !  where E = 0.5*rho*u^2 + P/(gamma-1)
  REAL*8 :: U(NX,3)      ! Integration variables
  REAL*8 :: UP(NX,3)     ! Stepped integration variables
  REAL*8 :: P(NX,3)      ! Primitives
  REAL*8 :: F(NX,3)      ! Fluxes

  ! Data bundles
  REAL*8 :: P0(8)       ! Initial conditions
  REAL*8 :: PARAMS(5)   ! Simulation parameters
  REAL*8 :: CVARS(4)    ! Simulation current values

  ! Utility variables
  REAL*8 :: IT = 0      ! Iteration number
  REAL*8 :: T = 0        ! Actual Time
  REAL*8 :: DT           ! Timestep size
  REAL*8 :: DX = L / NX  ! Grid spacing
  REAL*8 :: NDUMP = 0   ! Output number
  REAL*8 :: mark           ! Timing marker

  ! Budle simulation parameters and variables
  PARAMS(1) = CP 
  PARAMS(2) = ETA
  PARAMS(3) = GAM
  PARAMS(4) = L
  PARAMS(5) = NT

!###############################################################!

  !#################!
  !# PROGRAM START #!
  !#################!

  ! Initialize PGPLOT for real-time plotting
  IF (PLOTRT==1) THEN
    CALL PGBEGIN(0,'/XSERVE',2,2)
    !CALL PGENV(0.0,1.0,0.0,1.0,0,-2)
    CALL PGSVP(0.1,0.85,0.1,0.85)
    CALL PGSCH(1.75)
    CALL PGASK(.FALSE.)
    CALL PGSUBP(2,2)
  END IF

  ! Load initial conditions
  CALL INITCOND(ICS,NX,P,P0,NT)

  ! Initialize integration variables
  CALL GETVARS(NX,P,GAM,U)
  CALL GETCVARS(NDUMP,IT,T,DT,CVARS)

  ! Dump initial conditions
  CALL OUTPUT(DIR,ROOT,EXT,NX,P,SOLVER,PARAMS,CVARS)
  NDUMP = NDUMP + 1

  ! Plot ICs (only in real-time plotting)
  IF (PLOTRT==1) CALL PLOT(NX,P,PARAMS,CVARS,P0)

  ! Start timer
  CALL TIC(mark)

  ! Main Program Loop
  DO WHILE(T<=NT)

    CALL GETSTEP(P,NX,GAM,DX,CP,DT)    ! Get timestep

    ! Solve Euler Equations
    ! (for points 2, NX-1)
    SELECT CASE (SOLVER)

      CASE (1)
        ! Lax solver
        CALL LAX(NX,U,P,GAM,DT,DX,UP,F)

      CASE (2,3)
       ! MacCormack solver
        CALL MAC(NX,U,P,GAM,ETA,DT,DX,SOLVER,UP,F)

      CASE (4)
        ! Simple Godunov scheme
        CALL GOD(NX,U,P,GAM,DT,DX,UP,F)

      CASE (5)
        ! Godunov + HLL Riemann solver
        CALL HLL(NX,U,P,GAM,DT,DX,WSPD,UP,F)

      CASE (6)
        ! Godunov + HLLC Riemann solver
        CALL HLLC(NX,U,P,GAM,DT,DX,SOLVER,WSPD,AVG,UP,F)

      CASE (7)
        ! Godunov + HLLC + Domain reconstruction
        CALL HLLC(NX,U,P,GAM,DT,DX,SOLVER,WSPD,AVG,UP,F)

    END SELECT

    ! Calculate boundary grid points
    CALL BOUNDARY(NX,BCS,UP)

    ! Advance simulation
    CALL ADVANCE(NX,UP,DT,U,T,IT)
    ! Update primitives
    CALL GETPRIMS(NX,U,GAM,P)
    ! Update CVARS bundle
    CALL GETCVARS(NDUMP,IT,T,DT,CVARS)

    ! Show progress (in verbose mode)
    IF (VERB==1) THEN
      300 FORMAT('IT= ',I5.5,', T=',ES12.4,' (',F7.3,'%), DT=',ES12.4)
      WRITE(*,300), INT(IT),T,T/NT*100,DT
    END IF

    !Output data when apropriate
    IF (REAL(T)/NT > NDUMP/NOUT) THEN
      CALL OUTPUT(DIR,ROOT,EXT,NX,P,SOLVER,PARAMS,CVARS)
      NDUMP = NDUMP + 1
    END IF

    ! Plot (if real-time)
    IF (PLOTRT==1) CALL PLOT(NX,P,PARAMS,CVARS,P0)

    !IF (IT==10) STOP

  END DO

  ! Elapsed time
  CALL TOC(mark)

  ! Finalize PGPLOT on-screen facilty
  IF (PLOTRT==1) CALL PGEND()

  ! Hardcopy last simulation frame when asked
  IF (COPY == 1) THEN
    CALL PGBEGIN(0,'/PS',2,2)
    CALL PGENV(0.0,1.0,0.0,1.0,0,-2)
    CALL PGSCH(1.75)
    CALL PGASK(.FALSE.)
    CALL PGSUBP(2,2)
    CALL PLOT()
    CALL PGEND()
  END IF

END PROGRAM