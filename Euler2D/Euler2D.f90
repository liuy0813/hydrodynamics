!##########################################!
!#               EULER2D                  #!
!# ====================================== #!
!#              24 Mar 2009               #!
!#           by J.C. Toledo-Roy           #!
!##########################################!
!# This code solves 2D hidrodynamic Euler #!
!# equations with the HLL/HLLC Riemann    #!
!# solver. Several wavespeed estimation   #!
!# methods are available.                 #!
!#                                        #!
!# This program uses no global variables! #!
!##########################################!

program Euler2D

  ! Load Parameters
  use PARAMETERS

!###############################################################!

  implicit none

  !################!
  !# DECLARATIONS #!
  !################!

  ! Simulation Variables
  !  U(:,1)= rho      P(:,1)= u       F(:,1)=rho*u
  !  U(:,2)= rho*u    P(:,2)= rho     F(:,2)=rho*u^2+P
  !  U(:,3)= E        P(:,3)= P       F(:,3)=u(E+P)
  !  where E = 0.5*rho*u^2 + P/(gamma-1)
  real*8 :: U(NX,NY,4)      ! Integration variables
  real*8 :: UP(NX,NY,4)     ! Stepped integration variables
  real*8 :: P(NX,NY,4)      ! Primitives
  real*8 :: F(NX,NY,4)      ! Fluxes

  ! Data bundles
  real*8 :: CVARS(4)     ! Simulation current values <MODIFY>

  ! Utility variables
  real*8 :: IT = 0         ! Iteration number
  real*8 :: T = 0          ! Actual Time
  real*8 :: DT             ! Timestep size
  real*8 :: DX = LX / NX   ! Grid spacing in X
  real*8 :: DY = LY / NY   ! Grid spacing in Y
  real*8 :: NDUMP = 0      ! Number of outputs to write
  real*8 :: mark           ! Timing marker

!###############################################################!

  !#################!
  !# PROGRAM START #!
  !#################!

  ! Initialize PGPLOT for real-time plotting
!   if (PLOTRT==1) then
!     call PGBEGIN(0,'/XSERVE',2,2)
!     !call PGENV(0.0,1.0,0.0,1.0,0,-2)
!     call PGSVP(0.1,0.85,0.1,0.85)
!     call PGSCH(1.75)
!     call PGASK(.FALSE.)
!     call PGSUBP(2,2)
!   end if

  ! Load initial conditions
  call INITCOND(ICS,NX,P,P0,NT)

  ! Initialize integration variables
  call GETVARS(NX,P,GAM,U)
  call GETCVARS(NDUMP,IT,T,DT,CVARS)

  ! Dump initial conditions
  call OUTPUT(DIR,ROOT,EXT,NX,P,SOLVER,PARAMS,CVARS)
  NDUMP = NDUMP + 1

  ! Plot ICs (only in real-time plotting)
  if (PLOTRT==1) call PLOT(NX,P,PARAMS,CVARS,P0)

  ! Start timer
  call TIC(mark)

  ! Main Program Loop
  do while(T<=NT)

    call GETSTEP(P,NX,GAM,DX,CP,DT)    ! Get timestep

    ! Solve Euler Equations
    ! for non-boundary points
    select case (SOLVER)

      case (1)
        ! HLL Riemann solver
        call HLL(NX,U,P,GAM,DT,DX,WSPD,UP,F)

      case (2)
        ! HLLC Riemann solver
        call HLLC(NX,U,P,GAM,DT,DX,SOLVER,WSPD,AVG,UP,F)

      case default:
        print*, 'Invalid solver'
        stop

    end SELECT

    ! Calculate boundary grid points
    call BOUNDARY(NX,BCS,UP)

    ! Advance simulation
    call ADVANCE(NX,UP,DT,U,T,IT)
    ! Update primitives
    call GETPRIMS(NX,U,GAM,P)
    ! Update CVARS bundle
    call GETCVARS(NDUMP,IT,T,DT,CVARS)

    ! Show progress (in verbose mode)
    if (VERB==1) then
      300 FORMAT('IT= ',I5.5,', T=',ES12.4,' (',F7.3,'%), DT=',ES12.4)
      WRITE(*,300), int(IT),T,T/NT*100,DT
    end if

    !Output data when apropriate
    if (real(T)/NT > NDUMP/NOUT) then
      call OUTPUT(DIR,ROOT,EXT,NX,P,SOLVER,PARAMS,CVARS)
      NDUMP = NDUMP + 1
    end if

    ! Plot (if real-time)
    if (PLOTRT==1) call PLOT(NX,P,PARAMS,CVARS,P0)

    !if (IT==10) STOP

  end do

  ! Elapsed time
  call TOC(mark)

  ! Finalize PGPLOT on-screen facilty
!   if (PLOTRT==1) call PGEND()

  ! Hardcopy last simulation frame when asked
!   if (COPY == 1) then
!     call PGBEGIN(0,'/PS',2,2)
!     call PGENV(0.0,1.0,0.0,1.0,0,-2)
!     call PGSCH(1.75)
!     call PGASK(.FALSE.)
!     call PGSUBP(2,2)
!     call PLOT()
!     call PGEND()
!   end if

end PROGRAM