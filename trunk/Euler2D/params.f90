!#################################!
!#       PARAMETERS module       #!
!# Defines simulation parameters #!
!# for the Euler2D program       #!
!#################################!

module PARAMETERS

  !##############!
  !# PARAMETERS #!
  !##############!

  ! Define Simulation Parameters
  integer, parameter :: NX = 100        ! Number of grid points along X
  integer, parameter :: NY = 100        ! Number of grid points along Y
  real*8, parameter  :: LX = 1.0        ! Grid physical length along X
  real*8, parameter  :: LY = 1.0        ! Grid physical length along Y
  real*8, parameter  :: CP = 0.9        ! Courant Parameter
  real*8, parameter  :: GAM = 1.4       ! Heat capacity ratios

  integer, parameter :: NOUT = 20       ! Number of outputs
  integer, parameter :: VERB = 0        ! Print verbose progress

  ! PLOT sets the plotting options:
  !  0: do not plot
  !  1: plot each step to screen in real-time
  !  2: plot each step in real-time, hardcopy final step
  !  3: plot a hardcopy of each step
  integer, parameter :: PLOT = 0

  ! SOLVER specifies the solution algorithm:
  !  1: HLL Riemann solver (Godunov upwind scheme)
  !  2: HLLC Riemann solver (Godunov upwind scheme)
  integer, parameter :: SOLVER = 1

  ! WSPD sets the wavespeed estimation method
  !   1: PVRS scheme, 10.51 of Toro
  !   2: TRRS scheme, 10.53 of Toro
  !   3: TSRS scheme, 10.55 of Toro
  !   4: Adaptive Noniterative Riemann Solver;
  !      see section 9.5.2 of Toro
  integer, parameter :: WSPD = 4

  ! ICS allows for built-in initial conditions
  !  0: Generic custom ICs, set in sub INITCOND
  integer, parameter :: ICS = 1

  ! BCS sets boundary conditions type
  !  1: Free-flow (transmission)
  !  2: Reflection
  !  3: Periodic
  integer, parameter :: BCS = 1

  ! Final Integration time
  real*8, parameter :: NT = 1.0

  ! Output Filename Template
  ! Iteration number will be placed between root and extension
  ! DIR: data dump directory path (must exist beforehand)
  ! ROOT: root file name
  ! EXT: extension
  character(*), parameter :: DIR = './data/'
  character(*), parameter :: ROOT = 'Explosion'
  character(*), parameter :: EXT = 'dat'

! --------------------------------
! User-defined parameters END here
! --------------------------------

  ! Bundle simulation parameters
  real*8, parameter :: PARAMS(5) = (/LX,LY,CP,GAM,NT/)
  integer, parameter :: FLAGS(7) = (/SOLVER,WSPD,ICS,BCS,NOUT,VERB,PLOT/)
!  (/NX,NY,LX,LY,CP,GAM,NT/)

END MODULE
