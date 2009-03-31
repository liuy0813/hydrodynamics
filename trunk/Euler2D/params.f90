!#################################!
!#       PARAMETERS module       #!
!# Defines simulation parameters #!
!# for the EulerND program       #!
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
  integer, parameter :: PLOTRT = 0      ! Plot in real-time
  integer, parameter :: COPY = 0        ! Make a hardcopy

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
  integer, parameter :: ICS = 0

  ! BCS sets boundary conditions type
  !  1: Free-flow (transmission)
  !  2: Reflection
  !  3: Periodic
  integer, parameter :: BCS = 1

  ! Final Integration time
  ! Toro tests re-set their own NT; otherwise set here
  real*8 :: NT = 0.1

  ! Output Filename Template
  ! Iteration number will be placed between root and extension
  ! DIR: data dump directory path (must be created beforehand)
  ! ROOT: root file name
  ! EXT: extension
  character(*), parameter :: DIR = './data/'
  character(*), parameter :: ROOT = 'Riemann'
  character(*), parameter :: EXT = 'dat'

  ! Data bundles
  real*8 :: PARAMS(6)

  ! Bundle simulation parameters
  PARAMS(1) = NX
  PARAMS(2) = NY
  PARAMS(3) = LX
  PARAMS(4) = LY
  PARAMS(5) = CP
  PARAMS(6) = NT

END MODULE