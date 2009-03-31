!#################################!
!#       PARAMETERS module       #!
!# Defines simulation parameters #!
!# for the Euler1D program       #!
!#################################!

MODULE PARAMETERS

  !##############!
  !# PARAMETERS #!
  !##############!

  ! Define Simulation Parameters
  INTEGER, PARAMETER :: NX = 100        ! Number of grid points
  REAL*8, PARAMETER  :: L = 1.0         ! Grid physical length
  REAL*8, PARAMETER  :: CP = 0.9        ! Courant Parameter
  REAL*8, PARAMETER  :: GAM = 1.4       ! Heat capacity ratios
  INTEGER, PARAMETER :: NOUT = 20       ! Number of outputs
  INTEGER, PARAMETER :: VERB = 1        ! Print verbose progress
  INTEGER, PARAMETER :: PLOTRT = 1      ! Plot in real-time
  INTEGER, PARAMETER :: COPY = 0        ! Make a hardcopy

  ! SOLVER specifies the solution algorithm:
  !  1: Lax-Friedrichs method
  !  2: MacCormack method
  !  3: MacCormack method with flux-corrected viscosity
  !  4: Simple Godunov upwind scheme
  !  5: Godunov upwind scheme + HLL Riemann solver
  !  6: Godunov upwind scheme + HLLC Riemann solver
  !  7: Godunov + HLLC with Linear Data Reconstruction
  INTEGER, PARAMETER :: SOLVER = 6

  ! FOR METHODS 5, 6 and 7 ONLY
  ! WSPD sets the wavespeed estimation method:
  ! HLL (method 5) can use any (extra info ignored)
  ! HLLC (method 6) *must* use one of 4,5,6,7
  !  (Direct estimates, HLL only!)
  !   1: Davis direct, 10.37 of Toro
  !   2: Davis minmax, 10.38 of Toro
  !   3: Max single 'S+', 10.43 of Toro
  !  (Pressure-based estimates)
  !   4: PVRS scheme, 10.51 of Toro
  !   5: TRRS scheme, 10.53 of Toro
  !   6: TSRS scheme, 10.55 of Toro
  !   7: Adaptive Noniterative Riemann Solver;
  !      see section 9.5.2 of Toro
  INTEGER, PARAMETER :: WSPD = 7

  ! FOR METHOD 7 ONLY
  ! Averaging function to be used for data reconstruction
  !  1: simple arithmetic average
  !  2: minmod
  !  3: van Albada
  INTEGER, PARAMETER :: AVG = 1

  ! Artificial vicosity (for MacCormack schemes)
  REAL*8, PARAMETER :: ETA = 0.3

  ! ICS allows for built-in initial conditions
  !  0: Generic custom ICs, set in sub INITCOND
  !  1-5: Riemann problem tests of Toro
  !  6: Custom Riemann problem test
  INTEGER, PARAMETER :: ICS = 1

  ! BCS sets boundary conditions type
  !  1: Free-flow (transmission)
  !  2: Reflection
  INTEGER, PARAMETER :: BCS = 1

  ! Final Integration time
  ! Toro tests re-set their own NT; otherwise set here
  REAL*8 :: NT = 0.1

  ! Output Filename Template
  ! Iteration number will be placed between root and extension
  ! DIR: data dump directory path (must be created beforehand)
  ! ROOT: root file name
  ! EXT: extension
  CHARACTER(*), PARAMETER :: DIR = './data/'
  CHARACTER(*), PARAMETER :: ROOT = 'Riemann'
  CHARACTER(*), PARAMETER :: EXT = 'dat'

END MODULE
