!############################!
!#  Set Initial Conditions  #!
!############################!

! This subroutine receives the following arguments:
!
!   ICS: type of initial conditions:
!        0: Custom initial conditions - specify here
!        1: Spherical explosion
!   NX, NY: number of points along X and Y directions
!   LX, LY: grid physical size along X and Y directions
!
! and asigns and returns the following variables:
!
!   P: the 2-dimensional vector containing the
!      initialized primitives
!   NT: the simulation final time

!############################!

subroutine INITCOND(ICS,NX,NY,LX,LY,P,NT)

  integer, intent(in) :: ICS
  integer, intent(in) :: NX
  integer, intent(in) :: NY
  integer, intent(in) :: LX
  integer, intent(in) :: LY
  real*8, intent(out) :: P(NX,NY,3)
  real*8, intent(out) :: NT

  ! parameters for spherical explosion
  real*8 :: radius = 0.1
  real*8 :: xcenter = 0.5
  real*8 :: ycenter = 0.5
  real*8 :: prboost = 1000.0

  real*8 :: dx = LX / NX
  real*8 :: dy = LY / NY

  if (ICS==0) then
    ! Custom Initial Conditions
    ! The vectors of PRIMITIVES must be initialized here
  else
    ! Built-in test initial conditions
    select case (ICS)

    ! Spherical explosion
    case (1):
      do i=1,NX
        do j=1,NY
        
        end do
      end do
    

    case default:
      print*, 'No initial conditions loaded.'
      stop

    end select
  end if

end subroutine
