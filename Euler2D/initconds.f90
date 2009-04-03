!############################!
!#        INITCOND          #!
!#  Set initial conditions  #!
!############################!

! IN arguments:
!   ICS, int : type of initial conditions:
!   NX, int : number of grid points along X
!   NY, int : number of grid points along Y
!   LX, double : grid physical size along X
!   LY, double : grid physical size along Y
! OUT arguments:
!   P, double(NX,NY,4): containing the initialized primitives
!   NT, double: the simulation final time

!############################!

subroutine INITCOND(ICS,NX,NY,LX,LY,P,NT)

  implicit none

  integer, intent(in) :: ICS
  integer, intent(in) :: NX
  integer, intent(in) :: NY
  real*8, intent(in) :: LX
  real*8, intent(in) :: LY
  real*8, intent(out) :: P(NX,NY,4)
  real*8, intent(out) :: NT

  integer :: i, j
  real*8 :: x, y

  ! spherical explosion
  real*8 :: radius, xcen, ycen, dx, dy, d
  radius = 0.1
  xcen = 0.5
  ycen = 0.5
  dx = LX / NX
  dy = LY / NY

  select case (ICS)

  case (0)    ! Custom Initial Conditions
    ! The vectors of PRIMITIVES must be initialized here
    print*, 'No initial conditions loaded.'
    stop

  case (1)   ! Spherical explosion

    do i=1,NX
      do j=1,NY
        P(i,j,1) = 0.0      ! u
        P(i,j,2) = 0.0      ! v
        P(i,j,3) = 1.0      ! rho
        x = i*dx
        y = j*dy
        d = sqrt((x-xcen)**2+(y-ycen)**2)
        if (d <= radius) then
          P(i,j,4) = 1000.0     ! P
        else
          P(i,j,4) = 1.0     ! P
        end if
      end do
    end do

  case default
    print*, 'No initial conditions loaded.'
    stop

  end select


end subroutine
