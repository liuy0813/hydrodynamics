!#######################!
!#       GETSTEP       #!
!# Calculates timestep #!
!#######################!
!
! IN arguments:
!   NX, int : number of grid points along X
!   NY, int : number of grid points along Y
!   P, double(NX,NY,4) : vector of primitives
!   PARAMS, double(5) : vector with simulation parameters (physical)
! OUT arguments:
!   DT, double : calculated timestep
!
!############################!
subroutine GETSTEP(NX,NY,P,PARAMS,DT)

  implicit none

  integer, intent(in) :: NX
  integer, intent(in) :: NY
  real*8, intent(in) :: P(NX,NY,4)
  real*8, intent(in) :: PARAMS(5)
  real*8, intent(out) :: DT

  real*8 :: LX,LY,CP,GAM
  real*8 :: dx, dy, cs
  real*8 :: umax, uloc
  real*8 :: vmax, vloc
  integer :: i, j

  LX = PARAMS(1)
  LY = PARAMS(2)
  CP = PARAMS(3)
  GAM = PARAMS(4)
  dx = LX / NX
  dy = LY / NY

  ! The timestep is calculated following the
  ! Courant-Friedrichs-Lewis criterion:
  ! DT = DX / umax * CP
  ! where DX is the grid spacing, CP is the
  ! Courant parameter and umax is the maximum
  ! value of local flow speed + local sound speed.
  ! For 2D and 3D, the returned timestep is the
  ! minimum timestep along all dimensions.

  umax = 0
  vmax = 0
  do i=1,NX
    do j=1,NY
      cs = sqrt(GAM*P(i,j,4)/P(i,j,3))
      uloc = abs(P(i,j,1)) + cs
      vloc = abs(P(i,j,2)) + cs
      if (uloc > umax) then
        umax = uloc
      end if
      if (vloc > vmax) then
        vmax = vloc
      end if
    end do
  end do
  DT = min(dx/umax*CP,dy/vmax*CP)

end subroutine