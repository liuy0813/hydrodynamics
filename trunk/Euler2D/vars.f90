!# This file contains several subroutines to get
!# several simulation variables and parameters.

!#############################!
!#         GETVARS           #!
!# Calculates integration    #!
!# variables from primitives #!
!#############################!
!
! IN arguments:
!   NX, int : number of grid points along X
!   NY, int : number of grid points along Y
!   P, double(NX,NY,4) : vector of primitives
!   GAM, double : heat capacities ratio, Cp/Cv
! OUT arguments:
!   U, double(NX,NY,4) : integration variables derived from primitives
!
!############################!

subroutine GETVARS(NX,NY,P,GAM,U)

  integer, intent(in) :: NX
  integer, intent(in) :: NY
  real*8, intent(in) :: P(NX,NY,4)
  real*8, intent(in) :: GAM
  real*8, intent(out) :: U(NX,NY,4)

  U(:,:,1) = P(:,:,3)
  U(:,:,2) = P(:,:,3)*P(:,:,1)
  U(:,:,3) = P(:,:,3)*P(:,:,2)
  U(:,:,4) = 0.5*P(:,:,3)*(P(:,:,1)**2+P(:,:,2)**2) + P(:,:,4)/(GAM-1)

end subroutine

!##############################!
!#          GETPRIMS          #!
!# Calculates primitives from #!
!# integration variables      #!
!##############################!
!
! IN arguments:
!   NX, int : number of grid points along X
!   NY, int : number of grid points along Y
!   U, double(NX,NY,4) : vector of integration variables
!   GAM, double : heat capacities ratio, Cp/Cv
! OUT arguments:
!   P, double(NX,NY,4) : primitives derived from integration variables
!
!############################!
subroutine GETPRIMS(NX,NY,U,GAM,P)

  integer, intent(in) :: NX
  integer, intent(in) :: NY
  real*8, intent(in) :: U(NX,NY,4)
  real*8, intent(in) :: GAM
  real*8, intent(out) :: P(NX,NY,4)

  P(:,:,1) = U(:,:,2)/U(:,:,1)
  P(:,:,2) = U(:,:,3)/U(:,:,1)
  P(:,:,3) = U(:,:,1)
  P(:,:,4) = (GAM-1)*(U(:,:,4)-0.5*(U(:,:,2)**2+U(:,:,3)**2)/U(:,:,1))

end subroutine

!##########################!
!#        GETFLUX         #!
!# Calculates fluxes from #!
!# integration variables  #!
!##########################!
!
! IN arguments:
!   NX, int: number of grid points along X
!   NY, int: number of grid points along Y
!   P, double(NX,NY,4) : vector of primitives
!   GAM, double: heat capacities ratio, Cp/Cv
! OUT arguments:
!   F, double(NX,NY,4) : x-fluxes derived from primitives
!   G, double(NX,NY,4) : y-fluxes derived from primitives
!
!##########################!
subroutine GETFLUX(NX,NY,P,GAM,F,G)

  integer, intent(in) :: NX
  integer, intent(in) :: NY
  real*8, intent(in) :: P(NX,NY,4)
  real*8, intent(in) :: GAM
  real*8, intent(out) :: F(NX,NY,4)
  real*8, intent(out) :: G(NX,NY,4)

  F(:,:,1) = P(:,:,3)*P(:,:,1)
  F(:,:,2) = P(:,:,3)*(P(:,:,1)**2)+P(:,:,4)
  F(:,:,3) = P(:,:,3)*P(:,:,1)*P(:,:,2)
  F(:,:,4) = P(:,:,1)*(0.5*P(:,:,3)*(P(:,:,1)**2+&
             P(:,:,2)**2)+GAM/(GAM-1)*P(:,:,4))

  G(:,:,1) = P(:,:,3)*P(:,:,2)
  G(:,:,2) = P(:,:,3)*P(:,:,1)*P(:,:,2)
  G(:,:,3) = P(:,:,3)*P(:,:,2)**2+P(:,:,4)
  G(:,:,4) = P(:,:,2)*(0.5*P(:,:,3)*(P(:,:,1)**2+&
             P(:,:,2)**2)+GAM/(GAM-1)*P(:,:,4))

end subroutine

!############################!
!#        GETSTATE          #!
!# Updates the STATE bundle #!
!# with current sim state   #!
!############################!
!
! IN arguments:
!   NDUMP, int : number of last dumped output
!   IT, int : current iteration number
!   T, double : current simulation time
!   DT, double : current timestep
! OUT arguments:
!   STATE, double(4) : bundle with simulation current state
!
!##########################!

subroutine GETSTATE(NDUMP,IT,T,DT,STATE)

  integer, intent(in) :: NDUMP
  integer, intent(in) :: IT
  real*8, intent(in) :: T
  real*8, intent(in) :: DT
  real*8, intent(out) :: STATE(4)

  STATE(1) = NDUMP
  STATE(2) = IT
  STATE(3) = T
  STATE(4) = DT

end subroutine