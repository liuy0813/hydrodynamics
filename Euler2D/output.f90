!########################!
!#       OUTPUT         #!
!# Writes primitives to #!
!# output data file     #!
!########################!
!
! IN arguments:
!   DIR, char(*) : data output directory
!   ROOT, char(*) : root of data filename
!   EXT, char(*) : extension of data filename
!   NX, int : number of grid points along X
!   NY, int : number of grid points along Y
!   P, double(NX,NY,4) : vector of primitives
!   STATE, double(4) : vector with current simulation state
!   PARAMS, double(5) : vector with simulation parameters (physical)
!   FLAGS, integer(7) : vector with simulation parameters (flags)
! OUT arguments:
!   None
!
!############################!
subroutine OUTPUT(DIR,ROOT,EXT,NX,NY,P,STATE,PARAMS,FLAGS)

  implicit none

  character(*), intent(in) :: DIR, ROOT, EXT
  integer, intent(in) :: NX
  integer, intent(in) :: NY
  real*8, intent(in) :: P(NX,NY,4)
  real*8, intent(in) :: STATE(4)
  real*8, intent(in) :: PARAMS(5)
  integer, intent(in) :: FLAGS(7)

  character(50) :: fname
  integer :: time(8)
  integer :: i,j
  real*8 :: dx,dy

  ! De-bundle simulation parameters and variables
  integer :: NDUMP, IT
  real*8 :: T, DT
  real*8 :: LX, LY, CP, GAM, NT
  integer :: SOLVER, WSPD, ICS, BCS, NOUT, VERB, PLOT

  NDUMP = int(STATE(1))
  IT    = int(STATE(2))
  T     = STATE(3)
  DT    = STATE(4)
  LX  = PARAMS(1)
  LY  = PARAMS(2)
  CP  = PARAMS(3)
  GAM = PARAMS(4)
  NT  = PARAMS(5)
  SOLVER = FLAGS(1)
  WSPD   = FLAGS(2)
  ICS    = FLAGS(3)
  BCS    = FLAGS(4)
  NOUT   = FLAGS(5)
  VERB   = FLAGS(6)
  PLOT   = FLAGS(7)

  ! Build output filename
  100 format(A,A,A,I4.4,A,A)
  write(FNAME,100) DIR, ROOT, '.', NDUMP, '.', EXT

  open(unit=10,file=FNAME,status='REPLACE')

  101 format('Dumping output ',I4.4,' to ', A)
  write(*,101) NDUMP, FNAME

  ! Header
  call DATE_AND_TIME(VALUES=TIME)
  201 format('# Dumped on ',I2.2,'/',I2.2,'/',I4,1X,I2.2,':',I2.2,':',I2.2)
  write(10,201) TIME(3),TIME(2),TIME(1),TIME(5),TIME(6),TIME(7)
  202 format('# Using METHOD=',I1,' with WSPD=',I1)
  write(10,202) SOLVER, WSPD
  write(10,'(A)') '# Grid settings:'
  203 format('# NX=',I4,', NY=',I4,', LX=',F7.4,', LY=',F7.4)
  write(10,203) NX, NY, LX, LY
  write(10,'(A)') '# Physical parameters:'
  204 format('# CP=',F5.3,', GAM=',F4.2)
  write(10,204) CP,GAM
  write(10,'(A)') '# Other options:'
  205 format('# ICS=',I1,', BCS=',I1,', NOUT=',I3,', VERB=',I1,', PLOT='I1)
  write(10,205) ICS, BCS, NOUT, VERB, PLOT
  write(10,'(A)') '#'
  206 format('# IT=',I5.5,', T=',F7.4,' / ',F7.4,' (',F7.4,'%), DT=',ES12.4)
  write(10,206) IT,T,NT,T/NT*100,DT
  write(10,'(A)') '#'
  207 format('# Data Columns:  X  Y  U  V  RHO  P')
  write(10,207)
  write(10,'(A)') '#'

  ! Data - GNUPLOT format
  300 format(F5.3,1X,F5.3,4(1X,ES12.5))
  dx = LX / NX
  dy = LY / NY
  do i=1,NX
    do j=1,NY
      write(10,300) i*dx, j*dy, P(i,j,:)
    end do
    write(10,'(A)') ''
  end do

  close(unit=10)

end subroutine
