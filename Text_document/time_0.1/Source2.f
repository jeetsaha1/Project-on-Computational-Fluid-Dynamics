program channel_flow
implicit none

!========================
! Parameters & constants
!========================
integer, parameter :: dp = kind(1.0d0)

!========================
! Grid variables
!========================
integer :: nx, ny, n
integer :: i, j
real(dp) :: LX, LY, dxx, dyy

!========================
! Time control
!========================
real(dp) :: dt, time, FinalTime
integer :: itrn

!========================
! Physical parameters
!========================
real(dp) :: Re, Pr, Ri, Rd, omega, Error

!========================
! Geometry (step)
!========================
real(dp) :: xs, H1, H2
integer :: is, jH1, jH2

!========================
! Arrays
!========================
real(dp), allocatable :: x(:), y(:)
real(dp), allocatable :: dx(:), dy(:)
real(dp), allocatable :: u(:,:), v(:,:), T(:,:)
integer,  allocatable :: flag(:,:)

!========================
! Initialization
!========================
nx = 101
ny = 41
n  = 3

LX = 5.0_dp
LY = 1.0_dp

dt        = 0.001_dp
FinalTime = 0.01_dp
time      = 0.0_dp
itrn      = 0

Re    = 100.0_dp
Pr    = 6.2_dp
Ri    = 1.0_dp
Rd    = 1.0_dp / (Re * Pr)
omega = 0.7_dp
Error = 1.0e-6_dp

xs = 2.5_dp
H1 = 1.0_dp
H2 = 0.5_dp

!========================
! Allocate memory
!========================
allocate(x(nx), dx(nx))
allocate(y(ny), dy(ny))
allocate(u(nx,ny), v(nx,ny), T(nx,ny))
allocate(flag(nx,ny))

!========================
! Grid generation
!========================
dxx = LX / real(nx-1,dp)
dyy = LY / real(ny-1,dp)

x(1) = 0.0_dp
do i = 2, nx
   x(i) = x(i-1) + dxx
end do
dx(:) = dxx

y(1) = 0.0_dp
do j = 2, ny
   y(j) = y(j-1) + dyy
end do
dy(:) = dyy

!========================
! Map step to grid
!========================
is  = int(xs / dxx) + 1
jH1 = int(H1 / dyy) + 1
jH2 = int(H2 / dyy) + 1

!========================
! Fluid / solid mask
!========================
flag = 0

do i = 1, is-1
   do j = 1, jH1
      flag(i,j) = 1
   end do
end do

do i = is, nx
   do j = 1, jH2
      flag(i,j) = 1
   end do
end do

!========================
! Initial conditions
!========================
u = 0.0_dp
v = 0.0_dp
T = 0.0_dp

do i = 1, nx
   do j = 1, ny
      if (flag(i,j) == 1) then
         T(i,j) = 1.0_dp
      end if
   end do
end do

!========================
! Time loop
!========================
do while (time < FinalTime)

   itrn = itrn + 1
   time = time + dt

   ! ---- Boundary conditions (example) ----
   do i = 1, nx
      u(i,1)  = 0.0_dp
      v(i,1)  = 0.0_dp
      T(i,1)  = 1.0_dp
   end do

   ! ---- Solver placeholder ----
   ! (Add momentum, pressure correction, temperature here)

end do

!========================
! Output (simple test)
!========================
open(10,file='geometry.dat',status='replace')
write(10,*) 'variables = "x","y","flag"'
write(10,*) 'zone i=',nx,', j=',ny,', f=point'

do i = 1, nx
   do j = 1, ny
      write(10,*) x(i), y(j), flag(i,j)
   end do
end do
close(10)

print *, 'Program finished successfully'

end program channel_flow

