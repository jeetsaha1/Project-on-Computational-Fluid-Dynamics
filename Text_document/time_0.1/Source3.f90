
!
!
!                         **********************************************
!                         *                                            *
!                         *                                            *
!                         *                                            *
!**************************                                            *****************************
program channel_flow
    implicit none

    !----------------------------
    ! Integer variables
    !----------------------------
    integer :: i, j, k, l
    integer :: nx, ny, n, itrn
    integer :: j1, k1, i1, i2, i3
    integer :: nxp1, nxp2
    integer :: ii
    character(len=20) :: fname
    integer, parameter :: nprof = 7
    integer :: iprof(nprof)
    integer :: p
    integer :: iterp
    real(kind=8) :: bigp_time

    !----------------------------
    ! Real (double precision) scalars
    !----------------------------
    real(kind=8) :: Da, Re, dt, FinalTime, omega, Error, alpha
    real(kind=8) :: LX, LY, h, time, rhs, pcof, div
    real(kind=8) :: bigu, bigv, bigp, bigT
    real(kind=8) :: dxx, bigdiv, bigex, ddt, dyy
    real(kind=8) :: ua, LL, ub, LS, LNS
    real(kind=8) :: Rd, Pr, Ri

    !----------------------------
    ! 1D allocatable arrays
    !----------------------------
    real(kind=8), allocatable :: x(:), dx(:), y(:)
    real(kind=8), allocatable :: ys(:), yt(:), yb(:)
    real(kind=8), allocatable :: dy(:), dxr(:), dxl(:)
    real(kind=8), allocatable :: dyt(:), dyb(:)
    real(kind=8), allocatable :: DN(:), B(:)
    real(kind=8), allocatable :: C1(:), D1(:), uu(:)

    !----------------------------
    ! 2D allocatable arrays
    !----------------------------
    real(kind=8), allocatable :: u(:,:), u1(:,:)
    real(kind=8), allocatable :: v(:,:), v1(:,:)
    real(kind=8), allocatable :: cpp(:,:), cpp1(:,:)
    real(kind=8), allocatable :: vp(:,:), si(:,:)
    real(kind=8), allocatable :: T(:,:), T1(:,:)
    real(kind=8), allocatable :: AN(:,:), BN(:,:), CN(:,:)
    real(kind=8), allocatable :: DW(:,:), JN(:,:)
    real(kind=8), allocatable :: A(:,:), ff(:,:)

    !----------------------------
    ! 3D allocatable array
    !----------------------------
    real(kind=8), allocatable :: EN(:,:,:)

    !----------------------------
    ! Geometry / misc
    !----------------------------
    real(kind=8) :: area_fluid
    real(kind=8) :: xs, H1, H2

    integer :: is, jH1, jH2
    integer, allocatable :: flag(:,:)

    !----------------------------
    ! Initial values
    !----------------------------
    time = 0.0_8
    itrn = 0






!**********************************************************************!
!******************* INITIALIZATION OF SOME VARIABLES *****************!
!**********************************************************************!

nx = 501          ! number of grid points in x-direction
ny = 101          ! number of grid points in y-direction
n  = 3            ! system dimension (u, v, T)

LX = 5.0_8        ! channel length
LY = 1.0_8        ! channel height

Ri = 1.0_8        ! Richardson number
Pr = 6.2_8        ! Prandtl number
Re = 100.0_8      ! Reynolds number

dt        = 5.0d-4    ! time step
FinalTime = 0.02_8      ! final simulation time
omega     = 0.7_8      ! relaxation factor
Error     = 1.0E-5_8   ! convergence tolerance

ddt = 1.0_8 / dt
Rd  = 1.0_8 / (Re * Pr)

! Geometry of backward-facing step
xs = 2.5_8        ! step location
H1 = 1.0_8        ! inlet height
H2 = 0.5_8        ! outlet height after step



!**********************************************************************!
!**********************************************************************!







!**********************************************************************!
!************************** MEMORY ALLOCATION *************************!
!**********************************************************************!

allocate(AN(n,n))     ! coefficient matrix (left)
allocate(BN(n,n))     ! coefficient matrix (center)
allocate(CN(n,n))     ! coefficient matrix (right)
allocate(DN(n))       ! RHS vector

! Grid arrays
allocate(x(nx))
allocate(dx(nx))
allocate(y(ny))
allocate(ys(ny))
allocate(yt(ny))
allocate(yb(ny))
allocate(dy(ny))
allocate(dxl(nx))
allocate(dxr(nx))
allocate(dyt(ny))
allocate(dyb(ny))

! Flow variables
allocate(u(nx,ny))     ! x-velocity (n)
allocate(u1(nx,ny))    ! x-velocity (n+1)
allocate(v(nx,ny))     ! y-velocity (n)
allocate(v1(nx,ny))    ! y-velocity (n+1)
allocate(T(nx,ny))     ! temperature (n)
allocate(T1(nx,ny))    ! temperature (n+1)

! Pressure and stream function
allocate(cpp(nx,ny))   ! pressure correction
allocate(cpp1(nx,ny))  ! pressure correction increment
allocate(vp(nx,ny))    ! pressure
allocate(si(nx,ny))    ! stream function

! Vargas algorithm arrays
allocate(DW(n,nx))
allocate(JN(n,nx))
allocate(EN(n,n,nx))
allocate(A(n,n))
allocate(B(n))

! Auxiliary arrays
allocate(C1(nx))
allocate(D1(nx))
allocate(ff(n,n))

! Geometry mask
allocate(flag(nx,ny))
flag = 0

!**********************************************************************!
!**********************************************************************!









!**********************************************************************!
!*********************** CALCULATION OF dx ****************************!
!**********************************************************************!

dxx = LX / real(nx - 1, kind=8)

! Uniform grid spacing
do i = 1, nx
    dx(i) = dxx
end do

! x-coordinate
x(1) = 0.0d0
do i = 2, nx
    x(i) = x(i-1) + dx(i-1)
end do

! Right and left staggered spacings
dxr(1) = 0.5d0 * (dx(1) + dx(2))
do i = 2, nx-1
    dxr(i) = 0.5d0 * (dx(i) + dx(i+1))
    dxl(i) = 0.5d0 * (dx(i) + dx(i-1))
end do

dxl(nx) = 0.5 * (dx(nx) + dx(nx-1))
dxl(1)  = dxl(2)
dxr(nx) = dxr(nx-1)

! Write grid data
open(unit=10, file='Dx.dat',  status='replace')
open(unit=20, file='dxl.dat', status='replace')
open(unit=30, file='dxr.dat', status='replace')

do i = 1, nx
    write(10,*) i, dx(i),  x(i)
    write(20,*) i, dxl(i), x(i)
    write(30,*) i, dxr(i), x(i)
end do

close(10)
close(20)
close(30)

!**********************************************************************!
!**********************************************************************!




!**********************************************************************!
!*********************** CALCULATION OF dy ****************************!
!**********************************************************************!

dyy = LY / real(ny - 1, kind=8)

! Uniform grid spacing
do j = 1, ny
    dy(j) = dyy
end do

! y-coordinate
y(1) = 0.0d0
do j = 2, ny
    y(j) = y(j-1) + dy(j-1)
end do

! Cell face locations
yt(1) = 0.5d0 * (y(1) + y(2))
yb(ny) = 0.5d0 * (y(ny) + y(ny-1))

do j = 2, ny-1
    yt(j) = 0.5d0 * (y(j) + y(j+1))
    yb(j) = 0.5d0 * (y(j) + y(j-1))
end do

! Staggered grid spacing
dyt(1) = 0.5d0 * (dy(1) + dy(2))
do j = 2, ny-1
    dyt(j) = 0.5d0 * (dy(j) + dy(j+1))
    dyb(j) = 0.5d0 * (dy(j) + dy(j-1))
end do

dyb(ny) = 0.5d0 * (dy(ny) + dy(ny-1))
dyb(1)  = dyb(2)
dyt(ny) = dyt(ny-1)

! Write grid data
open(unit=40, file='Dy.dat',  status='replace')
open(unit=50, file='Dyt.dat', status='replace')
open(unit=60, file='Dyb.dat', status='replace')

do j = 1, ny
    write(40,*) j, dy(j),  y(j)
    write(50,*) j, dyt(j), y(j)
    write(60,*) j, dyb(j), y(j)
end do

close(40)
close(50)
close(60)

! Step geometry indices
is  = int(xs / dxx) + 1
jH1 = int(H1 / dyy) + 1
jH2 = int(H2 / dyy) + 1

!**********************************************************************!
!**********************************************************************!


!**********************************************************************!
!******************** CREATE FLUID / SOLID MASK ************************!
!**********************************************************************!

! Initialize entire domain as SOLID (0)
flag(:,:) = 0

!------------------- Fluid region BEFORE the step ---------------------!
do i = 1, is-1
    do j = 1, jH1
        flag(i,j) = 1     ! Fluid cell
    end do
end do

!------------------- Fluid region AFTER the step ----------------------!
do i = is, nx
    do j = 1, jH2
        flag(i,j) = 1     ! Fluid cell
    end do
end do

!**********************************************************************!
!**********************************************************************!




!**********************************************************************!
!***************************** GEOMETRY *******************************!
!**********************************************************************!

open(unit=24, file='geometry.dat', status='replace')

write(24,*) 'variables = "x", "y", "flag"'
write(24,*) 'zone i=', ny, ', j=', nx, ', f=point'

do i = 1, nx
    do j = 1, ny
        write(24,*) x(i), y(j), flag(i,j)
    end do
end do

close(24)

!**********************************************************************!
!**********************************************************************!






!**********************************************************************!
!************************* INITIAL CONDITIONS *************************!
!**********************************************************************!

do i = 1, nx
    do j = 1, ny

        ! Default initialization
        u(i,j)    = 0.0d0
        u1(i,j)   = 0.0d0
        v(i,j)    = 0.0d0
        v1(i,j)   = 0.0d0
        vp(i,j)   = 0.0d0
        cpp(i,j)  = 0.0d0
        cpp1(i,j) = 0.0d0

        if (flag(i,j) == 1) then
            ! Fluid region
            T(i,j)  = 1.0d0
            T1(i,j) = 1.0d0
        else
            ! Solid region
            T(i,j)  = 0.0d0
            T1(i,j) = 0.0d0
        end if

    end do
end do


!**********************************************************************!
!********************* END OF INITIAL CONDITIONS **********************!
!**********************************************************************!

! Copy initial fields to time-(n+1) arrays
do i = 1, nx
    do j = 1, ny
        u1(i,j) = u(i,j)
        v1(i,j) = v(i,j)
        T1(i,j) = T(i,j)
    end do
end do


















write(*,*) 'DEBUG: time = ', time
write(*,*) 'DEBUG: FinalTime = ', FinalTime


!**********************************************************************!
!*************************** TIME LOOP ********************************!
!**********************************************************************!

do while (time < FinalTime)

   itrn = itrn + 1

   if (mod(itrn,100) == 0) then
      write(*,'("Iter=",I6," Time=",F8.4)') itrn, time
end if

   !====================================================
   ! MAIN SOLVER LOOP (Y–SWEEP)
   !====================================================
   do k = 2, ny-1

      !------------------------------------------------
      ! Boundary initialization for Vargas algorithm
      !------------------------------------------------
      JN(:,1)  = (/ u(1,k), v(1,k), T(1,k) /)
      EN(:,:,1) = 0.0d0

      JN(:,nx) = (/ u(nx,k), v(nx,k), T(nx,k) /)
      EN(:,:,nx) = 0.0d0



      !================================================
      ! X–SWEEP: ASSEMBLE + FORWARD VARGAS SWEEP
      !================================================
      do i = 2, nx-1

         if (flag(i,k) == 0) then
            JN(:,i) = (/ u(i,k), v(i,k), T(i,k) /)
            EN(:,:,i) = 0.0d0
            cycle
        end if



         !================================================
         ! U–MOMENTUM
         !================================================
         AN(1,1) = -0.25d0 * dy(k)*(u(i,k)+u(i-1,k))  &
                   - (1.0d0/Re)*(dy(k)/dx(i))
         AN(1,2) = 0.0d0
         AN(1,3) = 0.0d0

         BN(1,1) = (dxr(i)*dy(k))/dt                                   &
                   + 0.25d0*dy(k)*(u(i+1,k)-u(i-1,k))                 &
                   + (1.0d0/Re)*( dy(k)*(1/dx(i+1)+1/dx(i))          &
                   + dxr(i)*(1/dyt(k)+1/dyb(k)) )

         BN(1,2) = 0.0d0
         BN(1,3) = 0.0d0

         CN(1,1) = 0.25d0*dy(k)*(u(i,k)+u(i+1,k))  &
                   - (1.0d0/Re)*(dy(k)/dx(i+1))
         CN(1,2) = 0.0d0
         CN(1,3) = 0.0d0

         DN(1) = u(i,k)*(dxr(i)*dy(k))/dt                              &
                 - dy(k)*(vp(i+1,k)-vp(i,k))                           &
                 + (1.0d0/Re)*dxr(i)*(u1(i,k+1)/dyt(k)+u1(i,k-1)/dyb(k))

         !================================================
         ! V–MOMENTUM
         !================================================
         AN(2,1) = 0.0d0
         AN(2,2) = - (1.0d0/Re)*(dyt(k)/dxl(i))
         AN(2,3) = 0.0d0

         BN(2,1) = 0.0d0
         BN(2,2) = (dx(i)*dyt(k))/dt                                   &
                   + (1.0d0/Re)*( dyt(k)*(1/dxr(i)+1/dxl(i)) )
         BN(2,3) = 0.0d0

         CN(2,1) = 0.0d0
         CN(2,2) = - (1.0d0/Re)*(dyt(k)/dxr(i))
         CN(2,3) = 0.0d0

         DN(2) = v(i,k)*(dx(i)*dyt(k))/dt                              &
                 + dx(i)*(vp(i,k)-vp(i,k+1))                           &
                 + (1.0d0/Re)*(v1(i,k+1)*(dx(i)/dy(k+1))               &
                 + v1(i,k-1)*(dx(i)/dy(k))                             &
                 + Ri*T1(i,k)*dx(i)*dy(k))

         !================================================
         ! TEMPERATURE
         !================================================
         AN(3,1) = 0.0d0
         AN(3,2) = 0.0d0
         AN(3,3) = -Rd*(dy(k)/dxl(i))

         BN(3,1) = 0.0d0
         BN(3,2) = 0.0d0
         BN(3,3) = (dx(i)*dy(k))/dt                                   &
                   + Rd*( dy(k)/dxr(i)+dy(k)/dxl(i)                  &
                   + dx(i)/dyt(k)+dx(i)/dyb(k) )

         CN(3,1) = 0.0d0
         CN(3,2) = 0.0d0
         CN(3,3) = -Rd*(dy(k)/dxr(i))

         DN(3) = T(i,k)*(dx(i)*dy(k))/dt                              &
                 + Rd*dx(i)*(T1(i,k+1)/dyt(k)+T1(i,k-1)/dyb(k))

         !================================================
         ! VARGAS FORWARD SWEEP
         !================================================
         A(:,:) = BN(:,:) - matmul(AN(:,:), EN(:,:,i-1))
         B(:)   = DN(:)   - matmul(AN(:,:), JN(:,i-1))
         
         
         do j1 = 1, 3
            if (abs(A(j1,j1)) < 1.0d-12) then
               A(j1,j1) = 1.0d0
               B(j1)    = JN(j1,i-1)
        end if
        end do

         call matinv(A,3,FF)

         JN(:,i)   = matmul(FF,B)
         EN(:,:,i) = matmul(FF,CN)

      end do   ! i-loop


      DW(:,nx) = JN(:,nx)
      !================================================
      ! BACKWARD SUBSTITUTION
      !================================================
      do ii = nx-1, 2, -1
         DW(:,ii) = JN(:,ii) - matmul(EN(:,:,ii), DW(:,ii+1))
      end do

      !================================================
      ! UPDATE SOLUTION
      !================================================
      do ii = 2, nx-1
         u1(ii,k) = DW(1,ii)
         v1(ii,k) = DW(2,ii)
         T1(ii,k) = DW(3,ii)
      end do

   end do   ! k-loop
   
      !----------------------------------------------------------------------!
!--------------------- BOUNDARY CONDITIONS ----------------------------!
!----------------------------------------------------------------------!

!====================== LOWER WALL (y = 1) ============================!
do i = 1, nx
    if (flag(i,1) == 0) cycle

    u(i,1)  = -u(i,2)
    u1(i,1) = -u1(i,2)

    v(i,1)  = 0.0
    v1(i,1) = 0.0

    T(i,1)  = 1.0
    T1(i,1) = 1.0

    cpp(i,1) = cpp(i,2)
end do


!================ UPPER WALL BEFORE STEP (i < is, j = jH1) =============!
do i = 1, is-1
    if (flag(i,jH1) == 0) cycle

    u(i,jH1)  = -u(i,jH1-1)
    u1(i,jH1) = -u1(i,jH1-1)

    v(i,jH1)  = 0.0
    v1(i,jH1) = 0.0

    T(i,jH1)  = 0.0
    T1(i,jH1) = 0.0

    cpp(i,jH1) = cpp(i,jH1-1)
end do


!================ UPPER WALL AFTER STEP (i = is, j = jH2) ==============!
do i = is, nx
    if (flag(i,jH2) == 0) cycle

    u(i,jH2)  = -u(i,jH2-1)
    u1(i,jH2) = -u1(i,jH2-1)

    v(i,jH2)  = 0.0
    v1(i,jH2) = 0.0

    T(i,jH2)  = 0.0
    T1(i,jH2) = 0.0

    cpp(i,jH2) = cpp(i,jH2-1)
end do


!==================== VERTICAL STEP FACE (x = is) =====================!
do j = jH2+1, jH1
    if (flag(is,j) == 0) cycle

    u(is,j)  = 0.0
    u1(is,j) = 0.0

    v(is,j)  = 0.0
    v1(is,j) = 0.0

    T(is,j)  = 0.0
    T1(is,j) = 0.0

    cpp(is,j) = cpp(is-1,j)
end do


!========================== INLET (x = 1) ==============================!
do j = 1, jH1
    if (flag(1,j) == 0) cycle

    u(1,j)  = 2.0 - u(2,j)
    u1(1,j) = 2.0 - u1(2,j)

    v(1,j)  = 0.0
    v1(1,j) = 0.0

    T(1,j)  = 0.0
    T1(1,j) = 0.0

    cpp(1,j) = cpp(2,j)
end do


!========================== OUTLET (x = nx) ============================!
do j = 1, ny
    if (flag(nx,j) == 0) cycle

    u(nx,j)  = u(nx-1,j)
    u1(nx,j) = u1(nx-1,j)

    v(nx,j)  = 0.0
    v1(nx,j) = 0.0

    T(nx,j)  = T(nx-1,j)
    T1(nx,j) = T1(nx-1,j)

    cpp(nx,j) = cpp(nx-1,j)
end do

!----------------------------------------------------------------------!
!-------------------- END OF BOUNDARY CONDITIONS ----------------------!
!----------------------------------------------------------------------!
   
   
!----------------------------------------------------------------------!
!------------------ PRESSURE POISSON ITERATION ------------------------!
!----------------------------------------------------------------------!
do iterp = 1, 300

   bigp = 0.0d0
   cpp1(:,:) = 0.0d0

   do j = 2, ny-1
      do i = 2, nx-1
         if (flag(i,j) == 0) cycle

         div = (u1(i,j)-u1(i-1,j))*dy(j) &
             + (v1(i,j)-v1(i,j-1))*dx(i)

         pcof = (2.0d0/3.0d0)*dt * ( &
                 dy(j)*(1/dxr(i)+1/dxl(i)) &
               + dx(i)*(1/dyt(j)+1/dyb(j)) )

         rhs  = (2.0d0/3.0d0)*dt * ( &
                 dy(j)*(cpp(i+1,j)/dxr(i)+cpp(i-1,j)/dxl(i)) &
               + dx(i)*(cpp(i,j+1)/dyt(j)+cpp(i,j-1)/dyb(j)) )

         cpp1(i,j) = omega * ((-div + rhs)/pcof)
         bigp = max(bigp, abs(cpp1(i,j)))
      end do
   end do

   cpp = cpp + cpp1

   ! Neumann BC
   do i = 1, nx
      cpp(i,1)  = cpp(i,2)
      cpp(i,ny) = cpp(i,ny-1)
   end do
   do j = 1, ny
      cpp(1,j)  = cpp(2,j)
      cpp(nx,j) = cpp(nx-1,j)
   end do

   if (bigp < Error) exit
   bigp_time = bigp
end do

!----------------------------------------------------------------------!

  vp(:,:) = vp(:,:) + cpp(:,:)







!----------------------------------------------------------------------!
!------------ vVELOCITY CORRERCTION ---------------------!
!----------------------------------------------------------------------!

do j = 2, ny-1
   do i = 2, nx-1
      if (flag(i,j) == 0) cycle

      if (flag(i+1,j) == 1) then
         u1(i,j) = u1(i,j) - (2.0d0*dt)/(3.0d0*dxr(i)) &
                   *(cpp(i+1,j)-cpp(i,j))
      end if

      if (flag(i,j+1) == 1) then
         v1(i,j) = v1(i,j) - (2.0d0*dt)/(3.0d0*dyt(j)) &
                   *(cpp(i,j+1)-cpp(i,j))
      end if
   end do
end do





!====================================================
   ! 5) CONVERGENCE & UPDATE
   !====================================================
bigu = 0.0d0 ; bigv = 0.0d0 ; bigT = 0.0d0 ; bigdiv = 0.0d0

   do j = 2, ny-1
      do i = 2, nx-1
         if (flag(i,j) == 0) cycle

         bigu = max(bigu, abs(u1(i,j)-u(i,j)))
         bigv = max(bigv, abs(v1(i,j)-v(i,j)))
         bigT = max(bigT, abs(T1(i,j)-T(i,j)))

         div = (u1(i,j)-u1(i-1,j))*dy(j)+(v1(i,j)-v1(i,j-1))*dx(i)
         bigdiv = max(bigdiv, abs(div))

         u(i,j) = u1(i,j)
         v(i,j) = v1(i,j)
         T(i,j) = T1(i,j)
      end do
   end do

    if (mod(itrn, 100) == 0) then
       write(*,'(A,F10.6,3X,A,ES10.3,3X,A,ES10.3,3X,A,ES10.3)') &
        'Time=', time, &
        'bigu=', bigu, &
        'bigv=', bigv, &
        'bigT=', bigT
        call flush(6)
    end if


    if (mod(itrn,100)==0) then
       write(*,'("CFL ˜ ",F8.4)') maxval(abs(u))*dt/minval(dx)
       end if


      if (mod(itrn,50) == 0) then
         write(*,'( &
         "Iter=",I6, &
         " Time=",F10.5, &
          " bigu=",ES10.3, &
          " bigv=",ES10.3, &
          " bigT=",ES10.3, &
          " bigp=",ES10.3, &
          " bigdiv=",ES10.3 )') &
          itrn, time, bigu, bigv, bigT, bigp_time, bigdiv
        end if



   time = time + dt

   if (mod(itrn,1) == 0) then
      write(*,'("Iter=",I6," Time=",F10.5," div=",ES10.3)') itrn,time,bigdiv
   end if

end do   ! TIME LOOP






!**********************************************************************!
! FINAL TIME CHECK
!**********************************************************************!

if (time >= FinalTime) then
   write(*,*) 'Final time reached: ', time
end if

!**********************************************************************!
! WRITE SOLUTION FIELDS
!**********************************************************************!

open(70,  file='u.dat')
open(80,  file='v.dat')
open(90,  file='p.dat')
open(100, file='T.dat')

write(70,*)  'variables = "x", "y", "u"'
write(80,*)  'variables = "x", "y", "v"'
write(90,*)  'variables = "x", "y", "p"'
write(100,*) 'variables = "x", "y", "T"'

write(70,*)  'zone i=', ny, ', j=', nx, ', f=point'
write(80,*)  'zone i=', ny, ', j=', nx, ', f=point'
write(90,*)  'zone i=', ny, ', j=', nx, ', f=point'
write(100,*) 'zone i=', ny, ', j=', nx, ', f=point'

do i = 1, nx
   do j = 1, ny
      write(70,*)  x(i), y(j), u(i,j)
      write(80,*)  x(i), y(j), v(i,j)
      write(90,*)  x(i), y(j), vp(i,j)
      write(100,*) x(i), y(j), T(i,j)
   end do

   write(70,*)
   write(80,*)
   write(90,*)
   write(100,*)
end do

close(70)
close(80)
close(90)
close(100)

!**********************************************************************!
! U–, V–, T– PROFILES AT FIXED X-LOCATIONS
!**********************************************************************!



iprof = (/ 350, 400, 450, 500, 550, 600, 650 /)

!------------------- Open files -------------------!
do p = 1, nprof

   write(fname,'("u",I4.4,".dat")') iprof(p)
   open(unit = 160 + (p-1)*3, file = fname, status='replace')

   write(fname,'("v",I4.4,".dat")') iprof(p)
   open(unit = 161 + (p-1)*3, file = fname, status='replace')

   write(fname,'("T",I4.4,".dat")') iprof(p)
   open(unit = 162 + (p-1)*3, file = fname, status='replace')

end do


!------------------- Write profiles -------------------!
do j = 1, ny
   do p = 1, nprof
      if (iprof(p) <= nx) then
         write(160 + (p-1)*3, *) y(j), u(iprof(p), j)
         write(161 + (p-1)*3, *) y(j), v(iprof(p), j)
         write(162 + (p-1)*3, *) y(j), T(iprof(p), j)
      end if
   end do
end do

!------------------- Close files -------------------!
do p = 1, nprof
   close(160 + (p-1)*3)
   close(161 + (p-1)*3)
   close(162 + (p-1)*3)
end do


!**********************************************************************!
!*************************** STREAMLINE ******************************!
!**********************************************************************!

! Initialize streamline at bottom
do i = 1, nx
   si(i,1) = 0.0d0
end do

! Compute streamline field
do i = 1, nx
   do j = 2, ny-1

      ! Skip solid cells
      if (flag(i,j) == 0) then
         si(i,j) = si(i,j-1)
         cycle
      end if

      si(i,j) = si(i,j-1) &
         - (dy(j)/4.0d0) * ( &
             (u(i,j+1)*dy(j) + u(i,j)*dy(j+1)) / dyt(j) &
           + (u(i,j-1)*dy(j) + u(i,j)*dy(j-1)) / dyb(j) &
           )

   end do
end do

!-------------------- Output streamline --------------------!
open(130, file='streamLine.dat', status='replace')

write(130,*) 'variables = "x", "y", "psi"'
write(130,*) 'zone i=', nx, ', j=', ny, ', f=point'

do i = 1, nx
   do j = 1, ny
      if (flag(i,j) == 0) cycle
      write(130,*) x(i), y(j), si(i,j)
   end do
   write(130,*) ' '
end do

close(130)

!**********************************************************************!
!-------------------------- AVERAGE VELOCITY --------------------------!
!**********************************************************************!

ua          = 0.0d0
area_fluid = 0.0d0

do i = 1, nx
   do j = 1, ny

      ! Skip solid cells
      if (flag(i,j) == 0) cycle

      ua          = ua + u(i,j) * dx(i) * dy(j)
      area_fluid = area_fluid + dx(i) * dy(j)

   end do
end do

if (area_fluid > 0.0d0) then
   ua = ua / area_fluid
else
   write(*,*) 'ERROR: zero fluid area'
   stop
end if

open(131, file='averagevelocity.dat', status='replace')
write(131,'(A,ES14.6)') 'average_velocity = ', ua
close(131)
!
write(*,*) 'Program finished. Press ENTER to exit.'
read(*,*)










!**********************************************************************!
!******************** CALCULATE INVERSE OF MATRIX *********************!
!**********************************************************************!
contains


subroutine matinv(M, n, Minv)

   implicit none
   integer, intent(in) :: n
   real(kind=8), intent(in)  :: M(3,3)
   real(kind=8), intent(out) :: Minv(3,3)

   real(kind=8) :: A, B, C, D, E, F, G, H, K, Det

   Minv(:,:) = 0.0d0

   if (n == 2) then

      Det = M(1,1)*M(2,2) - M(1,2)*M(2,1)

      if (abs(Det) < 1.0d-12) then
         write(*,*) 'matinv ERROR: singular 2x2 matrix'
         stop
      end if

      Minv(1,1) =  M(2,2) / Det
      Minv(1,2) = -M(1,2) / Det
      Minv(2,1) = -M(2,1) / Det
      Minv(2,2) =  M(1,1) / Det

   else if (n == 3) then

      A = M(2,2)*M(3,3) - M(2,3)*M(3,2)
      B = M(2,3)*M(3,1) - M(2,1)*M(3,3)
      C = M(2,1)*M(3,2) - M(2,2)*M(3,1)
      D = M(1,3)*M(3,2) - M(1,2)*M(3,3)
      E = M(1,1)*M(3,3) - M(1,3)*M(3,1)
      F = M(3,1)*M(1,2) - M(1,1)*M(3,2)
      G = M(1,2)*M(2,3) - M(1,3)*M(2,2)
      H = M(1,3)*M(2,1) - M(1,1)*M(2,3)
      K = M(1,1)*M(2,2) - M(1,2)*M(2,1)

      Det = M(1,1)*A + M(1,2)*B + M(1,3)*C

      if (abs(Det) < 1.0d-12) then
         write(*,*) 'matinv ERROR: singular 3x3 matrix'
         stop
      end if

      Minv(1,1) = A / Det
      Minv(1,2) = D / Det
      Minv(1,3) = G / Det
      Minv(2,1) = B / Det
      Minv(2,2) = E / Det
      Minv(2,3) = H / Det
      Minv(3,1) = C / Det
      Minv(3,2) = F / Det
      Minv(3,3) = K / Det

   else
      write(*,*) 'matinv ERROR: unsupported matrix size'
      stop
   end if

end subroutine matinv
end program channel_flow
