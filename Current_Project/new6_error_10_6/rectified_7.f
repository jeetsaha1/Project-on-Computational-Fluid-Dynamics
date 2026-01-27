

!
!
!                         **********************************************
!                         *                                            *
!                         *                                            *
!                         *                                            *
!**************************                                            *****************************
      program channel_flow
      implicit none

      integer::i,j,k,l,nx,ny,n,itrn,j1,k1,i1,i2,i3,nxp1,nxp2,ii

      real(kind=8)::Da,Re,dt,FinalTime,omega,Error,alpha,
     $  LX,LY,h,time,rhs,pcof,div,bigu,bigv,bigp,bigT,
     $  dxx,bigdiv,bigex,ddt,dyy,ua,LL,ub,LS,LNS,Rd,Pr,Ri

      real(kind=8),dimension(:),allocatable :: x,dx,y,ys,yt,dy,dxr,dxl,
     $ dyt,dyb,DN,B,C1,D1,uu,yb

      real(kind=8),dimension(:,:),allocatable :: u,u1,v,v1,cpp,cpp1,vp,
     $ AN,BN,DW,JN,CN,A,ff,si,T,T1

      real(kind=8),dimension(:,:,:),allocatable :: EN
      integer :: j_top, jmax
      real(kind=8) :: x_step, h_step
      real(kind=8) :: Tref, Sref
      integer i_step, j_step
      real(kind=8), dimension(:,:), allocatable :: Sgen, Sh, Sf
      real(kind=8), dimension(:), allocatable :: Nu_bottom, Nu_top
      real(kind=8) :: Nu_avg, Nu_sum, L_sum
      real(kind=8) :: Sgen_avg, Sgen_sum, Area_sum
      real(kind=8) :: Knff, chi, phi
      real(kind=8) :: U0, Lref, Sgen_avg_nd
      real(kind=8) :: Br
      logical, dimension(nx,ny) :: solid


      ! -------- Base fluid properties (non-dimensional ratios) ----------
      real(kind=8) :: rho_f, cp_f, k_f

      ! -------- Nanoparticle properties ----------
      real(kind=8) :: rho_p, cp_p, k_p

      ! -------- Effective nanofluid properties ----------
      real(kind=8) :: rho_nf, rhocp_nf, rhocp_f
      real(kind=8) :: Cp_ratio

!**********************************************************************!
!******************* INITIALIZATION OF SOME VARIABLES *****************!
!**********************************************************************!
      nx=405 ! no. of grid in x-direction
      ny=81 ! no. of grid in y-direction
      n=3 ! dimenssion
      LX=5.0


      Ri=0.0  !Richardson number
      LY=1.0
      Pr=6.2   !Prandtl number
      dt=0.001 ! time step
      FinalTime=10.0 ! Final Time
      omega=0.7
      Error=1.0E-06
      ddt=1.0/dt
      Re=30.0! Reynolds Number
      Rd=1.0/(Re*Pr)
      !T=300 ! Absolute Temperature
      !h=5E-8 ! Channel height
      x_step = 2.5      ! location of step (in x)
      h_step = 0.5*LY   ! step height (lower half blocked)


      phi  = 0.05      ! nanoparticle volume fraction (0 if pure fluid)
      chi  = 1.0E-04      ! viscous dissipation coefficient (set as per model)
      Knff = 1.0      ! effective thermal conductivity ratio

      U0   = 0.001     ! inlet velocity (matches u(1,j)=2.0-u(2,j))
      Lref = LY      ! channel height




      ! ---- Effective properties (simple models) ----
      Re = Re / (1.0 + 2.5*phi)          ! viscosity increase
      Pr = Pr * (1.0 + 2.5*phi)          ! thermal diffusivity change
      Br = 1.0/Re

      !Br = Br * (1.0 + 2.5*phi)          ! viscous entropy increases


      ! -------- Base fluid (water) ----------
      rho_f = 1.0
      cp_f  = 1.0
      k_f   = 1.0

      ! -------- Nanoparticle (Al2O3) ----------
      rho_p = 3970.0 / 1000.0    ! normalized
      cp_p  = 765.0  / 4180.0
      k_p   = 40.0   / 0.613


      ! -------- Effective density ----------
      rho_nf = (1.0-phi)*rho_f + phi*rho_p

! -------- Effective heat capacity ----------
      rhocp_f  = rho_f * cp_f
      rhocp_nf = (1.0-phi)*rhocp_f + phi * rho_p * cp_p

      Cp_ratio = rhocp_nf / rhocp_f

      ! -------- Maxwell thermal conductivity model ----------
      Knff = ( k_p + 2.0*k_f - 2.0*phi*(k_f - k_p) ) /
     &       ( k_p + 2.0*k_f + phi*(k_f - k_p) )



      Rd = Knff / (Re*Pr*Cp_ratio)

!**********************************************************************!
!**********************************************************************!







!**********************************************************************!
!************************** memory allocation *************************!
!**********************************************************************!
      allocate(AN(n,n)) ! coefficient of v(i-1,j) on (n+1) time step
      allocate(BN(n,n)) ! coefficient of v(i,j) on (n+1) time step
      allocate(CN(n,n)) ! coefficient of v(i+1,j) on (n+1) time step
      allocate(DN(n)) ! constant term / Remaining Term

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

      allocate(u(nx,ny)) ! Along x-direction Velocity component in n-th time step
      allocate(u1(nx,ny)) ! Along x-direction Velocity component in (n+1)-th time step
      allocate(v(nx,ny)) ! Along y-direction  Velocity component in n-th time step
      allocate(v1(nx,ny)) ! Along y-direction  Velocity component in (n+1)-th time step
      allocate(T(nx,ny))  !Temperature component in n-th time step
      allocate(T1(nx,ny))  !Temperature component in (n+1)-th time step


      allocate(cpp(nx,ny)) ! Pressure correction in n-th time step
      allocate(cpp1(nx,ny)) ! Pressure correction in (n+1)-th step
      allocate(vp(nx,ny)) ! Pressure
      allocate(si(nx,ny)) ! Stream Line

      allocate(DW(n,nx))
      allocate(JN(n,nx))
      allocate(EN(n,n,nx))
      allocate(A(n,n))
      allocate(B(n))

      allocate(c1(nx))
      allocate(d1(nx))

      allocate(ff(n,n))

      !allocate(solid(nx,ny))
      allocate(Sgen(nx,ny))
      allocate(Sh(nx,ny))
      allocate(Sf(nx,ny))


      allocate(Nu_bottom(nx))
      allocate(Nu_top(nx))


!**********************************************************************!
!**********************************************************************!







!**********************************************************************!
!************************* calculation of dx **************************!
!**********************************************************************!
      dxx=LX/(nx-1)
      do i=1,nx
        dx(i)=dxx
      end do
      x(1)=0.0 !0.02
      do i=2,nx
        x(i)=x(i-1)+dx(i-1)
      end do

      dxr(1)=0.5*(dx(1)+dx(2))
      do i=2,nx-1
        dxr(i)=0.5*(dx(i)+dx(i+1))
        dxl(i)=0.5*(dx(i)+dx(i-1))
      end do
      dxl(nx)=0.5*(dx(nx)+dx(nx-1))
      dxl(1)=dxl(2)
      dxr(nx)=dxr(nx-1)

      open(10,file='Dx.dat')
      open(20,file='dxl.dat')
      open(30,file='dxr.dat')
      do i=1,nx
        write(10,*)i,dx(i),x(i)
        write(20,*)i,dxr(i),x(i)
        write(30,*)i,dxl(i),x(i)
      end do
      close(10)
      close(20)
      close(30)
!**********************************************************************!
!**********************************************************************!



!**********************************************************************!
!************************** calculation of dy *************************!
!**********************************************************************!
      dyy=LY/(ny-1)

      do j=1,ny
       dy(j)=dyy
      end do

      y(1)=0.0
      do j=2,ny
        y(j)=y(j-1)+dy(j-1)
      end do
      yb(ny)=0.5*( y(ny)+y(ny-1) )
      yt(1)=0.5*( y(1)+y(2) )
      do j=2,ny-1
        yb(j)=0.5*( y(j)+y(j-1) )
        yt(j)=0.5*( y(j)+y(j+1) )
      end do

      dyt(1)=0.5*(dy(1)+dy(2))
      do j=2,ny-1
        dyt(j)=0.5*(dy(j)+dy(j+1))
        dyb(j)=0.5*(dy(j)+dy(j-1))
      end do
      dyb(ny)=0.5*(dy(ny)+dy(ny-1))
      dyb(1)=dyb(2)
      dyt(ny)=dyt(ny-1)

      open(40,file='Dy.dat')
      open(50,file='Dyt.dat')
      open(60,file='Dyb.dat')
      do j=1,ny
        write(40,*)j,dy(j),y(j)
        write(50,*)j,dyt(j),y(j)
        write(60,*)j,dyb(j),y(j)
      end do
      close(40)
      close(50)
      close(60)
!**********************************************************************!
!**********************************************************************!
!*****************************Geometry*********************************!



      do i = 1, nx
        if (x(i) .LE. x_step) i_step = i
      end do

      do j = 1, ny
        if (y(j) .LE. (LY - h_step)) j_step = j
      end do

C================= ZONE 1: BEFORE STEP ==================
      open(24,file='geometry.dat')
      write(24,*) 'variables = "x", "y"'

      write(24,*) 'zone i=',i_step,', j=',ny,', f=point'

      do 10 j = 1, ny
        do 20 i = 1, i_step
          write(24,*) x(i), y(j)
20      continue
        write(24,*)
10    continue


C================= ZONE 2: AFTER STEP ==================
      write(24,*) 'zone i=',nx-i_step+1,', j=',j_step,', f=point'

      do 30 j = 1, j_step
        do 40 i = i_step, nx
          write(24,*) x(i), y(j)
40      continue
        write(24,*)
30    continue

      close(24)




      do i = 1, nx
      do j = 1, ny
      solid(i,j) = .false.
      if (x(i) > x_step .and. y(j) > (LY - h_step)) then
      solid(i,j) = .true.
      end if
      end do
      end do



!**********************************************************************!
!************************* initial conditions *************************!
!**********************************************************************!
      do i=1,nx
        do j=1,ny
          u(i,j)=0.0
          u1(i,j)=0.0

          v(i,j)=0.0
          v1(i,j)=0.0

          T(i,j)=1.0
          T1(i,j)=1.0

          vp(i,j)= 0.0

          cpp(i,j)=0.0
          cpp1(i,j)=0.0
        end do
      end do
!**********************************************************************!
!********************* end of initial conditions **********************!
!**********************************************************************!
      x_step = 2.5
      h_step = 0.5*LY






      open(333,file='u.dat')
      open(334,file='v.dat')
      open(335,file='T.dat')

      do i=1,nx
        do j=1,ny
          u1(i,j)=u(i,j)
          v1(i,j)=v(i,j)
          T1(i,j)=T(i,j)
        end do
      end do

      close(335)
      close(334)
      close(333)


555   continue







!----------------------------------------------------------------------!
!--------------------- BOUNDARY CONDITIONS ----------------------------!
!----------------------------------------------------------------------!
!.......................... Lower wall ................................!
      do i=1,nx
        u(i,1)=0.0- u(i,2)
        u1(i,1)=0.0- u1(i,2)

        v(i,1)=0.0
        v1(i,1)=0.0

        T(i,1)=1.0
        T1(i,1)=1.0

        cpp(i,1)=cpp(i,2)
        cpp1(i,1)=cpp1(i,2)

C        vp(i,1)=vp(i,2)
      end do


!==================== UPPER WALL (VARIABLE HEIGHT) ====================!
      do i = 1, nx

        if (x(i) .le. x_step) then
          j_top = ny
        else
          j_top = j_step
        end if

        u(i,j_top)  = -u(i,j_top-1)
        u1(i,j_top) = -u1(i,j_top-1)

        v(i,j_top)  = 0.0
        v1(i,j_top) = 0.0

        T(i,j_top)  = 1.0
        T1(i,j_top) = 1.0

        cpp(i,j_top)  = cpp(i,j_top-1)
        cpp1(i,j_top) = cpp1(i,j_top-1)

      end do


!.............................. Inlet .................................!
      do j=1,ny
        if (y(j) .le. (LY-h_step)) then
        u(1,j)  = U0 *(LY/(LY-h_step)) - u(2,j)
        u1(1,j) = U0 *(LY/(LY-h_step)) - u1(2,j)
        else
        u(1,j)  = 0.0
        u1(1,j) = 0.0
        end if


        v(1,j)=0.0
        v1(1,j)=0.0

        T(1,j)=0.0
        T1(1,j)=0.0

        cpp(1,j)=cpp(2,j)
        cpp1(1,j)=cpp1(2,j)

C        vp(1,j)=vp(2,j)
        end do


!.............................. outlet ................................!
        do j=1,ny

        u(nx,j)=u(nx-1,j)
        u1(nx,j)=u1(nx-1,j)

        v(nx,j)=0.0
        v1(nx,j)=0.0

        T(nx,j)=T(nx-1,j)
        T1(nx,j)=T1(nx-1,j)

        cpp(nx,j)=cpp(nx-1,j)
        cpp1(nx,j)=cpp1(nx-1,j)

C        vp(nx,j)=vp(nx-1,j)
         end do


!----------------------------------------------------------------------!
!-------------------- End of Boundary Conditions ----------------------!
!----------------------------------------------------------------------!









!----------------------------------------------------------------------


!**********************************************************************!
!******** Steady mixed convection boundary layer flow equation ********!
!**********************************************************************!


!----------------------------------------------------------------------!
!-------------------------- Plain Channel -----------------------------!
!----------------------------------------------------------------------!

      do 1001 k=2,ny-1
!..................... Value On The Right Boundary ....................!
        DW(1,nx)=u(nx,k)
        DW(2,nx)=v(nx,k)
        DW(3,nx)=T(nx,k)
!...................... Boundary Condition For JN .....................!
        JN(1,1)=u(1,k)
        JN(2,1)=v(1,k)
        JN(3,1)=T(1,k)

        JN(1,nx)=u(nx,k)
        JN(2,nx)=v(nx,k)
        JN(3,nx)=T(nx,k)
!...................... Boundary Condition For EN .....................!
        do i=1,3
          do j=1,3
            EN(I,J,1)=0.0
            EN(I,J,nx)=0.0
          end do
        end do

!----------------------------------------------------------------------!

        do 1002 j=2,nx-1

!----------------------------------------------------------------------!
!............... U(ALONG X-DIRECTION) MOMENTUM EQUATION ...............!


          AN(1,1)=-(0.25*dy(k)*(u(j,k)+u(j-1,k)))-(1.0/Re)
     $      *(dy(k)/dx(j))


          AN(1,2)=0.0

          AN(1,3)=0.0


          BN(1,1)=((dxr(j)*dy(k))/dt)

     $      +0.25*dy(k)*(u(j+1,k)-u(j-1,k))



     $      +0.25*((dy(k+1)/dyt(k))*(v(j,k)*dx(j+1)+v(j+1,k)*dx(j)))

     $      -0.25*((dy(k-1)/dyb(k))*(v(j,k-1)*dx(j+1)

     $      +v(j+1,k-1)*dx(j)))



     $      +((1.0/Re)*((dy(k)*((1.0/dx(j+1))+(1.0/dx(j)))

     $      +(dxr(j)*((1.0/dyt(k))+(1.0/dyb(k)))))))



          BN(1,2)=0.0

          BN(1,3)=0.0


          CN(1,1)=(0.25*dy(k)*(u(j,k)+u(j+1,k)))

     $      -((1.0/Re)*(dy(k)/dx(j+1)))

          CN(1,2)=0.0

          CN(1,3)=0.0



          DN(1)=(u(j,k)*(dxr(j)*dy(k))/dt)-0.25

     $       *(u1(j,k+1)*(dy(k)/dyt(k))*(v(j,k)*dx(j+1)

     $       +v(j+1,k)*dx(j)))+0.25*(u1(j,k-1)*(dy(k)/dyb(k))

     $       *(v(j,k-1)*dx(j+1)+v(j+1,k-1)*dx(j)))

     $       -(dy(k)*(vp(j+1,k)-vp(j,k)))

     $       +((1.0/Re)*(dxr(j)*((u1(j,k+1)/dyt(k))

     $       +(u1(j,k-1)/dyb(k)))))



!............... V (ALONG Y-DIRECTION) MOMENTUM EQUATION ..............!


          AN(2,1)=0.0

          AN(2,2)=-((0.25*dx(j)*((u1(j-1,k+1)*dy(k))

     $       +(u1(j-1,k)*dy(k+1))))/dxl(j))

     $       -((1.0/Re)*(dyt(k)/dxl(j)))

          AN(3,1)=0.0



          BN(2,1)=0.0

          BN(2,2)=((dx(j)*dyt(k))/dt)

     $      +0.25*dx(j)*(v(j,k+1)-v(j,k-1))

     $      +((0.25*dx(j+1)*((u(j,k+1)*dy(k))+

     $       (u(j,k)*dy(k+1))))/dxr(j))

     $      -((0.25*dx(j-1)*((u(j-1,k+1)*dy(k))

     $      + (u(j-1,k)*dy(k+1))))/dxl(j))

     $      +((1.0/Re)*((dyt(k)*((1.0/dxr(j))

     $      +(1.0/dxl(j))))+((dx(j)*((1.0/dy(k))+(1.0/dy(k+1)))))))





          BN(2,3)=0.0



          CN(2,1)=0.0

          CN(2,2)=((0.25*dx(j)*((u(j,k+1)*dy(k))+

     $      (u(j,k)*dy(k+1))))/dxr(j))

     $     -((1.0/Re)*(dyt(k)/dxr(j)))




          CN(2,3)=0.0




          DN(2)=v(j,k)*((dx(j)*dyt(k))/dt)

     $      +dx(j)*(vp(j,k)-vp(j,k+1))-0.25*dx(j)

     $      *((v1(j,k+1)*(v(j,k)+v(j,k+1))))

     $      +0.25*dx(j)*(v1(j,k-1)*(v(j,k)+v(j,k-1)))

     $      +(1.0/Re)*(v1(j,k+1)*(dx(j)/dy(k+1))

     $      +v1(j,k-1)*(dx(j)/dy(k))+(Ri)*T1(j,k)*dx(j)*dy(k))



!..................... T(Temperature) EQUATION ........................!

          AN(3,1)=0.0

          AN(3,2)=0.0

          AN(3,3)=-(0.5*dy(k)*(dx(j)/dxl(j))*u(j-1,k))

     $      -(Rd*(dy(k)/dxl(j)))

          BN(3,1)=0.0

          BN(3,2)=0.0

          BN(3,3)=((dx(j)*dy(k))/dt)+0.5*dy(k)


     $      *(((dx(j+1)/dxr(j))*u(j,k))-((dx(j-1)/dxl(j))*u(j-1,k)))

     $      +0.5*dx(j)*(((dy(k+1)/dyt(k))*v(j,k))

     $      -((dy(k-1)/dyb(k))*v(j,k-1)))

     $      +(Rd*((dy(k)/dxr(j))+(dy(k)/dxl(j))+(dx(j)/dyt(k))

     $      +(dx(j)/dyb(k))))




          CN(3,1)=0.0

          CN(3,2)=0.0

          CN(3,3)=(0.5*dy(k)*(dx(j)/dxr(j))*u(j,k))-(Rd*(dy(k)/dxr(j)))



          DN(3)=(T(j,k)*((dx(j)*dy(k))/dt))-0.5*dx(j)

     $      *((T1(j,k+1)*v(j,k)*(dy(k)/dyt(k)))

     $      -(T1(j,k-1)*v(j,k-1)*(dy(k)/dyb(k))))

     $      +(Rd*dx(j)*((T1(j,k+1)/dyt(k))+(T1(j,k-1)/dyb(k))))




!----------------------------------------------------------------------!
!------------------------- VARGAS ALOGORITHM --------------------------!
!----------------------------------------------------------------------!


          do J1=1,3
            DO L=1,3
              B(J1)=0.0
              A(J1,L)=0.0
              DO K1=1,3
                B(J1)=B(J1)+AN(J1,K1)*JN(K1,J-1)
                A(J1,L)=A(J1,L)+AN(J1,K1)*EN(K1,L,J-1)
              end do
            end do
          end do

          DO J1=1,3
            B(J1)=DN(J1)-B(J1)
            DO K1=1,3
              A(J1,K1)=BN(J1,K1)-A(J1,K1)
            end do
          end do

          CALL matinv(A,3,FF)

          DO J1=1,3
            DO L=1,3
              JN(J1,J)=0.0
              EN(J1,L,J)=0.0
              DO K1=1,3
                JN(J1,J)=JN(J1,J)+FF(J1,K1)*B(K1)
                EN(J1,L,J)=EN(J1,L,J)+FF(J1,K1)*CN(K1,L)
              end do
            end do
          end do

1002    continue



        DO I=nx-1,2,-1
          DO J=1,3
            DW(J,I)=0.0
            DO K1=1,3
              DW(J,I)=DW(J,I)+EN(J,K1,I)*DW(K1,I+1)
            end do
          end do

          DO J=1,3
            DW(J,I)=JN(J,I)-DW(J,I)
          end do

        end do

        DO i=2,nx-1
          u1(i,k)=DW(1,i)
          v1(i,k)=DW(2,i)
          T1(i,k)=DW(3,i)
        end do



1001  CONTINUE
!----------------------------------------------------------------------!
!--------------------- end of VARGAS ALOGORITHM -----------------------!
!----------------------------------------------------------------------!


!----------------------------------------------------------------------!
!----------------- PRESSURE CORRECTION CALCULATION --------------------!
!----------------------------------------------------------------------!
5     continue

!----------------------------------------------------------------------!
!--------------------- BOUNDARY CONDITIONS ----------------------------!
!----------------------------------------------------------------------!
!.......................... Lower wall ................................!
      do i=1,nx
        u(i,1)=0.0- u(i,2)
        u1(i,1)=0.0- u1(i,2)

        v(i,1)=0.0
        v1(i,1)=0.0

        T(i,1)=1.0
        T1(i,1)=1.0

        cpp(i,1)=cpp(i,2)
        cpp1(i,1)=cpp1(i,2)

C        vp(i,1)=vp(i,2)
      end do


!==================== UPPER WALL (VARIABLE HEIGHT) ====================!
      do i = 1, nx

        if (x(i) .le. x_step) then
          j_top = ny
        else
          j_top = j_step
        end if

        u(i,j_top)  = -u(i,j_top-1)
        u1(i,j_top) = -u1(i,j_top-1)

        v(i,j_top)  = 0.0
        v1(i,j_top) = 0.0

        T(i,j_top)  = 1.0
        T1(i,j_top) = 1.0

        cpp(i,j_top)  = cpp(i,j_top-1)
        cpp1(i,j_top) = cpp1(i,j_top-1)

      end do





!.............................. Inlet .................................!
      do j=1,ny
        if (y(j) .le. (LY-h_step)) then
        u(1,j)  = U0 *(LY/(LY-h_step)) - u(2,j)
        u1(1,j) = U0 *(LY/(LY-h_step)) - u1(2,j)
        else
        u(1,j)  = 0.0
        u1(1,j) = 0.0
        end if


        v(1,j)=0.0
        v1(1,j)=0.0

        T(1,j)=0.0
        T1(1,j)=0.0

        cpp(1,j)=cpp(2,j)
        cpp1(1,j)=cpp1(2,j)

C        vp(1,j)=vp(2,j)

!.............................. outlet ................................!
        u(nx,j)=u(nx-1,j)
        u1(nx,j)=u1(nx-1,j)

        v(nx,j)=0.0
        v1(nx,j)=0.0

        T(nx,j)=T(nx-1,j)
        T1(nx,j)=T1(nx-1,j)

        cpp(nx,j)=cpp(nx-1,j)
        cpp1(nx,j)=cpp1(nx-1,j)

C        vp(nx,j)=vp(nx-1,j)
      end do
!----------------------------------------------------------------------!
!-------------------- End of Boundary Conditions ----------------------!


!**********************************************************************!
!**********************************************************************!
!**********************************************************************!
!**** PRESSURE CORRECTION EQUATION – ***************************!
!**********************************************************************!

      do i = 2, nx-1
         if (x(i) .le. x_step) then
          jmax = ny-1
        else
          jmax = j_step-1
        end if
        
        do j = 2, jmax

          div = (u1(i,j)-u1(i-1,j))*dy(j)
     $        + (v1(i,j)-v1(i,j-1))*dx(i)

          pcof = (2.0/3.0)*(dt*dy(j))*(1.0/dxr(i)+1.0/dxl(i))
     $         + (2.0/3.0)*(dt*dx(i))*(1.0/dyt(j)+1.0/dyb(j))

          rhs  = (2.0/3.0)*(dt*dy(j))
     $         * ( cpp(i+1,j)/dxr(i) + cpp1(i-1,j)/dxl(i) )
     $         + (2.0/3.0)*(dt*dx(i))
     $         * ( cpp(i,j+1)/dyt(j) + cpp1(i,j-1)/dyb(j) )

          cpp1(i,j) = omega * (-div + rhs) / pcof

        end do
      end do


!**********************************************************************!







!**********************************************************************!
!**********************************************************************!
!**** PRESSURE CORRECTION CONVERGENCE –  ************************!
!**********************************************************************!

      bigp = 0.0

      do i = 2, nx-1

        if (x(i) .le. x_step) then
          jmax = ny-1
        else
          jmax = j_step-1
        end if

        do j = 2, jmax
          bigp = max(bigp, abs(cpp(i,j)-cpp1(i,j)))
          cpp(i,j) = cpp1(i,j)
        end do

      end do




!****************** condition on pressure correction ******************!
      if (bigp .gt. Error) then
        write(*,*) 'bigp =', bigp
        goto 5
      end if


!**********************************************************************!
!**** PRESSURE & VELOCITY CORRECTION – ZONE 1 *************************!
!**********************************************************************!

      do i = 2, nx-1

        if (x(i) .le. x_step) then
          jmax = ny-1
        else
          jmax = j_step-1
        end if

        do j = 2, jmax
          vp(i,j) = vp(i,j) + cpp(i,j)

          u1(i,j) = u1(i,j)
     $     - (2.0*dt/(3.0*dxr(i))) * (cpp(i+1,j)-cpp(i,j))




          v1(i,j) = v1(i,j)
     $      - (2.0*dt/(3.0*dyt(j))) * (cpp(i,j+1)-cpp(i,j))
        end do

      end do




!**********************************************************************!




!**********************************************************************!
!**** CONVERGENCE CRITERIA (U, V, T) –  **************************!
!**********************************************************************!

      bigu = 0.0
      bigv = 0.0
      bigT = 0.0

      do i = 2, nx-1

        if (x(i) .le. x_step) then
          jmax = ny-1
        else
          jmax = j_step-1
        end if

        do j = 2, jmax

          if (abs(u(i,j)-u1(i,j)) .gt. bigu) then
            bigu = abs(u(i,j)-u1(i,j))
          end if
          u(i,j) = u1(i,j)

          if (abs(v(i,j)-v1(i,j)) .gt. bigv) then
            bigv = abs(v(i,j)-v1(i,j))
          end if
          v(i,j) = v1(i,j)

          if (abs(T(i,j)-T1(i,j)) .gt. bigT) then
            bigT = abs(T(i,j)-T1(i,j))
          end if
          T(i,j) = T1(i,j)

        end do
      end do

!**********************************************************************!





      itrn = itrn + 1
      time = time + dt

      if ((bigu .gt. Error) .or.
     &    (bigv .gt. Error) .or.
     &    (bigT .gt. Error)) then
        write(*,*) 'bigu=',bigu,'bigv=',bigv,'bigT=',bigT
        goto 555
      end if




!**********************************************************************!
!****** CONVERGENCE CRITERIA(upon divergence at a particular time)*****!
!**********************************************************************!

      bigdiv = 0.0

!==================== DIVERGENCE  ==============================!

      do i = 2, nx-1

        if (x(i) .le. x_step) then
          jmax = ny-1
        else
          jmax = j_step-1
        end if

        do j = 2, jmax
          div = (u1(i,j)-u1(i-1,j))*dy(j)
     $        + (v1(i,j)-v1(i,j-1))*dx(i)
          bigdiv = max(bigdiv, abs(div))
        end do

      end do




!==================== CONVERGENCE CHECK ================================!

      if (bigdiv .gt. 1.0E-4) then
        write(*,*) 'bigdiv = ', bigdiv
        goto 555
      end if

      if (time .lt. FinalTime) then
        write(*,*) 'Time = ', time
        goto 555
      end if

999   continue





!**********************************************************************!
!******************** ENTROPY GENERATION ******************************!
!**********************************************************************!

C-------------------- INITIALIZATION ----------------------------------!

      do i = 1, nx
        do j = 1, ny
          Sh(i,j)   = 0.0
          Sf(i,j)   = 0.0
          Sgen(i,j) = 0.0
        end do
      end do



C==================== ENTROPY (NANOFLUID) ===================!

      do i = 2, nx-1

        if (x(i) .le. x_step) then
          jmax = ny-1
        else
          jmax = j_step-1
        end if

        do j = 2, jmax

! -------- Thermal entropy (nanofluid, NON-DIMENSIONAL) ----------
          Sh(i,j) = (Knff/(Re*Pr*Cp_ratio)) * (
     &      ((T(i+1,j)-T(i,j))/dxr(i))**2
     &    + ((T(i,j)-T(i-1,j))/dxl(i))**2
     &    + ((T(i,j+1)-T(i,j))/dyt(j))**2
     &    + ((T(i,j)-T(i,j-1))/dyb(j))**2 )

! -------- Friction entropy (nanofluid, NON-DIMENSIONAL) ----------
          Sf(i,j) = chi * Br / ((1.0-phi)**2.5) * (
     &      2.0*((u(i,j+1)-u(i,j))/dyt(j))**2
     &    + 2.0*((v(i+1,j)-v(i,j))/dxr(i))**2
     &    + (((u(i+1,j)-u(i,j))/dx(i+1))
     &     + ((v(i,j+1)-v(i,j))/dy(j+1)))**2 )

! -------- Total entropy generation -----------------------------
          Sgen(i,j) = Sh(i,j) + Sf(i,j)

        end do
      end do





!**********************************************************************!
!************* AVERAGE ENTROPY GENERATION (WHOLE CHANNEL) *************!
!**********************************************************************!

      Sgen_sum = 0.0
      Area_sum = 0.0
      

      do i = 2, nx-1

         if (x(i) .le. x_step) then
            jmax = ny-1
            else
            jmax = j_step-1
         end if

      do j = 2, jmax
         Sgen_sum = Sgen_sum + Sgen(i,j)*dx(i)*dy(j)
         Area_sum = Area_sum + dx(i)*dy(j)
      end do

      end do



      Sgen_avg = Sgen_sum / Area_sum


!**********************************************************************!
!**************** WRITE AVERAGE ENTROPY *******************************!
!**********************************************************************!

      open(431,file='AvgEntropy.dat')
      write(431,*) 'Average non-dimensional entropy generation'
      write(431,*) Sgen_avg
      close(431)













!**********************************************************************!
!************* BOTTOM WALL NUSSELT NUMBER (LOCAL) *********************!
!**********************************************************************!

      do i = 1, nx
        Nu_bottom(i) = 0.0
      end do

      do i = 2, nx-1
        Nu_bottom(i) = -( T(i,2) - T(i,1) ) / dyb(2)
      end do


!**********************************************************************!
!************* AVERAGE NUSSELT NUMBER (BOTTOM WALL) *******************!
!**********************************************************************!



      Nu_sum = 0.0
      L_sum  = 0.0

      do i = 2, nx-1
          Nu_sum = Nu_sum + Nu_bottom(i)*dx(i)
          L_sum  = L_sum  + dx(i)
      end do

      Nu_avg = Nu_sum / L_sum





!**********************************************************************!
!**************** WRITE AVERAGE NUSSELT NUMBER ************************!
!**********************************************************************!

      open(421,file='AvgNusselt.dat')
      write(421,*) 'Average bottom wall Nusselt number'
      write(421,*) Nu_avg
      close(421)








      do i = 1, nx
         if (x(i) .LE. x_step) i_step = i
        end do

      do j = 1, ny
         if (y(j) .LE. (LY - h_step)) j_step = j
       end do















      open(70,file='u.dat')
      write(70,*) 'variables = "x","y","u"'

C==================== BEFORE STEP ====================

      write(70,*) 'zone i=', i_step, ', j=', ny, ', f=point'

      do j = 1, ny
         do i = 1, i_step
            write(70,*) x(i), y(j), u(i,j)
         end do
         write(70,*)
      end do

C==================== AFTER STEP =====================

      write(70,*) 'zone i=', nx-i_step+1, ', j=', j_step, ', f=point'

      do j = 1, j_step
         do i = i_step, nx
            write(70,*) x(i), y(j), u(i,j)
         end do
         write(70,*)
      end do

      close(70)













      open(80,file='v.dat')
      write(80,*) 'variables = "x","y","v"'

      write(80,*) 'zone i=', i_step, ', j=', ny, ', f=point'

      do j = 1, ny
         do i = 1, i_step
            write(80,*) x(i), y(j), v(i,j)
         end do
         write(80,*)
      end do

      write(80,*) 'zone i=', nx-i_step+1, ', j=', j_step, ', f=point'

      do j = 1, j_step
         do i = i_step, nx
            write(80,*) x(i), y(j), v(i,j)
         end do
         write(80,*)
      end do

      close(80)












      open(90,file='p.dat')
      write(90,*) 'variables = "x","y","p"'

      write(90,*) 'zone i=', i_step, ', j=', ny, ', f=point'

      do j = 1, ny
         do i = 1, i_step
            write(90,*) x(i), y(j), vp(i,j)
         end do
         write(90,*)
      end do

      write(90,*) 'zone i=', nx-i_step+1, ', j=', j_step, ', f=point'

      do j = 1, j_step
         do i = i_step, nx
            write(90,*) x(i), y(j), vp(i,j)
         end do
         write(90,*)
      end do

      close(90)







      open(100,file='T.dat')
      write(100,*) 'variables = "x","y","T"'

      write(100,*) 'zone i=', i_step, ', j=', ny, ', f=point'

      do j = 1, ny
         do i = 1, i_step
            write(100,*) x(i), y(j), T(i,j)
         end do
         write(100,*)
      end do

      write(100,*) 'zone i=', nx-i_step+1, ', j=', j_step, ', f=point'

      do j = 1, j_step
         do i = i_step, nx
            write(100,*) x(i), y(j), T(i,j)
         end do
         write(100,*)
      end do

      close(100)






!*************************u-velocity profile***************************!


      open(160,file='u350.dat')
      open(170,file='v350.dat')
      open(180,file='T350.dat')

      open(190,file='u400.dat')
      open(200,file='v400.dat')
      open(210,file='T400.dat')

      open(220,file='u450.dat')
      open(230,file='v450.dat')
      open(240,file='T450.dat')

      open(250,file='u500.dat')
      open(260,file='v500.dat')
      open(270,file='T500.dat')

      open(280,file='u550.dat')
      open(290,file='v550.dat')
      open(300,file='T550.dat')

      open(310,file='u600.dat')
      open(320,file='v600.dat')
      open(330,file='T600.dat')

      open(340,file='u650.dat')
      open(350,file='v650.dat')
      open(360,file='T650.dat')


      do j=1,ny
      write(160,*)y(j),u(350,j)
      write(170,*)y(j),v(350,j)
      write(180,*)y(j),T(350,j)

      write(190,*)y(j),u(400,j)
      write(200,*)y(j),v(400,j)
      write(210,*)y(j),T(400,j)

      write(220,*)y(j),u(450,j)
      write(230,*)y(j),v(450,j)
      write(240,*)y(j),T(450,j)

      write(250,*)y(j),u(500,j)
      write(260,*)y(j),v(500,j)
      write(270,*)y(j),T(500,j)

      write(280,*)y(j),u(550,j)
      write(290,*)y(j),v(550,j)
      write(300,*)y(j),T(550,j)

      write(310,*)y(j),u(600,j)
      write(320,*)y(j),v(600,j)
      write(330,*)y(j),T(600,j)

      write(340,*)y(j),u(650,j)
      write(350,*)y(j),v(650,j)
      write(360,*)y(j),T(650,j)


      end do

      close(160)
      close(170)
      close(180)

      close(190)
      close(200)
      close(210)

      close(220)
      close(230)
      close(240)

      close(250)
      close(260)
      close(270)

      close(280)
      close(290)
      close(300)

      close(310)
      close(320)
      close(330)

      close(340)
      close(350)
      close(360)


!**********************************************************************!
!*************************** stream line ******************************!
!**********************************************************************!
      do i=1,nx
        si(i,1)=0.0
      end do

C==================== STREAMLINE  ====================

      do i = 1, nx

        if (x(i) .le. x_step) then
          jmax = ny-1
        else
          jmax = j_step-1
        end if

        do j = 2, jmax

          si(i,j) = si(i,j-1)
     $    - (dy(j)/4.0)
     $    * ( (u(i,j+1)*dy(j) + u(i,j)*dy(j+1))/dyt(j)
     $      + (u(i,j-1)*dy(j) + u(i,j)*dy(j-1))/dyb(j) )

        end do
      end do



      open(130,file='streamLine.dat')
      write(130,*) 'variables = "x","y","psi"'

      write(130,*) 'zone i=', i_step, ', j=', ny, ', f=point'

      do j = 1, ny
         do i = 1, i_step
            write(130,*) x(i), y(j), si(i,j)
         end do
         write(130,*)
      end do

      write(130,*) 'zone i=', nx-i_step+1, ', j=', j_step, ', f=point'

      do j = 1, j_step
         do i = i_step, nx
            write(130,*) x(i), y(j), si(i,j)
         end do
         write(130,*)
      end do

      close(130)




!**********************************************************************!
!**********************************************************************!
!-------------------------- Average velocity --------------------------!



      ua = 0.0
      Area_sum = 0.0

      do i = 1, nx
           if (x(i) .le. x_step) then
              jmax = ny
           else
               jmax = j_step
           end if

           do j = 1, jmax
               ua = ua + u(i,j)*dx(i)*dy(j)
               Area_sum = Area_sum + dx(i)*dy(j)
        end do
      end do

      ua = ua / Area_sum


      open(131,file='averagevelocity.dat')
      write(131,*) 'averagevelocity=',ua
      close(131)


      end program channel_flow











!**********************************************************************!
!******************** CALCULATE INVERSE OF MATRIX *********************!
!**********************************************************************!
      subroutine matinv(M,n,Minv)

      implicit none
      integer,intent(in) :: n
      real(kind=8),dimension(3,3),intent(in)::M
      real(kind=8),dimension(3,3),intent(inout) ::Minv
      real(kind=8) :: A,B,C,D,E,F,G,H,K,Det

      if(n.eq.2)then
        Det=M(1,1)*M(2,2)-M(1,2)*M(2,1)
        if(Det.eq.0.0)then
           write(*,*)"matinv-ERROR: singular matrix"
           stop
        end if
        Minv(1,1)=M(2,2)/Det
        Minv(1,2)=-M(1,2)/Det
        Minv(2,1)=-M(2,1)/Det
        Minv(2,2)=M(1,1)/Det
      end if

      if(n.eq.3)then
        A=M(2,2)*M(3,3)-M(2,3)*M(3,2)
        B=M(2,3)*M(3,1)-M(2,1)*M(3,3)
        C=M(2,1)*M(3,2)-M(2,2)*M(3,1)
        D=M(1,3)*M(3,2)-M(1,2)*M(3,3)
        E=M(1,1)*M(3,3)-M(1,3)*M(3,1)
        F=M(3,1)*M(1,2)-M(1,1)*M(3,2)
        G=M(1,2)*M(2,3)-M(1,3)*M(2,2)
        H=M(1,3)*M(2,1)-M(1,1)*M(2,3)
        K=M(1,1)*M(2,2)-M(1,2)*M(2,1)
        Det=M(1,1)*(M(2,2)*M(3,3)-M(2,3)*M(3,2))
     $     +M(1,2)*(M(2,3)*M(3,1)-M(2,1)*M(3,3))
     $     +M(1,3)*(M(2,1)*M(3,2)-M(2,2)*M(3,1))

        if(Det.eq.0.0)then
          write(*,*)"matinv-ERROR: singular matrix"
          stop
        end if

        Minv(1,1)=A/Det
        Minv(1,2)=D/Det
        Minv(1,3)=G/Det
        Minv(2,1)=B/Det
        Minv(2,2)=E/Det
        Minv(2,3)=H/Det
        Minv(3,1)=C/Det
        Minv(3,2)=F/Det
        Minv(3,3)=K/Det
      end if

      end subroutine matinv
!***********************************************************************!
!****************************************************************

