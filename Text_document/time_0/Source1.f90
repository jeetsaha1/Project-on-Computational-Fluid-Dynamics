
!
!
!                         **********************************************
!                         *                                            *
!                         *                                            *
!                         *                                            *
!**************************                                            *****************************
      program channel_flow
      implicit none

      integer::i,j,k,l,nx,ny,n,itrn,j1,k1,i1,i2,i3,nxp1,nxp2

      real(kind=8)::Da,Re,dt,FinalTime,omega,Error,alpha,
     &  LX,LY,h,time,rhs,pcof,div,bigu,bigv,bigp,bigT,
     &  dxx,bigdiv,bigex,ddt,dyy,ua,LL,ub,LS,LNS,Rd,Pr,Ri

      real(kind=8),dimension(:),allocatable :: x,dx,y,ys,yt,dy,dxr,dxl,
     & dyt,dyb,DN,B,C1,D1,uu,yb

      real(kind=8),dimension(:,:),allocatable :: u,u1,v,v1,cpp,cpp1,vp,
     & AN,BN,DW,JN,CN,A,ff,si,T,T1

      real(kind=8),dimension(:,:,:),allocatable :: EN
      real(kind=8) :: area_fluid

      real(kind=8) :: xs, H1, H2
      integer :: is, jH1, jH2
      integer, dimension(:,:), allocatable :: flag

      time = 0.0
      itrn = 0


!**********************************************************************!
!******************* INITIALIZATION OF SOME VARIABLES *****************!
!**********************************************************************!
      nx=501 ! no. of grid in x-direction
      ny=101 ! no. of grid in y-direction
      n=3 ! dimenssion
      LX=5.0

      Ri=1.0  !Richardson number
      LY=1.0
      Pr=6.2   !Prandtl number
      dt=0.001 ! time step
      FinalTime=5.0 ! Final Time
      omega=0.1
      Error=1.0E-05
      ddt=1.0/dt
      Re=100.0! Reynolds Number
      Rd=1.0/(Re*Pr)
      !T=300 ! Absolute Temperature
      !h=5E-8 ! Channel height
      xs = 2.5
      H1 = 1.0
      H2 = 0.5


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
      allocate(flag(nx,ny))
      flag(:,:) = 0

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
      
      
      
      is  = INT(xs / dxx) + 1
      jH1 = INT(H1 / dyy) + 1
      jH2 = INT(H2 / dyy) + 1

      
!**********************************************************************!
!************************* Create a Fluid / Solid Mask **************************!
!**********************************************************************!



      do i = 1, is-1
         do j = 1, jH1
          flag(i,j) = 1                   ! Fuild region before step
        end do
      end do


      do i = is, nx
         do j = 1, jH2                    ! Fuild region after step
            flag(i,j) = 1
         end do
      end do



      
!**********************************************************************!
!**********************************************************************!
!*****************************Geometry*********************************!


      !open(24,file='geometry.dat')
      !write(24,*)'variables = "x", "y"'
      !write(24,*)'zone i=',ny,', j=',nx
      !do i=1,nx
       ! do j=1,ny
        !  write(24,*)x(i),y(j)
        !end do
      !end do
      !close(24)
      
      
      
      
      
      !open(24,file='geometry.dat')
      !write(24,*) 'variables = "x", "y"'
      !write(24,*) 'zone f=point'

      !do i = 1, nx
       ! do j = 1, ny
        !    if (flag(i,j) == 1) then
         !         write(24,*) x(i), y(j)
          !            end if
           !             end do
            !            end do
             !
              !          close(24)
              

       open(24,file='geometry.dat')
       write(24,*) 'variables = "x", "y", "flag"'
       write(24,*) 'zone i=',ny,', j=',nx,', f=point'

       do i = 1, nx
         do j = 1, ny
             write(24,*) x(i), y(j), flag(i,j)
         end do
       end do

       close(24)








!**********************************************************************!
!************************* initial conditions *************************!
!**********************************************************************!
      do i = 1, nx
         do j = 1, ny

         if (flag(i,j) == 0) then
         u(i,j)  = 0.0
         u1(i,j) = 0.0
         v(i,j)  = 0.0
         v1(i,j) = 0.0
         T(i,j)  = 0.0
         T1(i,j) = 0.0
         vp(i,j)   = 0.0
         cpp(i,j)  = 0.0
         cpp1(i,j) = 0.0
         else
         u(i,j)  = 0.0
         u1(i,j) = 0.0
          v(i,j)  = 0.0
          v1(i,j) = 0.0
          T(i,j)  = 1.0
          T1(i,j) = 1.0
          vp(i,j)   = 0.0
          cpp(i,j)  = 0.0
          cpp1(i,j) = 0.0
          end if

          end do
        end do

!**********************************************************************!
!********************* end of initial conditions **********************!
!**********************************************************************!
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


do while (time < FinalTime)

   555 continue
   ! solver body

end do








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
       ! cpp1(i,1)=cpp1(i,2)

C        vp(i,1)=vp(i,2)
      end do


!.......................... Upper wall ................................!
!      Upper wall BEFORE the step (i < is, j = jH1)
       do i = 1, is-1
          u(i,jH1)  = -u(i,jH1-1)
          u1(i,jH1) = -u1(i,jH1-1)

          v(i,jH1)  = 0.0
          v1(i,jH1) = 0.0

          T(i,jH1)  = 0.0      ! cold wall
          T1(i,jH1) = 0.0

          cpp(i,jH1)  = cpp(i,jH1-1)
          !cpp1(i,jH1) = cpp1(i,jH1-1)
        end do



       !. Upper wall AFTER the step (i = is, j = jH2)
       do i = is, nx
          u(i,jH2)  = -u(i,jH2-1)
          u1(i,jH2) = -u1(i,jH2-1)

          v(i,jH2)  = 0.0
          v1(i,jH2) = 0.0

          T(i,jH2)  = 0.0      ! cold wall
          T1(i,jH2) = 0.0

          cpp(i,jH2)  = cpp(i,jH2-1)
          !cpp1(i,jH2) = cpp1(i,jH2-1)
        end do


        !Vertical step face
        do j = jH2+1, jH1
           u(is,j)  = 0.0
           u1(is,j) = 0.0

           v(is,j)  = 0.0
           v1(is,j) = 0.0

           T(is,j)  = 0.0
           T1(is,j) = 0.0

           cpp(is,j)  = cpp(is-1,j)
         !  cpp1(is,j) = cpp1(is-1,j)
        end do




!.............................. Inlet .................................!
      do j=1,jH1
        u(1,j)=2.0-u(2,j)
        u1(1,j)=2.0-u1(2,j)

        v(1,j)=0.0
        v1(1,j)=0.0

        T(1,j)=0.0
        T1(1,j)=0.0

        cpp(1,j)=cpp(2,j)
        !cpp1(1,j)=cpp1(2,j)

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
        !cpp1(nx,j)=cpp1(nx-1,j)

!        vp(nx,j)=vp(nx-1,j)
      end do
!----------------------------------------------------------------------!
!-------------------- End of Boundary Conditions ----------------------!
!----------------------------------------------------------------------!

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
           !=====================================
           ! SKIP SOLID CELLS
           !=====================================
           if (flag(j,k) == 0) then
               DW(1,j) = 0.0
               DW(2,j) = 0.0
               DW(3,j) = 0.0

               JN(1,j) = 0.0
               JN(2,j) = 0.0
               JN(3,j) = 0.0
              cycle
           end if
!----------------------------------------------------------------------!
!............... U(ALONG X-DIRECTION) MOMENTUM EQUATION ...............!


          AN(1,1)=-(0.25*dy(k)*(u(j,k)+u(j-1,k)))-(1.0/Re)
     &      *(dy(k)/dx(j))


          AN(1,2)=0.0

          AN(1,3)=0.0


          BN(1,1)=((dxr(j)*dy(k))/dt)

     &      +0.25*dy(k)*(u(j+1,k)-u(j-1,k))



     &      +0.25*((dy(k+1)/dyt(k))*(v(j,k)*dx(j+1)+v(j+1,k)*dx(j)))

     &      -0.25*((dy(k-1)/dyb(k))*(v(j,k-1)*dx(j+1)

     &      +v(j+1,k-1)*dx(j)))



     &      +((1.0/Re)*((dy(k)*((1.0/dx(j+1))+(1.0/dx(j)))

     &      +(dxr(j)*((1.0/dyt(k))+(1.0/dyb(k)))))))



          BN(1,2)=0.0

          BN(1,3)=0.0


          CN(1,1)=(0.25*dy(k)*(u(j,k)+u(j+1,k)))

     &      -((1.0/Re)*(dy(k)/dx(j+1)))

          CN(1,2)=0.0

          CN(1,3)=0.0



          DN(1)=(u(j,k)*(dxr(j)*dy(k))/dt)-0.25

     &       *(u1(j,k+1)*(dy(k)/dyt(k))*(v(j,k)*dx(j+1)

     &       +v(j+1,k)*dx(j)))+0.25*(u1(j,k-1)*(dy(k)/dyb(k))

     &       *(v(j,k-1)*dx(j+1)+v(j+1,k-1)*dx(j)))

     &       -(dy(k)*(vp(j+1,k)-vp(j,k)))

     &       +((1.0/Re)*(dxr(j)*((u1(j,k+1)/dyt(k))

     &       +(u1(j,k-1)/dyb(k)))))



!............... V (ALONG Y-DIRECTION) MOMENTUM EQUATION ..............!


          AN(2,1)=0.0

          AN(2,2)=-((0.25*dx(j)*((u1(j-1,k+1)*dy(k))

     &       +(u1(j-1,k)*dy(k+1))))/dxl(j))

     &       -((1.0/Re)*(dyt(k)/dxl(j)))

          AN(3,1)=0.0



          BN(2,1)=0.0

          BN(2,2)=((dx(j)*dyt(k))/dt)

     &      +0.25*dx(j)*(v(j,k+1)-v(j,k-1))

     &      +((0.25*dx(j+1)*((u(j,k+1)*dy(k))+

     &       (u(j,k)*dy(k+1))))/dxr(j))

     &      -((0.25*dx(j-1)*((u(j-1,k+1)*dy(k))

     &      + (u(j-1,k)*dy(k+1))))/dxl(j))

     &      +((1.0/Re)*((dyt(k)*((1.0/dxr(j))

     &      +(1.0/dxl(j))))+((dx(j)*((1.0/dy(k))+(1.0/dy(k+1)))))))





          BN(2,3)=0.0



          CN(2,1)=0.0

          CN(2,2)=((0.25*dx(j)*((u(j,k+1)*dy(k))+

     &      (u(j,k)*dy(k+1))))/dxr(j))

     &     -((1.0/Re)*(dyt(k)/dxr(j)))




          CN(2,3)=0.0




          DN(2)=v(j,k)*((dx(j)*dyt(k))/dt)

     &      +dx(j)*(vp(j,k)-vp(j,k+1))-0.25*dx(j)

     &      *((v1(j,k+1)*(v(j,k)+v(j,k+1))))

     &      +0.25*dx(j)*(v1(j,k-1)*(v(j,k)+v(j,k-1)))

     &      +(1.0/Re)*(v1(j,k+1)*(dx(j)/dy(k+1))

     &      +v1(j,k-1)*(dx(j)/dy(k))+(Ri)*T1(j,k)*dx(j)*dy(k))



!..................... T(Temperature) EQUATION ........................!

          AN(3,1)=0.0

          AN(3,2)=0.0

          AN(3,3)=-(0.5*dy(k)*(dx(j)/dxl(j))*u(j-1,k))

     &      -(Rd*(dy(k)/dxl(j)))

          BN(3,1)=0.0

          BN(3,2)=0.0

          BN(3,3)=((dx(j)*dy(k))/dt)+0.5*dy(k)


     &      *(((dx(j+1)/dxr(j))*u(j,k))-((dx(j-1)/dxl(j))*u(j-1,k)))

     &      +0.5*dx(j)*(((dy(k+1)/dyt(k))*v(j,k))

     &      -((dy(k-1)/dyb(k))*v(j,k-1)))

     &      +(Rd*((dy(k)/dxr(j))+(dy(k)/dxl(j))+(dx(j)/dyt(k))

     &      +(dx(j)/dyb(k))))




          CN(3,1)=0.0

          CN(3,2)=0.0

          CN(3,3)=(0.5*dy(k)*(dx(j)/dxr(j))*u(j,k))-(Rd*(dy(k)/dxr(j)))



          DN(3)=(T(j,k)*((dx(j)*dy(k))/dt))-0.5*dx(j)

     &      *((T1(j,k+1)*v(j,k)*(dy(k)/dyt(k)))

     &      -(T1(j,k-1)*v(j,k-1)*(dy(k)/dyb(k))))

     &      +(Rd*dx(j)*((T1(j,k+1)/dyt(k))+(T1(j,k-1)/dyb(k))))




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
         if (flag(i,1) == 1) then
        u(i,1)=0.0- u(i,2)
        u1(i,1)=0.0- u1(i,2)

        v(i,1)=0.0
        v1(i,1)=0.0

        T(i,1)=1.0
        T1(i,1)=1.0

        cpp(i,1)=cpp(i,2)
       ! cpp1(i,1)=cpp1(i,2)
        end if
!        vp(i,1)=vp(i,2)
      end do


!.......................... Upper wall ................................!
      ! do i=1,nx
       ! u(i,ny)=0.0- u(i,ny-1)
        !u1(i,ny)=0.0- u1(i,ny-1)
         !
!        v(i,ny)=0.0
 !       v1(i,ny)=0.0
  !      v(i,ny-1)=0.0
   !     v1(i,ny-1)=0.0
    !
     !   T(i,ny)=1.0
      !  T1(i,ny)=1.0

       ! cpp(i,ny)=cpp(i,ny-1)
        !cpp1(i,ny)=cpp1(i,ny-1)

!        vp(i,ny)=vp(i,ny-1)
     ! end do
     
     !  UPPER WALL — BEFORE STEP
      !do i = 1, is-1
      !   u(i,jH1)  = -u(i,jH1-1)
      !   u1(i,jH1) = -u1(i,jH1-1)

      !   v(i,jH1)  = 0.0
      !   v1(i,jH1) = 0.0

      !   T(i,jH1)  = 0.0   ! cold wall
      !   T1(i,jH1) = 0.0
         
      !   cpp(i,jH1)  = cpp(i,jH1-1)
      !   cpp(i,jH1) = cpp(i,jH1-1)
       ! end do
        
        
       !UPPER WALL — AFTER STEP
      ! do i = is, nx
      !    u(i,jH2)  = -u(i,jH2-1)
      !    u1(i,jH2) = -u1(i,jH2-1)

      !    v(i,jH2)  = 0.0
      !    v1(i,jH2) = 0.0

      !    T(i,jH2)  = 0.0
      !    T1(i,jH2) = 0.0
          
      !    cpp(i,jH2)  = cpp(i,jH2-1)
      !    cpp(i,jH2) = cpp(i,jH2-1)
      !  end do


       !VERTICAL STEP FACE
       do j = jH2+1, jH1
          if (flag(is,j) == 0) cycle
          u(is,j)  = 0.0
          u1(is,j) = 0.0

           v(is,j)  = 0.0
           v1(is,j) = 0.0



           T(is,j)  = 0.0
           T1(is,j) = 0.0
           
           cpp(is,j)  = cpp(is-1,j)
          ! cpp(is,j) = cpp(is-1,j)
        end do

      ! UPPER WALL (apply only to fluid cells)
      do i = 1, nx
         if (flag(i,ny) == 0) cycle

            u(i,ny)  = -u(i,ny-1)
            u1(i,ny) = -u1(i,ny-1)

            v(i,ny)  = 0.0
            v1(i,ny) = 0.0

            T(i,ny)  = 0.0
            T1(i,ny) = 0.0

            cpp(i,ny)  = cpp(i,ny-1)
           ! cpp1(i,ny) = cpp1(i,ny-1)
        end do


!.............................. Inlet .................................!
      do j = 1, jH1
         u(1,j)  = 2.0 - u(2,j)
         u1(1,j) = 2.0 - u1(2,j)

         v(1,j)  = 0.0
         v1(1,j) = 0.0

         T(1,j)  = 0.0
         T1(1,j) = 0.0
         
        cpp(1,j)=cpp(2,j)
       ! cpp1(1,j)=cpp1(2,j)

        end do

!.............................. outlet ................................!
      do j = 1, ny
         if (flag(nx,j) == 0) cycle

         u(nx,j)  = u(nx-1,j)
         u1(nx,j) = u1(nx-1,j)

         v(nx,j)  = v(nx-1,j)
         v1(nx,j) = v1(nx-1,j)


         T(nx,j)  = T(nx-1,j)
         T1(nx,j) = T1(nx-1,j)

         cpp(nx,j)  = cpp(nx-1,j)
         !cpp1(nx,j) = cpp1(nx-1,j)
        end do

-----------------------------------------------------------!
!-------------------- End of Boundary Conditions ----------------------!
!----------------------------------------------------------------------!


!**********************************************************************!
!**********************************************************************!
      cpp1(:,:) = 0.0
      do  j = 2, ny-1
      do  i = 2, nx-1

         if (flag(i,j) == 0) then
            cpp(i,j) = 0.0
            cycle
        end if

        div = (u1(i,j)-u1(i-1,j))*dy(j) + (v1(i,j)-v1(i,j-1))*dx(i)


         pcof = (2.0/3.0)*(dt*dy(j))*(1.0/dxr(i)+1.0/dxl(i)) &
              + (2.0/3.0)*(dt*dx(i))*(1.0/dyt(j)+1.0/dyb(j))

        rhs = (2.0/3.0)*(dt*dy(j))*(cpp(i+1,j)/dxr(i) + cpp(i-1,j)/dxl(i)) &
            + (2.0/3.0)*(dt*dx(i))*(cpp(i,j+1)/dyt(j) + cpp(i,j-1)/dyb(j))


        cpp1(i,j) = omega * ((-div + rhs) / pcof)



      end do
      end do



!**********************************************************************!

!**********************************************************************!
      bigp=0.0
      do j=2,ny-1
        do i=2,nx-1
           if (flag(i,j) == 0) then
              cpp(i,j) = 0.0
              cycle
           end if
              bigp = max(bigp, abs(cpp1(i,j)))
          cpp(i,j) = cpp(i,j) + cpp1(i,j)
        end do
      end do


!****************** condition on pressure correction ******************!
      if(bigp.gt.Error)then
        write(*,*)'bigp=',bigp
        go to 5
      end if

      ! Pressure BCs
      do i = 1, nx
         cpp(i,1)  = cpp(i,2)
         cpp(i,ny) = cpp(i,ny-1)
        end do

      do j = 1, ny
         cpp(1,j)  = cpp(2,j)
         cpp(nx,j) = cpp(nx-1,j)
        end do



!**********************************************************************!
!*********  CORRECTED PRESSURE & CORRECTED VELOCITIES *****************!
!**********************************************************************!

      do j = 2, ny-1
         do i = 2, nx-1

            if (flag(i,j) == 0) cycle   ! ? no end if needed

               vp(i,j) = vp(i,j) + cpp(i,j)

               if (flag(i+1,j) == 1) then
               u1(i,j) = u1(i,j) - ((2.0*dt)/(3.0*dxr(i))) * (cpp(i+1,j) - cpp(i,j))
               end if
               
               if (flag(i,j+1) == 1) then
               v1(i,j) = v1(i,j) - ((2.0*dt)/(3.0*dyt(j))) * (cpp(i,j+1) - cpp(i,j))
               end if
         end do
      end do

!**********************************************************************!




!**********************************************************************!
!CONVERGENCE CRITERIA OF VELOCITIES AND TEMPERATURE AT PARTICULAR TIME!
!**********************************************************************!
      bigu = 0.0
      bigv = 0.0
      bigT = 0.0

      do j = 2, ny-1
        do i = 2, nx-1

            if (flag(i,j) == 0) cycle   ! skip solid cells

                bigu = max(bigu, abs(u1(i,j) - u(i,j)))
                bigv = max(bigv, abs(v1(i,j) - v(i,j)))
                bigT = max(bigT, abs(T1(i,j) - T(i,j)))

        end do
      end do

      ! Update solution AFTER convergence check
      do j = 2, ny-1
         do i = 2, nx-1
            if (flag(i,j) == 0) cycle
               u(i,j) = u1(i,j)
               v(i,j) = v1(i,j)
               T(i,j) = T1(i,j)
        end do
       end do

       itrn = itrn + 1
       time = time + dt

       if ((bigu.gt.Error) .or. (bigv.gt.Error) .or. (bigT.gt.Error)) then
          write(*,*) 'bigu=',bigu,' bigv=',bigv,' bigT=',bigT
          stop
        end if

!**********************************************************************!

!**********************************************************************!
!****** CONVERGENCE CRITERIA(upon divergence at a particular time)*****!
!**********************************************************************!
      bigdiv=0.0
      do i=2,nx-1
        do j=2,ny-1
          if (flag(i,j) == 0) cycle
          div=(u1(i,j)-u1(i-1,j))*dy(j)+(v1(i,j)-v1(i,j-1))*dx(i)

          !if(div.gt.bigdiv)bigdiv=div
          bigdiv = max(bigdiv, abs(div))
        end do
      end do

      if (bigdiv.gt.1.0E-4) then
         write(*,*) 'DIVERGENCE ERROR: bigdiv =', bigdiv
         stop
      end if



      if (time.ge.FinalTime) then
         write(*,*) 'Final time reached:', time
         goto 999
        end if




999   continue

      open(70,file='u.dat')
      write(70,*)'variables = "x", "y", "z"'
      write(70,*)'zone i=',ny,', j=',nx,', f=point'

      open(80,file='v.dat')
      write(80,*)'variables = "x", "y", "z"'
      write(80,*)'zone i=',ny,', j=',nx,', f=point'

      open(90,file='p.dat')
      write(90,*)'variables = "x", "y", "z"'
      write(90,*)'zone i=',ny,', j=',nx,', f=point'

      open(100,file='T.dat')
      write(100,*)'variables = "x", "y", "z"'
      write(100,*)'zone i=',ny,', j=',nx,', f=point'


      do i=1,nx
        do j=1,ny
          write(70,*)x(i),y(j),u(i,j)
          write(80,*)x(i),y(j),v(i,j)
          write(90,*)x(i),y(j),vp(i,j)
          write(100,*)x(i),y(j),T(i,j)
        end do
        write(70,*)' '
        write(80,*)' '
        write(90,*)' '
        write(100,*)' '
      end do

      close(70)
      close(80)
      close(90)
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

      do i=1,nx
        do j=2,ny-1
         if (flag(i,j) == 0) then
            si(i,j) = si(i,j-1)   ! keep previous value
            cycle
         end if
         
         si(i,j)=si(i,j-1)-
     &  (dy(j)/4.0)*((u(i,j+1)*dy(j)+u(i,j)*dy(j+1))/dyt(j)
     &  +(u(i,j-1)*dy(j)+u(i,j)*dy(j-1))/dyb(j))
        end do
      end do

      open(130,file='streamLine.dat')
      write(130,*)'variable = "x", "y", "z"'
      write(130,*) 'zone i=',nx-1,', j=',ny,', f=point'


      do i=2,nx
        do j=1,ny
          if (flag(i,j) == 0) cycle
          write(130,*)x(i),y(j),si(i,j)
        end do
        write(130,*)' '
      end do
      close(130)
!**********************************************************************!
!**********************************************************************!
!-------------------------- Average velocity --------------------------!



      ua = 0.0
      area_fluid = 0.0

      do i = 1, nx
         do j = 1, ny

         if (flag(i,j) == 0) cycle   ! skip solid cells

            ua = ua + u(i,j) * dx(i) * dy(j)
            area_fluid = area_fluid + dx(i) * dy(j)

         end do
       end do

       ua = ua / area_fluid

       open(131,file='averagevelocity.dat')
       write(131,*) 'average_velocity = ', ua
       close(131)









!**********************************************************************!
!******************** CALCULATE INVERSE OF MATRIX *********************!
!**********************************************************************!
      contains
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
     &     +M(1,2)*(M(2,3)*M(3,1)-M(2,1)*M(3,3))
     &     +M(1,3)*(M(2,1)*M(3,2)-M(2,2)*M(3,1))

        if(abs(Det) < 1.0d-12)then
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
      end program channel_flow

!***********************************************************************!
!****************************************************************

