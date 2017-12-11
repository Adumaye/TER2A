program main
  implicit none



  real*8, dimension(3,100)::U,Unext
  real*8::dt,dx,tf
  integer:: i,j,Nx,Nt

  real*8::rhog,rhod,pg,pd,vg,vd
  real*8::mVg,mVd,Eg,Ed

  !! INITIALISATION
  Nx=100
  Nt=10
  tf=1.
  dt=tf/Nt
  dx=1./Nx

  rhog=1.
  pg=1.
  vg=0.

  rhod=0.125
  pd=0.1
  vd=0.

  mVg=rhog*(Vg)
  mVd=rhod*(Vd)
  Eg=1/2*rhog*vg**2
  Ed=1/2*rhod*vd**2

  do i=0,Nx
     if (i<Nx/2) then
        U(1,i)=rhog
        U(2,i)=mVg
        U(3,i)=Eg

     else
        U(1,i)=rhod
        U(2,i)=mVd
        U(3,i)=Ed
     end if
  end do


  do i=0,Nx+1
     print*, i*dx, U(2,i), U(1,i), U(3,i)
  end do

  !!Boucle en temps et en espace

  do j=0,Nt
     do i=1,Nx
        Unext(:,i)=U(:,i)-dt/dx*(F(i,U(:,:))-F(i-1,U(:,:)))
        !print*, F(i,U(:,:)),F(i-1,U(:,:))
     end do
     U=Unext
  end do

  do i=0,Nx+1
     ! print*, i*dx, U(2,i)/U(1,i)
  end do
contains

  Function F(i,U)
    Integer,intent(in)::i
    real*8,dimension(3)::F
    real*8,dimension(3,100),intent(in)::U
    real*8::b
    !F=(1./2.)*(U(:,i+1)+U(:,i))-(b/2.)*(U(:,i+1)-U(:,i))
    if (U(2,i)<0) then
       F=U(2,i+1)/U(1,i+1)
    else
       F=U(2,i)/U(1,i)
    end if
  end Function F


end program main

