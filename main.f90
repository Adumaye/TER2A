program main
  use variable
  use fonction
  implicit none


  real*8, dimension(:,:),allocatable::U,Unext
  real*8::dt,tf,t,h
  integer:: i,j,Nx,Nt

  real*8::rhog,rhod,pg,pd,vg,vd,mVg,mVd,Eg,Ed
  real*8,dimension(3)::Flux_i

  !!! INITIALISATION  --> /!\ MODIFIE ENERGIE !!
  h=get_dx()
  Nx=int(1./h)+1
  Tf=1.

  Allocate(U(3,0:Nx+1),Unext(3,0:Nx+1))
  rhog=1.
  pg=1.
  vg=0.

  rhod=0.125
  pd=0.1
  vd=0.

  mVg=rhog*(Vg)
  mVd=rhod*(Vd)
  Eg=1/2*vg**2 + pg*1./(get_gamma()-1.)
  Ed=1/2*vd**2 +pd*1./(get_gamma()-1.)

  do i=0,Nx+1,1
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

  !!! BOUCLE EN TEMPS
  t=0.
  Nt=0
  do while (t<tf)
    Unext=0.
    !!! PARCOURE DES INTERFACES
    do i=0,Nx
      FLux_i=Flux(U(:,i),U(:,i+1))
      Unext(:,i)=Unext(:,i)-Flux_i
      Unext(:,i+1)=Unext(:,i+1)+Flux_i
    end do

    dt=get_dt(U)
    U(:,1:Nx)=U(:,1:Nx)+dt/dx*Unext(:,1:Nx)
    do i=0,Nx+1
      call write(i*h,U(:,i),i)
      print*,i*h
    end do
    t=t+dt
    Nt=Nt+1
  end do


  deallocate(U,Unext)

end program main
