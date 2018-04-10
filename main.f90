program main
  use variable
  use fonction
  implicit none


  real*8, dimension(:,:,:),allocatable::U,Unext
  real*8::dt,tf,t,h_X,h_Y
  integer:: i,j,Nx,Ny,Nt

  real*8::rhog,rhod,pg,pd,vg_X,vg_Y,vd_X,vd_Y,mVg_X,mVg_Y,mVd_X,mVd_Y,Eg,Ed
  real*8,dimension(4)::Flux_i,Flux_j

  !!! INITIALISATION  --> /!\ MODIFIE ENERGIE !!
  h_X=get_dx()
  h_Y=get_dy()
  Nx=int(1./h_X)+1
  Ny=int(1./h_Y)+1
  Tf=1.d0

  Allocate(U(4,0:Nx+1,0:Ny+1),Unext(4,0:Nx+1,0:Ny+1))
  rhog=1.d0
  pg=1.d0
  vg_X=0.d0
  vg_Y=0.d0

  rhod=0.125
  pd=0.1
  vd_X=0.d0
  vd_Y=0.d0

  mVg_X=rhog*(Vg_X)
  mVd_X=rhod*(Vd_X)
  mVg_Y=rhog*Vg_Y
  mVd_Y=rhod*Vd_Y
  Eg=1.d0/2.d0*(Vg_X**2 + Vg_Y**2) + pg*1.d0/(get_gamma()-1.d0)
  Ed=1.d0/2.d0*(Vd_X**2 + Vd_Y**2) +pd*1.d0/(get_gamma()-1.d0)

!! CI en X
  ! do i=0,Nx+1
  !   if (i<Nx/2) then
  !     U(1,i,:)=rhog
  !     U(2,i,:)=mVg_X
  !     U(3,i,:)=mVg_Y
  !     U(4,i,:)=Eg
  !   else
  !     U(1,i,:)=rhod
  !     U(2,i,:)=mVd_X
  !     U(3,i,:)=mVd_Y
  !     U(4,i,:)=Ed
  !   end if
  ! end do

!! CI en Y
do j=0,Ny+1
  if (j<Ny/2) then
    U(1,:,j)=rhog
    U(2,:,j)=mVg_X
    U(3,:,j)=mVg_Y
    U(4,:,j)=Eg
  else
    U(1,:,j)=rhod
    U(2,:,j)=mVd_X
    U(3,:,j)=mVd_Y
    U(4,:,j)=Ed
  end if
end do

  !!! BOUCLE EN TEMPS
  t=0.d0
  Nt=0
  do while (t<tf)
    if (t==0) then
      call system ('mkdir Data/')
    end if
    Unext=0.d0
    !!! PARCOURE DES INTERFACES
    do j=0,Ny
      do i=0,Nx
        Flux_i=Flux_X(U(:,i,j),U(:,i+1,j))
        Flux_j=Flux_Y(U(:,i+1,j),U(:,i+1,j+1))
        Unext(:,i,j)=Unext(:,i,j)-Flux_i
        Unext(:,i+1,j)=Unext(:,i+1,j)+Flux_i
        Unext(:,i+1,j)=Unext(:,i+1,j)-Flux_j
        Unext(:,i+1,j+1)=Unext(:,i+1,j+1)+Flux_j
      end do
    end do

    dt=get_dt(U)
    U(:,1:Nx,1:Ny)=U(:,1:Nx,1:Ny)+dt/dx*Unext(:,1:Nx,1:Ny)
    do j=1,Ny
      do i=1,Nx
        call write(i*h_x,j*h_y,U(:,i,j),Nt)
      end do
    end do
    t=t+dt
    Nt=Nt+1
    print*,t
  end do


  deallocate(U,Unext)

end program main
