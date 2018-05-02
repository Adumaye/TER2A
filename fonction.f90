module fonction  !!! Ce module calculera les flux au interfaces
  use variable


contains
  ! function get_F_X(U)
  !   real*8,dimension(4),intent(in)::U
  !   real*8,dimension(4)::get_F_X
  !   get_F_X(1)=U(2)
  !   get_F_X(2)= U(2)*U(2)/U(1) + get_P(U)
  !   get_F_X(3)=U(2)*U(3)/U(1)
  !   get_F_X(4)=(U(4) + get_P(U))*U(2)/U(1)
  ! end function
  !
  ! function get_F_Y(U)
  !   real*8,dimension(4),intent(in)::U
  !   real*8,dimension(4)::get_F_Y
  !   get_F_Y(1)=U(3)
  !   get_F_Y(2)=U(2)*U(3)/U(1)
  !   get_F_Y(3)= U(3)*U(3)/U(1) + get_P(U)
  !   get_F_Y(4)=(U(4) + get_P(U))*U(3)/U(1)
  ! end function

  function get_F_X(U)
    real*8,dimension(4),intent(in)::U
    real*8,dimension(4)::get_F_X
    get_F_X(:)=2.d0*U(:)
  end function

  function get_F_Y(U)
    real*8,dimension(4),intent(in)::U
    real*8,dimension(4)::get_F_Y
    get_F_Y(:)=0.d0*U(:)
  end function


  function Flux_X(Ui,Ui_1)
    real*8,dimension(4),intent(in)::Ui,Ui_1
    real*8,dimension(4)::Flux_X
    Flux_X=0.5d0*(get_F_X(Ui_1)+get_F_X(Ui))-get_b_X(Ui,Ui_1)*0.5d0*(Ui_1-Ui)
  end function

  function Flux_Y(Uj,Uj_1)
    real*8,dimension(4),intent(in)::Uj,Uj_1
    real*8,dimension(4)::Flux_Y
    Flux_Y=0.5d0*(get_F_Y(Uj_1)+get_F_Y(Uj))-get_b_Y(Uj,Uj_1)*0.5d0*(Uj_1-Uj)
  end function



  function get_b_X(Ui,Ui_1)
    real*8,dimension(4),intent(in)::Ui,Ui_1
    real*8::get_b_X,ci,ci_1,vxi,vxi_1
    vxi=Ui(2)/Ui(1)
    vxi_1=Ui_1(2)/Ui_1(1)
    ci=get_c(Ui)
    ci_1=get_c(Ui_1)
    get_b_X=abs(vxi)

    if (abs(vxi+ci)>get_b_X) then
      get_b_X=abs(vxi+ci)
    end if
    if (abs(vxi-ci)>get_b_X) then
      get_b_X=abs(vxi-ci)
    end if
    if (abs(vxi_1)>get_b_X) then
      get_b_X=abs(vxi_1)
    end if
    if (abs(vxi_1+ci_1)>get_b_X) then
      get_b_X=abs(vxi_1+ci_1)
    end if
    if (abs(vxi_1-ci_1)>get_b_X) then
      get_b_X=abs(vxi_1-ci_1)
    end if
  end function

  function get_b_Y(Uj,Uj_1)
    real*8,dimension(4),intent(in)::Uj,Uj_1
    real*8::get_b_Y,cj,cj_1,vyj,vyj_1
    vyj=Uj(3)/Uj(1)
    vyj_1=Uj_1(3)/Uj_1(1)
    cj=get_c(Uj)
    cj_1=get_c(Uj_1)
    get_b_Y=abs(vyj)

    if (abs(vyj+cj)>get_b_Y) then
      get_b_Y=abs(vyj+cj)
    end if
    if (abs(vyj-cj)>get_b_Y) then
      get_b_Y=abs(vyj-cj)
    end if
    if (abs(vyj_1)>get_b_Y) then
      get_b_Y=abs(vyj_1)
    end if
    if (abs(vyj_1+cj_1)>get_b_Y) then
      get_b_Y=abs(vyj_1+cj_1)
    end if
    if (abs(vyj_1-cj_1)>get_b_Y) then
      get_b_Y=abs(vyj_1-cj_1)
    end if
  end function


  function get_dt(U)
    real*8,dimension(:,:,:),intent(in)::U
    real*8::a,b,max,get_dt
    integer::i,j
    max=0.d0
    do j=0,int(1.d0/get_dy())
      do i=0,int(1.d0/get_dx())
        b=get_b_X(U(:,i,j),U(:,i+1,j))
        if (b>max) then
          max=b
        end if
        a=get_b_Y(U(:,i+1,j),U(:,i+1,j+1))
        if (a>max) then
          max=a
        end if
      end do
    end do
    get_dt=0.9d0 *( 1.d0/(2.d0*max* (1.d0/get_dx() + 1.d0/get_dy()) ) )
  end function


  subroutine write(x,y,U,i)
    real*8,intent(in)::x,y
    real*8,dimension(4),intent(in)::U
    integer,intent(in)::i
    character*14 :: name
    write(name,'(a7,I3.3,a4)') "Data/t_",i,'.vtk'
    if (x==0 .and. y==0) then
      open(1,file=name,form="formatted")
      write(1,*)  "# vtk DataFile Version 3.0 "
      write(1,*)  "2D Unstructured Grid"
      write(1,*)  "ASCII"
      write(1,*)  "DATASET UNSTRUCTURED_GRID"
      write(1,*)  "POINTS ", (int(1./get_dx())+1)*(int(1./get_dy())+1) ," float "
      write(1,*)  'CELL_DATA ', get_dx()*get_dy()
      write(1,*)  'LOOKUP_TABLE default'

      ! Write(1,'(A26)')"# vtk DataFile Version 3.0"
      ! Write(1,'(A4)')"cell"
      ! Write(1,'(A5)')"ASCII"
      ! Write(1,'(A25)')"DATASET STRUCTURED_POINTS"
      ! Write(1,'(A11,3(I5,1X))')"DIMENSIONS ",(int(1./get_dx())+1),(int(1./get_dy())+1),1
      ! Write(1,'(A12)')"ORIGIN 0 0 0"
      ! Write(1,'(A7,3(ES11.4,1X))')"SPACING",get_dx(),get_dy(),1.D0
      ! Write(1,'(A10,I6)')"POINT_DATA",(int(1./get_dx())+1)*(int(1./get_dy())+1)
      ! Write(1,'(A7,1X,A4,1X,A5)')"SCALARS","cell","float"
      ! Write(1,'(A12,1X,A7)')"LOOKUP_TABLE","default"

    else
      open(1,file=name, form="formatted",position="append")
    end if
  !  Write(1,*) U(1)

   write(1,*) x,y, U(1), U(2)/U(1), U(3)/U(1), get_P(U)
    close(1)
  end subroutine write

end module
