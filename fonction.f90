module fonction  !!! Ce module calculera les flux au interfaces
  use variable


contains
  function get_F(U)
    real*8,dimension(4),intent(in)::U
    real*8,dimension(4)::get_F
    get_F(1)=U(2)
    get_F(2)= U(2)*U(2)/U(1) + get_P(U)
    get_F(3)=U(2)*U(3)/U(1)
    get_F(4)=(1./(Gamma-1.)*get_P(U)/U(1)+ get_V(U)**2/2. + get_P(U))*U(2)/U(1)
  end function



  function get_b_X(Ui,Ui_1)
    real*8,dimension(4),intent(in)::Ui,Ui_1
    real*8::get_b_X,ci,ci_1,vxi,vxi_1
    vxi=Ui(2)/Ui(1)
    vxi_1=Ui_1(2)/Ui_1(1)
    ci=get_c(Ui)
    ci_1=get_c(Ui_1)
    get_b_X=vxi

    if (vxi+ci>get_b_X) then
      get_b_X=vxi+ci
    end if
    if (vxi-ci>get_b_X) then
      get_b_X=vxi-ci
    end if
    if (vxi_1>get_b_X) then
      get_b_X=vxi_1
    end if
    if (vxi_1+ci_1>get_b_X) then
      get_b_X=vxi_1+ci_1
    end if
    if (vxi_1-ci_1>get_b_X) then
      get_b_X=vxi_1-ci_1
    end if
  end function

  function get_b_Y(Uj,Uj_1)
    real*8,dimension(4),intent(in)::Uj,Uj_1
    real*8::get_b_Y,cj,cj_1,vyj,vyj_1
    vyj=Uj(2)/Uj(1)
    vyj_1=Uj_1(2)/Uj_1(1)
    cj=get_c(Uj)
    cj_1=get_c(Uj_1)
    get_b_Y=vyj

    if (vyj+cj>get_b_Y) then
      get_b_Y=vyj+cj
    end if
    if (vyj-cj>get_b_Y) then
      get_b_Y=vyj-cj
    end if
    if (vyj_1>get_b_Y) then
      get_b_Y=vyj_1
    end if
    if (vyj_1+cj_1>get_b_Y) then
      get_b_Y=vyj_1+cj_1
    end if
    if (vyj_1-cj_1>get_b_Y) then
      get_b_Y=vyj_1-cj_1
    end if
  end function

  function Flux_X(Ui,Ui_1)
    real*8,dimension(4),intent(in)::Ui,Ui_1
    real*8,dimension(4)::Flux_X
    Flux_X=0.5d0*(get_F(Ui_1)+get_F(Ui))-get_b_X(Ui,Ui_1)*0.5d0*(Ui_1-Ui)
  end function

  function Flux_Y(Uj,Uj_1)
    real*8,dimension(4),intent(in)::Uj,Uj_1
    real*8,dimension(4)::Flux_Y
    Flux_Y=0.5d0*(get_F(Uj_1)+get_F(Uj))-get_b_Y(Uj,Uj_1)*0.5d0*(Uj_1-Uj)
  end function

  function get_dt(U)
    real*8,dimension(:,:,:),intent(in)::U
    real*8::a,b,max,get_dt
    integer::i,j
    max=0.d0
    do j=0,int(abs(1.d0/get_dy()))
      do i=0,int(abs(1.d0/get_dx()))
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
    get_dt=0.9d0 *get_dx()/(2.d0*max)
  end function


  subroutine write(x,y,U,i)
    real*8,intent(in)::x,y
    real*8,dimension(4),intent(in)::U
    integer,intent(in)::i
    character*10 :: name
    write(name,'(a7,I3.3)') "Data/t_",i
    if (x==0 .and. y==0) then
      open(1,file=name,form="formatted")
    else
      open(1,file=name, form="formatted",position="append")
    end if
    write(1,*) x,y, U(1), U(2)/U(1),get_P(U),get_T(U)
    close(1)
  end subroutine write

end module
