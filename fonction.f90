module fonction  !!! Ce module calculera les flux au interfaces
use variable



contains
function get_F(U)
  real*8,dimension(3),intent(in)::U
  real*8,dimension(3)::get_F
  get_F(1)=U(2)
  get_F(2)=U(2)*U(2)/(U(1)) + get_P(U)
  get_F(3)=1./(Gamma-1.)*get_P(U)/U(1)+(U(2)/U(1))**2/2.
end function



function get_b(Ui,Ui_1)
  real*8,dimension(3),intent(in)::Ui,Ui_1
  real*8::get_b,ci,ci_1,vxi,vxi_1
  vxi=Ui(2)/Ui(1)
  vxi_1=Ui_1(2)/Ui_1(1)
  ci=get_c(Ui)
  ci_1=get_c(Ui_1)
  get_b=vxi

  if (vxi+ci>get_b) then
     get_b=vxi+ci
  end if
  if (vxi-ci>get_b) then
    get_b=vxi-ci
  end if
  if (vxi_1>get_b) then
    get_b=vxi_1
  end if
  if (vxi_1+ci_1>get_b) then
     get_b=vxi_1+ci_1
  end if
  if (vxi_1-ci_1>get_b) then
    get_b=vxi_1-ci_1
  end if
end function



function Flux(Ui,Ui_1)
    real*8,dimension(3),intent(in)::Ui,Ui_1
    real*8,dimension(3)::Flux
    Flux=0.5d0*(get_F(Ui_1)+get_F(Ui))-get_b(Ui,Ui_1)*0.5d0*(Ui_1-Ui)
end function

function get_dt(U)
  real*8,dimension(:,:),intent(in)::U
  real*8::b,max,get_dt
  integer::i
  max=0.d0
  do i=0,int(abs(1.d0/get_dx()))
    b=get_b(U(:,i),U(:,i+1))
    if (b>max) then
      max=b
    end if
  end do
  get_dt=0.9d0 *get_dx()/(2.d0*max)
end function


subroutine write(x,U,i)
  real*8,intent(in)::x
  real*8,dimension(3),intent(in)::U
  integer,intent(in)::i
  character*10 :: name
  write(name,'(a7,I3.3)') "Data/t_",i
  if (x==0) then
    open(1,file=name,form="formatted")
  else
    open(1,file=name, form="formatted",position="append")
  end if
  write(1,*) x, U(1), U(2)/U(1),get_P(U),get_T(U)
  close(1)
end subroutine write

end module
