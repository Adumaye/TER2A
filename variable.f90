module variable  !!! Ce module sert a récupérer toute les données néccéssaire au problème

  !PARAMETRE FIXE DU PROBLEME
  real*8,parameter:: Kb=1.38d0
  real*8,parameter:: N=400.d0
  real*8,parameter:: Gamma=1.4d0
  real*8,parameter:: dx=10.d0**(-2)
  real*8,parameter:: dy=10.d0**(-2)


contains
! function get_V(U)
!   real*8,dimension(4),intent(in)::U
!   real*8:: get_V
!   get_V= sqrt( (U(2)/U(1))**2 + (U(3)/U(1))**2 )
! end function

! function get_T(U)
!   real*8,dimension(4),intent(in)::U
!   real*8::get_T
!   get_T=get_P(U)/(N*Kb)
! end function

function get_P(U)
  real*8,dimension(4),intent(in)::U
  real*8::get_P
  ! get_P=(U(4)- get_V(U)**2 /2.d0)*U(1)*(Gamma-1.d0)
  get_P=(2.d0*U(4)-(U(2)**2+U(3)**2)/U(1)) /5.d0
end function

function get_c(U)
  real*8,dimension(4),intent(in)::U
  real*8::get_c
  get_c=sqrt(Gamma*get_P(U)/U(1))
end function

function get_dx()
  real*8::get_dx
  get_dx=dx
end function

function get_dy()
  real*8::get_dy
  get_dy=dy
end function

! function get_gamma()
!   real*8::get_gamma
!   get_gamma=Gamma
! end function
end module
