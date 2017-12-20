module variable  !!! Ce module sert a récupérer toute les données néccéssaire au problème

  !PARAMETRE FIXE DU PROBLEME
  real*8,parameter:: Kb=1.38!*10**(-23)
  real*8,parameter:: N=400
  real*8,parameter:: Gamma=1.4
  real*8,parameter:: dx=10**(-2)


contains

function get_T(U)
  real*8,dimension(3),intent(in)::U
  real*8::get_T
  get_T=U(3)/(N*Kb*1.5)
end function

function get_P(U)
  real*8,dimension(3),intent(in)::U
  real*8::get_P
  get_P=U(3)/(1.5*dx) !on prend dy et dz egaux a 1.
end function

function get_c(U)
  real*8,dimension(3),intent(in)::U
  real*8::get_c
  get_c=sqrt(Gamma*get_P(U)/U(1))
end function

function get_dx()
  real*8::get_dx
  get_dx=dx
end function
end module
