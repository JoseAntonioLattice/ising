program main

  use iso_fortran_env, only: dp => real64
  implicit none

  integer, parameter :: L = 8, Nm = 10**3, Nterm = 500, Nt = 100, Nskip = 10
  real(dp) :: r
  integer :: i, isweeps,iskip, it
  integer :: spin(L,L), ip(L), im(L)
  real(dp) :: T(Nt), beta(Nt)
  real(dp), parameter :: Tmin = 0.1_dp, Tmax = 4.0_dp, DT = Tmax - Tmin
  
  ip = [(i+1, i = 1, L)] ; ip(L) = 1
  im = [(i-1, i = 1, L)] ; im(1) = L
  T = [(Tmin + DT*i/(Nt - 1), i = 0, Nt-1)]
  beta = 1/T

  spin = 1

  do it = 1, Nt
     do isweeps = 1, Nterm
        call sweeps(spin,beta(it))
     end do

     do isweeps = 1, Nterm
        do iskip = 1, Nskip
           call sweeps(spin,beta(it))
        end do
     end do    
  end do
contains

  subroutine sweeps(spin,beta)
    integer, intent(inout) :: spin(L,L)
    real(dp), intent(in) ::  beta
    integer :: i, j

    do i = 1, L
       do j = 1, L
          call metropolis(spin,[i,j],beta)
       end do
    end do
  end subroutine sweeps

  subroutine metropolis(spin,x,beta)
    integer, intent(inout) :: spin(L,L)
    real(dp), intent(in) :: beta
    integer, intent(in) :: x(2)
    integer :: DH
    real(dp) :: r

    DH = DE(spin,x)

    if( DH <= 0 )then
       spin(x(1),x(2)) = -spin(x(1),x(2))
    else
       call random_number(r)
       if( r <= exp(-DH*beta)) spin(x(1),x(2)) = -spin(x(1),x(2))
    end if
    
  end subroutine metropolis

  function DE(spin,x)
    integer, intent(in) :: spin(L,L)
    integer, intent(in) :: x(2)
    integer :: DE, i,j

    i = x(1)
    j = x(2)

    DE = 2*spin(i,j)*(spin(ip(i),j) + spin(i,ip(j))+spin(im(i),j) + spin(i,im(j)))
    
  end function DE
  
end program main
