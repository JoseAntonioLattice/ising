program main

  use iso_fortran_env, only: dp => real64
  implicit none

  integer, parameter :: L = 60, Nm = 1000, Nterm = 500, Nt = 10, Nskip = 10
  real(dp) :: r
  integer :: i, isweeps,iskip, it
  integer :: spin(L,L), ip(L), im(L)
  real(dp) :: E(Nm), M(Nm)
  real(dp) :: T(Nt), beta(Nt)
  real(dp), parameter :: Tmin = 0.1_dp, Tmax = 4.0_dp, DT = Tmax - Tmin
  
  ip = [(i+1, i = 1, L)] ; ip(L) = 1
  im = [(i-1, i = 1, L)] ; im(1) = L
  T = [(Tmin + DT*i/(Nt - 1), i = 0, Nt-1)]
  beta = 1/T

  spin = 1

  call random_init(.true.,.false.)
  
  open(unit = 10, file = "datosF.dat")
  do it = 1, Nt
     do isweeps = 1, Nterm
        call sweeps(spin,beta(it))
     end do

     do isweeps = 1, Nm
        do iskip = 1, Nskip
           call sweeps(spin,beta(it))
        end do
        E(isweeps) = energy_density(spin)
        M(isweeps) = 1.0_dp*abs(sum(spin))/L**2
     end do
     write(10,*) T(it), avr(E),stderr(E), avr(M), stderr(M)
  end do
contains

  function avr(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: avr

    avr = sum(x)/size(x)
  end function avr

  function var(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: var, avg

    avg = avr(x)
    var = sum((x-avg)**2)/(size(x) - 1)
    
  end function var

  function stderr(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: stderr

    stderr = sqrt(var(x)/size(x))
  end function stderr
  
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

  function energy_density(spin)
    integer, intent(in) :: spin(L,L)
    real(dp) :: energy_density
    integer :: E, i, j

    E = 0
    do i = 1, L
       do j = 1, L
          E = E + spin(i,j) * (spin(ip(i),j) + spin(i,ip(j)))
       end do
    end do
    energy_density = -real(E,dp)/L**2
        
  end function energy_density
  
end program main
