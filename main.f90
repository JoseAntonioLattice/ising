program main

  use iso_fortran_env, only: dp => real64
  implicit none

  integer, parameter :: L = 16, Nm = 10000, Nterm = 2000, Nt = 20, Nskip = 100
  real(dp) :: r
  integer :: i, isweeps,iskip, it
  integer :: spin(L,L), ip(L), im(L), En
  real(dp) :: E(Nm), M(Nm)
  real(dp) :: T(Nt), beta(Nt)
  real(dp), parameter :: Tmin = 0.01_dp, Tmax = 1.0_dp, DT = Tmax - Tmin
  real(dp) :: q = 1.1_dp
  
  
  ip = [(i+1, i = 1, L)] ; ip(L) = 1
  im = [(i-1, i = 1, L)] ; im(1) = L
  T = [(Tmin + DT*i/(Nt - 1), i = 0, Nt-1)]
  beta = 1/T

  spin = 1

  
  En = energy(spin)
  open(unit = 10, file = "q=1.1.dat")
  do it = 1, Nt
     do isweeps = 1, Nterm
        call sweeps(spin,beta(it),q,En)
     end do

     do isweeps = 1, Nm
        do iskip = 1, Nskip
           call sweeps(spin,beta(it),q,En)
        end do
        E(isweeps) = energy_density(spin)
        M(isweeps) = 1.0_dp*abs(sum(spin))/L**2
     end do
     write(10,*) T(it), avr(E),stderr(E), jackknife(1.0_dp*E,20), avr(M), stderr(M), jackknife(1.0_dp*M,20)
     write(*,*) T(it), avr(E),stderr(E), jackknife(1.0_dp*E,20), avr(M), stderr(M), jackknife(1.0_dp*M,20)
     FLUSH(10)
     !write(*,*) spin
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

  function jackknife(x,bins)
    real(dp) :: jackknife, x(:)
    integer, intent(in) :: bins
    integer :: MM, NN, i
    real(dp) :: xbar, sum_x
    real(dp) :: x_m(bins)


    NN = size(x)
    MM = NN/bins

    xbar = avr(x)
    sum_x = sum(x)
    x_m = 1.0_dp/(NN-MM) * [(sum_x - sum(x(MM*(i-1)+1:MM*i)),i=1,bins)]

    jackknife = sqrt( real(bins - 1,dp)/bins * sum( (x_m - xbar)**2) )
  end function jackknife
  
  subroutine sweeps(spin,beta,q,E)
    integer, intent(inout) :: spin(L,L)
    real(dp), intent(in) ::  beta
    integer, intent(inout) :: E
    integer :: i, j
    real(dp) :: q
    
    do i = 1, L
       do j = 1, L
          call metropolis(spin,[i,j],beta,q,E)
       end do
    end do
  end subroutine sweeps

  subroutine metropolis(spin,x,beta,q, E)
    integer, intent(inout) :: spin(L,L)
    real(dp), intent(in) :: beta, q
    integer, intent(in) :: x(2)
    integer, intent(inout) :: E
    integer :: DH, H
    real(dp) :: r, p

    !DH = DE(spin,x)
    !if( DH <= 0 )then
    !   spin(x(1),x(2)) = -spin(x(1),x(2))
    !else
    !   call random_number(r)
    !   if( r <= exp(-DH*beta)) spin(x(1),x(2)) = -spin(x(1),x(2))
    !end if
    
    H = E!energy(spin) 
    DH = DE(spin,x)
    !p = min(1.0_dp, (expq(-beta*(H+DH),q)/expq(-beta*H,q))**q)
    p = min(1.0_dp,(expq(1.0_dp*DH/((1.0_dp-q)*H-1/beta),q)**q))
    call random_number(r)
    if( r <= p ) then
       spin(x(1),x(2)) = -spin(x(1),x(2))
       E = E + DH
    end if
  end subroutine metropolis

  function DE(spin,x)
    integer, intent(in) :: spin(L,L)
    integer, intent(in) :: x(2)
    integer :: DE, i,j

    i = x(1)
    j = x(2)

    DE = 2*spin(i,j)*(spin(ip(i),j) + spin(i,ip(j)) + spin(im(i),j) + spin(i,im(j)))
    
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
    energy_density = 2.0_dp-real(E,dp)/L**2
        
  end function energy_density

    function energy(spin) result(E)
    integer, intent(in) :: spin(L,L)
    integer :: E, i, j

    E = 0
    do i = 1, L
       do j = 1, L
          E = E + spin(i,j) * (spin(ip(i),j) + spin(i,ip(j)))
       end do
    end do
    E = 2*L**2-E
  end function energy
  
  function expq(x,q)
    real(dp), intent(in) :: x, q
    real(dp) :: expq, arg
    arg = 1.0_dp - q
    expq = (1.0_dp + arg*x)**(1/arg)
    
  end function expq

  function escort_prob(x,q)
    real(dp), intent(in) :: x, q 
    real(dp) :: escort_prob
    escort_prob = (expq(x,q))**q
  end function escort_prob
  
end program main
