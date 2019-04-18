!***************************************************************************
subroutine fugacities( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    x, & ! vector of mole fractions
    type_k, & ! vector of binary interaction types
    coef_ab,& ! vector of a, b coefficients
    fuga,   & ! vector of fugacity coefficients (output results)
    Gex,    & ! value of excess Gibbs energy (output results)
    delta   & ! matrix of binary interaction parameters
    )
!{
  implicit none
  integer :: n
  real(8) :: P,T,Pc(n),Tc(n),w(n),x(n),fuga(n),Gex
  integer :: type_k(n)
  real(8) :: coef_ab(n)
  real(8) :: delta(n,n)
  real(8) :: lnphi(n),V

  call fugacities_general(P,T,n,Pc,Tc,w,x,type_k,coef_ab,V,lnphi,fuga,Gex,delta)
!}
end subroutine fugacities

subroutine ln_fugacities( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    x, & ! vector of mass fractions
    type_k, & ! vector of binary interaction types
    coef_ab,& ! vector of a, b coefficients
    lnphi   & ! vector of ln(fugacity coefficients) (output results)
    )
!{
  implicit none
  integer :: n
  real(8) :: P,T,Pc(n),Tc(n),w(n),x(n),fuga(n),Gex
  integer :: type_k(n)
  real(8) :: coef_ab(n)
  real(8) :: delta(n,n)
  real(8) :: lnphi(n), V

  call fugacities_general(P,T,n,Pc,Tc,w,x,type_k,coef_ab,V,lnphi,fuga,Gex,delta)
!}
end subroutine ln_fugacities

subroutine molar_volume( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    x, & ! vector of mass fractions
    type_k, & ! vector of binary interaction types
    coef_ab,& ! vector of a, b coefficients
    V       &
    )
!{
  implicit none
  integer :: n
  real(8) :: P,T,Pc(n),Tc(n),w(n),x(n),fuga(n),Gex
  integer :: type_k(n)
  real(8) :: coef_ab(n)
  real(8) :: delta(n,n)
  real(8) :: lnphi(n), V

  call fugacities_general(P,T,n,Pc,Tc,w,x,type_k,coef_ab,V,lnphi,fuga,Gex,delta)
!}
end subroutine molar_volume

subroutine molar_volume2( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    kij, & ! matrix of BIPs
    x, & ! vector of mass fractions
    V  & ! output molar volume
    )
!{
  implicit none
  integer :: n
  real(8) :: P,T,Pc(n),Tc(n),w(n),x(n)  
  real(8) :: kij(n,n)
  real(8) :: V

  integer :: i,j
  real(8) :: a(n),b(n),am,bm
  character :: state = 'm'

  ! PR coefficients of pure components

  call calculate_a_b2(T,n,Pc,Tc,w,a,b)

  ! PR coefficients of the mixture
  bm = 0d0
  am = 0d0

  do i=1,n
    bm      = bm + x(i)*b(i)
    am      = am + x(i)**2*a(i)

    do j=i+1,n
      am    = am + 2.*x(i)*x(j)*(1.-kij(i,j))*dsqrt(a(i)*a(j))
    enddo
  enddo  

  call PR_vol(P,T,am,bm,V,state)
!}
end subroutine molar_volume2

subroutine fugacities_general( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    x, & ! vector of mass fractions
    type_k, & ! vector of binary interaction types
    coef_ab,& ! vector of a, b coefficients
    V,      & ! value of volume
    lnphi,  & ! vector of log(phi)
    fuga,   & ! vector of fugacity coefficients (output results)
    Gex,    & ! value of excess Gibbs energy (output results)
    delta   & ! matrix of binary interaction parameters
    )
!{
  implicit none
  integer :: n
  real(8) :: P,T,Pc(n),Tc(n),w(n),x(n),fuga(n),Gex
  integer :: type_k(n)
  real(8) :: coef_ab(n)

  real(8) :: a(n), b(n), delta(n,n) ! delta is binary interaction parameter (k_ij)
  real(8) :: am, bm, Tr, kappa, alpha
  real(8) :: lnphi(n), V, Z 
  real(8) :: lnphip(n),vv(n),zi

  integer :: i,j
  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0
  character :: state = 'm'

  ! PR coefficients of pure components

  call calculate_a_b(T,n,Pc,Tc,w,coef_ab,a,b)
  !call calculate_kij(0,T,n,Pc,Tc,w,type_k,delta)
  call calculate_kij_from_table(0,T,n,delta)

  ! PR coefficients of the mixture
  bm = 0d0
  am = 0d0

  do i=1,n
    bm      = bm + x(i)*b(i)
    am      = am + x(i)**2*a(i)

    do j=i+1,n
      am    = am + 2.*x(i)*x(j)*(1.-delta(i,j))*dsqrt(a(i)*a(j))
    enddo
  enddo  

  call PR_vol(P,T,am,bm,V,state)
  !print *, "x1,P,V=",x(1),P,V
  do i=1,n
    call PR_vol(P,T,a(i),b(i),vv(i),state)
  enddo

  Z = P*V/R_gas/T
  Gex = 0.0d0
  do i=1,n
  !{
    lnphi(i)=-b(i)/bm
    do j=1,n
      lnphi(i) = lnphi(i) +2d0*x(j)*dsqrt(a(i)*a(j))*(1.-delta(i,j))/am
    enddo

    lnphi(i) = - lnphi(i) * sqr2*.25d0 * am/bm/R_gas/T &
                 *dlog( (V+bm*(1d0+sqr2))/(V+bm*(1d0-sqr2)) ) &
             - dlog(Z-bm/V*Z) +b(i)/bm*(Z-1d0)

    fuga(i) = dexp(lnphi(i))

    zi = P*vv(i)/R_gas/T
    lnphip(i) = - sqr2/4. * a(i)/b(i)/R_gas/T &
                 *dlog( (vv(i)+b(i)*(1.+sqr2))/(vv(i)+b(i)*(1.-sqr2)) ) &
             - dlog(zi-b(i)/vv(i)*zi) +(zi-1.)

    Gex = Gex + x(i)*lnphi(i) - x(i)*lnphip(i)
    if(x(i)>0) Gex = Gex + x(i)*dlog(x(i))
    !print *, "fuga(", i, ")=", fuga(i)

  !}
  enddo
!}
end subroutine fugacities_general

subroutine fugacities2( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    kij, & ! matrix of BIPs
    x, & ! vector of mole fractions
    fuga, & ! vector of fugacity coefficients (output)
    Gex  & ! value of excess Gibbs energy (output)
    )
!{
  implicit none
  integer :: n
  real(8) :: P,T,Pc(n),Tc(n),w(n),kij(n,n),x(n),fuga(n),Gex

  real(8) :: a(n),b(n)
  real(8) :: am,bm,Tr,kappa,alpha
  real(8) :: lnphi(n),V,Z 
  real(8) :: lnphip(n),vv(n),zi

  integer :: i,j
  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0
  character :: state = 'm'

  ! PR coefficients of pure components

  call calculate_a_b2(T,n,Pc,Tc,w,a,b)

  ! PR coefficients of the mixture
  bm = 0d0
  am = 0d0

  do i=1,n
    bm      = bm + x(i)*b(i)
    am      = am + x(i)**2*a(i)

    do j=i+1,n
      am    = am + 2.*x(i)*x(j)*(1.-kij(i,j))*dsqrt(a(i)*a(j))
    enddo
  enddo  

  call PR_vol(P,T,am,bm,V,state)
  do i=1,n
    call PR_vol(P,T,a(i),b(i),vv(i),state)
  enddo

  Z = P*V/R_gas/T
  Gex = 0.0d0
  do i=1,n
  !{
    lnphi(i)=-b(i)/bm
    do j=1,n
      lnphi(i) = lnphi(i) +2d0*x(j)*dsqrt(a(i)*a(j))*(1.-kij(i,j))/am
    enddo

    lnphi(i) = - lnphi(i) * sqr2*.25d0 * am/bm/R_gas/T &
                 *dlog( (V+bm*(1d0+sqr2))/(V+bm*(1d0-sqr2)) ) &
             - dlog(Z-bm/V*Z) +b(i)/bm*(Z-1d0)

    fuga(i) = dexp(lnphi(i))

    zi = P*vv(i)/R_gas/T
    lnphip(i) = - sqr2/4. * a(i)/b(i)/R_gas/T &
                 *dlog( (vv(i)+b(i)*(1.+sqr2))/(vv(i)+b(i)*(1.-sqr2)) ) &
             - dlog(zi-b(i)/vv(i)*zi) +(zi-1.)

    Gex = Gex + x(i)*lnphi(i) - x(i)*lnphip(i)
    if(x(i)>0) Gex = Gex + x(i)*dlog(x(i))
  !}
  enddo
!}
end subroutine fugacities2

subroutine fugacities3( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    kij, & ! matrix of BIPs
    x, & ! vector of mole fractions
    lnphi  & ! vector of fugacity coefficients (output)
    )
!{
  implicit none
  integer :: n
  real(8) :: P,T,Pc(n),Tc(n),w(n),kij(n,n),x(n),lnphi(n)

  real(8) :: a(n),b(n)
  real(8) :: am,bm,Tr,kappa,alpha
  real(8) :: V,Z

  integer :: i,j
  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0
  character :: state = 'm'

  ! PR coefficients of pure components

  call calculate_a_b2(T,n,Pc,Tc,w,a,b)

  ! PR coefficients of the mixture
  bm = 0d0
  am = 0d0

  do i=1,n
    bm      = bm + x(i)*b(i)
    am      = am + x(i)**2*a(i)

    do j=i+1,n
      am    = am + 2.*x(i)*x(j)*(1.-kij(i,j))*dsqrt(a(i)*a(j))
    enddo
  enddo  

  call PR_vol(P,T,am,bm,V,state)  

  Z = P*V/R_gas/T
  do i=1,n
  !{
    lnphi(i)=-b(i)/bm
    do j=1,n
      lnphi(i) = lnphi(i) +2d0*x(j)*dsqrt(a(i)*a(j))*(1.-kij(i,j))/am
    enddo

    lnphi(i) = - lnphi(i) * sqr2*.25d0 * am/bm/R_gas/T &
                 *dlog( (V+bm*(1d0+sqr2))/(V+bm*(1d0-sqr2)) ) &
             - dlog(Z-bm/V*Z) +b(i)/bm*(Z-1d0)    
  !}
  enddo
!}
end subroutine fugacities3
!***************************************************************************
subroutine pressure_PR_EoS( &
    rho, & ! density (Unit: kg/m^3)
    T,   & ! temperature (Unit: K)
    n,   & ! number of species
    Pc,  & ! vector of critical pressures
    Tc,  & ! vector of critical temperatures
    w,   & ! vector of acentric factors
    M,   & ! vector of molecular weights
    x,   & ! vector of mass fractions
    type_k, & ! vector of binary interaction types
    coef_ab,& ! vector of a, b coefficients
    P    & ! pressure (Unit: Pa) (OUTPUT)
    )
!{
  implicit none
  integer :: n
  real(8) :: rho,P,T,Pc(n),Tc(n),w(n),M(n),x(n)
  integer :: type_k(n)
  real(8) :: coef_ab(n)

  real(8) :: a(n), b(n), delta(n,n)
  real(8) :: am, bm, Tr, kappa, alpha, V

  integer :: i,j
  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0
  character :: state = 'm'

  ! PR coefficients of pure components
  call calculate_a_b(T,n,Pc,Tc,w,coef_ab,a,b)
  call calculate_kij(0,T,n,Pc,Tc,w,type_k,delta)

  V = 0.0d0
  do i=1,n
    V = V + x(i) * M(i)
  enddo
  V = V*1d-3/rho
  !print *,"V,rho=",V,rho

  ! PR coefficients of the mixture
  bm = 0d0
  am = 0d0

  do i=1,n
    bm      = bm + x(i)*b(i)
    am      = am + x(i)**2*a(i)

    do j=i+1,n
      am    = am + 2.*x(i)*x(j)*(1.-delta(i,j))*dsqrt(a(i)*a(j))
    enddo
  enddo  

  P = R_gas*T/(V-bm) - am/(V*V+2d0*bm*V-bm*bm)
!}
end subroutine pressure_PR_EoS
!***************************************************************************
subroutine density_PR_EoS( &
    P,   & ! pressure (Unit: Pa) 
    T,   & ! temperature (Unit: K)
    n,   & ! number of species
    Pc,  & ! vector of critical pressures
    Tc,  & ! vector of critical temperatures
    w,   & ! vector of acentric factors
    M,   & ! vector of molecular weights
    x,   & ! vector of mass fractions
    type_k, & ! vector of binary interaction types
    coef_ab,& ! vector of a, b coefficients
    rho  & ! density (Unit: kg/m^3) (OUTPUT)
    )
!{
  implicit none
  integer :: n
  real(8) :: rho,P,T,Pc(n),Tc(n),w(n),M(n),x(n)
  integer :: type_k(n)
  real(8) :: coef_ab(n)

  real(8) :: a(n), b(n), delta(n,n)
  real(8) :: am, bm, Tr, kappa, alpha, V

  integer :: i,j
  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0
  character :: state = 'm'

  ! PR coefficients of pure components
  call calculate_a_b(T,n,Pc,Tc,w,coef_ab,a,b)
  !call calculate_kij(0,T,n,Pc,Tc,w,type_k,delta)
  call calculate_kij_from_table(0,T,n,delta)

  rho = 0.0d0
  do i=1,n
    rho = rho + x(i) * M(i)
  enddo

  ! PR coefficients of the mixture
  bm = 0d0
  am = 0d0

  do i=1,n
    bm = bm + x(i)*b(i)
    am = am + x(i)**2*a(i)

    do j=i+1,n
      am = am + 2.*x(i)*x(j)*(1.-delta(i,j))*dsqrt(a(i)*a(j))
    enddo
  enddo  

  !P = R_gas*T/(V-bm) - am/(V*V+2d0*bm*V-bm*bm)
  call PR_vol(P,T,am,bm,V,state)

  rho = rho*1d-3/V
  !print *,"V,rho=",V,rho

!}
end subroutine density_PR_EoS

!***************************************************************************
!***************************************************************************
subroutine density_PR_EoS2( &
    P,   & ! pressure (Unit: Pa) 
    T,   & ! temperature (Unit: K)
    x,   & ! vector of mass fractions
    n,   & ! number of species
    Pc,  & ! vector of critical pressures
    Tc,  & ! vector of critical temperatures
    w,   & ! vector of acentric factors
    M,   & ! vector of molecular weights
    kij, & ! matrix of BIPs
    rho  & ! density (Unit: kg/m^3) (OUTPUT)
    )
!{
  implicit none
  integer :: n
  real(8) :: P,T,x(n),Pc(n),Tc(n),w(n),M(n),kij(n,n)
  real(8) :: rho

  real(8) :: a(n), b(n)
  real(8) :: am, bm, Tr, kappa, alpha, V

  integer :: i,j
  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0
  character :: state = 'm'

  ! PR coefficients of pure components
  call calculate_a_b2(T,n,Pc,Tc,w,a,b)

  rho = 0.0d0
  do i=1,n
    rho = rho + x(i) * M(i)
  enddo

  ! PR coefficients of the mixture
  bm = 0d0
  am = 0d0

  do i=1,n
    bm = bm + x(i)*b(i)
    am = am + x(i)**2*a(i)

    do j=i+1,n
      am = am + 2.*x(i)*x(j)*(1.-kij(i,j))*dsqrt(a(i)*a(j))
    enddo
  enddo  

  !P = R_gas*T/(V-bm) - am/(V*V+2d0*bm*V-bm*bm)
  call PR_vol(P,T,am,bm,V,state)

  rho = rho*1d-3/V
  !print *,"V,rho=",V,rho

!}
end subroutine density_PR_EoS2

!***************************************************************************
subroutine findEquilibrium_Fugacity( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    type_k, & ! vector of binary interaction types
    coef_ab,& ! vector of a, b coefficients
    x_a, & ! vector of mass fractions: alpha phase
    x_b  & ! vector of mass fractions: beta phase
    )
!{
  implicit none
  integer :: n
  integer :: type_k(n)
  real(8) :: coef_ab(n), delta(n,n), f1, f2, xa1_old, xb1_old, dx
  real(8) :: P,T,Pc(n),Tc(n),w(n),x_a(n),x_b(n),fuga_a(n),fuga_b(n)
  real(8) :: G_a, G_b, s
  integer :: i

  ! Gex Method for equilibrium
  integer, parameter :: m=2000
  real(8) :: x1(n),x2(n),xo1(n),xo2(n),G(0:m),x(0:m,n),G1,xr,dGdx, tmp_d2Gdx2
  integer :: j,j1,i1,counter

  x1=(/1. , 0. /)
  x2=(/0. , 1. /)
                
  dx = 1d0/dfloat(m)
  do i=0,m 
    xr = dfloat(i)/dfloat(m)
    !xr = xr**3
    x(i,:) = xr*x1(:)+(1d0-xr)*x2(:)
  enddo
 do i=0,m
   call fugacities(P,T,n,Pc,Tc,w,x(i,:),type_k,coef_ab,fuga_a,G(i),delta)
   !print *, x(i,1), G(i)!, tmp_d2Gdx2
 enddo
     
  dGdx = 1d10
  !print *, 'd2Gdx2 =', dGdx
  do i=1,m-1
    !tmp_d2Gdx2 = ( (G(i+1)-G(i))/(x(i+1,1)-x(i,1)) - (G(i)-G(i-1))/(x(i,1)-x(i-1,1)) ) &
    !           / ( x(i+1,1) - x(i-1,1) ) * 2.0d0
    tmp_d2Gdx2 = (G(i+1)-2d0*G(i)+G(i-1))/dx**2 
    !print *, x(i,1), G(i), tmp_d2Gdx2
    if (dGdx > tmp_d2Gdx2 ) then
      dGdx = tmp_d2Gdx2
      i1 = i
    endif
  enddo
  !print *, 'd2Gdx2 =', dGdx
  if (dGdx .gt. 0) then
    x_a(1) = 0.0d0
    x_a(1) = 1.0/x_a(1)
    x_a(2) = x_a(1)
    x_b(1) = x_a(1)
    x_b(2) = x_a(1)
    return  
  endif
       
  i = i1-1
  j = i1+1 
          
  do while ( G(j) + (G(j+1)-G(j)) / (x(j+1,1)-x(j,1)) * (x(i,1)-x(j,1)) >=G(i) .and. j<m)
    j = j+1
    do while ( G(i) + (G(i-1)-G(i)) / (x(i-1,1)-x(i,1)) * (x(j,1)-x(i,1)) >=G(j) .and. i>1)
      i = i-1
    enddo
  enddo
      
  x2=x(i,:)
  x1=x(j,:)
  x_a = x1 
  x_b = x2
         
  !open(unit=888,file='func_G.dat')
  !do i=0,m
  !  write(888,*) x(i,1), G(i)
  !enddo
  !print*, 'i=',i,'j=',j
  !print*, 'x2=',x(i,1),'x1',x(j,1)
  !!stop
       
  call fugacities(P,T,n,Pc,Tc,w,x_a,type_k,coef_ab,fuga_a,G_a,delta)
  call fugacities(P,T,n,Pc,Tc,w,x_b,type_k,coef_ab,fuga_b,G_b,delta)
  print *, "LLE Residue @ alpha =", x_a(1)*fuga_a(1) -x_b(1)*fuga_b(1)
  print *, "LLE Residue @ beta  =", x_a(2)*fuga_a(2) -x_b(2)*fuga_b(2)
!}
end subroutine findEquilibrium_Fugacity

!***************************************************************************
! function PR_vol
! Calculates volume from Peng-Robinson EOS at given pressure and temperature
! 
! Input:  
!      P: pressure in Pa
!      T: temperature in K
!      a,b: coefficients of PR EOS
!      optional:
!      state: 'L' for liquid and 'V' or 'G' for vapor
!          if not specifoed, the state with minimum Gibbs energy
!
! Output:  
!      V: volume in m^3/mol
!***************************************************************************
subroutine PR_vol(P,T,a,b,PRV,state)
!{
  implicit none
  !character, intent(in), optional :: state
  character :: state
  real(8) :: P,T,a,b,PRV, &
             a1,a2,a3,a4, &
             V,Vmin,G,Gmin,Z,VV,VL
  complex(8),dimension(3) :: vol
  integer :: i,j
  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0

  ! Make up the coefficient vector for the cubic volume equation
  ! P*V^3 + (P*b-R*T)*V^2 + (a-3*P*b^2-2*R*T*b)*V + (P*b^3+R*T*b^2-a*b)=0

  a1 = P                       ! V^3 term
  a2 = P*b-R_gas*T             ! V^2 term
  a3 = a-3.d0*P*b*b-2.d0*R_gas*T*b ! V term
  a4 = P*b*b*b+R_gas*T*b*b-a*b ! constant term

  call cubic_solve(a1,a2,a3,a4,vol)

  !  print *,vol(1), vol(2), vol(3)

  if (imag(vol(2)) .ne. 0d0) then
  !{
    ! only one real root
    ! if (present(state)) print*,'Only one root; no option for state.'

    PRV = real(vol(1))
    ! V=Vmin
    !Gmin = -R_gas*T*dlog(v-b)+R_gas*T*b/(v-b)&
    !     + 0.25*a*sqr2/b*dlog((v+b-sqr2*b)/(v+b+sqr2*b)) &
    !     - a*v/(v**2+2*b*v-b**2)
  !}
  else !if (present(state)) then
  !{
    select case (state)
      case ('L','l')
        PRV = dmin1(real(vol(1)),real(vol(2)),real(vol(3)))
      case ('V','v','G','g')
        PRV = dmax1(real(vol(1)),real(vol(2)),real(vol(3)))
      case default
        !print*, vol(1), vol(2), vol(3)
        !print*,'Incorrect state indicator'
        !stop
    !end select
  !}
  !else
  !{
    ! decide based on minimum Gibbs energy
      VV = dmax1(real(vol(1)),real(vol(2)),real(vol(3)))
      VL = dmin1(real(vol(1)),real(vol(2)),real(vol(3)))
  
      vol(1) = VV
      vol(2) = VL
      
      Gmin = 1.d20
      do j=1,2
      !{
        v = vol(j)
        if (v>b) then
          !G = -R_gas*T*dlog(v-b) + R_gas*T*b/(v-b) &
          !  + 0.25*a*sqr2/b*dlog((v+b-sqr2*b)/(v+b+sqr2*b)) &
          !  - a*v/(v**2+2.*b*v-b**2)
          Z = P*V/R_gas/T
          G = R_gas*T*( Z-1. - dlog(Z-Z*b/V) &
                       -0.25*sqr2*a/b/R_gas/T &
                        *dlog((v+b+sqr2*b)/(v+b-sqr2*b)) &
                      )
  
          if (G<Gmin) then
            Gmin = G
            Vmin = V
          endif
        endif
      !}
      enddo
  
      PRV = Vmin

    end select
  !}
  endif

!}
end subroutine PR_vol

subroutine PR_phase_(P,T,a,b,phase)
!{
  implicit none
  real(8) :: P,T,a,b, &
             a1,a2,a3,a4, &
             V,Vmin,G,Gmin,Z,VV,VL
  complex(8),dimension(3) :: vol
  integer :: i,j,phase
  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0

  ! Make up the coefficient vector for the cubic volume equation
  ! P*V^3 + (P*b-R*T)*V^2 + (a-3*P*b^2-2*R*T*b)*V + (P*b^3+R*T*b^2-a*b)=0

  a1 = P                       ! V^3 term
  a2 = P*b-R_gas*T             ! V^2 term
  a3 = a-3.d0*P*b*b-2.d0*R_gas*T*b ! V term
  a4 = P*b*b*b+R_gas*T*b*b-a*b ! constant term

  call cubic_solve(a1,a2,a3,a4,vol)

  if (imag(vol(2)) .ne. 0d0) then
  !{
    ! only one real root
    ! if (present(state)) print*,'Only one root; no option for state.'

    phase = -1
    ! V=Vmin
    !Gmin = -R_gas*T*dlog(v-b)+R_gas*T*b/(v-b)&
    !     + 0.25*a*sqr2/b*dlog((v+b-sqr2*b)/(v+b+sqr2*b)) &
    !     - a*v/(v**2+2*b*v-b**2)
  !}
  else !if (present(state)) then
  !{
    ! decide based on minimum Gibbs energy
      VV = dmax1(real(vol(1)),real(vol(2)),real(vol(3)))
      VL = dmin1(real(vol(1)),real(vol(2)),real(vol(3)))
  
      vol(1) = VV
      vol(2) = VL
      
      Gmin = 1.d20
      do j=1,2
      !{
        v = vol(j)
        if (v>b) then
          !G = -R_gas*T*dlog(v-b) + R_gas*T*b/(v-b) &
          !  + 0.25*a*sqr2/b*dlog((v+b-sqr2*b)/(v+b+sqr2*b)) &
          !  - a*v/(v**2+2.*b*v-b**2)
          Z = P*V/R_gas/T
          G = R_gas*T*( Z-1. - dlog(Z-Z*b/V) &
                       -0.25*sqr2*a/b/R_gas/T &
                        *dlog((v+b+sqr2*b)/(v+b-sqr2*b)) &
                      )
  
          if (G<Gmin) then
            Gmin = G
            Vmin = V
            if (j .eq. 1) then
              phase = 0
            else
              phase = 1
            endif
          endif
        endif
      !}
      enddo
  
  !}
  endif

!}
end subroutine PR_phase_

subroutine PR_phase( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    x, & ! vector of mass fractions
    type_k, & ! vector of binary interaction types
    coef_ab,& ! vector of a, b coefficients
    phase)   ! value of volume
!{
  implicit none
  integer :: i, j, n
  real(8) :: P,T,Pc(n),Tc(n),w(n),x(n)
  integer :: type_k(n), phase
  real(8) :: coef_ab(n)

  real(8) :: a(n), b(n), delta(n,n) ! delta is BIP
  real(8) :: am, bm, Tr, kappa, alpha

  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0

  call calculate_a_b(T,n,Pc,Tc,w,coef_ab,a,b)
  call calculate_kij(0,T,n,Pc,Tc,w,type_k,delta)

  ! PR coefficients of the mixture
  bm = 0d0
  am = 0d0

  do i=1,n
    bm      = bm + x(i)*b(i)
    am      = am + x(i)**2*a(i)

    do j=i+1,n
      am    = am + 2.*x(i)*x(j)*(1.-delta(i,j))*dsqrt(a(i)*a(j))
    enddo
  enddo  

  call PR_phase_(P,T,am,bm,phase)
!}
end subroutine PR_phase


!***************************************************************************
subroutine cubic_solve(a,b,c,d,roots)
!{
  implicit none
  complex(8),dimension(3) :: roots
  real(8) :: a,b,c,d
  integer :: nroot

  ! ----------------------------------------------------------------------
  ! Solve a cubic equation where a, b, c, and d are real.
  !   a*x**3 + b*x**2 + c*x + d = 0
  !
  ! Variables used:
  !   a, b, c, d  ... coefficients (input)
  !   y1, y2, y3  ... three transformed solutions
  !   y2r, y2i    ... real and imaginary parts of a pair of complex roots
  !   x(i)        ... three (generally) complex solutions (output)
  !   nroot       ... number of roots
  !
  ! Formula used are given in Tuma, "Engineering Mathematics Handbook", p7
  !   (McGraw Hill, 1978).
  !   Step 0: If a is 0. use the quadrati! formula to avoid dividing by 0.
  !   Step 1: Calculate p and q
  !           p = ( 3*c/a - (b/a)**2 ) / 3
  !           q = ( 2*(b/a)**3 - 9*b*c/a/a + 27*d/a ) / 27
  !   Step 2: Calculate discriminant D
  !           D = (p/3)**3 + (q/2)**2
  !   Step 3: Depending on the sign of D, we follow different strategy.
  !           If D<0, three distinct real roots.
  !           If D=0, three real roots of which at least two are equal.
  !           If D>0, one real and two complex roots.
  !   Step 3a: For D>0 and D=0,
  !           Calculate u and v
  !           u = cubic_root(-q/2 + sqrt(D))
  !           v = cubic_root(-q/2 - sqrt(D))
  !           Find the three transformed roots
  !           y1 = u + v
  !           y2 = -(u+v)/2 + i (u-v)*sqrt(3)/2
  !           y3 = -(u+v)/2 - i (u-v)*sqrt(3)/2
  !   Step 3b Alternately, for D<0, a trigonometri! formulation is more convenient
  !           y1 =  2 * sqrt(|p|/3) * cos(phi/3)
  !           y2 = -2 * sqrt(|p|/3) * cos((phi+pi)/3)
  !           y3 = -2 * sqrt(|p|/3) * cos((phi-pi)/3)
  !           where phi = acos(-q/2/sqrt(|p|**3/27))
  !                 pi  = 3.141592654...
  !   Step 4  Finally, find the three roots
  !           x = y - b/a/3
  !
  ! ----------------------------------------------------------------------
  ! Instructor: Nam Sun Wang
  ! ----------------------------------------------------------------------

  ! Declare variables
  complex(8) :: x(3)
  real(8) :: pi=3.1415926535897932384626433832795d0
  real(8) :: DD, p, q, phi, temp1, temp2, y1,y2,y3, u, v, y2r, y2i

  ! Step 0: If a is 0 use the quadratic formula. -------------------------
  if(a .eq. 0.) then
  !{

    if(b .eq. 0.)then 
    !{
      if(c .eq. 0.)then
        ! We have a non-equation; therefore, we have no valid solution
        nroot = 0
      else
        ! We have a linear equation with 1 root.
        nroot = 1
        x(1) = dcmplx(-d/c, 0.)
      endif
    !}
    else
    !{
      ! We have a true quadratic equation.  Apply the quadratic formula to find two roots.
      nroot = 2
      DD = c*c-4.d0*b*d
      if(DD .ge. 0.)then
        x(1) = dcmplx((-c+dsqrt(DD))/2.d0/b, 0.)
        x(2) = dcmplx((-c-dsqrt(DD))/2.d0/b, 0.)
      else
        x(1) = dcmplx(-c/2.d0/b, +dsqrt(-DD)/2.d0/b)
        x(2) = dcmplx(-c/2.d0/b, -dsqrt(-DD)/2.d0/b)
      endif
    !}
    endif

  !}
  else
  !{

    ! Cubic equation with 3 roots
    nroot = 3

    ! Step 1: Calculate p and q --------------------------------------------
    p  = c/a - b*b/(a*a*3.d0)
    q  = (2.d0*b*b*b/(a*a*a) - 9.d0*b*c/(a*a) + 27.d0*d/a) / 27.d0

    ! Step 2: Calculate DD (discriminant) ----------------------------------
    DD = p*p*p/27.d0 + q*q*.25d0

    ! Step 3: Branch to different algorithms based on DD -------------------
    if(DD .lt. 0.)then
      ! Step 3b:
      ! 3 real unequal roots -- use the trigonometric formulation
      phi = dacos(-q*0.5d0/dsqrt(dabs(p*p*p)/27.d0))
      temp1=2.d0*dsqrt(dabs(p)/3.d0)
      y1 =  temp1*dcos(phi/3.d0)
      y2 = -temp1*dcos((phi+pi)/3.d0)
      y3 = -temp1*dcos((phi-pi)/3.d0)
    else
      ! Step 3a:
      ! 1 real root & 2 conjugate complex roots OR 3 real roots (some are equal)
      temp1 = -q/2.d0 + dsqrt(DD)
      temp2 = -q/2.d0 - dsqrt(DD)
      u = dabs(temp1)**(1./3.d0)
      v = dabs(temp2)**(1./3.d0)
      if(temp1 .lt. 0.) u=-u
      if(temp2 .lt. 0.) v=-v
      y1  = u + v
      y2r = -(u+v)/2.d0
      y2i =  (u-v)*dsqrt(3.d0)/2.d0
    endif

    ! Step 4: Final transformation -----------------------------------------
    temp1 = b/a/3.d0
    y1 = y1-temp1
    y2 = y2-temp1
    y3 = y3-temp1
    y2r=y2r-temp1

    ! Assign answers -------------------------------------------------------
    if(DD .lt. 0.)then
      x(1) = dcmplx( y1,  0.)
      x(2) = dcmplx( y2,  0.)
      x(3) = dcmplx( y3,  0.)
    elseif(DD .eq. 0.)then
      x(1) = dcmplx( y1,  0.)
      x(2) = dcmplx(y2r,  0.)
      x(3) = dcmplx(y2r,  0.)
    else
      x(1) = dcmplx( y1,  0.)
      x(2) = dcmplx(y2r, y2i)
      x(3) = dcmplx(y2r,-y2i)
    endif

  !}
  endif

  roots = x

!}
end subroutine cubic_solve

!************************************************************************
! INTERNAL CALLS
subroutine calculate_a_b( &
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    coef_ab, & ! vector of a,b coefficients
    a, & ! vector of a (OUTPUT)
    b  & ! vector of b (OUTPUT)
    )
!{
  implicit none
  integer :: n
  real(8) :: T,Pc(n),Tc(n),w(n),a(n),b(n),dadT(n)
  real(8) :: coef_ab(n)

  call calculate_a_b_and_dadT(T,n,Pc,Tc,w,coef_ab,a,b,dadT)
!}
end subroutine calculate_a_b
subroutine calculate_a_b_and_dadT( &
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    coef_ab, & ! vector of a,b coefficients
    a, & ! vector of a (OUTPUT)
    b, & ! vector of b (OUTPUT)
    dadT & ! vector of dadT (OUTPUT)
    )
!{
  implicit none
  integer :: n, i
  real(8) :: T,Pc(n),Tc(n),w(n),a(n),b(n),dadT(n)
  real(8) :: coef_ab(n)

  real(8) :: Tr, kappa, alpha, alpha2, cc
  real(8) :: R_gas = 8.3144621d0

  do i=1,n
    b(i)  =  0.0777960739* R_gas*Tc(i) / Pc(i)

    Tr = T/Tc(i)

    if (coef_ab(i) .lt. 0) then
      if (w(i) .le. 0.491d0) then
        kappa = (0.374640d0+1.54226d0*w(i)-0.269920d0*w(i)**2)
      else
        kappa = (0.379642d0+1.48503d0*w(i)-0.164423d0*w(i)**2+0.0166667d0*w(i)**3)
      endif

      alpha = 1d0+kappa*(1d0-dsqrt(Tr))
      alpha2= alpha**2

      cc = 0.457235529d0*(R_gas*Tc(i))**2/Pc(i)

      a(i)  =  cc*alpha2

      dadT(i) = -cc*kappa/(dsqrt(Tc(i)*T))*alpha

!     print *, "NON-ASPHALTENE BEGIN"
!     print *, "a=",a(i),"b=",b(i),"c=",kappa,"R2T2/P=",(R_gas*Tc(i))**2/Pc(i)
!     print *, "alpha=",alpha
!     print *, "NON-ASPHALTENE END"
    else
      ! Method of Sabbagh et al. 2006
      ! Construct asphaltene species using a monomer data
      ! I will NOT use this part of code currently
      ! ----------------------------------------------------------------
      ! NOT USED NOT USED NOT USED NOT USED NOT USED NOT USED NOT USED 
      ! BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN
      kappa = (0.3796+1.485*w(i)-0.1644*w(i)**2+0.01667*w(i)**3)
      alpha = dabs( (1d0+kappa*(1d0-Tr)**3)**.3333333333333333d0 )


      a(i)  = (R_gas*Tc(i))**2/Pc(i)*alpha

      a(i)  = a(i)*coef_ab(i)**2
      b(i)  = b(i)*coef_ab(i)

      !print *, "a=",a(i),"b=",b(i),"c=",kappa,"R2T2/P=",(R_gas*Tc(i))**2/Pc(i)
      !print *, "alpha=",alpha
      ! NOT USED NOT USED NOT USED NOT USED NOT USED NOT USED NOT USED 
      ! END END END END END END END
      ! -----------------------------------------------------------------
    endif
  end do
!}
end subroutine calculate_a_b_and_dadT
!************************************************************************
!************************************************************************
! INTERNAL CALLS
subroutine calculate_a_b2( &
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    a, & ! vector of a (OUTPUT)
    b  & ! vector of b (OUTPUT)
    )
!{
  implicit none
  integer :: n,i
  real(8) :: T,Pc(n),Tc(n),w(n),a(n),b(n)

  real(8) :: Tr, kappa, alpha, alpha2, cc
  real(8) :: R_gas = 8.3144621d0

  do i=1,n
    b(i)  =  0.0777960739* R_gas*Tc(i) / Pc(i)

    Tr = T/Tc(i)
    
    if (w(i) .le. 0.491d0) then
       kappa = (0.374640d0+1.54226d0*w(i)-0.269920d0*w(i)**2)
    else
       kappa = (0.379642d0+1.48503d0*w(i)-0.164423d0*w(i)**2+0.0166667d0*w(i)**3)
    endif

    alpha = 1d0+kappa*(1d0-dsqrt(Tr))
    alpha2= alpha**2

    cc = 0.457235529d0*(R_gas*Tc(i))**2/Pc(i)

    a(i)  =  cc*alpha2
  end do
!}
end subroutine calculate_a_b2
subroutine calculate_a_b_and_dadT2( &
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    a, & ! vector of a (OUTPUT)
    b, & ! vector of b (OUTPUT)
    dadT & ! vector of dadT (OUTPUT)
    )
!{
  implicit none
  integer :: n, i
  real(8) :: T,Pc(n),Tc(n),w(n),a(n),b(n),dadT(n)

  real(8) :: Tr, kappa, alpha, alpha2, cc
  real(8) :: R_gas = 8.3144621d0

  do i=1,n
    b(i)  =  0.0777960739* R_gas*Tc(i) / Pc(i)

    Tr = T/Tc(i)

    if (w(i) .le. 0.491d0) then
       kappa = (0.374640d0+1.54226d0*w(i)-0.269920d0*w(i)**2)
    else
       kappa = (0.379642d0+1.48503d0*w(i)-0.164423d0*w(i)**2+0.0166667d0*w(i)**3)
    endif

    alpha = 1d0+kappa*(1d0-dsqrt(Tr))
    alpha2= alpha**2

    cc = 0.457235529d0*(R_gas*Tc(i))**2/Pc(i)

    a(i)  =  cc*alpha2

    dadT(i) = -cc*kappa/(dsqrt(Tc(i)*T))*alpha
  end do
!}
end subroutine calculate_a_b_and_dadT2
!************************************************************************
subroutine calculate_kij( &
    bSet, & ! >=1 then set every thing, <=0 then just get values restored
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    type_k, & ! vector of binary interaction types
    kij & ! matrix of binary interation parameters
    )
!   coef_ab,& ! vector of a,b coefficients
!{
  implicit none
  integer :: n, bSet
  integer :: type_k(n)
  real(8) :: Pc(n), Tc(n), w(n), coef_ab(n)
  real(8) :: T, kij(n,n), Fij(n,n), dFdT(n,n)

  coef_ab = -1

  call calculate_kij_and_dkijdT(bSet,T,n,Pc,Tc,w,coef_ab,type_k,kij,Fij,dFdT)
!}
end subroutine calculate_kij
subroutine calculate_kij_and_dkijdT( &
    bSet, & ! >=1 then set every thing, <=0 then just get values restored
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    coef_ab,& ! vector of a,b coefficients
    type_k, & ! vector of binary interaction types
    kij, & ! matrix of binary interation parameters
    Fij, & ! maxtrix of F(i,j)
    dFdT & ! matrix of d(Fij)/dT
    )
!{
  implicit none
  integer :: n, i, j, bSet
  integer :: type_k(n)
  real(8) :: Pc(n), Tc(n), w(n), coef_ab(n)
  real(8) :: T, T0, kij(n,n), ai(n), bi(n), Fij(n,n), dFdT(n,n)

  ! PPR78 static data and code
  ! PPR78 PPR78 PPR78 PPR78 PPR78 PPR78 PPR78 PPR78 PPR78
  integer :: ni, nj, k, l
  real(8) :: sumG, sqrta_bi, sqrta_bj, r1_298, tmp1, tmp2
  logical :: first_time=.true.
  integer :: n_species, n_pure, NS
  real(8) :: A(22,22), B(22,22)!TOTAL 22 GROUPS in current PPR78 method

  !GROUP INFO n_species x 23, T_up(*,*)
  real(8), allocatable, dimension(:,:) :: G, T_up, T_up2
  real(8), allocatable, dimension(:,:) :: sBIP,rBIP
  real(8), allocatable, dimension(:,:) :: kij_stored  ! stored value for retrieval
  !last column Dipole Moment (NOT USED HERE)
  save first_time,n_species,n_pure,A,B,G,r1_298,T_up,kij_stored,sBIP,rBIP, T_up2

  ! FIRST TIME PROCESSING FIRST TIME PROCESSING FIRST TIME PROCESSING
  ! BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN
  !{
  if ( first_time ) then
    first_time = .false.
    r1_298 = 1.d0 / 298.15d0

    !if (type_k(1) .lt. 0) then 
    !{
      ! READ PPR78 DATA to be STATIC in memory
      allocate( kij_stored(n,n) )

      ! 0) read A, B data 
      open(unit = 111, file = 'PPR78AB.dat+', &
        access = 'sequential', status = 'unknown')
      do i=1,22
        read(111,*) (A(i,j),j=1,22)
        !print*, (A(i,j),j=1,21)
      enddo
      do i=1,22
        read(111,*) (B(i,j),j=1,22)
        !print*, (B(i,j),j=1,21)
      enddo
      close(111)

      ! 1) read a header file for species' database
      !      n_species (it is database size, used species No. is n <= n_species)
      open(unit = 111, file = 'species.dat+', &
        access = 'sequential', status = 'unknown')
      read(111,*) n_species, n_pure
      close(111)
      !print*,'first time loading: n_species =',n_species
      NS = n_species*3 + n_pure

      ! 2) read group data for all species
      allocate( G(1:n_species*3+n_pure,1:23) );!times 3 means P.N.A. + water + n-decane + ...
      open(unit = 111, file = 'groups.dat+', &
        access = 'sequential', status = 'unknown')
      do i=1,n_species*3+n_pure
        read(111,*) (G(i,j),j=1,23)
        !print*, (G(i,j),j=1,22)
        
        ! calculate the fractions
        sumG = 0.0d0
        do j=1,22
          sumG = sumG + G(i,j)
        enddo
        sumG = 1.0d0 / sumG
        do j=1,22
          G(i,j) = G(i,j)*sumG 
        enddo
        !print*, (G(i,j),j=1,22)

      enddo
      close(111)

      ! 3) read T_up
      allocate( T_up(1:NS,1:NS) )
      open(unit = 111, file = 'BIP_T_up.dat+', &
        access = 'sequential', status = 'unknown')
      do i=1,NS
        read(111,*) (T_up(i,j),j=1,NS)
        !print*, (T_up(i,j),j=1,NS)
      enddo
      close(111)

      allocate( T_up2(1:NS,1:NS) )
      open(unit = 111, file = 'BIP_T_up2.dat+', &
        access = 'sequential', status = 'unknown')
      do i=1,NS
        read(111,*) (T_up2(i,j),j=1,NS)
        !print*, (T_up2(i,j),j=1,NS)
      enddo
      close(111)

      ! 4) read BIP sign
      allocate( sBIP(1:NS,1:NS) )
      open(unit = 111, file = 'BIP_sign.dat+', &
        access = 'sequential', status = 'unknown')
      !print*, 'sBIP:'
      do i=1,NS
        read(111,*) (sBIP(i,j),j=1,NS)
        !print*, (sBIP(i,j),j=1,NS)
      enddo
      close(111)

      ! 5) read BIP value
      allocate( rBIP(1:NS,1:NS) )
      open(unit = 111, file = 'BIP_value.dat+', &
        access = 'sequential', status = 'unknown')
      !print*, 'rBIP:'
      do i=1,NS
        read(111,*) (rBIP(i,j),j=1,NS)
        !print*, (rBIP(i,j),j=1,NS)
      enddo
      close(111)
    !}
    !endif

  endif
  !}
  ! END END END END END END END END END END END END END END END END
  ! FIRST TIME PROCESSING FIRST TIME PROCESSING FIRST TIME PROCESSING

  !if ( .not. first_time ) then
  !  print*, G(100,10)
  !endif

  if (bSet<=0) then
    do i=1,n
      do j=1,n
        kij(i,j) = kij_stored(i,j)
      end do
    end do

    return
  end if

  do i=1,n
    kij(i,i)=0.0
    kij_stored(i,i)=kij(i,i)

    do j=i+1,n
    !{
      ! ONE         1: water            2: toluene
      if ( (type_k(i).eq.1 .AND. type_k(j).eq.2)   &
       .OR.(type_k(j).eq.1 .AND. type_k(i).eq.2) ) then
        kij(i,j) = 0.1935 - 2.1344/(4.*(T-581.15)**2+12.3108)
        dFdT(i,j) = 2.1344*8*(T-581.15)/(4.*(T-581.15)**2+12.3108)**2

        if(T>581.15) then
            kij(i,j)=0.0183
            dFdT(i,j)=0
        endif

        !kij(i,j)=0.02012
        kij(j,i) = kij(i,j)
        dFdT(j,i) = dFdT(i,j)

      ! TWO         1: water            3: n-decane
      elseif ( (type_k(i).eq.1 .AND. type_k(j).eq.3)   &
           .OR.(type_k(j).eq.1 .AND. type_k(i).eq.3) ) then
        kij(i,j) = 0.3237 - 34.9097/(4.*(T-631.9)**2+82.6681)
        dFdT(i,j) = 34.9097*8*(T-631.9)/(4.*(T-631.9)**2+82.6681)**2

        if(T>631.9) then
            kij(i,j)=-0.099
            dFdT(i,j)=0
        endif

        kij(j,i) = kij(i,j)
        dFdT(j,i) = dFdT(i,j)

      ! THREE       2: toluene          3: n-decane
      elseif ( (type_k(i).eq.2 .AND. type_k(j).eq.3)   &
           .OR.(type_k(j).eq.2 .AND. type_k(i).eq.3) ) then
        kij(i,j) = -0.01
        dFdT(i,j)=0
        kij(j,i) = kij(i,j)
        dFdT(j,i) = dFdT(i,j)

      ! FOUR        4: saturates        5: asphaltenes
      elseif ( (type_k(i).eq.4 .AND. type_k(j).eq.5)   &
           .OR.(type_k(j).eq.4 .AND. type_k(i).eq.5) ) then
        kij(i,j) = -0.00256
        dFdT(i,j)=0
        kij(j,i) = kij(i,j)
        dFdT(j,i) = dFdT(i,j)

      ! FIVE        1: water            6: PC 1
      elseif ( (type_k(i).eq.1 .AND. type_k(j).eq.6)   &
           .OR.(type_k(j).eq.1 .AND. type_k(i).eq.6) ) then
        kij(i,j) = 0.08
        dFdT(i,j)=0
        kij(j,i) = kij(i,j)
        dFdT(j,i) = dFdT(i,j)

      ! SIX         1: water            7: PC 2
      elseif ( (type_k(i).eq.1 .AND. type_k(j).eq.7)   &
           .OR.(type_k(j).eq.1 .AND. type_k(i).eq.7) ) then
        kij(i,j) = 0.12
        dFdT(i,j)=0
        kij(j,i) = kij(i,j)
        dFdT(j,i) = dFdT(i,j)

      ! SEVEN       1: water            8: PC 3
      elseif ( (type_k(i).eq.1 .AND. type_k(j).eq.8)   &
           .OR.(type_k(j).eq.1 .AND. type_k(i).eq.8) ) then
        kij(i,j) = 0.16
        dFdT(i,j)=0
        kij(j,i) = kij(i,j)
        dFdT(j,i) = dFdT(i,j)

      ! EIGHT       1: water            9: PC 4
      elseif ( (type_k(i).eq.1 .AND. type_k(j).eq.9)   &
           .OR.(type_k(j).eq.1 .AND. type_k(i).eq.9) ) then
        kij(i,j) = 0.18
        dFdT(i,j)=0
        kij(j,i) = kij(i,j)
        dFdT(j,i) = dFdT(i,j)

      ! NINE        1: water           10: PC 5
      elseif ( (type_k(i).eq.1 .AND. type_k(j).eq.10)   &
           .OR.(type_k(j).eq.1 .AND. type_k(i).eq.10) ) then
        kij(i,j) = 0.20
        dFdT(i,j)=0
        kij(j,i) = kij(i,j)
        dFdT(j,i) = dFdT(i,j)

      ! TEN         1: water           11: PC 6
      elseif ( (type_k(i).eq.1 .AND. type_k(j).eq.11)   &
           .OR.(type_k(j).eq.1 .AND. type_k(i).eq.11) ) then
        kij(i,j) = 0.25
        dFdT(i,j)=0
        kij(j,i) = kij(i,j)
        dFdT(j,i) = dFdT(i,j)

      ! ELEVEN      1: water           12: PC 7
      elseif ( (type_k(i).eq.1 .AND. type_k(j).eq.12)   &
           .OR.(type_k(j).eq.1 .AND. type_k(i).eq.12) ) then
        kij(i,j) = 0.30
        dFdT(i,j)=0
        kij(j,i) = kij(i,j)
        dFdT(j,i) = dFdT(i,j)

      ! PPR78 PPR78 PPR78 PPR78 PPR78 PPR78 PPR78 PPR78 PPR78 
      !{_________________________________________________________
      !elseif ( (type_k(i).lt.0 .AND. type_k(j).lt.0) ) then
      elseif ( (type_k(i).lt.0 ) ) then
        ni = - type_k(i)
        nj = - type_k(j)

        T0 = T
        if (T_up(ni,nj)>0 .AND. T>T_up(ni,nj)) T0=T_up(ni,nj)
        call calculate_a_b(T0,n,Pc,Tc,w,coef_ab,ai,bi)

        sqrta_bi = dsqrt(ai(i))/bi(i)
        sqrta_bj = dsqrt(ai(j))/bi(j)

        sumG = 0.0d0
        dFdT(i,j) = 0.0d0
        do k=1,22
          do l=1,22
            if ( dabs(A(k,l)) > 1.0d-5 ) then
              sumG = sumG + (G(ni,k)-G(nj,k))*(G(ni,l)-G(nj,l)) &
                          *  A(k,l)*(298.15d0/T0)**(B(k,l)/A(k,l)-1.d0)
              !sumG = sumG + (G(ni,k)-G(nj,k))*(G(ni,l)-G(nj,l))*A(k,l)
              dFdT(i,j) = dFdT(i,j) &
                          + (G(ni,k)-G(nj,k))*(G(ni,l)-G(nj,l)) &
              * (A(k,l)-B(k,l))*r1_298*(298.15d0/T0)**(B(k,l)/A(k,l))
            endif
          enddo 
        enddo  

        Fij(i,j) = sumG*1d6
        Fij(j,i) = Fij(i,j)
        dFdT(i,j) = dFdT(i,j)*1d6
        dFdT(j,i) = dFdT(i,j)

        kij(i,j) = (-0.5d0*Fij(i,j) &
                  -(sqrta_bi-sqrta_bj)*(sqrta_bi-sqrta_bj)) &
                 *   0.5d0/(sqrta_bi*sqrta_bj)
        kij(j,i) = kij(i,j)

        if (T_up(ni,nj)>0 .AND. T>T_up(ni,nj)) then
          tmp1 = kij(i,j)

          T0 = T_up(ni,nj)-1
          call calculate_a_b(T0,n,Pc,Tc,w,coef_ab,ai,bi)
 
          sqrta_bi = dsqrt(ai(i))/bi(i)
          sqrta_bj = dsqrt(ai(j))/bi(j)
 
          sumG = 0.0d0
          dFdT(i,j) = 0.0d0
          do k=1,22
            do l=1,22
              if ( dabs(A(k,l)) > 1.0d-5 ) then
                sumG = sumG + (G(ni,k)-G(nj,k))*(G(ni,l)-G(nj,l)) &
                            *  A(k,l)*(298.15d0/T0)**(B(k,l)/A(k,l)-1.d0)
                !sumG = sumG + (G(ni,k)-G(nj,k))*(G(ni,l)-G(nj,l))*A(k,l)
                dFdT(i,j) = dFdT(i,j) &
                            + (G(ni,k)-G(nj,k))*(G(ni,l)-G(nj,l)) &
                * (A(k,l)-B(k,l))*r1_298*(298.15d0/T0)**(B(k,l)/A(k,l))
              endif
            enddo 
          enddo  
 
          Fij(i,j) = sumG*1d6
          Fij(j,i) = Fij(i,j)
          dFdT(i,j) = dFdT(i,j)*1d6
          dFdT(j,i) = dFdT(i,j)
 
          tmp2 = (-0.5d0*Fij(i,j) &
                    -(sqrta_bi-sqrta_bj)*(sqrta_bi-sqrta_bj)) &
                   *   0.5d0/(sqrta_bi*sqrta_bj)

          T0 = T_up(ni,nj)

          if (T_up2(ni,nj)>0 .AND. T>T_up2(ni,nj)) then
            kij(i,j) = (tmp1-tmp2)*(T_up2(ni,nj)-T0) + tmp1
          else
            kij(i,j) = (tmp1-tmp2)*(T-T0) + tmp1
          endif

          kij(j,i) = kij(i,j)
        end if

        !print *, ni,nj,sBIP(ni,nj),rBIP(ni,nj)
        if (sBIP(ni,nj) .gt. 0) then
          kij(i,j) = rBIP(ni,nj)
          kij(j,i) = kij(i,j)
          print *, 'kij', kij(i,j)
        end if

      !}_________________________________________________________
      ! DEFAULT: zero, which means no binary interactions
      else
        kij(i,j) = 0.0d0
        kij(j,i) = kij(i,j)
        Fij(i,j) = 0d0
        Fij(j,i) = Fij(i,j)
        dFdT(i,j) = 0d0
        dFdT(j,i) = dFdT(i,j)
      endif

      kij_stored(i,j)=kij(i,j)
      kij_stored(j,i)=kij(i,j)
    !}
    enddo
  enddo
!}
end subroutine calculate_kij_and_dkijdT
!************************************************************************
!************************************************************************
subroutine calculate_kij_from_table( &
  bSet, & ! >=1 then set every thing, <=0 then just get values restored
  T, &    ! temperature (Unit: K)
  n, &    ! number of species  
  kij &   ! matrix of binary interation parameters
  )
!{
  implicit none
  integer :: n, bSet
  real(8) :: T, kij(n,n)

  integer :: nT, i, j, k, idT, idx
  real(8) :: Ta, Tb, dT, T1
  logical :: first_time=.true.

  real(8), allocatable, dimension(:,:) :: kij_T
  real(8), allocatable, dimension(:,:) :: kij_stored
  
  save first_time,kij_T,Ta,Tb,nT,kij_stored

  ! FIRST TIME PROCESSING FIRST TIME PROCESSING FIRST TIME PROCESSING
  ! BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN
  !{
  if ( first_time ) then
    first_time = .false.

    open(unit = 111, file = 'constant/BIP_nT.dat+', &
         access = 'sequential', status = 'unknown')
    read(111,*) Ta, Tb, nT
    close(111)

    ! read kij vs T data
    allocate( kij_T(1:nT,1:n*n) )
    allocate( kij_stored(n,n) )
    open(unit = 111, file = 'constant/kij_T.dat+', &
         access = 'sequential', status = 'unknown')
    do k=1,nT
       read(111,*) (kij_T(k,j),j=1,n*n)          
    enddo    
    close(111)
    !}
    !endif
  endif
  !}
  ! END END END END END END END END END END END END END END END END
  ! FIRST TIME PROCESSING FIRST TIME PROCESSING FIRST TIME PROCESSING

  if (bSet<=0) then
    do i=1,n
      do j=1,n
        kij(i,j) = kij_stored(i,j)
      end do
    end do

    return
  end if

  if (T .GT. Ta .AND. T .LT. Tb) then
     dT = (Tb - Ta)/(nT - 1)
     idT = floor((T - Ta)/dT)
     T1 = Ta + idT*dT
     idT = idT + 1  

     do i=1,n
        do j=1,n
           idx = i*n + j
           kij(i,j) = kij_T(idT,idx) + (kij_T(idT+1,idx) - kij_T(idT,idx))*(T - T1)/dT 
        end do
     end do
  else
     if (T .LE. Ta) then
        do i=1,n
           do j=1,n
              idx = i*n + j
              kij(i,j) = kij_T(1,idx) 
           end do
        end do
     else
        do i=1,n
           do j=1,n
              idx = i*n + j
              kij(i,j) = kij_T(nT,idx) 
           end do
        end do
     endif
  endif
  
!}
end subroutine calculate_kij_from_table
!************************************************************************
!***************************************************************************
subroutine fugacities_n_its_derivatives( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    x, & ! vector of mass fractions
    type_k,     & ! vector of binary interaction types
    coef_ab,    & ! vector of a, b coefficients
    lnphi,      & ! vector of fugacity coefficients (output results)
    dlnphi_dxj, & ! matrix of dlnphi_i/dx_j         (output results)
    V           & ! molar volume, scalar            (output results)
    )
!{
  implicit none
  integer :: n
  real(8) :: P,T,Pc(n),Tc(n),w(n),x(n),lnphi(n),dlnphi_dxj(n,n),V,&
             am, dam_dxi(n), d2am_dxidxj(n,n), dV_dxi(n), coef_ab(n)
  integer :: type_k(n)

  call fugacities_n_its_derivatives_general( &
         P,T,n,Pc,Tc,w,x,type_k,coef_ab,lnphi,dlnphi_dxj,V,&
         am, dV_dxi, dam_dxi, d2am_dxidxj);

!}
end subroutine fugacities_n_its_derivatives

subroutine fugacities_n_its_derivatives_general( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    x, & ! vector of mass fractions
    type_k,     & ! vector of binary interaction types
    coef_ab,    & ! vector of a, b coefficients
    lnphi,      & ! vector of fugacity coefficients (output results)
    dlnphi_dxj, & ! matrix of dlnphi_i/dx_j         (output results)
    V,          & ! molar volume, scalar            (output results)
    am, dV_dxi, dam_dxi, d2am_dxidxj &
    )
!{
  implicit none
  integer :: n
  real(8) :: P,T,Pc(n),Tc(n),w(n),x(n)
  integer :: type_k(n)
  real(8) :: coef_ab(n)

  real(8) :: a(n), b(n), delta(n,n) ! delta is binary interaction parameter (k_ij)
  real(8) :: am, bm, Tr, kappa, alpha, &
             dam_dxi(n), d2am_dxidxj(n,n), dV_dxi(n)
  real(8) :: lnphi(n), dlnphi_dxj(n,n), V, Z, Vs, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8

  integer :: i,j
  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0
  character :: state = 'm'

  ! PR coefficients of pure components
  call calculate_a_b(T,n,Pc,Tc,w,coef_ab,a,b)
  !call calculate_kij(0,T,n,Pc,Tc,w,type_k,delta)
  call calculate_kij_from_table(0,T,n,delta)

  ! PR coefficients of the mixture
  am = 0d0
  bm = 0d0
  dam_dxi = 0d0
  d2am_dxidxj = 0d0
  do i=1,n
    bm      = bm + x(i)*b(i)
    do j=1,n
      tmp1             = 2d0*(1d0-delta(i,j))*dsqrt(a(i)*a(j))
      am               = am + .5*x(i)*x(j)*tmp1
      dam_dxi(i)       = dam_dxi(i) + x(j)*tmp1
      d2am_dxidxj(i,j) = tmp1                
    enddo
  enddo  

  call PR_vol(P,T,am,bm,V,state)

  Vs= V*V + 2d0*V*bm - bm*bm
  tmp1 = 2d0*am*(V-bm)
  tmp2 = ((Vs/(V-bm))**2)*R_gas*T
  tmp3 = 1d0 / (2d0*am*(V+bm) - tmp2)
  do i=1,n
    dV_dxi(i) = ( Vs*dam_dxi(i) - (tmp1+tmp2)*b(i) ) * tmp3
  end do

  Z = P*V/R_gas/T
  tmp1 = 1d0/(V-bm)
  tmp2 = 1d0/(dsqrt(8d0)*bm*R_gas*T)
  tmp3 = dlog((V+(1.+sqr2)*bm)/(V+(1.-sqr2)*bm))
  tmp5 = 1.0/(bm*bm)
  tmp6 = P/(bm*R_gas*T)
  tmp7 = tmp2*dsqrt(8d0)
  tmp8 = 1.0/(am*am)

  do i=1,n
  !{
    tmp4 = dam_dxi(i)/am-b(i)/bm

    lnphi(i) = - am*tmp2*tmp3*tmp4 - dlog(Z*(1.-bm/V)) + b(i)/bm*(Z-1.)

    do j=1,n
      dlnphi_dxj(i,j) = &
        - tmp1*(dV_dxi(j)-b(j)) - b(i)*b(j)*tmp5*(Z-1.) + b(i)*tmp6*dV_dxi(j) &
        - (dam_dxi(j) - am*b(j)/bm)*tmp2*tmp3*tmp4  &
        - am*tmp7*tmp4/Vs*( V*b(j) - bm*dV_dxi(j) ) &
        - am*tmp2*tmp3*( - dam_dxi(i)*dam_dxi(j)*tmp8 + d2am_dxidxj(i,j)/am + b(i)*b(j)*tmp5 )
    enddo
  !}
  enddo
!}
end subroutine fugacities_n_its_derivatives_general

subroutine fugacities_n_its_derivatives2( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    x, & ! vector of mole fractions
    kij, & ! matrix of BIPs
    lnphi, & ! vector of fugacity coefficients (output)
    dlnphi_dxj & ! matrix of dlnphi_i/dx_j (output)
    )
!{
  implicit none
  integer :: n
  real(8) :: P,T,Pc(n),Tc(n),w(n),x(n),kij(n,n), lnphi(n),dlnphi_dxj(n,n)

  integer :: i,j
  real(8) :: a(n), b(n)
  real(8) :: am, bm, dam_dxi(n), d2am_dxidxj(n,n), dV_dxi(n)
  real(8) :: V, Z, Vs, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8
  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0
  character :: state = 'm'
  
  ! PR coefficients of pure components
  call calculate_a_b2(T,n,Pc,Tc,w,a,b)  

  ! PR coefficients of the mixture
  am = 0d0
  bm = 0d0
  dam_dxi = 0d0
  d2am_dxidxj = 0d0
  do i=1,n
    bm      = bm + x(i)*b(i)
    do j=1,n
      tmp1             = 2d0*(1d0-kij(i,j))*dsqrt(a(i)*a(j))
      am               = am + .5*x(i)*x(j)*tmp1
      dam_dxi(i)       = dam_dxi(i) + x(j)*tmp1
      d2am_dxidxj(i,j) = tmp1        
    enddo
  enddo  

  call PR_vol(P,T,am,bm,V,state)

  Vs= V*V + 2d0*V*bm - bm*bm
  tmp1 = 2d0*am*(V-bm)
  tmp2 = ((Vs/(V-bm))**2)*R_gas*T
  tmp3 = 1d0 / (2d0*am*(V+bm) - tmp2)
  do i=1,n
    dV_dxi(i) = ( Vs*dam_dxi(i) - (tmp1+tmp2)*b(i) ) * tmp3
  end do

  Z = P*V/R_gas/T
  tmp1 = 1d0/(V-bm)
  tmp2 = 1d0/(dsqrt(8d0)*bm*R_gas*T)
  tmp3 = dlog((V+(1.+sqr2)*bm)/(V+(1.-sqr2)*bm))
  tmp5 = 1.0/(bm*bm)
  tmp6 = P/(bm*R_gas*T)
  tmp7 = tmp2*dsqrt(8d0)
  tmp8 = 1.0/(am*am)

  do i=1,n
  !{
    tmp4 = dam_dxi(i)/am-b(i)/bm

    lnphi(i) = - am*tmp2*tmp3*tmp4 - dlog(Z*(1.-bm/V)) + b(i)/bm*(Z-1.)

    do j=1,n
      dlnphi_dxj(i,j) = &
        - tmp1*(dV_dxi(j)-b(j)) - b(i)*b(j)*tmp5*(Z-1.) + b(i)*tmp6*dV_dxi(j) &
        - (dam_dxi(j) - am*b(j)/bm)*tmp2*tmp3*tmp4  &
        - am*tmp7*tmp4/Vs*( V*b(j) - bm*dV_dxi(j) ) &
        - am*tmp2*tmp3*( - dam_dxi(i)*dam_dxi(j)*tmp8 + d2am_dxidxj(i,j)/am + b(i)*b(j)*tmp5 )
    enddo
  !}
  enddo
!}
end subroutine fugacities_n_its_derivatives2

subroutine fugacities_n_its_derivatives3( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    x, & ! vector of mole fractions
    kij, & ! matrix of BIPs
    lnphi, & ! vector of fugacity coefficients (output)
    dlnphi_dxj, & ! matrix of dlnphi_i/dx_j (output)
    V  & ! molar volume (output)
    )
!{
  implicit none
  integer :: n
  real(8) :: P,T,Pc(n),Tc(n),w(n),x(n),kij(n,n), lnphi(n),dlnphi_dxj(n,n),V

  integer :: i,j
  real(8) :: a(n), b(n)
  real(8) :: am, bm, dam_dxi(n), d2am_dxidxj(n,n), dV_dxi(n)
  real(8) :: Z, Vs, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8
  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0
  character :: state = 'm'
  
  ! PR coefficients of pure components
  call calculate_a_b2(T,n,Pc,Tc,w,a,b)  

  ! PR coefficients of the mixture
  am = 0d0
  bm = 0d0
  dam_dxi = 0d0
  d2am_dxidxj = 0d0
  do i=1,n
    bm      = bm + x(i)*b(i)
    do j=1,n
      tmp1             = 2d0*(1d0-kij(i,j))*dsqrt(a(i)*a(j))
      am               = am + .5*x(i)*x(j)*tmp1
      dam_dxi(i)       = dam_dxi(i) + x(j)*tmp1
      d2am_dxidxj(i,j) = tmp1        
    enddo
  enddo  

  call PR_vol(P,T,am,bm,V,state)

  Vs= V*V + 2d0*V*bm - bm*bm
  tmp1 = 2d0*am*(V-bm)
  tmp2 = ((Vs/(V-bm))**2)*R_gas*T
  tmp3 = 1d0 / (2d0*am*(V+bm) - tmp2)
  do i=1,n
    dV_dxi(i) = ( Vs*dam_dxi(i) - (tmp1+tmp2)*b(i) ) * tmp3
  end do

  Z = P*V/R_gas/T
  tmp1 = 1d0/(V-bm)
  tmp2 = 1d0/(dsqrt(8d0)*bm*R_gas*T)
  tmp3 = dlog((V+(1.+sqr2)*bm)/(V+(1.-sqr2)*bm))
  tmp5 = 1.0/(bm*bm)
  tmp6 = P/(bm*R_gas*T)
  tmp7 = tmp2*dsqrt(8d0)
  tmp8 = 1.0/(am*am)

  do i=1,n
  !{
    tmp4 = dam_dxi(i)/am-b(i)/bm

    lnphi(i) = - am*tmp2*tmp3*tmp4 - dlog(Z*(1.-bm/V)) + b(i)/bm*(Z-1.)

    do j=1,n
      dlnphi_dxj(i,j) = &
        - tmp1*(dV_dxi(j)-b(j)) - b(i)*b(j)*tmp5*(Z-1.) + b(i)*tmp6*dV_dxi(j) &
        - (dam_dxi(j) - am*b(j)/bm)*tmp2*tmp3*tmp4  &
        - am*tmp7*tmp4/Vs*( V*b(j) - bm*dV_dxi(j) ) &
        - am*tmp2*tmp3*( - dam_dxi(i)*dam_dxi(j)*tmp8 + d2am_dxidxj(i,j)/am + b(i)*b(j)*tmp5 )
    enddo
  !}
  enddo
!}
end subroutine fugacities_n_its_derivatives3
!***************************************************************************

!***************************************************************************
! water species index is assumed to be n !!!
subroutine findEquilibrium_search( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    k, & ! major oil species' index
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    type_k, & ! vector of binary interaction types
    coef_ab,& ! vector of a, b coefficients
    x_a, & ! vector of mass fractions: alpha phase
    x_b, & ! vector of mass fractions: beta phase
    s_min& ! minimum evaluation value
    )
!{
  implicit none
  integer :: n, k
  integer :: type_k(n)
  real(8) :: coef_ab(n), delta(n,n), dx, xk
  real(8) :: P,T,Pc(n),Tc(n),w(n),x_a(n),x_b(n),fuga_a(n),fuga_b(n),x_bo(n)
  real(8) :: G_a, G_b, s_min
  integer :: i, j, opposite_phase, im

  ! Gex Method for equilibrium
  integer, parameter :: m=1000, m2=100
  real(8) :: x1(n),x2(n),xr,s,xa(0:m,n),xa2(0:m2,n)

  ! check whether x_a is the oil-rich phase or the opposite one
  opposite_phase = 0
  if (x_a(n)>x_b(n)) opposite_phase = 1

  ! initialize the phases
  ! copy x_a & x_b to x1(oil phase) & x2(water phase)
  x1=x_a
  x2=x_b
  if (opposite_phase .eq. 1) then
    x1=x_b
    x2=x_a
  endif

  ! phase grid size
  s=0
  do j=1,n-1
    if (j .ne. k) s=s+x1(j)
  enddo
  dx= (1.d0-s)/dfloat(m)

  ! generating phase points xa(m,n) for oil phases
  do i=0,m 
    xa(i,:) = x1
    ! oil-rich phase
    xr = dfloat(i)*dx
    xa(i,k) = xr
    xa(i,n) = 1.d0-s-xr
  enddo                                                                  
  
  call findEquilibrium_search_core(P,T,k,n,Pc,Tc,w,type_k,coef_ab,m,xa,im,x1,x2)
  !print *,'S'
  
!  if (im .le. 2) return                                                                  
!  xk = xa(im-2,k)                                                                        
!  dx =(xa(im+2,k)-xk)/dfloat(m2)                                                         
!  do i=0,m2                                                                              
!    xa2(i,:) = x1                                                                        
!    xr = xk + dx*dfloat(i)                                                               
!    xa2(i,k) = xr                                                                        
!    xa2(i,n) = 1.d0-s-xr                                                                 
!  enddo                                                                                  
!                                                                                         
!  call findEquilibrium_search_core(P,T,k,n,Pc,Tc,w,type_k,coef_ab,m2,xa2,im,x1,x2)       
                                                                                          
!  if (im .le. 2) return                                                                  
!  xk = xa2(im-2,k)                                                                       
!  dx =(xa2(im+2,k)-xk)/dfloat(m2)                                                        
!  do i=0,m2                                                                              
!    xa2(i,:) = x1                                                                        
!    xr = xk + dx*dfloat(i)                                                               
!    xa2(i,k) = xr                                                                        
!    xa2(i,n) = 1.d0-s-xr                                                                 
!  enddo                                                                                  
!                                                                                         
!  call findEquilibrium_search_core(P,T,k,n,Pc,Tc,w,type_k,coef_ab,m2,xa2,im,x1,x2)       
                                                                                         
  call fugacities(P,T,n,Pc,Tc,w,x1,type_k,coef_ab,fuga_a,G_a,delta)
  call fugacities(P,T,n,Pc,Tc,w,x2,type_k,coef_ab,fuga_b,G_b,delta)

  !print *,'dxf', fuga_a*x1 - fuga_b*x2

  ! evaluate the equilibrium criteria
  s = 0
  do j=1,n
    s = max(s , dabs(x1(j)*fuga_a(j)-x2(j)*fuga_b(j)))
  enddo
  s_min = s

  if (s>1d-3) then ! obtained through testing
    x1=x2
  endif

  ! reset the phase
  if (x1(n) .lt. x2(n)) then
    x_a = x1    
    x_b = x2    
    if (opposite_phase .eq. 1) then
      x_a = x2    
      x_b = x1    
    endif
  else
    x_a = x2    
    x_b = x1    
    if (opposite_phase .eq. 1) then
      x_a = x1    
      x_b = x2    
    endif
  endif
                                                                                          
!}
end subroutine findEquilibrium_search
!***************************************************************************

!***************************************************************************
! water species index is assumed to be n !!!
subroutine findEquilibrium_search_core( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    k, & ! major oil species' index
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    type_k, & ! vector of binary interaction types
    coef_ab,& ! vector of a, b coefficients
    m,   & ! first index size of xa
    xa,  & ! matrix of mass fractions: alpha phase
    im,  & ! index of solution
    x_a, & ! vector of mass fractions: alpha phase - output
    x_b  & ! vector of mass fractions: beta phase  - output
    )
!{
  implicit none
  integer :: n,k
  integer :: type_k(n)
  real(8) :: coef_ab(n), delta(n,n)
  real(8) :: P,T,Pc(n),Tc(n),w(n),x_a(n),x_b(n),fuga_a(n),fuga_b(n),x_bo(n)
  real(8) :: G_a, G_b
  integer :: i,j,i1,i_min,total_converged_points,im

  integer :: m
  real(8) :: x1(n),x2(n),Ga(0:m),xa(0:m,n),fa(0:m,n),vs(0:m),&
             xr,s,min_s,xb(0:m,n),s0, x_a0(n), x_b0(n)

  ! BEGIN tmp calculations (will stop here for debugging)
  x_a(1)=0.999; x_a(2)=0.001;
  x_b(1)=0.001; x_b(2)=0.999;
  do i=0, m
    call fugacities(P,T,n,Pc,Tc,w,x_a,type_k,coef_ab,fuga_a,G_a,delta)
    call fugacities(P,T,n,Pc,Tc,w,x_b,type_k,coef_ab,fuga_b,G_b,delta)
    x1 = x_a; x2 = x_b
    x_a(1) = x2(1)*fuga_b(1)/fuga_a(1)
    x_a(2) = x2(2)*fuga_b(2)/fuga_a(2)
    x_b(1) = x1(1)*fuga_a(1)/fuga_b(1)
    x_b(2) = x1(2)*fuga_a(2)/fuga_b(2)
    x_a0 = x_a; x_b0 = x_b;
    s = x_a(1) + x_a(2)
    x_a(1) = x_a(1) / s
    x_a(2) = x_a(2) / s
    s = x_b(1) + x_b(2)
    x_b(1) = x_b(1) / s
    x_b(2) = x_b(2) / s
    x_a = (x_a + x1)*0.5
    x_b = (x_b + x2)*0.5

    if (dabs(x_a(2)-x_b(2))<1d-3) then
      x_b(1) = 0.001
      x_b(2) = 0.999
    endif

    print *, 'X', x_a0, x_b0, x_a, x_b
    s = 0
    do j=1,2
      s = max(s, dabs(x_a(j) - x1(j)))
      s = max(s, dabs(x_b(j) - x2(j)))
    enddo
    if (s<1d-6) exit
  enddo
  return
  ! END tmp calculations (will stop here for debugging)

  ! calculating fugacities
  min_s = 1d100
  total_converged_points = 0
  do i=0, m
    if (xa(i,k) .le. 1d-4) then
      vs(i) = -1
      cycle
    endif

    call fugacities(P,T,n,Pc,Tc,w,xa(i,:),type_k,coef_ab,fa(i,:),Ga(i),delta)

    ! initial water phase
    do j=1,n-1
      x_b(j) = 1d-4 ! this value must <= 1e-4 (tested value)
    enddo
    x_b(n) = 1 - 1d-4*dfloat(n-1)

    ! converge the water phase
    do i1=1,100
      call fugacities(P,T,n,Pc,Tc,w,x_b,type_k,coef_ab,fuga_b,G_b,delta)

      s = 0
      x_bo = x_b
      !do j=1,n-1
      !  x_b(j) = xa(i,j)*fa(i,j)/fuga_b(j)
      !  s = s + x_b(j)
      !enddo
      !x_b(n) = 1 - s
      ! NEW ITERATION CONVERGE BETTER THAN FORMER ONE
      do j=1,n
        x_b(j) = xa(i,j)*fa(i,j)/fuga_b(j)
        s = s + x_b(j)
      enddo
      do j=1,n
        x_b(j) = x_b(j)/s
      enddo

      s = 0
      do j=1,n
        s = max(s, dabs(x_b(j)-x_bo(j)))
      enddo

      if (s .le. 1d-8) exit
    enddo
    !print *,'i1',i1

    if (s .gt. 1d-8) then 
      vs(i) = -1
      cycle
    endif

    xb(i,:) = x_b

    ! re-calcuate fugacity for the water phase
    call fugacities(P,T,n,Pc,Tc,w,x_b,type_k,coef_ab,fuga_b,G_b,delta)

    ! evaluate the equilibrium criteria
    s = 0
    do j=1,n
      s = max(s, dabs(xa(i,j)*fa(i,j)-x_b(j)*fuga_b(j)))
    enddo
    s0= s
    s = s / dabs(xa(i,k)-x_b(k))**1.9 ! this coef is tested (CANNOT BE TOO BIG)

    !if (s<10) then
      print *,'S',i,xa(i,:),x_b,fa(i,:),fuga_b,s,s0
    !  print *,'S',xa(i,k),s,s0
    !endif
    !print *,'S',xa(i,n),s
    vs(i) = s

    ! find the minimum
    if (s .lt. min_s) then
      min_s = s
      i_min = i
      total_converged_points = total_converged_points + 1
    endif
  enddo                                                                                   

  ! IMPORTANT CONDITION: because near UCST, there may not be any converged
  ! points
  if (total_converged_points .eq. 0) then
    x_a=xa(0,:)
    x_b=xa(0,:)
    im = 0
  else
    !print *, 'min_s', min_s
    x_a=xa(i_min,:)
    x_b=xb(i_min,:)
    im = i_min
  endif

!}
end subroutine findEquilibrium_search_core

!***************************************************************************
! water species index is assumed to be n !!!
subroutine findEquilibrium_new( &
    bPrint, & ! 1: print debugging ; 0: do not print
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    k, & ! major oil species' index
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    type_k, & ! vector of binary interaction types
    coef_ab,& ! vector of a, b coefficients
    x_a, & ! vector of mass fractions: alpha phase
    x_b, & ! vector of mass fractions: beta phase
    s_min& ! minimum evaluation value
    )
!{
  implicit none
  integer :: n, k, bPrint
  integer :: type_k(n)
  real(8) :: coef_ab(n), delta(n,n)
  real(8) :: P,T,Pc(n),Tc(n),w(n),x_a(n),x_b(n),fuga_a(n),fuga_b(n),&
             x_bo(n),x_ao(n),tmp(n)
  real(8) :: G_a, G_b, s_min
  integer :: i, j, opposite_phase, i1, n_backward

  ! Gex Method for equilibrium
  real(8) :: x10(n),x1(n),x2(n),s,s1,fixed

!  ! check whether x_a is the oil-rich phase or the opposite one
!  opposite_phase = 0
!  if (x_a(n)>x_b(n)) opposite_phase = 1

  ! initialize the phases
  ! copy x_a & x_b to x1(oil phase) & x2(water phase)
  x1=x_a
  x2=x_b
!  if (opposite_phase .eq. 1) then
!    x1=x_b
!    x2=x_a
!  endif
!  x10 = x1

  ! phase grid size
  s=0
  do j=1,n-1
    if (j .ne. k) s=s+x1(j)
  enddo
  fixed = 1-s

  x1(k) = fixed-0e-3
  x1(n) = 0e-3

  ! initial water phase
  do j=1,n-1
    x2(j) = 0d-3 
  enddo
  x2(n) = 1 - 0d-3*dfloat(n-1)

  ! converge the water phase
  n_backward = 0
  do i=1,1000
    call fugacities(P,T,n,Pc,Tc,w,x1,type_k,coef_ab,fuga_a,G_a,delta)
    call fugacities(P,T,n,Pc,Tc,w,x2,type_k,coef_ab,fuga_b,G_b,delta)

    if (bPrint>0) then
      print *,'x1',x1
      print *,'x2',x2
      print *,'f1',fuga_a
      print *,'f2',fuga_b
      !print *,'sum',sum(x1),sum(x2)
    endif

    x_ao = x1
    x_bo = x2

    !print *,'k',k
    s1=0
    do j=1,n-1
      if (j .ne. k) s1=s1 + x1(j)*fuga_a(j)/fuga_b(j)
    enddo

    x1(k) = (fixed*fuga_a(n)-(1-s1)*fuga_b(n)) &
                            /                  &
            (fuga_a(n)-fuga_b(n)*fuga_a(k)/fuga_b(k))
    x1(n) = (fixed*fuga_a(k)-(1-s1)*fuga_b(k)) &
                            /                  &
            (fuga_a(k)-fuga_b(k)*fuga_a(n)/fuga_b(n))
    x2(n) =-(fixed*fuga_a(k)-(1-s1)*fuga_b(k)) &
                            /                  &
            (fuga_b(k)-fuga_a(k)*fuga_b(n)/fuga_a(n))
    x2(k) =-(fixed*fuga_a(n)-(1-s1)*fuga_b(n)) &
                            /                  &
            (fuga_b(n)-fuga_a(n)*fuga_b(k)/fuga_a(k))

    if (bPrint>0) then
      print *,'Xa', x1
      print *,'Xb', x2
    endif

    x1(n) = x1(n)*0.5d0 + x_ao(n)*0.5d0
    x1(k) = fixed - x1(n)

    if (bPrint>0) then
      print *,'X1', x1
      print *,'X2', x2
    endif

    if (x1(k)<0) x1(k)=-x1(k)
    if (x1(n)<0) x1(n)=-x1(n)
    if (x2(k)<0) x2(k)=-x2(k)
    if (x2(n)<0) x2(n)=-x2(n)

    ! scale x1
    s = x1(k) + x1(n)
    x1(k) = x1(k)/s*fixed
    x1(n) = fixed-x1(k)

    do j=1,n-1
      if (j .ne. k) then
        x2(j) = x_ao(j)*fuga_a(j)/fuga_b(j)
      endif 
    enddo

    ! scale x2
    s = 0
    do j=1,n
      s = s + x2(j)
    enddo
    do j=1,n
      x2(j) = x2(j)/s
    enddo

    ! re-init if converges for 3 times
    if (n_backward < 3) then
      if (dabs(x2(n)-x1(n))<1d-3) then
        s = 0                     !
        do j=1,n-1                !
          s = s + x2(j)           !
        enddo                     !
                                   
        !do j=1,n-1                
        !  x2(j) = x2(j)/s*1d-4    
        !enddo                     
        !x2(n) = 1-1d-4            
                                   
        do j=1,n-1                !
          x2(j) = 0               !
        enddo                     !
        x2(n) = 1                 !

        n_backward = n_backward + 1
       endif
    endif

    if (x2(n) < x1(n)) then
      tmp = x2
      x2 = x1
      x1 = tmp
      do j=1,n-1
        if (j.ne.k) x1(j)=x10(j)
      enddo
    endif

    s = 0
    do j=1,n
      s = max(s, dabs(x2(j)-x_bo(j)))
    enddo
    s = max(s, dabs(x1(k)-x_ao(k)))
    s = max(s, dabs(x1(n)-x_ao(n)))

    if (s .le. 1d-6) exit

  enddo

  ! evaluate the equilibrium criteria
  call fugacities(P,T,n,Pc,Tc,w,x1,type_k,coef_ab,fuga_a,G_a,delta)
  call fugacities(P,T,n,Pc,Tc,w,x2,type_k,coef_ab,fuga_b,G_b,delta)
  !print *, 'df',x1*fuga_a-x2*fuga_b
  s = 0
  do j=1,n
    s = max(s , dabs(x1(j)*fuga_a(j)-x2(j)*fuga_b(j)))
  enddo
  s_min = s
  !print *,'times',i,s
  !print *,sum(x1),sum(x2)

  ! reset the phase
  x_a = x1    
  x_b = x2    
!  if (opposite_phase .eq. 1) then
!    x_a = x2    
!    x_b = x1    
!  endif
                                                                                          
!}
end subroutine findEquilibrium_new
!***************************************************************************

subroutine findEquilibrium_fix_water_b( &
    bPrint, & ! 1: print debugging ; 0: do not print
    bInitial, & ! 1: using a fresh value as I.C.; 0: using existing value as I.C.
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    kij, & ! matrix of BIPs
    x_a, & ! vector of mass fractions: alpha phase
    x_b, & ! vector of mass fractions: beta phase
    s_min& ! minimum evaluation value
    )
!{
  implicit none
  integer :: n, bPrint, bInitial, bCloseToMiscible
  real(8) :: P,T,Pc(n),Tc(n),w(n),kij(n,n),x_a(n),x_b(n),s_min
  real(8) :: fuga_a(n),fuga_b(n),x_bo(n),x_ao(n),tmp(n),K(n)
  real(8) :: G_a,G_b
  integer :: i,j,opposite_phase,i1

  ! Gex Method for equilibrium
  real(8) :: x1(n),x2(n),s1,s2,fixed1,fixed2

  if (n<3) then
    print *, 'Error: species number must >=3'
    print *, 'n=', n
    return
  endif

  ! initialize the phases
  x1=x_a
  x2=x_b

  ! phase grid size
  s1 = 0
  do j=3,n-1
    s1=s1+x1(j)
  enddo
  fixed1 = 1-s1
  fixed2 = 1-x2(n)
  !print *, 'fixed1', fixed1
  !print *, 'fixed2', fixed2

  bCloseToMiscible = 0
  !if (dabs(x1(n)-x2(n)) .lt. 1e-2) bCloseToMiscible = 1

  ! set up x1 & x2
  if (bInitial==1 .or. bCloseToMiscible==1) then 
    x1(n) = x1(n) - 0.1
    if (x1(n)<0) then 
      x1(n) = 0
      do j=1,2
        x1(j) = fixed1*0.5
      enddo
    else
      do j=1,2
        x1(j) = x1(j) + 0.05
      enddo
    endif

    !do j=1,n-2
    !  x2(j) = 1d0/dble(n-2)
    !enddo
    !x2(n-1) = 0d0
    !do j=1,n-1
    !  x2(j) = fixed2/dble(n-1)
    !enddo

  endif

  !print *, x1, x2

  ! converge the water phase
  do i=1,100000
    x_ao = x1
    x_bo = x2

    call fugacities2(P,T,n,Pc,Tc,w,kij,x1,fuga_a,G_a)
    call fugacities2(P,T,n,Pc,Tc,w,kij,x2,fuga_b,G_b)

    K = fuga_a/fuga_b

    !if (bPrint>0) then
      !print *,'x1',x1
      !print *,'x2',x2
      !print *,'f1',fuga_a
      !print *,'f2',fuga_b
      !print *,'sum',sum(x1),sum(x2)
    !endif

    s1=1.0
    s2=1.0
    do j=3,n-1
      s1 = s1 - x1(j)
      s2 = s2 - x1(j)*K(j)
    enddo
    s1 = s1 - x2(n)/K(n)
    s2 = s2 - x2(n)

    x1(2) = (s1*K(1)-s2)/(K(1)-K(2))
    x1(1) = s1 - x1(2)

    if (x1(2)<0) x1(2) = -x1(2)
    if (x1(1)<0) x1(1) = -x1(1)

    do j=1,n-1
      x2(j) = x_ao(j)*K(j)
    enddo
    x1(n) = x_bo(n)/K(n)

    ! scale x1
    s1 = x1(1) + x1(2) + x1(n)

    s2 = 0
    do j=1,n-1
      s2 = s2 + x2(j)
    enddo

    j=1; x1(j) = x1(j)/s1*fixed1
    j=2; x1(j) = x1(j)/s1*fixed1
    j=n; x1(j) = x1(j)/s1*fixed1

    do j=1,n-1
      x2(j) = x2(j)/s2*fixed2
    enddo

    !do j=1,n
    !  x1(j) = (x1(j) + x_ao(j))*0.5d0
    !  x2(j) = (x2(j) + x_bo(j))*0.5d0
    !enddo

    if (bPrint>0) then
      print *,'Xa', x1
      print *,'Xb', x2
    endif

    s1 = 0
    do j=1,n
      s1 = max(s1, dabs(x2(j)-x_bo(j)))
      s1 = max(s1, dabs(x1(j)-x_ao(j)))
    enddo

    !print *, 'S1', s1

    if (s1 .le. 1d-8) exit

  enddo

  ! evaluate the equilibrium criteria
  call fugacities2(P,T,n,Pc,Tc,w,kij,x1,fuga_a,G_a)
  call fugacities2(P,T,n,Pc,Tc,w,kij,x2,fuga_b,G_b)
  !print *, 'steps:', i
  !print *, x1*fuga_a
  !print *, x2*fuga_b
  !print *, 'df',x1*fuga_a-x2*fuga_b
  s1 = 0
  do j=1,n
    s1 = max(s1, dabs(x1(j)*fuga_a(j)-x2(j)*fuga_b(j)))
  enddo
  s_min = s1
  !print *,'times',i,s
  !print *,sum(x1),sum(x2)

  x_a = x1    
  x_b = x2    
                                                                                          
!}
end subroutine findEquilibrium_fix_water_b
!***************************************************************************

!***************************************************************************
! water species index is assumed to be n !!!
subroutine findEquilibrium_new2( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    k, & ! major oil species' index
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    kij, & ! matrix of BIPs
    x_a, & ! vector of mass fractions: alpha phase
    x_b, & ! vector of mass fractions: beta phase
    s_min  & ! minimum evaluation value
    )
!{
  implicit none
  integer :: n, k, bPrint
  real(8) :: P,T,Pc(n),Tc(n),w(n),kij(n,n),x_a(n),x_b(n),s_min
  real(8) :: fuga_a(n),fuga_b(n),x_bo(n),x_ao(n),tmp(n)
  real(8) :: G_a,G_b
  integer :: i,j,opposite_phase,i1,n_backward

  ! Gex Method for equilibrium
  real(8) :: x10(n),x1(n),x2(n),s,s1,fixed

  bPrint = 0 ! 1: print debugging ; 0: do not print

  ! initialize the phases
  ! copy x_a & x_b to x1(oil phase) & x2(water phase)
  x1=x_a
  x2=x_b

  ! phase grid size
  s=0
  do j=1,n-1
    if (j .ne. k) s=s+x1(j)
  enddo
  fixed = 1-s

  x1(k) = fixed-0e-3
  x1(n) = 0e-3

  ! initial water phase
  do j=1,n-1
    x2(j) = 1d-6 
  enddo
  x2(n) = 1 - 1d-6*dfloat(n-1)

  ! print *,'x1: ',x1
  ! print *,'x2: ',x2
  ! print *,'Ts: ',T

  ! converge the water phase
  n_backward = 0
  do i=1,100000
    call fugacities2(P,T,n,Pc,Tc,w,kij,x1,fuga_a,G_a)
    call fugacities2(P,T,n,Pc,Tc,w,kij,x2,fuga_b,G_b)

    if (bPrint>0) then
      print *,'x1: ',x1
      print *,'x2: ',x2
      print *,'f1: ',fuga_a
      print *,'f2: ',fuga_b
      print *,'sum: ',sum(x1),' ',sum(x2)
    endif

    x_ao = x1
    x_bo = x2

    !print *,'k',k
    s1=0
    do j=1,n-1
      if (j .ne. k) s1=s1 + x1(j)*fuga_a(j)/fuga_b(j)
    enddo

    x1(k) = (fixed*fuga_a(n)-(1-s1)*fuga_b(n)) &
                            /                  &
            (fuga_a(n)-fuga_b(n)*fuga_a(k)/fuga_b(k))

    if (x1(k)<0) x1(k)=-x1(k)

    do j=1,n-1
      x2(j) = fuga_a(j)*x1(j)/fuga_b(j)
    enddo

    ! x1(n), x2(n)
    x1(n) = 1 
    x2(n) = 1 
    do j=1,n-1
      x1(n) = x1(n) - x1(j)
      x2(n) = x2(n) - x2(j)
    end do

    if (x1(n)<0) x1(n)=-x1(n)
    if (x2(n)<0) x2(n)=-x2(n)

    ! scale x1
    s = x1(k) + x1(n)
    x1(k) = x1(k)/s*fixed
    x1(n) = fixed-x1(k)

    ! scale x2
    s = 0
    do j=1,n
      s = s + x2(j)
    enddo
    do j=1,n
      x2(j) = x2(j)/s
    enddo

    ! re-init if converges for 3 times
    if (n_backward < 10) then
      if (dabs(x2(n)-x1(n))<1d-3) then
        s = 0                     !
        do j=1,n-1                !
          s = s + x2(j)           !
        enddo                     !
                                   
        !do j=1,n-1                
        !  x2(j) = x2(j)/s*1d-4    
        !enddo                     
        !x2(n) = 1-1d-4            
                                   
        do j=1,n-1                !
          x2(j) = 0               !
        enddo                     !
        x2(n) = 1                 !

        n_backward = n_backward + 1
       endif
    endif

    s = 0
    do j=1,n
      s = max(s, dabs(x2(j)-x_bo(j)))
    enddo
    s = max(s, dabs(x1(k)-x_ao(k)))
    s = max(s, dabs(x1(n)-x_ao(n)))

    if (s .le. 1d-6) then
      !print *,'X1',x1
      !print *,'X2',x2
      exit
    endif

  enddo

  ! evaluate the equilibrium criteria
  call fugacities2(P,T,n,Pc,Tc,w,kij,x1,fuga_a,G_a)
  call fugacities2(P,T,n,Pc,Tc,w,kij,x2,fuga_b,G_b)
  !print *, 'df',x1*fuga_a-x2*fuga_b
  s = 0
  do j=1,n
    s = max(s , dabs(x1(j)*fuga_a(j)-x2(j)*fuga_b(j)))
  enddo
  s_min = s

  ! print *,'iters: ',i
  ! print *,'s: ',s
  ! print *,'x1: ',x1
  ! print *,'x2: ',x2
  ! print *,'sum(x1): ',sum(x1),' sum(x2): ',sum(x2)

  ! reset the phase
  x_a = x1    
  x_b = x2    
                                                                                          
!}
end subroutine findEquilibrium_new2
!***************************************************************************

!***************************************************************************
subroutine phase_stability_iteration2( &
    stable, & ! OUTPUT: 1 - stable; 0 - unstable
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    kij, & ! matrix of BIPs
    z, & ! vector of original phase mole fractions
    h, & ! vector of original phase chemical potentials
    Yi,& ! vector of inital guess
    Ye,& ! vector of final results
    y  & ! vector of final mole fractions
    )
!{
  implicit none
  ! input and output variables
  integer :: n, stable
  real(8) :: P,T,Pc(n),Tc(n),w(n),kij(n,n),z(n),h(n),Yi(n),Ye(n),y(n)

  ! local variables
  ! used
  integer :: i, j
  real(8) :: lnphi_y(n), sumY, y1(n), y2(n), max_dy, max_dya, dy
  real(8) :: sumY1, sumY2, one, g, beta, r, dgdY(n), dgdY_max

  one = 1d0 + 1d-8

  Ye = Yi

  sumY = 0d0
  do i=1,n
    sumY = sumY + Ye(i)
  end do
  do i=1,n
    y(i) = Ye(i)/sumY 
    y1(i) = y(i)
    y2(i) = y(i)
  end do

  sumY1 = sumY
  sumY2 = sumY

  do j=1,1008

    !print *, 'phase stability iteration # ', j

    call fugacities3(P,T,n,Pc,Tc,w,kij,y,lnphi_y)

    g=1.0d0
    beta=0.0d0
    dgdY_max=0.0d0
    do i=1,n
       dgdY(i) = dlog( Ye(i) ) + lnphi_y(i) - h(i)
       beta = beta + ( Ye(i) - z(i) )*dgdY(i)
       g = g + Ye(i)*( dgdY(i) - 1.0d0 )
       if (dabs(dgdY(i)) > dgdY_max) then
          dgdY_max = dabs(dgdY(i))
       end if
    end do    

    do i=1,n
      Ye(i) = dexp( h(i) - lnphi_y(i) )
    end do
    
    sumY = 0d0
    do i=1,n
      sumY = sumY + Ye(i)
    end do
    do i=1,n
      y(i) = Ye(i)/sumY 
    end do

    !print *, 'dgdY max ', dgdY_max
    if (j > 20) then
       if (dgdY_max < 1.0d-8) then
          exit
       end if

       r = 2.0d0*g/beta
       !print *, 'r ', r
       if (r > 0.8d0) then
          Ye = z
          sumY = 1.0d0       
          exit
       end if
    end if

    ! max_dy  = 0d0
    ! max_dya = 0d0
    ! do i=1,n

    !   dy = dabs(y(i) - y1(i))
    !   if (dy > max_dy) max_dy = dy

    !   ! when oscillation occurs
    !   ! two values are got alternately
    !   dy = dabs(y(i) - y2(i))
    !   if (dy > max_dya) max_dya = dy

    ! end do    

    ! if (j>20) then
    !   if (max_dy < 1d-9 .or. max_dya < 1d-9) then
    !     exit
    !   end if
    ! end if    

    ! y2 = y1
    ! y1 = y
    ! sumY2 = sumY1
    ! sumY1 = sumY

  end do

  stable = 1
  if (sumY > one) stable = 0

  !print *, 'stability test'
  !print *, '# iterations ', j
  !print *, j, stable,  sumY1, sumY2
  !print *, 'y        z'
  !do i=1,n
  !   print *, y(i), '  ', z(i)
  !end do
  !print *
!}
end subroutine phase_stability_iteration2


! test for phase stability and provides initial K values for LLE calculation
subroutine phase_stability2( &
    stable, & ! OUTPUT: 1 - stable; 0 - unstable
    lnK, & ! initial lnK values for LLE (corresponding to test phase)
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    kij, & ! matrix of BIPs
    z  & ! vector of inital mole fractions
    )
!{
  implicit none
  ! input and output variables
  integer :: n, stable
  real(8) :: P,T,Pc(n),Tc(n),w(n),kij(n,n),z(n),lnK(n)

  ! local variables
  ! used
  integer :: i, j, stablej 
  real(8) :: lnphi_z(n), h(n), Ye(n), y(n), K(n), Yi(n,4)
  real(8) :: lnphi_y(n), sumY, dGmin, dG

  !print *, 'z', z

  call fugacities3(P,T,n,Pc,Tc,w,kij,z,lnphi_z)

  do i=1,n
    h(i) = lnphi_z(i)
    if (z(i) .gt. 1e-8) then
      h(i) = h(i) + dlog( z(i) )
    end if

    K(i) = Pc(i)/P * dexp ( 5.42*(1-Tc(i)/T) )
  end do
  
  do i=1,n
     Yi(i,1) = K(i)*z(i)
     Yi(i,2) = z(i)/K(i)
     Yi(i,3) = dexp(h(i))
     Yi(i,4) = (Yi(i,1) + Yi(i,2) + Yi(i,3))/3.0d0
  end do
                              
  dGmin=1.0d0
  stable=1
  do j=1,4
     call phase_stability_iteration2(stablej,P,T,n,Pc,Tc,w,kij,z,h,Yi(:,j),Ye,y)
     if (stablej == 0) then
        call fugacities3(P,T,n,Pc,Tc,w,kij,y,lnphi_y)
        dG=0.0d0
        do i=1,n
           dG = dG + y(i)*(lnphi_y(i) + dlog( y(i) ) - h(i))
        end do
        if (dG < dGmin) then
           lnK = lnphi_z - lnphi_y
        end if
     end if
     stable = stable*stablej
  end do

  if (stable==1) then
     lnK = 0.0d0
  end if

  ! print *, 'stable', stable
  ! print *, 'lnK test phase'
  ! do i=1,n
  !    print *, lnK(i)
  ! end do
  ! print *
!}
end subroutine phase_stability2
!***************************************************************************

!***************************************************************************
! Latest LLE solver based on Rachford-Rice method
subroutine species_LLE4( &
    P, &! pressure (Pa)
    T, &! temperature (K)
    n, &! # of species
    Pc,& ! vector of Pc (Pa)
    Tc,& ! vector of Tc (K)
    w, & ! vector of acentric factor
    kij, &! matrix of BIPs
    xL,&! molar fractions of left  side of interface
    xR,&! molar fractions of right side of interface
    c1, c0, &! original mixing fractions: c1*xL+(1-c1)*xR resulting ==> c0*x1 + (1-c0)*x2
    x1, &! equilibrium molar fractions for left side (OUTPUT)
    x2, &! equilibrium molar fractions for right side (OUTPUT)
    n_miscible &! >0: miscible; <=0: immiscible (OUTPUT)
    )
!{
  implicit none
  ! INPUT & OUTPUT
  integer :: n, n_miscible
  real(8) :: P,T,Pc(n),Tc(n),w(n),kij(n,n),xL(n),xR(n),c1,c0,x1(n),x2(n)
  ! local variables
  real(8) :: lnphi1(n), lnphi2(n), lnK(n)
  real(8) :: beta  ! mole fraction of phase 2, here R or water phase
  real(8) :: lnK_old(n)  ! values of the last step
  integer :: i, j, maxLoop=100000, stable
  real(8) :: tol = 1d-7, K(n), xC(n)
  real(8) :: x1j_tmp, err_lnK(n), err_norm_lnK        


  do j=1,n                 
    xC(j) = xL(j)*c1+xR(j)*(1-c1)
  end do

  call phase_stability2(stable,lnK,P,T,n,Pc,Tc,w,kij,xC)

  if (stable == 1) then
     n_miscible = 1
     c0 = 0.0d0
     x1 = xC
     x2 = xC
  else
     n_miscible = -1
     do i=1,maxLoop
! Successive Substitution iterations
        !print *, 'SSI # ', i
        lnK_old = lnK
        K = dexp(lnK)
        call calc_beta_RR(beta,n,xC,K)
        do j=1,n
           x1(j) = xC(j)/(1.0d0 + beta*(K(j) - 1.0d0))
           x2(j) = K(j)*x1(j)
        end do
        c0 = 1.0d0 - beta

        call fugacities3(P,T,n,Pc,Tc,w,kij,x1,lnphi1)
        call fugacities3(P,T,n,Pc,Tc,w,kij,x2,lnphi2)

        lnK = lnphi1 - lnphi2
        err_lnK = dabs(lnK - lnK_old)
        err_norm_lnK = 0.0d0
        do j=1,n
           if (err_lnK(j) > err_norm_lnK) then
              err_norm_lnK = err_lnK(j)
           end if
        end do
        !print *, 'err_norm_lnK ', err_norm_lnK 
        !print *
        if (err_norm_lnK < tol) then
           exit
        end if
     end do

     ! check oil-rich and water-rich
     if (x1(n)>x2(n)) then
        do j=1,n
           x1j_tmp = x1(j)
           x1(j) = x2(j)
           x2(j) = x1j_tmp          
        end do
        c0 = 1.0d0 - c0
     end if     
  end if  
!}
end subroutine species_LLE4


!
subroutine calc_beta_RR( &
     beta, &! OUTPUT mole fraction of phase 2
     n, &! # of species
     xC, &! total mole numbers vector
     K &! K-values vector
     )
!{
  implicit none
  ! INPUT & OUTPUT
  integer :: n
  real(8) :: xC(n),K(n),beta
  ! local variables
  integer :: i, j, MAX_BRENT_ITERS=1008
  real(8) :: tol=1.0d-7, a, b, fa, fb, s1, s2, Km1(n), tmp
  real(8) :: b_old, fb_old, b_tmp, fb_tmp

  Km1 = K - 1.0d0

  a = 0.0d0
  fa = 0.0d0
  do j=1,n     
     fa = fa + xC(j)*Km1(j)/(1.0d0 + a*Km1(j))
  end do

  b = 1.0d0
  fb = 0.0d0
  do j=1,n     
     fb = fb + xC(j)*Km1(j)/(1.0d0 + b*Km1(j))
  end do

  if (dabs(fa) < dabs(fb)) then
     tmp = a
     a = b
     b = tmp
     tmp = fa
     fa = fb
     fb = tmp
  end if

  b_old = a
  fb_old = fa

  do i=1,MAX_BRENT_ITERS

     b_tmp = b
     fb_tmp = fb

     if (dabs(fb - fb_old) < tol) then
        b = 0.5d0*(b + a)
     else
        s1 = b - fb*(b - b_old)/(fb - fb_old)
        s2 = 0.5d0*(b + a)
        if (dabs(b - s1) < dabs(b - s2)) then
           b = s1
        else
           b = s2
        end if
     end if

     fb = 0.0d0
     do j=1,n     
        fb = fb + xC(j)*Km1(j)/(1.0d0 + b*Km1(j))
     end do

     b_old = b_tmp
     fb_old = fb_tmp
     
     if (fb*fa > 0) then
        a = b_old
        fa = fb_old
     end if

     if (dabs(fa) < dabs(fb)) then
        tmp = a
        a = b
        b = tmp
        tmp = fa
        fa = fb
        fb = tmp
     end if

     if (dabs(fb) < tol) then
        exit
     end if
  end do

  beta = b
  !print *, 'Phase fraction (beta) calculation'
  !print *, 'beta ', beta
  !print *, '# Brent iterations ', i  
  !print *, 'error f(beta) ', dabs(fb)
!}
end subroutine calc_beta_RR
!***************************************************************************

!***************************************************************************
! binary LLE solver 
subroutine binaryLLE( &
    P, &! pressure (Pa)
    T, &! temperature (K)
    n, &! # of species
    Pc,& ! vector of Pc (Pa)
    Tc,& ! vector of Tc (K)
    w, & ! vector of acentric factor
    kij, & ! matrix of BIPs    
    x1, & ! equilibrium molar fractions for oil phase (OUTPUT)
    x2, & ! equilibrium molar fractions for water phase (OUTPUT)
    n_miscible & ! >0: miscible; <=0: immiscible (OUTPUT)
    )
!{
  implicit none
  ! INPUT & OUTPUT
  integer :: n, n_miscible
  real(8) :: P,T,Pc(n),Tc(n),w(n),kij(n,n),x1(n),x2(n)
  ! local variables
  real(8) :: lnphi1(n), lnphi2(n), lnK(n), K(n), lnK_old(n)    
  integer :: i, j, maxLoop=100008
  real(8) :: tol=1.0d-7
  real(8) :: x1j_tmp, err_lnK(n), err_norm_lnK        


  n_miscible = -1
  x1(2) = 1.0d-6
  x1(1) = 1.0d0 - x1(2)
  x2(1) = 1.0d-6
  x2(2) = 1.0d0 - x2(1)  

  call fugacities3(P,T,n,Pc,Tc,w,kij,x1,lnphi1)
  call fugacities3(P,T,n,Pc,Tc,w,kij,x2,lnphi2)

  lnK = lnphi1 - lnphi2

  do i=1,maxLoop
     ! Successive Substitution iterations
     !print *, 'SSI # ', i
     lnK_old = lnK
     K = dexp(lnK)
     
     x1(1) = ( 1.0d0 - K(2) )/( K(1) - K(2) )
     x2(1) = K(1)*x1(1)          
     x1(2) = 1.0d0 - x1(1)
     x2(2) = 1.0d0 - x2(1)

     call fugacities3(P,T,n,Pc,Tc,w,kij,x1,lnphi1)
     call fugacities3(P,T,n,Pc,Tc,w,kij,x2,lnphi2)

     lnK = lnphi1 - lnphi2
     err_lnK = dabs(lnK - lnK_old)
     err_norm_lnK = 0.0d0
     do j=1,n
        if (err_lnK(j) > err_norm_lnK) then
           err_norm_lnK = err_lnK(j)
        end if
     end do
     !print *, 'err_norm_lnK ', err_norm_lnK 
     !print *, 'x1 ', x1
     !print *, 'x2 ', x2
     if (err_norm_lnK < tol) then
        exit
     end if
  end do

  if (dabs( x1(n) - x2(n) ) < 1.0d-4) then
     n_miscible = 1
  end if

  ! check oil-rich and water-rich
  if (x1(n)>x2(n)) then
     do j=1,n
        x1j_tmp = x1(j)
        x1(j) = x2(j)
        x2(j) = x1j_tmp          
     end do     
  end if

  ! print *, '------binaryLLE result-------'
  ! print *, 'err_norm_lnK ', err_norm_lnK 
  ! print *, 'x1 ', x1
  ! print *, 'x2 ', x2
  ! print *, 'n_miscible ', n_miscible
  ! print *, '------binaryLLE result-------'  
!}
end subroutine binaryLLE
!***************************************************************************

!***************************************************************************
subroutine calc_dkijdT( &
    n,    & ! number of species
    a,    & ! vector of ai
    b,    & ! vector of bi
    dadT, & ! vector of d(ai)/dT
    F,    & ! matrix of F(i,j)
    dFdT, & ! matrix of d(Fij)/dT
    dkdT  & ! matrix of d(kij)/dT
    )
!{
  implicit none
  integer :: n, i, j
  real(8) :: a(n),b(n),dadT(n),F(n,n),dFdT(n,n),dkdT(n,n),sqrta_bi, sqrta_bj

  do i=1,n
    dkdT(i,i) = 0d0
    do j=i+1,n
        sqrta_bi = dsqrt(a(i))/b(i)
        sqrta_bj = dsqrt(a(j))/b(j)

        dkdT(i,j) = -(0.5d0*dFdT(i,j)+(sqrta_bi-sqrta_bj) &
                        *(sqrta_bi/a(i)*dadT(i)-sqrta_bj/a(j)*dadT(j))) & !
                   * 0.5d0/(sqrta_bi*sqrta_bj)                          & !
                    +(0.5d0*F(i,j)+(sqrta_bi-sqrta_bj)*(sqrta_bi-sqrta_bj)) &!
                        *(dadT(i)*a(j)+a(i)*dadT(j))                        &!
                   * .25d0/(sqrta_bi*sqrta_bj*a(i)*a(j))             

        dkdT(j,i) = dkdT(i,j)
    enddo
  enddo

!}
end subroutine calc_dkijdT
!***************************************************************************
subroutine calc_am_bm_and_damdT( &
    n,    & ! number of species
    x,    & ! vector of molar fractions
    a,    & ! vector of ai
    b,    & ! vector of bi
    dadT, & ! vector of d(ai)/dT
    kij,  & ! matrix of BIP
    dkdT, & ! matrix of d(kij)/dT
    am, bm, ami, &
    damdT,& ! value  of d(am)/dT
    damidT& ! vector of d(ami)/dT
    )
!{
  implicit none
  integer :: n, i, j
  real(8) :: x(n),a(n),b(n),dadT(n),kij(n,n),dkdT(n,n),damidT(n),&
             am, bm, ami(n), damdT,sqrt_aij

  am = 0d0
  bm = 0d0
  damdT = 0d0
  do i=1,n
    ami(i) = 0d0
    damidT(i) = 0d0

    do j=1,n
      sqrt_aij = dsqrt(a(i)*a(j))

      ami(i) = ami(i) + x(j)*(1-kij(i,j))*sqrt_aij

      damidT(i) = damidT(i) + x(j)*( - dkdT(i,j)*sqrt_aij            &
                    + .5d0*(1d0-kij(i,j))*(dadT(i)*a(j)+a(i)*dadT(j))/sqrt_aij )
      
    enddo

    am = am + x(i)*ami(i)
    bm = bm + x(i)*b(i)
    damdT = damdT + x(i)*damidT(i)
  enddo
!}
end subroutine calc_am_bm_and_damdT
!***************************************************************************
subroutine calc_dVdT( &
    T, V, V_star, am, bm, damdT, dVdT &
    )
!{
  implicit none
  real(8) :: T, V, V_star, am, bm, damdT, dVdT, r1_V_bm, r1_V_star
  real(8) :: R_gas = 8.3144621d0

  r1_V_bm = 1d0/(V - bm)
  r1_V_star = 1d0/V_star

  dVdT = (R_gas*r1_V_bm - damdT*r1_V_star) &
       / (R_gas*T*r1_V_bm*r1_V_bm - 2*am*(V+bm)*r1_V_star*r1_V_star)
!}
end subroutine calc_dVdT
!***************************************************************************

subroutine get_dkdT( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    x, & ! vector of mass fractions
    type_k, & ! vector of binary interaction types
    coef_ab,& ! vector of a, b coefficients
    kij,    & ! matrix of BIP
    dkdT    & ! matrix of d(BIP)/dT
    )
!{
  implicit none
  integer :: n
  real(8) :: P,T,Pc(n),Tc(n),w(n),x(n),lnphi(n),dlnphidT(n)
  integer :: type_k(n)
  real(8) :: coef_ab(n)

  real(8) :: a(n), b(n), dadT(n), kij(n,n), Fij(n,n), dFdT(n,n), dkdT(n,n), &
             am, bm, ami(n), damdT, damidT(n), V, V_star, dVdT, &
             C1, dC1dT, C2, dC2dT, C3(n), dC3dT(n), C(n), dCdT(n)

  integer :: i
  real(8) :: R_gas = 8.3144621d0
  character :: state = 'm'

  call calculate_a_b_and_dadT(T,n,Pc,Tc,w,coef_ab,a,b,dadT)

  call calculate_kij_and_dkijdT(1,T,n,Pc,Tc,w,coef_ab,type_k,kij,Fij,dFdT)

  call calc_dkijdT(n,a,b,dadT,Fij,dFdT,dkdT)
!}
end subroutine get_dkdT

subroutine get_dadT( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    x, & ! vector of mass fractions
    type_k, & ! vector of binary interaction types
    coef_ab,& ! vector of a, b coefficients
    a,      & ! vector of a(n)
    dadT,   & ! vector of da/dT
    ami, damidT, & ! vector of ami(n) and dami/dT
    am,  damdT   & ! value of damdT
    )
!{
  implicit none
  integer :: n
  real(8) :: P,T,Pc(n),Tc(n),w(n),x(n),lnphi(n),dlnphidT(n)
  integer :: type_k(n)
  real(8) :: coef_ab(n)

  real(8) :: a(n), b(n), dadT(n), kij(n,n), Fij(n,n), dFdT(n,n), dkdT(n,n), &
             am, bm, ami(n), damdT, damidT(n), V, V_star, dVdT, &
             C1, dC1dT, C2, dC2dT, C3(n), dC3dT(n), C(n), dCdT(n)

  integer :: i
  real(8) :: R_gas = 8.3144621d0
  character :: state = 'm'

  call calculate_a_b_and_dadT(T,n,Pc,Tc,w,coef_ab,a,b,dadT)

  call calculate_kij_and_dkijdT(1,T,n,Pc,Tc,w,coef_ab,type_k,kij,Fij,dFdT)

  call calc_dkijdT(n,a,b,dadT,Fij,dFdT,dkdT)

  call calc_am_bm_and_damdT(n,x,a,b,dadT,kij,dkdT,am,bm,ami,damdT,damidT)
!}
end subroutine get_dadT

subroutine get_dVdT( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    x, & ! vector of mass fractions
    type_k, & ! vector of binary interaction types
    coef_ab,& ! vector of a, b coefficients
    V,      & ! value of a(n)
    dVdT    & ! value of da/dT
    )
!{
  implicit none
  integer :: n
  real(8) :: P,T,Pc(n),Tc(n),w(n),x(n),lnphi(n),dlnphidT(n)
  integer :: type_k(n)
  real(8) :: coef_ab(n)

  real(8) :: a(n), b(n), dadT(n), kij(n,n), Fij(n,n), dFdT(n,n), dkdT(n,n), &
             am, bm, ami(n), damdT, damidT(n), V, V_star, dVdT, &
             C1, dC1dT, C2, dC2dT, C3(n), dC3dT(n), C(n), dCdT(n)

  integer :: i
  real(8) :: R_gas = 8.3144621d0
  character :: state = 'm'

  call calculate_a_b_and_dadT(T,n,Pc,Tc,w,coef_ab,a,b,dadT)

  call calculate_kij_and_dkijdT(1,T,n,Pc,Tc,w,coef_ab,type_k,kij,Fij,dFdT)

  call calc_dkijdT(n,a,b,dadT,Fij,dFdT,dkdT)

  call calc_am_bm_and_damdT(n,x,a,b,dadT,kij,dkdT,am,bm,ami,damdT,damidT)

  call PR_vol(P,T,am,bm,V,state)

  V_star = V*V + 2d0*bm*V - bm*bm

  call calc_dVdT(T,V,V_star,am,bm,damdT,dVdT)
!}
end subroutine get_dVdT

