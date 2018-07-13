! CALLING ORDER
! 
! 1) V, CP, CV, Hpar, Vpar
!   subroutine thermo_properties(P,T,n,Pc,Tc,w,x,Tb,SG,H_8,type_k, &
!                             V,CP,CV,Hpar,Vpar,dVdT,G,lnphi,G_only)
! 2) cond_m, vis_m
!   subroutine vis_n_cond(P,T,n,Pc,Tc,Vc,x,Cv_m,V,cond_m,vis_m)

!***************************************************************************
! subroutine vis_n_cond(P,T,n,Pc,Tc,Vc,w,M,k,dm0,x,Cv_m,V,cond_m,vis_m)
! Calculates thermal conductivity of mixture and viscosity
! from the Chung 1988 Ind. Eng. Chem. Res. paper.
! Input:  
!      P: pressure in Pa
!      T: temperature in K
!      n: number of speices
!      Pc: unit Pa
!      Tc: unit K
!      Vc: unit cm^3/mol
!      w: acentric factor
!      M: molecular weight
!      k: association parameter
!      dm0: dipole moment, will be normalized as dm=131.3dm0/dsqrt(VcTc)
!      REMOVED                  !ek: epsilon/k (k is Boltzmann constant), potential energy
!      REMOVED                  !si: sigma, potential distance
!      x: mole fraction of the mixture
!      Cv_m: specific heat at constant volume   in J/mol.K
!      V: volume in m^3/mol (computed from PR-EOS)
!
! Output:  
!      cond_m: thermal conductivity of mixture in (W/m.K)
!      vis_m : viscosity of mixture in Pa.sec (NOT Poise) (1P=100cP=0.1Pa.s)
!***************************************************************************

subroutine vis_n_cond(P,T,n,Pc,Tc,Vc,w,M,k,dm,x,Cv_m,V,cond_m,vis_m)
!{
  implicit none

  integer :: n
  real(8) :: P,T,Pc(n),Tc(n),Vc(n), w(n), M(n), k(n), dm(n), &  
             x(n),Cv_m,V,cond_m,vis_m, ek(n),si(n), &
             sigma3(n,n), eps_k(n,n), omega(n,n), kappa(n,n), MW(n,n) 
             
  real(8) :: beta_i(n), beta_ij(n,n), tw(n)

  real(8), dimension(7) :: a1,b1,c1,d1, B
  real(8), dimension(10):: a2,b2,c2,d2, E
  real(8) :: si_m3, Ts_m, ek_m, MW_m, w_m, dm_m4, k_m, Fc_m, Vc_m, Tc_m, vis_mo, &
        collision, Tr, alpha, beta, Z, Psi, y, G1, G2, q, dmr_m, etas, etass
  integer :: i,j
  real(8) :: R_gas = 8.3144621d0

  do i=1,n
    ek(i) = Tc(i)/1.2593         !Eq.(9-4.7)
    si(i) = 0.809*Vc(i)**(1./3.)
    beta_i(i)  = 0.7862 - 0.7109*w(i) + 1.3168*w(i)**2

    ! WATER
    if (dabs(Tc(i)-647).lt.0.5d0 .AND. dabs(Pc(i)-220.6d5).lt.0.5d5) then
      beta_i(i) = 1.0d0/.78d0       
      !print *, "vis_n_cond gets water", i
    endif
  enddo

  do i=1,n
    do j=1,n
      sigma3(i,j) = dsqrt( si(i) * si(j) )**3
      eps_k(i,j)  = dsqrt( ek(i) * ek(j) )
      omega(i,j)  = 0.5*( w (i) + w (j) )
      kappa(i,j)  = dsqrt( k (i) * k (j) )
      MW(i,j)     = 2.0*( M (i) * M (j) )/( M (i) + M (j) )
      beta_ij(i,j)= 0.5*( beta_i(i) + beta_i(j) )
    enddo
  enddo

  !a1 = (/ 2.4166, -0.5092, 6.6107,14.5425, 0.7927, -5.8634, 81.1710/)! Chung et al. 1988
  !a1 = (/ 2.4166, -0.5092, 6.6107,14.5425, 0.7927, -5.8634, 91.0890/)! Homer & Kayar 1994
  a1 = (/  2.4166, -0.5092, 6.6107,14.5425, 0.7927, -5.8634, 88.9077/)! Fit by Ping He
  b1 = (/  0.7482, -1.5094, 5.6207,-8.9139, 0.8202, 12.8005,114.1580/) 
  c1 = (/ -0.9186,-49.9912,64.7599,-5.6379,-0.6937,  9.5893,-60.8410/)
  d1 = (/121.7210, 69.9834,27.0389,74.3435, 6.3173,-65.5292,466.7750/)

  a2 = (/ 6.3240d+0, 1.2102d-3, 5.2835d+0, 6.6226d+0, 1.9745d+1, &
         -1.8999d+0, 2.4275d+1, 7.9716d-1,-2.3816d-1, 6.8629d-2/)
  b2 = (/ 5.0412d+1,-1.1536d-3, 2.5421d+2, 3.8096d+1, 7.6303d+0, &
         -1.2537d+1, 3.4495d+0, 1.1176d+0, 6.7695d-2, 3.4793d-1/)
  c2 = (/-5.1680d+1,-6.2571d-3,-1.6848d+2,-8.4641d+0,-1.4354d+1, &
          4.9853d+0,-1.1291d+1, 1.2348d-2,-8.1630d-1, 5.9260d-1/)
  d2 = (/ 1.1890d+3, 3.7283d-2, 3.8983d+3, 3.1418d+1, 3.1527d+1, &
         -1.8151d+1, 6.9347d+1,-4.1166d+0, 4.0253d+0,-7.2663d-1/)

  si_m3 = 0d0
  Ts_m  = 0d0
  ek_m  = 0d0
  MW_m  = 0d0
  w_m   = 0d0
  dm_m4 = 0d0
  k_m   = 0d0
  !Cv_m  = 0d0
  beta   = 0d0

  do i = 1, n
    do j = 1, n
      si_m3 = si_m3 + x(i)*x(j)*sigma3(i,j)
      ek_m  = ek_m  + x(i)*x(j)*sigma3(i,j)*eps_k(i,j)
      MW_m  = MW_m  + x(i)*x(j)*sigma3(i,j)**(2./3.)*eps_k(i,j)*dsqrt(MW(i,j))
      w_m   = w_m   + x(i)*x(j)*sigma3(i,j)*omega(i,j)
      !dm_m4 = dm_m4 + x(i)*x(j)/sigma3(i,j)*(dm(i)*dm(j))**2
      dm_m4 = dm_m4 + x(i)*x(j)/(eps_k(i,j)*sigma3(i,j))*(dm(i)*dm(j))**2
      k_m   = k_m   + x(i)*x(j)*kappa(i,j)

      beta   = beta   + x(i)*x(j)*sigma3(i,j)*beta_ij(i,j)
    enddo
    !Cv_m = Cv_m + x(i)*Cv(i)
  enddo

  !si_m = si_m**(1./3.)
  ek_m  = ek_m/si_m3
  Ts_m  = T/ek_m
  MW_m  = ( MW_m/ek_m/si_m3**(2./3.) )**2
  w_m   = w_m/si_m3
  !dm_m4 = dm_m4*si_m3
  dm_m4 = dm_m4*si_m3*ek_m
  beta   = beta/si_m3
  !print *, 'dm_m4', dm_m4

  Vc_m = si_m3/0.809**3 ! in cm^3/mol
  Tc_m = 1.2593*ek_m
  !print *, Vc_m, Tc_m

  ! reduced molecular dipole
  dmr_m = 131.3*dm_m4**0.25/dsqrt(Vc_m*Tc_m)
  !print *, 'dmr_m', dmr_m
  
  Fc_m = 1. - 0.2756*w_m + 0.059035*dmr_m**4 + k_m

  ! collision integral based on Lennard-Jones potensial
  collision = 1.16145 * (Ts_m)**(-0.14874)   &
            + 0.52487 *dexp(-0.77320*Ts_m)   &
            + 2.16178 *dexp(-2.43787*Ts_m)   &
            - 6.345d-4*Ts_m**0.14874*dsin(18.0323*Ts_m**(-0.76830)-7.27371)
  !print *, collision, -6.345d-4*Ts_m**0.14874*dsin(18.0323*Ts_m**(-0.76830)-7.27371)
  !print *, 'collision', collision

  ! conductivity (unit N.s/m^2)
  ! viscosity of the mixture at the limit of zero pressure        ! H. Ping:
  vis_mo = 26.69*Fc_m*dsqrt(MW_m*T)/si_m3**(2./3.)/collision*1d-6  ! unit Poise, 0.1 Pa s
  ! vis_mo has the unit kg/m/s
  !print *, 'vis_mo=', vis_mo

  Tr = T/Tc_m
  
  alpha = Cv_m/R_gas-1.5
  Z     = 2.0+10.5*Tr**2
  
  Psi   = 1.+alpha*( (0.215+0.28288*alpha-1.061*beta+0.26665*Z) &
                    /(0.6366+beta*Z+1.061*alpha*beta) )
  
  B = a1 + b1*w_m + c1*dmr_m**4 + d1*k_m
  
  !********************************************
  ! originally, V calculated as a mole averaged volume
  y  = Vc_m/(V*1d6)/6.  !V converted to cm^3/mol 
  !********************************************
  G1 = (1.-0.5*y)/(1.-y)**3
  G2 = ( B(1)/y*(1.-dexp(-B(4)*y)) + B(2)*G1*dexp(B(5)*y) + B(3)*G1 ) &
     / ( B(1)*B(4)+B(2)+B(3) )
  q  = 0.1272*dsqrt(Tc_m/Mw_m)/Vc_m**(2./3.)

  ! in W/m.K
  cond_m = 3117.9*vis_mo*Psi/MW_m*(1./G2+B(6)*y) + q*B(7)*y**2*dsqrt(Tr)*G2

  ! viscosity
  E = a2 + b2*w_m + c2*dmr_m**4 + d2*k_m

  G2 = ( E(1)/y*(1.-dexp(-E(4)*y)) + E(2)*G1*dexp(E(5)*y) + E(3)*G1 ) &
     / ( E(1)*E(4)+E(2)+E(3) )
  
  etass = E(7)*y**2*G2*dexp(E(8)+E(9)/Ts_m+E(10)/Ts_m**2)
  
  ! viscosity of the mixture in Poise
  vis_m  = vis_mo*(1./G2+E(6)*y) &                         ! vis_k
         + 36.344d-6*dsqrt(MW_m*Tc_m)/Vc_m**(2./3.)*etass   ! vis_p

  ! convert from Poise to Pa*sec
  vis_m = vis_m * 0.1

!}
end subroutine vis_n_cond

!***************************************************************************
! subroutine thermo(p,T,i)
! Calculates volume and specific heats from the Peng-Robinson EOS
! Input:  
!      P: pressure in Pa
!      T: temperature in K
!      n: # of species
!      Pc:
!      Tc:
!      w:
!      MW:
!      x: mole fraction of the mixture
!      Tb:
!      SG:
!      H_8:
!      type_k:
!
!
! Output:  
!      V: volume in m^3/mol
!      Cp: specific heat at constant pressure in J/mol.K
!      Cv: specific heat at constant volume   in J/mol.K
!      Vpar(i): partial molar volume of species i in m^3/mol
!      Hpar(i): partial molar enthalpy of species i in J/mol
!      dVdT: derivative of volume wrt temperature and constant P in m^3/mol.K
!
! Note: For volume calculation a single phase is assumed, i.e. no phase equilibrium
!
!                                                            R      da/dT|V   2
!                                                          (--- - -----------)
!                      T         V+b(1+sqrt(2))  d^2a|      V-b   V^2+2bV-b^2
! Cp-Cp(ideal gas)=---------.log(--------------).----| -T.---------------------- -R
!                  2sqrt(2)b     V+b(1-sqrt(2))  dT^2|V    -RT        2a(V+b)  
!                                                        ------- + -------------
!                                                        (V-b)^2   (V2+2bV-b2)^2
!                      T         V+b(1+sqrt(2))  d^2a|
! Cv-Cv(ideal gas)=---------.log(--------------).----|
!                  2sqrt(2)b     V+b(1-sqrt(2))  dT^2|V
!
!***************************************************************************
subroutine thermo_properties(P,T,n,Pc,Tc,w,MW,x,Tb,SG,H_8,type_k, &
                             V,CP,CV,CP_IG,Hpar,Hdep,Vpar,dVdT,G,lnphi,am,bm,G_only)
!{
  implicit none

  integer :: n
  integer, save :: check
  complex(8),dimension(3) :: vol
  integer :: i,j
  real(8) :: P,T,Pc(n),Tc(n),MW(n),w(n),x(n),dadT(n),d2adT2(n), &
    am, dadTm, d2adT2m, bm, V, a1,a2,a3,a4 , Tr, kappa, alpha,CP, CV, MWave,Vmin,G,Gmin
  real(8) :: Tb(n), SG(n), H_8(n), CP_IG(n), H_IG(n), coef_ab(n)
  real(8) :: apar(n),bpar(n), Hpar(n), Hdep(n), Vterm, Vpar(n), dapardT(n),dVdT,lnphi(n),Z,lnphip(n),vv(n),zi
  real(8) :: term1, term2, term3, term2a, term2b, term2c, term2d, term3a, term3b, term3c, term3d, deri1, deri2, deri3
  real(8) :: ac(n), a(n), b(n), delta(n,n), delta_d(n,n), delta_d2(n,n)
  integer :: type_k(n)
  integer :: G_only
  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0
  character :: state = 'm'

  coef_ab = -1

  ! PR coefficients of pure components and their derivatives
  do i=1,n
  !{
    Tr = T/Tc(i)
    if (w(i) .le. 0.491d0) then
      kappa = (0.374640d0+1.54226d0*w(i)-0.269920d0*w(i)**2)
    else
      kappa = (0.379642d0+1.48503d0*w(i)-0.164423d0*w(i)**2+0.0166667d0*w(i)**3)
    endif
    alpha = (1d0+kappa*(1d0-dsqrt(Tr)))**2

    ac(i) = 0.457235529d0*(R_gas*Tc(i))**2/Pc(i) 
    b(i)  =  0.0777960739* R_gas*Tc(i) / Pc(i)

    a(i)      =  ac(i)*alpha
    dadT(i)   = -ac(i)*kappa/dsqrt(T*Tc(i))*alpha
    d2adT2(i) =  ac(i)*kappa/dsqrt(T*Tc(i))/T/2.*(1.+kappa)
  !}
  enddo

  !call calculate_kij(0,T,n,Pc,Tc,w,type_k,delta)
  call calculate_kij_from_table(0,T,n,delta)

  ! PR coefficients of the mixture and its derivatives

  bm = 0d0
  am = 0d0
  dadTm = 0d0
  d2adT2m = 0d0

  do i=1,n
  !{
    bm      = bm + x(i)*b(i)
    am      = am + x(i)**2*a(i)
    dadTm   = dadTm + x(i)**2*dadT(i) 
    d2adT2m = d2adT2m + x(i)**2*d2adT2(i)
    delta_d (i,i) = 0.0d0
    delta_d2(i,i) = 0.0d0

    do j=i+1,n

      delta_d (i,j) = 0.0d0
      delta_d (j,i) = 0.0d0
      delta_d2(i,j) = 0.0d0
      delta_d2(j,i) = 0.0d0

      am      = am + 2.d0*x(i)*x(j)*(1.-delta(i,j))*dsqrt(a(i)*a(j))

      dadTm   = dadTm + x(i)*x(j)*(1.d0-delta(i,j))  &
                * ( dsqrt(a(i)/a(j))*dadT(j)+dsqrt(a(j)/a(i))*dadT(i) )  &
              + 2.d0*x(i)*x(j)* dsqrt(a(i)*a(j))*(-1.d0)*delta_d(i,j)

      d2adT2m = d2adT2m + x(i)*x(j) * (1.d0-delta(i,j)) &
                * ( dsqrt(a(i)/a(j)) * ( d2adT2(j)-0.5d0/a(j)*dadt(j)**2d0 ) + &
                    dsqrt(a(j)/a(i)) * ( d2adT2(i)-0.5d0/a(i)*dadt(i)**2d0 ) + &
                    1.d0/dsqrt(a(i)*a(j))*dadT(i)*dadT(j) )  &
              + 2.d0*x(i)*x(j)* ( (dsqrt(a(i)/a(j))*dadT(j)+dsqrt(a(j)/a(i))*dadT(i)) &
                               *(-1.)*delta_d(i,j) &
                               + dsqrt(a(i)*a(j))*(-1.d0)*delta_d2(i,j))
    enddo
  !}
  enddo  

  !print*, "*1"
  call PR_vol(P,T,am,bm,V,state)
  !print*, "*2"
  do i=1,n
    call PR_vol(P,T,a(i),b(i),vv(i),state)
  enddo
  Vterm = (V**2+2.*V*bm-bm**2)

  ! P.HE: checked to be correct!
  ! computed using direct differential operations on EoS 
  dVdT = (R_gas*Vterm-dadTm*(V-bm)) &
       / (P*(3.*V**2+2.*bm*V-3*bm**2)-R_gas*T*2.*(V+bm)+am)

  ! fugacities
  Z = P*V/R_gas/T
  G = 0.
  do i=1,n
  !{
    lnphi(i)=-b(i)/bm
    do j=1,n
      lnphi(i) = lnphi(i) +2.*x(j)*dsqrt(a(i)*a(j))*(1.-delta(i,j))/am
    enddo

    lnphi(i) = -lnphi(i) * sqr2/4. * am/bm/R_gas/T &
               *dlog( (V+bm*(1.+sqr2))/(V+bm*(1.-sqr2)) ) &
             - dlog(Z-bm/V*Z) + b(i)/bm*(Z-1.)

    Zi = P*vv(i)/R_gas/T
    lnphip(i) = -sqr2/4. * a(i)/b(i)/R_gas/T &
                *dlog( (vv(i)+b(i)*(1.+sqr2))/(vv(i)+b(i)*(1.-sqr2)) ) &
              - dlog(Zi-b(i)/vv(i)*Zi) + (Zi-1.)
    G = G + x(i)*lnphi(i) - x(i)*lnphip(i)
    !if(x(i)>0) G = G + x(i)*dlog(x(i))
    if(x(i)>1e-12) G = G + x(i)*dlog(x(i))
  !}
  enddo  

  ! if only G is required for LLE, dont calculate anything else
  if (G_only .gt. 0) return

  ! specific heats residuals
  CV = T/2./sqr2/bm*dlog( (V+bm*(1.+sqr2))/(V+bm*(1.-sqr2)) ) * d2adT2m
  CP = - R_gas + CV &
     - T*( R_gas  /(V-bm)    -     dadTm   / Vterm )**2  &
        /(-R_gas*T/(V-bm)**2 + 2.*am*(V+bm)/ Vterm**2 )
  ! adding the ideal gass specific heats
  do i=1,n
    call IG_CP_H(Tb(i),Tc(i),SG(i),H_8(i),MW(i),T,H_IG(i),CP_IG(i))
    cp = cp + x(i)* CP_IG(i)
    ! Ideal Gas Cp-Cv=R (Mayer's relation)
    cv = cv + x(i)*(CP_IG(i)-R_gas) 
  enddo

  ! calculating partial molar volumes and enthalpies
  do i=1, n
  !{
    ! Equation (15)
    bpar(i) = 1.*(b(i) - bm)
    apar(i) = 0d0
    dapardT(i) = 0d0
    do j=1, n
      apar(i) = apar(i) + x(j)*(1.-delta(i,j))*dsqrt(a(j))
      dapardT(i) = dapardT(i) + x(j)*(1.-delta(i,j))/2. &
                   *( dsqrt(a(j)/a(i))*dadT(i)+dsqrt(a(i)/a(j))*dadT(j) ) &
                 + x(j)* dsqrt(a(j)*a(i))*(-1.)*delta_d(i,j)
    enddo
    apar(i) = 2.*( -am + dsqrt(a(i))*apar(i) )

    dapardT(i) = 2.*( -dadTm + dapardT(i) )

    Vpar(i) = ( R_gas*T*bpar(i)/(V-bm)**2 - apar(i)/Vterm &
               + am*2.*(V-bm)*bpar(i)/Vterm**2 ) &
            / ( R_gas*T/(V-bm)**2 - am*2.*(V+bm)/Vterm**2 )

    Vpar(i) = Vpar(i) +  V

    ! Equation (16)
    Hdep(i) = P*Vpar(i) - R_gas*T + 1./2./sqr2/bm * ( -am+T*dadTm ) &     
                           * dlog((V+(1.+sqr2)*bm)/(V+(1.-sqr2)*bm)) &     
                           *( 1. - bpar(i)/bm                       &     
                             +(-apar(i)+T*dapardT(i))/(-am +T*dadTm)&     
                             - 2.*sqr2*(bm*Vpar(i)-V*b   (i))/Vterm &     
                               /dlog((V+(1.+sqr2)*bm)/(V+(1.-sqr2)*bm)) )  
!   A mistake is here for Guang's original code
!   bpar(i) should be b(i) in the following line
!                            - 2.*sqr2*(bm*Vpar(i)-V*bpar(i))/Vterm &     
!   correction has been made in the above Hpar(i) equation
!   CORRECT ONE IN THE FOLLOWING
!                            - 2.*sqr2*(bm*Vpar(i)-V*b   (i))/Vterm &     

    ! Compare Guang's Eq. with Ashwin's Eq.
!   Vpar(i) = H_IG(i) + P*Vpar(i) - R_gas*T &
!           + (am-T*dadTm)*(Vpar(i)-V*(bpar(i)+bm)/bm)/Vterm &
!           - 1./2./sqr2/bm*dlog((V+(1.+sqr2)*bm)/(V+(1.-sqr2)*bm)) &     
!           *(apar(i)+2.*am-T*(dapardT(i)+2.*dadTm)-(bpar(i)+bm)/bm*(am-T*dadTm))

    ! adding partial enthalpy of the ideal gas                            
    Hpar(i) = Hdep(i) + H_IG(i)                                           !
  !}
  enddo

!}
END subroutine thermo_properties
!***************************************************************************

!***************************************************************************
subroutine calc_v_Cp_h(P,T,x,n,Pc,Tc,w,MW,Tb,SG,H_8,kij, &
                       V,Cp,h)
!{
  implicit none

  integer :: n
  integer :: i,j
  real(8) :: P,T,x(n),Pc(n),Tc(n),w(n),MW(n),Tb(n), SG(n), H_8(n), kij(n,n)
  real(8) :: V, Cp, h
  real(8) :: dadT(n),d2adT2(n), am, dadTm, d2adT2m, bm, Tr, kappa, alpha, Vterm, h_IG_m, Cp_IG_m, term1
  real(8) :: CP_IG(n), H_IG(n)  
  real(8) :: ac(n), a(n), b(n)
  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0
  character :: state = 'm'  

  ! PR coefficients of pure components and their derivatives
  do i=1,n
  !{
    Tr = T/Tc(i)
    if (w(i) .le. 0.491d0) then
      kappa = (0.374640d0+1.54226d0*w(i)-0.269920d0*w(i)**2)
    else
      kappa = (0.379642d0+1.48503d0*w(i)-0.164423d0*w(i)**2+0.0166667d0*w(i)**3)
    endif
    alpha = (1d0+kappa*(1d0-dsqrt(Tr)))**2

    ac(i) = 0.457235529d0*(R_gas*Tc(i))**2/Pc(i) 
    b(i)  =  0.0777960739* R_gas*Tc(i) / Pc(i)

    a(i)      =  ac(i)*alpha
    dadT(i)   = -ac(i)*kappa/dsqrt(T*Tc(i))*alpha
    d2adT2(i) =  ac(i)*kappa/dsqrt(T*Tc(i))/T/2.*(1.+kappa)
  !}
  enddo  

  ! PR coefficients of the mixture and its derivatives

  bm = 0d0
  am = 0d0
  dadTm = 0d0
  d2adT2m = 0d0

  do i=1,n
  !{
    bm      = bm + x(i)*b(i)
    am      = am + x(i)**2*a(i)
    dadTm   = dadTm + x(i)**2*dadT(i) 
    d2adT2m = d2adT2m + x(i)**2*d2adT2(i)
    !kij_d (i,i) = 0.0d0
    !kij_d2(i,i) = 0.0d0

    do j=i+1,n

      !kij_d (i,j) = 0.0d0
      !kij_d (j,i) = 0.0d0
      !kij_d2(i,j) = 0.0d0
      !kij_d2(j,i) = 0.0d0

      am      = am + 2.d0*x(i)*x(j)*(1.-kij(i,j))*dsqrt(a(i)*a(j))

      dadTm   = dadTm + x(i)*x(j)*(1.d0-kij(i,j))  &
                * ( dsqrt(a(i)/a(j))*dadT(j)+dsqrt(a(j)/a(i))*dadT(i) )
              !+ 2.d0*x(i)*x(j)* dsqrt(a(i)*a(j))*(-1.d0)*kij_d(i,j)

      d2adT2m = d2adT2m + x(i)*x(j) * (1.d0-kij(i,j)) &
                * ( dsqrt(a(i)/a(j)) * ( d2adT2(j)-0.5d0/a(j)*dadT(j)**2d0 ) + &
                    dsqrt(a(j)/a(i)) * ( d2adT2(i)-0.5d0/a(i)*dadT(i)**2d0 ) + &
                    1.d0/dsqrt(a(i)*a(j))*dadT(i)*dadT(j) ) !&
              !+ 2.d0*x(i)*x(j)* ( (dsqrt(a(i)/a(j))*dadT(j)+dsqrt(a(j)/a(i))*dadT(i)) &
                               !*(-1.)*kij_d(i,j) &
                               !+ dsqrt(a(i)*a(j))*(-1.d0)*kij_d2(i,j))
    enddo
  !}
  enddo  

  call PR_vol(P,T,am,bm,V,state)  
  
  Vterm = (V**2+2.*V*bm-bm**2)
  
  Cp_IG_m = 0d0
  h_IG_m = 0d0
  ! adding the ideal gass specific heats
  do i=1,n
    call IG_CP_H(Tb(i),Tc(i),SG(i),H_8(i),MW(i),T,H_IG(i),CP_IG(i))
    Cp_IG_m = cp_IG_m + x(i)* CP_IG(i)
    h_IG_m = h_IG_m + x(i)* H_IG(i)
  enddo

  term1 = 1./2./sqr2/bm*dlog( (V+bm*(1.+sqr2))/(V+bm*(1.-sqr2)) )

  h = h_IG_m + P*V - R_gas*T + term1*(T*dadTm - am)
  ! specific heats residuals
  
  Cp = Cp_IG_m - R_gas
  Cp = Cp + T*term1*d2adT2m
  Cp = Cp - T*( R_gas  /(V-bm)    -     dadTm   / Vterm )**2  &
        /(-R_gas*T/(V-bm)**2 + 2.*am*(V+bm)/ Vterm**2 )
!}
END subroutine calc_v_Cp_h
!***************************************************************************

!***************************************************************************
subroutine calc_v_h(P,T,x,n,Pc,Tc,w,MW,Tb,SG,H_8,kij, &
                    V,h)
!{
  implicit none

  integer :: n
  integer :: i,j
  real(8) :: P,T,x(n),Pc(n),Tc(n),w(n),MW(n),Tb(n),SG(n),H_8(n),kij(n,n)
  real(8) :: V, Cp, h
  real(8) :: dadT(n), am, dadTm, bm, Tr, kappa, alpha, Vterm, h_IG_m, term1
  real(8) :: CP_IG(n), H_IG(n)  
  real(8) :: ac(n), a(n), b(n)
  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0
  character :: state = 'm'

  ! PR coefficients of pure components and their derivatives
  do i=1,n
  !{
    Tr = T/Tc(i)
    if (w(i) .le. 0.491d0) then
      kappa = (0.374640d0+1.54226d0*w(i)-0.269920d0*w(i)**2)
    else
      kappa = (0.379642d0+1.48503d0*w(i)-0.164423d0*w(i)**2+0.0166667d0*w(i)**3)
    endif
    alpha = (1d0+kappa*(1d0-dsqrt(Tr)))**2

    ac(i) = 0.457235529d0*(R_gas*Tc(i))**2/Pc(i) 
    b(i)  =  0.0777960739* R_gas*Tc(i) / Pc(i)

    a(i)      =  ac(i)*alpha
    dadT(i)   = -ac(i)*kappa/dsqrt(T*Tc(i))*alpha
  !}
  enddo

  ! PR coefficients of the mixture and its derivatives

  bm = 0d0
  am = 0d0
  dadTm = 0d0

  do i=1,n
  !{
    bm      = bm + x(i)*b(i)
    am      = am + x(i)**2*a(i)
    dadTm   = dadTm + x(i)**2*dadT(i)
    !kij_d (i,i) = 0.0d0

    do j=i+1,n

      !kij_d (i,j) = 0.0d0
      !kij_d (j,i) = 0.0d0      

      am      = am + 2.d0*x(i)*x(j)*(1.-kij(i,j))*dsqrt(a(i)*a(j))

      dadTm   = dadTm + x(i)*x(j)*(1.d0-kij(i,j))  &
                * ( dsqrt(a(i)/a(j))*dadT(j)+dsqrt(a(j)/a(i))*dadT(i) )  !&
              !+ 2.d0*x(i)*x(j)* dsqrt(a(i)*a(j))*(-1.d0)*kij_d(i,j)
    enddo
  !}
  enddo  

  call PR_vol(P,T,am,bm,V,state)  
  
  Vterm = (V**2+2.*V*bm-bm**2)
    
  h_IG_m = 0d0
  ! adding the ideal gass specific heats
  do i=1,n
    call IG_CP_H(Tb(i),Tc(i),SG(i),H_8(i),MW(i),T,H_IG(i),CP_IG(i))
    h_IG_m = h_IG_m + x(i)* H_IG(i)
  enddo

  term1 = 1./2./sqr2/bm*dlog( (V+bm*(1.+sqr2))/(V+bm*(1.-sqr2)) )

  h = h_IG_m + P*V - R_gas*T + term1*(T*dadTm - am)
  
END subroutine calc_v_h
!***************************************************************************

!***************************************************************************
subroutine calc_v_CvIG_Cp_hpar(P,T,x,n,Pc,Tc,w,MW,Tb,SG,H_8,kij, &
                               V,CvIG,Cp,hpar)
!{
  implicit none

  integer :: n
  integer :: i,j
  real(8) :: P,T,x(n),Pc(n),Tc(n),w(n),MW(n),Tb(n),SG(n),H_8(n),kij(n,n)
  real(8) :: V,CvIG,Cp,hpar(n)
  real(8) :: dadT(n),d2adT2(n),am,dadTm,d2adT2m,bm,Tr,kappa,alpha,Vterm,Cp_IG_m,term1
  real(8) :: apar(n),bpar(n),Hdep(n),Vpar(n),dapardT(n)
  real(8) :: CP_IG(n),H_IG(n)
  real(8) :: ac(n),a(n),b(n)
  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0
  character :: state = 'm'  

  ! PR coefficients of pure components and their derivatives
  do i=1,n
  !{
    Tr = T/Tc(i)
    if (w(i) .le. 0.491d0) then
      kappa = (0.374640d0+1.54226d0*w(i)-0.269920d0*w(i)**2)
    else
      kappa = (0.379642d0+1.48503d0*w(i)-0.164423d0*w(i)**2+0.0166667d0*w(i)**3)
    endif
    alpha = (1d0+kappa*(1d0-dsqrt(Tr)))**2

    ac(i) = 0.457235529d0*(R_gas*Tc(i))**2/Pc(i) 
    b(i)  =  0.0777960739* R_gas*Tc(i) / Pc(i)

    a(i)      =  ac(i)*alpha
    dadT(i)   = -ac(i)*kappa/dsqrt(T*Tc(i))*alpha
    d2adT2(i) =  ac(i)*kappa/dsqrt(T*Tc(i))/T/2.*(1.+kappa)
  !}
  enddo

  ! PR coefficients of the mixture and its derivatives

  bm = 0d0
  am = 0d0
  dadTm = 0d0
  d2adT2m = 0d0

  do i=1,n
  !{
    bm      = bm + x(i)*b(i)
    am      = am + x(i)**2*a(i)
    dadTm   = dadTm + x(i)**2*dadT(i) 
    d2adT2m = d2adT2m + x(i)**2*d2adT2(i)
    !kij_d (i,i) = 0.0d0
    !kij_d2(i,i) = 0.0d0

    do j=i+1,n

      !kij_d (i,j) = 0.0d0
      !kij_d (j,i) = 0.0d0
      !kij_d2(i,j) = 0.0d0
      !kij_d2(j,i) = 0.0d0

      am      = am + 2.d0*x(i)*x(j)*(1.-kij(i,j))*dsqrt(a(i)*a(j))

      dadTm   = dadTm + x(i)*x(j)*(1.d0-kij(i,j))  &
                * ( dsqrt(a(i)/a(j))*dadT(j)+dsqrt(a(j)/a(i))*dadT(i) )  !&
              !+ 2.d0*x(i)*x(j)* dsqrt(a(i)*a(j))*(-1.d0)*kij_d(i,j)

      d2adT2m = d2adT2m + x(i)*x(j) * (1.d0-kij(i,j)) &
                * ( dsqrt(a(i)/a(j)) * ( d2adT2(j)-0.5d0/a(j)*dadT(j)**2d0 ) + &
                    dsqrt(a(j)/a(i)) * ( d2adT2(i)-0.5d0/a(i)*dadT(i)**2d0 ) + &
                    1.d0/dsqrt(a(i)*a(j))*dadT(i)*dadT(j) )  !&
              !+ 2.d0*x(i)*x(j)* ( (dsqrt(a(i)/a(j))*dadT(j)+dsqrt(a(j)/a(i))*dadT(i)) &
                               !*(-1.)*kij_d(i,j) &
                               !+ dsqrt(a(i)*a(j))*(-1.d0)*kij_d2(i,j))
    enddo
  !}
  enddo  

  call PR_vol(P,T,am,bm,V,state)  
  
  Vterm = (V**2+2.*V*bm-bm**2)
  
  Cp_IG_m = 0d0
  ! adding the ideal gass specific heats
  do i=1,n
    call IG_CP_H(Tb(i),Tc(i),SG(i),H_8(i),MW(i),T,H_IG(i),CP_IG(i))
    Cp_IG_m = cp_IG_m + x(i)* CP_IG(i)
  enddo

  CvIG = Cp_IG_m - R_gas;

  term1 = 1./2./sqr2/bm*dlog( (V+bm*(1.+sqr2))/(V+bm*(1.-sqr2)) )    
  Cp = Cp_IG_m - R_gas
  Cp = Cp + T*term1*d2adT2m
  Cp = Cp - T*( R_gas  /(V-bm)    -     dadTm   / Vterm )**2  &
        /(-R_gas*T/(V-bm)**2 + 2.*am*(V+bm)/ Vterm**2 )

  ! calculating partial molar volumes and enthalpies
  do i=1, n
  !{
    ! Equation (15)
    bpar(i) = 1.*(b(i) - bm)
    apar(i) = 0d0
    dapardT(i) = 0d0
    do j=1, n
      apar(i) = apar(i) + x(j)*(1.-kij(i,j))*dsqrt(a(j))
      dapardT(i) = dapardT(i) + x(j)*(1.-kij(i,j))/2. &
                   *( dsqrt(a(j)/a(i))*dadT(i)+dsqrt(a(i)/a(j))*dadT(j) ) !&
                 !+ x(j)* dsqrt(a(j)*a(i))*(-1.)*kij_d(i,j)
    enddo
    apar(i) = 2.*( -am + dsqrt(a(i))*apar(i) )

    dapardT(i) = 2.*( -dadTm + dapardT(i) )

    Vpar(i) = ( R_gas*T*bpar(i)/(V-bm)**2 - apar(i)/Vterm &
               + am*2.*(V-bm)*bpar(i)/Vterm**2 ) &
            / ( R_gas*T/(V-bm)**2 - am*2.*(V+bm)/Vterm**2 )

    Vpar(i) = Vpar(i) +  V

    ! Equation (16)
    Hdep(i) = P*Vpar(i) - R_gas*T + 1./2./sqr2/bm * ( -am+T*dadTm ) &     
                           * dlog((V+(1.+sqr2)*bm)/(V+(1.-sqr2)*bm)) &     
                           *( 1. - bpar(i)/bm                       &     
                             +(-apar(i)+T*dapardT(i))/(-am +T*dadTm)&     
                             - 2.*sqr2*(bm*Vpar(i)-V*b   (i))/Vterm &     
                               /dlog((V+(1.+sqr2)*bm)/(V+(1.-sqr2)*bm)) )  
!   A mistake is here for Guang's original code
!   bpar(i) should be b(i) in the following line
!                            - 2.*sqr2*(bm*Vpar(i)-V*bpar(i))/Vterm &     
!   correction has been made in the above Hpar(i) equation
!   CORRECT ONE IN THE FOLLOWING
!                            - 2.*sqr2*(bm*Vpar(i)-V*b   (i))/Vterm &     

    ! Compare Guang's Eq. with Ashwin's Eq.
!   Vpar(i) = H_IG(i) + P*Vpar(i) - R_gas*T &
!           + (am-T*dadTm)*(Vpar(i)-V*(bpar(i)+bm)/bm)/Vterm &
!           - 1./2./sqr2/bm*dlog((V+(1.+sqr2)*bm)/(V+(1.-sqr2)*bm)) &     
!           *(apar(i)+2.*am-T*(dapardT(i)+2.*dadTm)-(bpar(i)+bm)/bm*(am-T*dadTm))

    ! adding partial enthalpy of the ideal gas                            
    hpar(i) = Hdep(i) + H_IG(i)                                           !
  !}
  enddo
!}
END subroutine calc_v_CvIG_Cp_hpar
!***************************************************************************

!***************************************************************************
subroutine calc_hpar(P,T,x,n,Pc,Tc,w,MW,Tb,SG,H_8,kij, &
                                 hpar)
!{
  implicit none

  integer :: n
  integer :: i,j
  real(8) :: P,T,x(n),Pc(n),Tc(n),w(n),MW(n),Tb(n),SG(n),H_8(n),kij(n,n)
  real(8) :: hpar(n)
  real(8) :: V,dadT(n),am,dadTm,bm,Tr,kappa,alpha,Vterm
  real(8) :: apar(n),bpar(n),Hdep(n),Vpar(n),dapardT(n)
  real(8) :: CP_IG(n),H_IG(n)
  real(8) :: ac(n),a(n),b(n)
  real(8) :: R_gas = 8.3144621d0, sqr2 = 1.414213562373095d0
  character :: state = 'm'

  ! PR coefficients of pure components and their derivatives
  do i=1,n
  !{
    Tr = T/Tc(i)
    if (w(i) .le. 0.491d0) then
      kappa = (0.374640d0+1.54226d0*w(i)-0.269920d0*w(i)**2)
    else
      kappa = (0.379642d0+1.48503d0*w(i)-0.164423d0*w(i)**2+0.0166667d0*w(i)**3)
    endif
    alpha = (1d0+kappa*(1d0-dsqrt(Tr)))**2

    ac(i) = 0.457235529d0*(R_gas*Tc(i))**2/Pc(i) 
    b(i)  =  0.0777960739* R_gas*Tc(i) / Pc(i)

    a(i)      =  ac(i)*alpha
    dadT(i)   = -ac(i)*kappa/dsqrt(T*Tc(i))*alpha
  !}
  enddo

  ! PR coefficients of the mixture and its derivatives

  bm = 0d0
  am = 0d0
  dadTm = 0d0

  do i=1,n
  !{
    bm      = bm + x(i)*b(i)
    am      = am + x(i)**2*a(i)
    dadTm   = dadTm + x(i)**2*dadT(i)
    !kij_d (i,i) = 0.0d0
    !kij_d2(i,i) = 0.0d0

    do j=i+1,n

      !kij_d (i,j) = 0.0d0
      !kij_d (j,i) = 0.0d0
      !kij_d2(i,j) = 0.0d0
      !kij_d2(j,i) = 0.0d0

      am      = am + 2.d0*x(i)*x(j)*(1.-kij(i,j))*dsqrt(a(i)*a(j))

      dadTm   = dadTm + x(i)*x(j)*(1.d0-kij(i,j))  &
                * ( dsqrt(a(i)/a(j))*dadT(j)+dsqrt(a(j)/a(i))*dadT(i) )  !&
              !+ 2.d0*x(i)*x(j)* dsqrt(a(i)*a(j))*(-1.d0)*kij_d(i,j)      
    enddo
  !}
  enddo  

  call PR_vol(P,T,am,bm,V,state)  
  
  Vterm = (V**2+2.*V*bm-bm**2)
  
  ! adding the ideal gass specific heats
  do i=1,n
    call IG_CP_H(Tb(i),Tc(i),SG(i),H_8(i),MW(i),T,H_IG(i),CP_IG(i))
  enddo

  ! calculating partial molar volumes and enthalpies
  do i=1, n
  !{
    ! Equation (15)
    bpar(i) = 1.*(b(i) - bm)
    apar(i) = 0d0
    dapardT(i) = 0d0
    do j=1, n
      apar(i) = apar(i) + x(j)*(1.-kij(i,j))*dsqrt(a(j))
      dapardT(i) = dapardT(i) + x(j)*(1.-kij(i,j))/2. &
                   *( dsqrt(a(j)/a(i))*dadT(i)+dsqrt(a(i)/a(j))*dadT(j) ) !&
                 !+ x(j)* dsqrt(a(j)*a(i))*(-1.)*kij_d(i,j)
    enddo
    apar(i) = 2.*( -am + dsqrt(a(i))*apar(i) )

    dapardT(i) = 2.*( -dadTm + dapardT(i) )

    Vpar(i) = ( R_gas*T*bpar(i)/(V-bm)**2 - apar(i)/Vterm &
               + am*2.*(V-bm)*bpar(i)/Vterm**2 ) &
            / ( R_gas*T/(V-bm)**2 - am*2.*(V+bm)/Vterm**2 )

    Vpar(i) = Vpar(i) +  V

    ! Equation (16)
    Hdep(i) = P*Vpar(i) - R_gas*T + 1./2./sqr2/bm * ( -am+T*dadTm ) &     
                           * dlog((V+(1.+sqr2)*bm)/(V+(1.-sqr2)*bm)) &     
                           *( 1. - bpar(i)/bm                       &     
                             +(-apar(i)+T*dapardT(i))/(-am +T*dadTm)&     
                             - 2.*sqr2*(bm*Vpar(i)-V*b   (i))/Vterm &     
                               /dlog((V+(1.+sqr2)*bm)/(V+(1.-sqr2)*bm)) )  
!   A mistake is here for Guang's original code
!   bpar(i) should be b(i) in the following line
!                            - 2.*sqr2*(bm*Vpar(i)-V*bpar(i))/Vterm &     
!   correction has been made in the above Hpar(i) equation
!   CORRECT ONE IN THE FOLLOWING
!                            - 2.*sqr2*(bm*Vpar(i)-V*b   (i))/Vterm &     

    ! Compare Guang's Eq. with Ashwin's Eq.
!   Vpar(i) = H_IG(i) + P*Vpar(i) - R_gas*T &
!           + (am-T*dadTm)*(Vpar(i)-V*(bpar(i)+bm)/bm)/Vterm &
!           - 1./2./sqr2/bm*dlog((V+(1.+sqr2)*bm)/(V+(1.-sqr2)*bm)) &     
!           *(apar(i)+2.*am-T*(dapardT(i)+2.*dadTm)-(bpar(i)+bm)/bm*(am-T*dadTm))

    ! adding partial enthalpy of the ideal gas                            
    hpar(i) = Hdep(i) + H_IG(i)                                           !
  !}
  enddo
!}
END subroutine calc_hpar
!***************************************************************************

!***************************************************************************
subroutine calc_v_CvIG(P,T,x,n,Pc,Tc,w,MW,Tb,SG,H_8,kij, &
                       V,CvIG)
!{
  implicit none

  integer :: n
  integer :: i,j
  real(8) :: P,T,x(n),Pc(n),Tc(n),w(n),MW(n),Tb(n),SG(n),H_8(n),kij(n,n)
  real(8) :: V,CvIG
  real(8) :: am,bm,Tr,kappa,alpha,Cp_IG_m
  real(8) :: CP_IG(n),H_IG(n)
  real(8) :: ac(n),a(n),b(n)
  real(8) :: R_gas = 8.3144621d0
  character :: state = 'm'

  ! PR coefficients of pure components and their derivatives
  do i=1,n
  !{
    Tr = T/Tc(i)
    if (w(i) .le. 0.491d0) then
      kappa = (0.374640d0+1.54226d0*w(i)-0.269920d0*w(i)**2)
    else
      kappa = (0.379642d0+1.48503d0*w(i)-0.164423d0*w(i)**2+0.0166667d0*w(i)**3)
    endif
    alpha = (1d0+kappa*(1d0-dsqrt(Tr)))**2

    ac(i) = 0.457235529d0*(R_gas*Tc(i))**2/Pc(i) 
    b(i)  =  0.0777960739* R_gas*Tc(i) / Pc(i)
    a(i)      =  ac(i)*alpha    
  !}
  enddo

  ! PR coefficients of the mixture and its derivatives

  bm = 0d0
  am = 0d0

  do i=1,n
  !{
    bm      = bm + x(i)*b(i)
    am      = am + x(i)**2*a(i)    
    do j=i+1,n
      am = am + 2.d0*x(i)*x(j)*(1.-kij(i,j))*dsqrt(a(i)*a(j))
    enddo
  !}
  enddo  

  call PR_vol(P,T,am,bm,V,state)  
  
  Cp_IG_m = 0d0
  ! adding the ideal gass specific heats
  do i=1,n
    call IG_CP_H(Tb(i),Tc(i),SG(i),H_8(i),MW(i),T,H_IG(i),CP_IG(i))
    Cp_IG_m = cp_IG_m + x(i)* CP_IG(i)
  enddo

  CvIG = Cp_IG_m - R_gas;  
!}
END subroutine calc_v_CvIG
!***************************************************************************

!function [H Cp H_Btu_lb Cp_Btu_lbF] = Kesler_Lee_Enthalpy_Star2(K, SG, Tc, H_8, T)
subroutine KL_IG_Cp_H(K,SG,Tc,H_8,T,H,Cp)
!{

! Equation (2) from Kesler & Lee, Hydrocarbon Processing, 55, 153-158, 1976
! INPUTS:
!           K   - Watson factor
!           SG  - Specific gravity (unitless @ 60F/60F) (Watson & Nelson, 1933)
!           Tc  - critical temperature (K) 
!           *****PLEASE NOTICE UNIT BTU/LB for H_8*****************************
!           H_8 - specific enthalpy @ Tr=0.8 (Btu/lb) 
!           *******************************************************************
!           T   - temperature (K) 
! OUTPUTS:
!           H          - specific enthalpy ( J/kg ) 
!           Cp         - specific heat capacity ( J/kg/K ) 

  implicit none
  real(8) :: K, SG, Tc, H_8, T, H, Cp, Tb
  real(8) :: CF, A, B, C, T_R, Tc_R, Tc_8_R, Cp_Btu_lbF, H_Btu_lb

! [R] = [K]*1.8;
! [F] = [R] - R;
! Fahrenheit to Rankine conversion constant
  
! unit conversion
  real(8) :: Btu = 1054.35026444d0 ! Joules/Btu
  real(8) :: lb  = 0.45359237d0    ! kg/lb
! F is the ratio for (delta K or delta C) / (delta F or delta R)
  real(8) :: F   = 1.8d0           ! Rankine/Kelvin 

  ! TEST
  ! Because: K = ((Tb+273.15d0)*1.8d0)**(1d0/3d0) / SG
  !Tb = (K*SG)**3d0 / 1.8d0

! variables in Eq. (2)
  CF = ((12.8/K-1)*(10.0/K-1)*100)**2;

  A = -0.32646 + 0.02678*K - CF*(0.084773-0.080809*SG);
  B = -(1.3892-1.2122*K+0.03803*K*K)*1d-4 + CF*(2.1773-2.0826*SG)*1d-4;
  C = -1.5393*1d-7 - CF*(0.78649-0.70423*SG)*1d-7;

! convert units
  T_R    = T *F;
  Tc_R   = Tc*F;
  Tc_8_R = Tc_R*0.8;  ! when Tr=0.8, Rankine value

! Cp
  Cp_Btu_lbF = A + B*T_R + C*T_R**2;
  Cp = Cp_Btu_lbF * Btu/lb*F;

! H
  H_Btu_lb = H_8 + A*(T_R - Tc_8_R) + 0.5d0*B*(T_R**2-Tc_8_R**2) + C*(T_R**3-Tc_8_R**3)/3.0d0;
  H = H_Btu_lb * Btu/lb;
!}
end subroutine KL_IG_Cp_H
 
subroutine IG_CP_H(Tb,Tc,SG,H_8,MW,T,H_IG,CP_IG)
!{
  implicit none
  real(8) :: Tb, Tc, SG, H_8, MW, T, H_IG, CP_IG, K
  real(8) :: A, B, C, D, E, F, G, T0, Dh
  integer :: idx
  real(8) :: coef(21,7)

  coef(1,:)=(/ 40.44511222D0,-0.262900074D0,2.45691231080674D-03,-4.90044507215163D-06,&
               4.63676556231701D-09,-2.15267561854545D-12,3.93327779392278D-16 /)
  coef(2,:)=(/53.47635091D0,-2.80909078497078D-2,1.7824939030726D-3,-3.52700975226758D-6,&
               3.22424229545744D-09,-1.45384913059128D-12,2.59807241264174D-16 /)
  coef(3,:)=(/ 52.73685516D0,-0.162378987D0,2.48350138530583D-03,-4.92933490889021D-06,&
               4.58866189384468D-09,-2.10205538802517D-12,3.80289149242926D-16 /)
  coef(4,:)=(/ 70.27885986D0,-0.178325529D0,2.60627390635009D-03,-4.99750782601509D-06,&
               4.53564774591324D-09,-2.04016382538111D-12,3.64226668188094D-16 /)
  coef(5,:)=(/ 36.52910783D0,-0.135119316D0,2.87344676574249D-03,-6.0083758893473D-06,&
               5.75305864655736D-09,-2.6806911852124D-12,4.90321870829123D-16 /)
  coef(6,:)=(/ 130.7909704D0,-0.590792605D0,4.04954778235944D-03,-7.37792578873729D-06,&
               6.61311825478702D-09,-2.98062036121945D-12,5.38658747325557D-16 /)
  coef(7,:)=(/ 36.2895121D0,-4.88082866754714D-2,3.50007722206282D-3,-7.551151686372D-6,&
               7.3039265025383D-09,-3.41889840909161D-12,6.26845145434836D-16 /)
  coef(8,:)=(/ 198.863437D0,-1.379460452D0,8.20469099708419D-03,-1.55175065496762D-05,&
               1.44566788156306D-08,-6.68260308297633D-12,1.22142183464829D-15 /)
  coef(9,:)=(/ 158.2805373D0,-1.013463244D0,6.85674164636535D-03,-1.28923832848463D-05,&
               1.15960009709803D-08,-5.08805790785147D-12,8.73941281749261D-16 /)
  coef(10,:)=(/ 93.51996032D0,1.129297835D0,-1.2265022385122D-03,3.5493091753117D-06,&
               -5.26174667748648D-09,3.28898900916197D-12,-7.34597063439489D-16 /)
  coef(11,:)=(/ 204.7694141D0,-1.8175202D0,1.00282408676624D-02,-1.94088416638732D-05,&
                1.88308452648452D-08,-9.06858405534736D-12,1.71783037243269D-15 /)
  coef(12,:)=(/ 109.9748292D0,-0.212145607D0,2.57329714581114D-03,-4.55670597006053D-06,&
                3.9006593589516D-09,-1.6968229081491D-12,2.99868150818351D-16 /)
  coef(13,:)=(/ 131.0837487D0,-0.172105686D0,2.84860630464068D-03,-5.04955481955849D-06,&
                4.29704104117114D-09,-1.86228387137944D-12,3.28979281975546D-16 /)
! CO2 - index -24
  coef(14,:)=(/ 23.50610388D0,3.80655731822383D-02,7.40233050458845D-05,-2.22713256071084D-07,&
               2.34374551216989D-10,-1.14647749938175D-13,2.16814726839333D-17/)
! ethanol - index -25
  coef(15,:)=(/ 53.77704843D0,-0.207301605D0,1.42227077769417D-03,-2.62982313956358D-06,&
              2.39648619777309D-09,-1.08912257461417D-12,1.9639329537972D-16/)
  !CO   carbon monoxide                 # -26
  coef(16,:)=(/ 28.50457641D0,1.02017609728895D-02,-6.15947212697459D-05,1.61354408202049D-07,&
             -1.78138017055272D-10,9.02011087812847D-14,-1.73591160415148D-17/)
  !NO2   nitrogen dioxide               # -27
  coef(17,:)=(/ 35.36589116D0,-4.41934679029063D-02,2.76561358668253D-04,-4.74332672039955D-07,&
              3.98251343504695D-10,-1.67881344119746D-13,2.83948907318438D-17/)
  !N2    nitrogen                       # -28
  coef(18,:)=(/ 28.71677138D0,7.34582726095694D-03,-4.54759070492041D-05,1.16406082454094D-07,&
             -1.22458456787057D-10,5.90449193986651D-14,-1.08748042093876D-17/)
  !O2    oxygen                         # -29
  coef(19,:)=(/ 29.79023995D0,-9.48853768266291D-03,2.85799016034241D-05,9.87286057866708D-09,&
             -5.66510613516015D-11,4.30015831906746D-14,-1.0218904651183D-17/)
  !C9H12 1,2,4-trimethylbenzene         # -30
  coef(20,:)=(/ 81.78321263D0,-0.29658314D0,3.09922456937352D-03,-5.88597175944603D-06,&
              5.34239989224643D-09,-2.40491097086861D-12,4.29578195599291D-16/)
  !C12H26        dodecane               # -31
  coef(21,:)=(/ 88.38144398D0,0.588033701D0,5.22064388340535D-04,-9.55903233901003D-07,&
              4.63376217143547D-10,-5.57235289191563D-14,-8.16098726300127D-18/)

  if (dabs(H_8+1d0) .lt. 1d-7) then ! H_8=-1 : water 
    A =  143.05d0                   ! FROM Hill-Peterson propulsion
    B = -183.54d0
    C =  82.751d0
    D = -3.6989d0

    T0 = 373.15 
    !     h of vaporization + h of heat from 0C - 100C
    !     pV1-pV0 is neglected because liquid water has small thermodynamic potential
    Dh = 40.68d3 + 7.5327d3

    H_IG = Dh + A*(T-T0) + B*0.252982212813470d0*(T**1.25-T0**1.25) &
         + C*0.066666666666667d0*(T**1.5-T0**1.5) + D*.005*(T**2-T0**2) 
    CP_IG= A + B*0.316227766016838d0*T**0.25 + C*.1*T**0.5 + D*.01*T
  elseif (dabs(H_8+2d0) .lt. 1d-7) then ! H_8=-2 : n-decane
    A =  3.1780d+1
    B =  7.4489d-1
    C = -1.0945d-4
    D = -2.2670d-7
    E =  9.3458d-11
    T0= 0.
    Dh= 0.
    H_IG = A*(T-T0) + B*(T**2-T0**2)/2. + C*(T**3-T0**3)/3. &  
                    + D*(T**4-T0**4)/4. + E*(T**5-T0**5)/5. + Dh
    CP_IG = A + B*T + C*T**2 + D*T**3 + E*T**4
  elseif (dabs(H_8+3d0) .lt. 1d-7) then ! H_8=-3 : methane (CH4)
    A =  3.4942D+1
    B = -3.9957D-2
    C =  1.9184D-4
    D = -1.5300D-7
    E =  3.9321D-11
    T0= 0.
    Dh= 0.
    H_IG = A*(T-T0) + B*(T**2-T0**2)/2. + C*(T**3-T0**3)/3. &  
                    + D*(T**4-T0**4)/4. + E*(T**5-T0**5)/5. + Dh
    CP_IG = A + B*T + C*T**2 + D*T**3 + E*T**4
  elseif (dabs(H_8+4d0) .lt. 1d-7) then ! H_8=-4 : toluene 
    A =-2.4097D+01
    B = 5.2187D-01
    C =-2.9827D-04
    D = 6.122D-08          
    E = 1.2576D-12
    T0= 0.
    Dh= 0.
    H_IG = A*(T-T0) + B*(T**2-T0**2)/2. + C*(T**3-T0**3)/3. &  
                    + D*(T**4-T0**4)/4. + E*(T**5-T0**5)/5. + Dh
    CP_IG = A + B*T + C*T**2 + D*T**3 + E*T**4
  elseif (dabs(H_8+5d0) .lt. 1d-7) then ! H_8=-5 : triacontane n-C30 alkane
    A =-2.4382D+01
    B = 2.8538D+00 
    C =-1.6082D-03
    D = 3.458D-07
    E = 0.0000D+00
    T0= 0.
    Dh= 0.
    H_IG = A*(T-T0) + B*(T**2-T0**2)/2. + C*(T**3-T0**3)/3. &  
                    + D*(T**4-T0**4)/4. + E*(T**5-T0**5)/5. + Dh
    CP_IG = A + B*T + C*T**2 + D*T**3 + E*T**4
  elseif (dabs(H_8+6d0) .lt. 1d-7) then ! H_8=-6 : n-pentacontane n-C50 alkane
    A = 8.224925043D2
    B =-4.234244214D0
    C = 0.032859313D0
    D =-6.79166D-05 
    E = 6.80157D-08 
    F =-3.33832D-11  
    G = 6.40795D-15
    T0= 0.
    Dh= 0.
    H_IG = A*(T-T0) + B*(T**2-T0**2)/2. + C*(T**3-T0**3)/3. &  
                    + D*(T**4-T0**4)/4. + E*(T**5-T0**5)/5. &
                    + F*(T**6-T0**6)/6. + G*(T**7-T0**7)/7. + Dh
    CP_IG = A + B*T + C*T**2 + D*T**3 + E*T**4 + F*T**5 + G*T**6
!********************************************************************************
! NEW long chain aromatics (1-ring & 2-ring systems)
  elseif (dabs(H_8+7d0) .lt. 1d-7) then ! H_8=-7 : 1-decylbenzene: Benzene-nC10
    A = 173.1732818896  
    B =-9.03973016640795D-02 
    C = 3.39474666972966D-03  
    D =-6.03330011877701D-06 
    E = 5.09415660904854D-09  
    F =-2.19792944862381D-12 
    G = 3.88510308397148D-16
    T0= 0.
    Dh= 0.
    H_IG = A*(T-T0) + B*(T**2-T0**2)/2. + C*(T**3-T0**3)/3. &  
                    + D*(T**4-T0**4)/4. + E*(T**5-T0**5)/5. &
                    + F*(T**6-T0**6)/6. + G*(T**7-T0**7)/7. + Dh
    CP_IG = A + B*T + C*T**2 + D*T**3 + E*T**4 + F*T**5 + G*T**6
  elseif (dabs(H_8+8d0) .lt. 1d-7) then ! H_8=-8 : 1-triacontylbenzene: Benzene-nC30
    A = 582.6646591834  
    B =-2.4805514657 
    C = 1.82052412634452D-02  
    D =-3.42722335146829D-05 
    E = 3.16278252274334D-08  
    F =-1.45785011219947D-11 
    G = 2.67363395228979D-15
    T0= 0.
    Dh= 0.
    H_IG = A*(T-T0) + B*(T**2-T0**2)/2. + C*(T**3-T0**3)/3. &  
                    + D*(T**4-T0**4)/4. + E*(T**5-T0**5)/5. &
                    + F*(T**6-T0**6)/6. + G*(T**7-T0**7)/7. + Dh
    CP_IG = A + B*T + C*T**2 + D*T**3 + E*T**4 + F*T**5 + G*T**6
  elseif (dabs(H_8+9d0) .lt. 1d-7) then ! H_8=-9 : 1-decylnaphthalene: Naph-nC10
    A = 322.5198462402  
    B =-1.7827759055 
    C = 1.12957787829348D-02  
    D =-2.15433848269945D-05 
    E = 2.03387912013662D-08  
    F =-9.58492180095023D-12 
    G = 1.79279170289216D-15
    T0= 0.
    Dh= 0.
    H_IG = A*(T-T0) + B*(T**2-T0**2)/2. + C*(T**3-T0**3)/3. &  
                    + D*(T**4-T0**4)/4. + E*(T**5-T0**5)/5. &
                    + F*(T**6-T0**6)/6. + G*(T**7-T0**7)/7. + Dh
    CP_IG = A + B*T + C*T**2 + D*T**3 + E*T**4 + F*T**5 + G*T**6
  elseif (dabs(H_8+ 10) .lt. 1d-7) then ! H_8=-10: 1-dodecylnaphthalene: Naph-nC12
    A = 356.0951868704  
    B =-1.9151590534 
    C = 1.23195580369635D-02  
    D =-2.35051512108933D-05 
    E = 2.21882658842386D-08  
    F =-1.04589728472607D-11 
    G = 1.95734500347755D-15
    T0= 0.
    Dh= 0.
    H_IG = A*(T-T0) + B*(T**2-T0**2)/2. + C*(T**3-T0**3)/3. &  
                    + D*(T**4-T0**4)/4. + E*(T**5-T0**5)/5. &
                    + F*(T**6-T0**6)/6. + G*(T**7-T0**7)/7. + Dh
    CP_IG = A + B*T + C*T**2 + D*T**3 + E*T**4 + F*T**5 + G*T**6
!********************************************************************************
  elseif (H_8.lt.-10.5d0 .AND. H_8.gt.-31.5d0) then ! new data
    idx = - int(H_8) - 10
    A = coef(idx, 1)
    B = coef(idx, 2)
    C = coef(idx, 3)
    D = coef(idx, 4)
    E = coef(idx, 5)
    F = coef(idx, 6)
    G = coef(idx, 7)
    T0= 0.
    Dh= 0.
    H_IG = A*(T-T0) + B*(T**2-T0**2)/2. + C*(T**3-T0**3)/3. &  
                    + D*(T**4-T0**4)/4. + E*(T**5-T0**5)/5. &
                    + F*(T**6-T0**6)/6. + G*(T**7-T0**7)/7. + Dh
    CP_IG = A + B*T + C*T**2 + D*T**3 + E*T**4 + F*T**5 + G*T**6
!********************************************************************************
  else
    K = ((Tb+273.15d0)*1.8d0)**(1d0/3d0) / SG
    !print *, 'K=', K, 'SG=', SG
    call KL_IG_Cp_H(K,SG,Tc,H_8,T,H_IG,CP_IG)
    ! convert the unit from 'per Kg' to 'per mole'
    H_IG = H_IG * MW * 1d-3
    CP_IG = CP_IG * MW * 1d-3
  endif
!}
end subroutine IG_CP_H

!************************************************************************
! Calculate mass diffusivities:
!   (m2/sec , intger , integer , kelvin , ratio , kg/m3) 
! equations are from Kutney's PhD thesis "2005"
!************************************************************************
!    sigma,           & ! vector of potential distances
!    epsilon_over_k,  & ! vecotr of potential energies
subroutine TLSM_diffusion_ij( &
    rho,             & ! density
    T,               & ! temperature
    n,               & ! # of species
    Pc,              & ! vector of Pc (Pa)
    Tc,              & ! vector of Tc (K)
    Vc,              & ! vector of Vc (cm^3/mol)
    MW,              & ! vector of molecular weights
    x,               & ! vector of molar fractions
    D,               & ! vector of mass diffusivity
    Dij)               ! matrix of binary mass diffusivity 
!{
  implicit none
  real(8), intent(in) :: T ! (K) Temperature distribution in space
  real(8), intent(in) :: rho ! (kg/m3) density ( position )
  integer, intent(in) :: n ! number of points in space , number of species
  real(8), intent(in), dimension(n):: x ! Concentration distribution ( position , species )
  real(8), intent(in), dimension(n):: Tc! (K) 
  real(8), intent(in), dimension(n):: Pc! (Pa)
  real(8), intent(in), dimension(n):: Vc! (cm^3/mol)

  real(8), intent(in), dimension(n):: MW ! (g/mol) molecular masses

  real(8), intent(out), dimension(n) :: D ! (cm2/sec) diffusion coefficients 

  real(8), dimension(n):: sigma ! (cm) molecular diameter, 
    !reference: Kurtney PhD thesis "2005". All are the known values, except water (p. 376)
  real(8), dimension(n):: eps_k ! (K) well-depth petential, 
    !reference: Kurtney PhD thesis "2005". All are the known values, except water (p.376)
  real(8), dimension(n,n) :: Dij          ! (cm2/sec) diffusion coefficients and then converted to m2/sec
                                          ! (solute , in solvent)

  real(8) :: rho_in_cm ! (g/cm3) density
  real(8) :: xi, xj, sum_xj    ! normalized molar fractions
  real(8), parameter :: R_gas  = 8.3144621d0   ! (J/mol.K) , Universal gas constant

  integer::  i, j ! loop counters i:solute_counter, j:solvent_counter

  real(8):: M2    ! molecular mass for solvent
  real(8):: M_1_2 ! mixture molecular mass (calculated from equation 204 )
  real(8):: sigma_1_2 ! Lennard-Jones soft-sphere diameter (calculated from equation 203)
  real(8):: eps_k_1_2 ! TLSM-Lennard-Jones 
    !energy well-depth parameter (calculated from equation 202)
  real(8):: sigma_2_TLSM ! TLSM Lennard-Jones soft-sphere diameter 
    !(calculated from equation 146)
  real(8):: intermediate_variable
  real(8), parameter:: N_A = 6.0221415e23 ! (molecules/mol) Avogadro's number
  real(8), parameter :: two2onethird = 2d0**(1d0/3d0) , &
                        onethird = 1d0/3d0

  ! Liu et al. 1998 data
  do i=1,n
    ! Silva et al. 1996 correlation
    ! THIS CORRELATION DOES NOT WORK WHEN n-alkane carbon # increases > 30
    !sigma(i) = 1d-8*( 0.17791+11.779*(Tc(i)/Pc(i)*1d5) &
    !                 -0.049029*(Tc(i)/Pc(i)*1d5)**2 )**onethird
    sigma(i) = 1d-8*0.809*Vc(i)**onethird
    eps_k(i) = 0.774*Tc(i)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Water is different because of strong H-bond !!!!!!!!!!!!!!!!!!!!!!!!
    if (dabs(Tc(i)-647).lt.0.5d0 .AND. dabs(Pc(i)-220.6d5).lt.0.5d5) then
      sigma(i) = 1.53091d-8
      eps_k(i) = 3788.51
      !print *, "TLSM gets water", i
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !if (i>50) then
    !  print *, "i,sigma", i, sigma(i)
    !endif
  enddo

  rho_in_cm = rho / 1000.

  ! (cm) ( water , toluene , thiophene ) from Silva et al 1998
  !sigma = (/ 1.53e-8, 5.42e-8, 4.86e-8 /)     
  !sigma = (/ 1.53e-8, 6.71e-8, 7.00e-8 /)

  ! (K) ( water , toluene , thiophene ) from Silva et al 1998
  !epsilon_over_k = (/ 3.789e3, 4.58e2, 4.48e2 /)  
  !epsilon_over_k = (/ 3.789e3, 4.35e2, 6.73e2 /)

  ! ( water , toluene , thiophene )
  !M_mass = (/ 18e0, 92.14e0, 84.14e0 /)  
  !M_mass = (/ 18e0, 142.3e0, 170.3e0 /)
  
  do i = 1, n ! loop over all species as a solute
    Dij(i,i) = 0.0;
    do j = i+1, n ! loop over all species as a solvent
      ! prepare xi & xj
      !{
      xi = x(i)
      xj = x(j)
      if ( x(i) < 1.d-10 ) xi = 1.d-10
      if ( x(j) < 1.d-10 ) xj = 1.d-10
      xi = xi/(xi+xj)
      xj = 1-xi
      !!!!M2 = MW(i)
      !!!!if (xj>xi) M2 = MW(j)
      !xi = 0.5
      !xj = 0.5
      !}

      !{
      ! Reference: Kurtney PhD thesis "2005":
      M_1_2 = xi*MW(i) + xj*MW(j) ! equation 204
        
      sigma_1_2 = xi*sigma(i) + xj*sigma(j) ! equation 203
        
      eps_k_1_2 = eps_k(i)**(xi) * eps_k(j)**(xj) ! equation 202
        
      ! "sigma_2_TLSM" is calculated in three steps: equation 146
      intermediate_variable = dsqrt(T/eps_k_1_2)
      intermediate_variable = (1.0 + 1.2 * intermediate_variable)**(1.0/3.0)
      sigma_2_TLSM = two2onethird * ( sigma_1_2 ** 2 ) / intermediate_variable
        
      ! "Diff" is calculated in six steps: equation 201
      intermediate_variable = 1.2588*M_1_2 - rho_in_cm*N_A*(dsqrt(sigma_2_TLSM)**3) 
      !intermediate_variable =  1.2588*M2    - rho_in_cm*N_A*(dsqrt(sigma_2_TLSM)**3) 
      intermediate_variable = (-0.75*rho_in_cm*N_A*(dsqrt(sigma_2_TLSM)**3))&
                            / intermediate_variable
      intermediate_variable = intermediate_variable - 0.27862*eps_k_1_2/T
        
      intermediate_variable = dsqrt(R_gas*T/M_1_2) * dexp( intermediate_variable )
      intermediate_variable = 669.138 * M_1_2 * intermediate_variable
      !intermediate_variable  = 669.138 * M2    * intermediate_variable
        
      Dij(i,j) = intermediate_variable / (rho_in_cm*N_A*sigma_2_TLSM) ! unit: cm^2/sec

      Dij(i,j) = (1d-4) * Dij(i,j) ! convert to m2/sec
        
      !}

      Dij(j,i) = Dij(i,j)
    end do
  end do

  ! HACKING
  !do i = 1, n ! loop over all species as a solute
  !  do j = 1, n ! loop over all species as a solvent
  !    Dij(i,j) = 5e-8
  !  end do
  !end do
  !Dij(1,2) =2.8d-8; Dij(2,1) = Dij(1,2)  
  !Dij(1,3) =2.5d-8; Dij(3,1) = Dij(1,3)  
  !Dij(2,3) =1.7d-8; Dij(3,2) = Dij(2,3)  
  !Dij(1,4) =5.3d-8; Dij(4,1) = Dij(1,4)  
  !Dij(2,4) =5.1d-8; Dij(4,2) = Dij(2,4)  
  !Dij(3,4) =4.9d-8; Dij(4,3) = Dij(3,4)  

  !Dij(1,3) = Dij(1,3)*.5; Dij(3,1) = Dij(1,3) 
  !Dij(2,3) = Dij(2,3)*.5; Dij(3,2) = Dij(2,3) 
  !Dij(3,4) = Dij(3,4)*.5; Dij(4,3) = Dij(3,4) 

  do i = 1, n ! loop over all species as a solute
    D(i) = 0.0d0

    ! S = sum_{i<>j}(xj/Dij)
    sum_xj = 0;
    do j = 1, n ! loop over all species as a solvent
      if ((Dij(i,j)>1d-16) .and. (i.ne.j)) then
        ! make sure diffusivity always yields a non-zero value
        if ( x(j) < 1.d-8 ) xj = 1.d-8
        D(i) = D(i) + xj/Dij(i,j)
        sum_xj = sum_xj + xj
      !else
      !  print *, 'i,j,Dij:', i, j, Dij(i,j)
      end if
    end do

    ! Di = (1-xi)/S
    if (D(i)>1e-16) then
      !D(i) = (1-x(i))/D(i)
      D(i) = sum_xj/D(i)
    else
      !print *, 'i,Di:', i, D(i)
      D(i) = 0.0d0
    end if
  end do

  !return
!}
end subroutine TLSM_diffusion_ij

! trace diffusion coef. will be used in Wesselingh & Krishna 1990 model
! to calc. the multicomponent diffusion coef.
subroutine TLSM_diffusion_trace_ij( &
    rho,             & ! density
    T,               & ! temperature
    n,               & ! # of species
    Pc,              & ! vector of Pc (Pa)
    Tc,              & ! vector of Tc (K)
    Vc,              & ! vector of Vc (cm^3/mol)
    MW,              & ! vector of molecular weights
    x,               & ! vector of molar fractions
    Dij)               ! matrix of binary mass diffusivity 
!{
  implicit none
  real(8), intent(in) :: T ! (K) Temperature distribution in space
  real(8), intent(in) :: rho ! (kg/m3) density ( position )
  integer, intent(in) :: n ! number of points in space , number of species
  real(8), intent(in), dimension(n):: x ! Concentration distribution ( position , species )
  real(8), intent(in), dimension(n):: Tc! (K) 
  real(8), intent(in), dimension(n):: Pc! (Pa)
  real(8), intent(in), dimension(n):: Vc! (cm^3/mol)

  real(8), intent(in), dimension(n):: MW ! (g/mol) molecular masses

  real(8), dimension(n):: sigma ! (cm) molecular diameter, 
    !reference: Kurtney PhD thesis "2005". All are the known values, except water (p. 376)
  real(8), dimension(n):: eps_k ! (K) well-depth petential, 
    !reference: Kurtney PhD thesis "2005". All are the known values, except water (p.376)
  real(8), dimension(n,n) :: Dij          ! (cm2/sec) diffusion coefficients and then converted to m2/sec
                                          ! (solute , in solvent)

  real(8) :: rho_in_cm ! (g/cm3) density
  real(8) :: xi, xj, sum_xj    ! normalized molar fractions
  real(8), parameter :: R_gas  = 8.3144621d0   ! (J/mol.K) , Universal gas constant

  integer::  i, j ! loop counters i:solute_counter, j:solvent_counter

  real(8):: M2    ! molecular mass for solvent
  real(8):: M_1_2 ! mixture molecular mass (calculated from equation 204 )
  real(8):: sigma_1_2 ! Lennard-Jones soft-sphere diameter (calculated from equation 203)
  real(8):: eps_k_1_2 ! TLSM-Lennard-Jones 
    !energy well-depth parameter (calculated from equation 202)
  real(8):: sigma_2_TLSM ! TLSM Lennard-Jones soft-sphere diameter 
    !(calculated from equation 146)
  real(8):: intermediate_variable
  real(8), parameter:: N_A = 6.0221415e23 ! (molecules/mol) Avogadro's number
  real(8), parameter :: two2onethird = 2d0**(1d0/3d0) , &
                        onethird = 1d0/3d0

  ! Liu et al. 1998 data
  do i=1,n
    ! Silva et al. 1996 correlation
    ! THIS CORRELATION DOES NOT WORK WHEN n-alkane carbon # increases > 30
    !sigma(i) = 1d-8*( 0.17791+11.779*(Tc(i)/Pc(i)*1d5) &
    !                 -0.049029*(Tc(i)/Pc(i)*1d5)**2 )**onethird
    sigma(i) = 1d-8*0.809*Vc(i)**onethird
    eps_k(i) = 0.774*Tc(i)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Water is different because of strong H-bond !!!!!!!!!!!!!!!!!!!!!!!!
    if (dabs(Tc(i)-647).lt.0.5d0 .AND. dabs(Pc(i)-220.6d5).lt.0.5d5) then
      sigma(i) = 1.53091d-8
      eps_k(i) = 3788.51
      !print *, "TLSM gets water", i
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !if (i>50) then
    !  print *, "i,sigma", i, sigma(i)
    !endif
  enddo

  rho_in_cm = rho / 1000.

  ! (cm) ( water , toluene , thiophene ) from Silva et al 1998
  !sigma = (/ 1.53e-8, 5.42e-8, 4.86e-8 /)     
  !sigma = (/ 1.53e-8, 6.71e-8, 7.00e-8 /)

  ! (K) ( water , toluene , thiophene ) from Silva et al 1998
  !epsilon_over_k = (/ 3.789e3, 4.58e2, 4.48e2 /)  
  !epsilon_over_k = (/ 3.789e3, 4.35e2, 6.73e2 /)

  ! ( water , toluene , thiophene )
  !M_mass = (/ 18e0, 92.14e0, 84.14e0 /)  
  !M_mass = (/ 18e0, 142.3e0, 170.3e0 /)
  
  do i = 1, n ! loop over all species as a solute
    Dij(i,i) = 0.0;
    do j = 1, n ! loop over all species as a solvent
      if (i.ne.j) then

        ! prepare xi & xj
        xi = .0001d0 !
        xj = .9999d0 !
 
        !{
        ! Reference: Kurtney PhD thesis "2005":
        M_1_2 = xi*MW(i) + xj*MW(j) ! equation 204
          
        sigma_1_2 = xi*sigma(i) + xj*sigma(j) ! equation 203
          
        eps_k_1_2 = eps_k(i)**(xi) * eps_k(j)**(xj) ! equation 202
          
        ! "sigma_2_TLSM" is calculated in three steps: equation 146
        intermediate_variable = dsqrt(T/eps_k_1_2)
        intermediate_variable = (1.0 + 1.2 * intermediate_variable)**(1.0/3.0)
        sigma_2_TLSM = two2onethird * ( sigma_1_2 ** 2 ) / intermediate_variable
          
        ! "Diff" is calculated in six steps: equation 201
        intermediate_variable = 1.2588*M_1_2 - rho_in_cm*N_A*(dsqrt(sigma_2_TLSM)**3) 
        !intermediate_variable =  1.2588*M2    - rho_in_cm*N_A*(dsqrt(sigma_2_TLSM)**3) 
        intermediate_variable = (-0.75*rho_in_cm*N_A*(dsqrt(sigma_2_TLSM)**3))&
                              / intermediate_variable
        intermediate_variable = intermediate_variable - 0.27862*eps_k_1_2/T
          
        intermediate_variable = dsqrt(R_gas*T/M_1_2) * dexp( intermediate_variable )
        intermediate_variable = 669.138 * M_1_2 * intermediate_variable
        !intermediate_variable  = 669.138 * M2    * intermediate_variable
          
        Dij(i,j) = intermediate_variable / (rho_in_cm*N_A*sigma_2_TLSM) ! unit: cm^2/sec
 
        Dij(i,j) = (1d-4) * Dij(i,j) ! convert to m2/sec
      !}
      endif
    end do
  end do
!}
end subroutine TLSM_diffusion_trace_ij

subroutine TLSM_diffusion_and_Wesselingh_Krishna_model( &
    rho,             & ! density
    T,               & ! temperature
    n,               & ! # of species
    Pc,              & ! vector of Pc (Pa)
    Tc,              & ! vector of Tc (K)
    Vc,              & ! vector of Vc (cm^3/mol)
    MW,              & ! vector of molecular weights
    x,               & ! vector of molar fractions
    Dij)               ! matrix of binary mass diffusivity 
!{
  implicit none
  integer :: i, j
  real(8), intent(in) :: rho ! (kg/m3) density ( position )
  real(8), intent(in) :: T ! (K) Temperature distribution in space
  integer, intent(in) :: n ! number of points in space , number of species
  real(8), intent(in), dimension(n):: Pc! (Pa)
  real(8), intent(in), dimension(n):: Tc! (K) 
  real(8), intent(in), dimension(n):: Vc! (cm^3/mol)
  real(8), intent(in), dimension(n):: MW ! (g/mol) molecular masses
  real(8), intent(in), dimension(n):: x ! Concentration distribution ( position , species )
  real(8), dimension(n,n) :: Dij,&!(cm2/sec)diffusion coefficients and converted to m2/sec
                             Dij_trace

  call TLSM_diffusion_trace_ij(rho,T,n,Pc,Tc,Vc,MW,x,Dij_trace)

  do i=1,n
    do j=1,n
      Dij(i,j) = Dij_trace(i,j)**((1+x(j)-x(i))*0.5) &
               * Dij_trace(j,i)**((1+x(i)-x(j))*0.5) 
    enddo
  enddo
!}
end subroutine TLSM_diffusion_and_Wesselingh_Krishna_model

! Liu et al., Ind. Eng. Chem. Res. 1997, 36, 246-252
subroutine TLSM_diffusion_trace_new( &
    P,               & ! Pressure (Pa)
    T,               & ! temperature (K)
    n,               & ! # of species
    Pc,              & ! vector of Pc (Pa)
    Tc,              & ! vector of Tc (K)
    Vc,              & ! vector of Vc (cm^3/mol)
    w,               & ! vector of acentric factors    
    MW,              & ! vector of molecular weights (g/mol)
    kij,             & ! matrix of BIPs  
    Dij)               ! matrix of binary mass diffusivity
!{
  implicit none
  real(8), intent(in) :: P ! (Pa) pressure
  real(8), intent(in) :: T ! (K) Temperature
  integer, intent(in) :: n ! number of points in space , number of species
  real(8), intent(in), dimension(n):: Tc! (K) 
  real(8), intent(in), dimension(n):: Pc! (Pa)
  real(8), intent(in), dimension(n):: Vc! (cm^3/mol)
  real(8), intent(in), dimension(n):: w!  acentric factor  
  real(8), intent(in), dimension(n):: MW ! (g/mol) molecular masses
  real(8), intent(in), dimension(n,n):: kij! BIPs

  real(8), dimension(n):: sigma ! (cm) molecular diameter, 
  real(8), dimension(n):: sigma_eff ! (cm) molecular diameter, 
    !reference: Kurtney PhD thesis "2005". All are the known values, except water (p. 376)
  real(8), dimension(n):: eps_k ! (K) well-depth petential, 
    !reference: Kurtney PhD thesis "2005". All are the known values, except water (p.376)
  real(8), dimension(n,n) :: Dij          ! (cm2/sec) diffusion coefficients and then converted to m2/sec
                                          ! (solute , in solvent)

  real(8), parameter :: R_gas  = 8.3144621d0   ! (J/mol.K) , Universal gas constant
  real(8) :: x(n), V, &
             rho1, rho1_star, & ! number density of solvent
             T12

  integer::  i, j ! loop counters i:solute_counter, j:solvent_counter

  real(8):: M2    ! molecular mass for solvent
  real(8):: M_1_2 ! mixture molecular mass (calculated from equation 204 )
  real(8):: sigma_1_2 ! Lennard-Jones soft-sphere diameter (calculated from equation 203)
  real(8):: sigma_eff_1_2 
  real(8):: eps_k_1_2 ! TLSM-Lennard-Jones 
    !energy well-depth parameter (calculated from equation 202)
  real(8):: sigma_2_TLSM ! TLSM Lennard-Jones soft-sphere diameter 
    !(calculated from equation 146)
  real(8):: intermediate_variable
  real(8), parameter:: N_A = 6.0221415e23 ! (molecules/mol) Avogadro's number
  real(8), parameter :: two2onethird = 2d0**(1d0/3d0) , &
                        onethird = 1d0/3d0

  ! Liu et al. 1998 data
  do i=1,n
    ! Silva et al. 1996 correlation
    ! THIS CORRELATION DOES NOT WORK WHEN n-alkane carbon # increases > 30
    !sigma(i) = 1d-8*( 0.17791+11.779*(Tc(i)/Pc(i)*1d5) &
    !                 -0.049029*(Tc(i)/Pc(i)*1d5)**2 )**onethird
    sigma(i) = 0.809d-8*Vc(i)**onethird
    eps_k(i) = 0.774d0*Tc(i)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Water is different because of strong H-bond !!!!!!!!!!!!!!!!!!!!!!!!
    if (dabs(Tc(i)-647).lt.0.5d0 .AND. dabs(Pc(i)-220.6d5).lt.0.5d5) then
      sigma(i) = 1.53091d-8
      eps_k(i) = 3788.51d0
      !print *, "TLSM gets water", i
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! equation (7)
    sigma_eff(i) = sigma(i)*2d0**(1d0/6d0) &
                 * (1d0+dsqrt(1.3229d0*T/eps_k(i)))**(-1d0/6d0)

    !if (i>50) then
    !  print *, "i,sigma", i, sigma(i)
    !endif
  enddo

  do i = 1, n ! loop over all species as a solute
    Dij(i,i) = 0.0;
    do j = 1, n ! loop over all species as a solvent
      if (i.ne.j) then

        x = 0d0
        x(i) = 1d-4
        x(j) = .9999d0
        
        call molar_volume2(P,T,n,Pc,Tc,w,kij,x,V)
        
        rho1 = N_A/(V*1d6) ! 1/cm^3

        rho1_star = rho1*sigma_eff(j)**3 ! equation (2)

        !{
        ! Reference: Magalhaes et al. Ind. Eng. Chem. Res. 49, 7697-7700, 2010
        M_1_2 = MW(i)*MW(j)/(MW(i)+MW(j)) ! equation (4)
          
        sigma_1_2 = (sigma(i)+sigma(j))*0.5d0 ! equation (6)

        eps_k_1_2 = dsqrt(eps_k(i)*eps_k(j)*(sigma(i)*sigma(j))**3) &
                  / (sigma_1_2**3)            ! equation (5)

        ! equation (7)
        sigma_eff_1_2 = sigma_1_2*2d0**(1d0/6d0) &
                     * (1d0+dsqrt(1.3229d0*T/eps_k_1_2))**(-1d0/6d0)
          
        T12 = T/eps_k_1_2 ! equation (3)
          
        Dij(i,j) = 21.16d0/(rho1*sigma_eff_1_2*sigma_eff_1_2) &
                 * dsqrt(5d2*R_gas*T/M_1_2) &
                 * dexp(-0.75d0*rho1_star/(1.2588d0-rho1_star)-0.27862d0/T12)
 
        Dij(i,j) = (1d-4) * Dij(i,j) ! convert to m2/sec
      !}
      endif
    end do
  end do
!}
end subroutine TLSM_diffusion_trace_new

subroutine new_TLSM_diffusion_Krishna_model( &
    P,               & ! pressure
    T,               & ! temperature
    n,               & ! # of species
    Pc,              & ! vector of Pc (Pa)
    Tc,              & ! vector of Tc (K)
    Vc,              & ! vector of Vc (cm^3/mol)
    w,               & ! vector of acentric factor    
    MW,              & ! vector of molecular weights
    kij,             & ! matrix of BIPs
    x,               & ! vector of molar fractions
    Dij)               ! matrix of binary mass diffusivity 
!{
  implicit none
  integer :: i, j
  real(8), intent(in) :: P ! (Pa) pressure
  real(8), intent(in) :: T ! (K) Temperature
  integer, intent(in) :: n ! number of points in space , number of species
  real(8), intent(in), dimension(n):: Pc! (Pa)
  real(8), intent(in), dimension(n):: Tc! (K) 
  real(8), intent(in), dimension(n):: Vc! (cm^3/mol)
  real(8), intent(in), dimension(n):: w!  acentric factors  
  real(8), intent(in), dimension(n):: MW ! (g/mol) molecular masses
  real(8), intent(in), dimension(n,n):: kij! BIPs
  real(8), intent(in), dimension(n):: x!  molar fractions
  !(cm2/sec)diffusion coefficients and converted to m2/sec
  real(8), dimension(n,n) :: Dij, Dij_trace

  call TLSM_diffusion_trace_new(P,T,n,Pc,Tc,Vc,w,MW,kij,Dij_trace)

  do i=1,n
    do j=1,n
      Dij(i,j) = Dij_trace(i,j)**((1+x(j)-x(i))*0.5) &
               * Dij_trace(j,i)**((1+x(i)-x(j))*0.5) 
    enddo
  enddo
!}
end subroutine new_TLSM_diffusion_Krishna_model

!***************************************************************************
subroutine findEquilibrium_fix_HC_ratios( &
    P, & ! pressure (Unit: Pa)
    T, & ! temperature (Unit: K)
    n, & ! number of species
    Pc,& ! vector of critical pressures
    Tc,& ! vector of critical temperatures
    w, & ! vector of acentric factors
    type_k, & ! vector of binary interaction types
    x1i, & ! vector of mass fractions: alpha phase initial guess (input)
    x2i, & ! vector of mass fractions: beta phase  initial guess (input)
    x_a, & ! vector of mass fractions: alpha phase
    x_b, & ! vector of mass fractions: beta phase
    n_miscible & ! >0: miscible; <=0: immiscible
    )
!{
  implicit none
  integer :: n, n_miscible
  integer :: type_k(n)
  real(8) :: coef_ab(n), delta(n,n), f1, f2, xa1_old, xb1_old, dx
  real(8) :: P,T,Pc(n),Tc(n),w(n),x_a(n),x_b(n),fuga_a(n),fuga_b(n),x1i(n),x2i(n)
  real(8) :: G_a, G_b
  integer :: i

  ! Gex Method for equilibrium
  integer, parameter :: m=500
  real(8) :: x1(n),x2(n),xo1(n),xo2(n),G(0:m),x(0:m,n),G1,xr,dGdx, tmp_d2Gdx2
  integer :: j,j1,i1,counter

  coef_ab = -1

  dx = 0.0d0
  do i=1,n-1
    x1(i) = x1i(i)
    x2(i) = 0.0
    dx = dx + x1(i)
  end do
  x1(n) = 0.0
  x2(n) = 1.0

  do i=1,n-1
    x1(i) = x1(i)/dx
  end do

  dx = 1d0/dfloat(m)                                                                        !
  do i=0,m                                                                                  !
    xr = dfloat(i)/dfloat(m)                                                                !
    !xr = xr**3                                                                             !
    x(i,:) = xr*x2(:)+(1d0-xr)*x1(:)                                                        !
  enddo                                                                                     !
  do i=0,m                                                                                   !
    call fugacities(P,T,n,Pc,Tc,w,x(i,:),type_k,coef_ab,fuga_a,G(i),delta)                   !
  enddo                                                                                      !
  !open(unit=888,file='func_G.dat')                                                         !
  !do i=0,m                                                                                 !
  !  write(888,*) x(i,n), G(i)                                                              !
  !enddo                                                                                    !
                                                                                            !
  dGdx = 1d10                                                                               !
  !print *, 'd2Gdx2 =', dGdx                                                                 !
  do i=1,m-1                                                                                !
    tmp_d2Gdx2 = (G(i+1)-2d0*G(i)+G(i-1))/dx**2                                             !
    if (dGdx > tmp_d2Gdx2 ) then                                                            !
      dGdx = tmp_d2Gdx2                                                                     !
      i1 = i                                                                                !
    endif                                                                                   !
  enddo                                                                                     !

  if (dGdx .gt. 0) then                                                                    !!
    n_miscible = 1
    return                                                                                 !!
  endif                                                                                    !!
  n_miscible = -1
                                                                                            !
  i = i1-1                                                                                 !!
  j = i1+1                                                                                 !!
                                                                                           !!
  do while ( G(j) + (G(j+1)-G(j)) / (x(j+1,n)-x(j,n)) * (x(i,n)-x(j,n)) >=G(i) .and. j<m)  !!
    j = j+1                                                                                !!
    do while ( G(i) + (G(i-1)-G(i)) / (x(i-1,n)-x(i,n)) * (x(j,n)-x(i,n)) >=G(j) .and. i>1)!!
      i = i-1                                                                              !!
    enddo                                                                                  !!
  enddo                                                                                    !!
                                                                                           !!
  x1=x(i,:)                                                                                !!
  x2=x(j,:)                                                                                !!
  x_a = x1                                                                                 !!
  x_b = x2                                                                                 !!
                                                                                            !
  !call fugacities(P,T,n,Pc,Tc,w,x_a,type_k,coef_ab,fuga_a,G_a,delta)                        !
  !call fugacities(P,T,n,Pc,Tc,w,x_b,type_k,coef_ab,fuga_b,G_b,delta)                        !
  !do i=0,n                                                                                 !!
  !  print*, 'i, xiA*fiA-xiB*fiB:', i, x_a(i)*fuga_a(i)-x_b(i)*fuga_b(i)
  !enddo                                                                                    !!
                                                                                            !
!}
end subroutine findEquilibrium_fix_HC_ratios

subroutine species_LLE( &
    P, &! pressure (Pa)
    T, &! temperature (K)
    n, &! # of species
    Pc,& ! vector of Pc (Pa)
    Tc,& ! vector of Tc (K)
    Vc,& ! vector of Vc (cm^3/mol)
    w, & ! vector of acentric factor
    MW,& ! molecular weight (g/mol)
    type_k, &! type of kij
    xL,&! molar fractions of left  side of interface
    xR,&! molar fractions of right side of interface
    c0, &! original mixing fractions: c0*xL + (1-c0)*xR resulting ==> 0.5*x1 + 0.5*x2
    x1, &! equilibrium molar fractions for left side    ! OUTPUT vector
    x2, &! equilibrium molar fractions for right side   ! OUTPUT vector
    n_miscible &! >0: miscible; <=0: immiscible        ! OUTPUT scalar
    )
!{
  implicit none
  ! INPUT & OUTPUT
  integer :: n, n_miscible
  real(8) :: P,T,Pc(n),Tc(n),Vc(n),w(n),MW(n),xL(n),xR(n),x1(n),x2(n)
  integer :: type_k(n)
  ! local variables
  real(8) :: lnphi1(n), lnphi2(n), & ! ln(fuga_coef)
             phi1(n), phi2(n),     &
             x1_old(n), x2_old(n)    ! values of the last step
  integer :: i, j, jm, x1m, maxLoop=100000
  real(8) :: tol = 1d-7, K(n), A(n), b(n), sumA, sumB
  real(8) :: dx1(n), dx2(n), max_dx, sum_x1M, sum_x2M
  real(8) :: c0
  real(8) :: am, bm ! NOT USED
  real(8) :: nC0

  ! unused local variables for subroutine calls only
  real(8) :: Tb(n), SG(n), H_8(n), rho, CP, CV, CP_IG(n), Hpar(n), Hpur(n), Vpar(n), dVdT, G

  ! setup the initial input
  sumA=0d0                   
  do j=1,n-1                 
    sumA = sumA + xL(j)      
  end do                     
                             
  do j=1,n-1                 
    x1(j) = xL(j) / sumA * (1d0-1d-4)
    x2(j) = 0.0              
  end do                     
  x1(n) = 1d-4
  x2(n) = 1.0                

  do i=1,maxLoop
  !{
    ! save the old values
    x1_old = x1
    x2_old = x2

    ! find max x1: x1m=x1(jm)=max(x1)
    x1m = x1(1)
    jm = 1
    do j=1,n-1
      if (x1m<x1(j)) then
        x1m = x1(j)
        jm = j
      end if
    end do

    !print *, 'jm', jm

    ! calculate initial fugacity
    call thermo_properties(P,T,n,Pc,Tc,w,MW,x1,Tb,SG,H_8,type_k, &
                rho,CP,CV,CP_IG,Hpar,Hpur,Vpar,dVdT,G,lnphi1,am,bm,1)
    call thermo_properties(P,T,n,Pc,Tc,w,MW,x2,Tb,SG,H_8,type_k, &
                rho,CP,CV,CP_IG,Hpar,Hpur,Vpar,dVdT,G,lnphi2,am,bm,1)
    !print *, 'ln(phi)', lnphi1, lnphi2

    phi1 = dexp(lnphi1)
    phi2 = dexp(lnphi2)

    do j=1,n
      K(j) = phi1(j)/phi2(j)
    end do

    do j=1,n
      A(j) = (xL(j) -xR(j) )*(1+K(jm))                  &
           /((xL(jm)-xR(jm))*(1+K(j)))
      b(j) = 2.*(xR(j)*(xL(jm)-xR(jm)) - xR(jm)*(xL(j)-xR(j))) &
           / ((xL(jm)-xR(jm))*(1+K(j)))
    end do

    sumA = 0.0
    sumB = 0.0
    do j=1,n-1
      sumA = sumA + A(j)*(phi1(n)-K(j)*phi2(n))
      sumB = sumB + b(j)*(K(j)-1.)
    end do

    ! x1(jm)
    x1(jm) = (phi1(n)-phi2(n) + sumB) / sumA

    ! x1(j<>jm, n)
    do j=1,n-1
      if (j .ne. jm) x1(j) = A(j)*x1(jm) + b(j)
    end do

    ! x2(j<>n)
    do j=1,n-1
      x2(j) = x1(j)*K(j)
    end do

    ! x1(n), x2(n)
    x1(n) = 1 
    x2(n) = 1 
    do j=1,n-1
      x1(n) = x1(n) - x1(j)
      x2(n) = x2(n) - x2(j)
    end do

    ! negative corrections
    do j=1,n
      if (x1(j) .lt. 0) x1(j) = -x1(j)
      if (x2(j) .lt. 0) x2(j) = -x2(j)
    end do

    ! normalize every species
    sum_x1M = 0.0d0
    sum_x2M = 0.0d0
    do j=1,n
      sum_x1M = sum_x1M + x1(j)
      sum_x2M = sum_x2M + x2(j)
    end do
    !print *, sum_x1M, sum_x2M
    do j=1,n
      x1(j) = x1(j)/sum_x1M
      x2(j) = x2(j)/sum_x2M
    end do

    ! check the step size
    dx1 = dabs(x1 - x1_old)
    dx2 = dabs(x2 - x2_old)
    max_dx = -1.0d0
    do j=1,n
      if (dx1(j)>max_dx) max_dx = dx1(j)
      if (dx2(j)>max_dx) max_dx = dx2(j)
    end do
    !print *, i, x1(1), x1(n), x2(1), x2(n) 

!    print *, 'dx', i,max_dx
!    print *, 'X1', i,x1
!    print *, 'X2', i,x2
!    !print *, 'C1', i,(C-x2)/(x1-x2)
!    print *, 'C', i, (0.5*(x1+x2)-xR)/(xL-xR)

    if (max_dx < tol) then
      !print *, 'it takes', i, 'loops'
      !print *, 'fuga diff', x1*dexp(lnphi1) - x2*dexp(lnphi2)
      exit
    end if

  !}
  end do
  
  ! check oil-rich and water-rich
  if (x1(n)>x2(n)) then
    do j=1,n
      max_dx = x1(j)
      x1(j) = x2(j)
      x2(j) = max_dx
    end do
  end if

  !print *, x1, x2

  c0 = 10d0
  do j=1,n
    nc0 = (0.5*(x1(j)+x2(j))-xR(j))/(xL(j)-xR(j))
    if (nc0<c0) c0=nc0
  end do
  !print *, 'c0', c0
  !print *, (0.5*(x1+x2)-xR)/(xL-xR)

  n_miscible = -1
  if ((dabs(x1(n) - x2(n)) < 1d-3) .OR. isnan(x1(1))) then
    n_miscible = 1
    print *, 'single phase', x1(n), x2(n)
    return
  end if

!}
end subroutine species_LLE

subroutine species_LLE2( &
    P, &! pressure (Pa)
    T, &! temperature (K)
    n, &! # of species
    Pc,& ! vector of Pc (Pa)
    Tc,& ! vector of Tc (K)
    Vc,& ! vector of Vc (cm^3/mol)
    w, & ! vector of acentric factor
    MW,& ! molecular weight (g/mol)
    type_k, &! type of kij
    xL,&! molar fractions of left  side of interface
    xR,&! molar fractions of right side of interface
    c0, &! original mixing fractions: .5*xL+.5*xR resulting ==> c0*x1 + (1-c0)*x2
    x1, &! equilibrium molar fractions for left side    ! OUTPUT vector
    x2, &! equilibrium molar fractions for right side   ! OUTPUT vector
    n_miscible &! >0: miscible; <=0: immiscible        ! OUTPUT scalar
    )
!{
  implicit none
  ! INPUT & OUTPUT
  integer :: n, n_miscible
  real(8) :: P,T,Pc(n),Tc(n),Vc(n),w(n),MW(n),xL(n),xR(n),x1(n),x2(n)
  integer :: type_k(n)
  ! local variables
  real(8) :: lnphi1(n), lnphi2(n), & ! ln(fuga_coef)
             phi1(n), phi2(n),     &
             x1_old(n), x2_old(n)    ! values of the last step
  integer :: i, j, jm, x1m, maxLoop=100000
  real(8) :: tol = 1d-7, K(n), A(n), sumA, sumB, xC(n)
  real(8) :: dx1(n), dx2(n), max_dx, sum_x1M, sum_x2M
  real(8) :: c0
  real(8) :: am, bm ! NOT USED
  real(8) :: nC0

  ! unused local variables for subroutine calls only
  real(8) :: Tb(n), SG(n), H_8(n), rho, CP, CV, CP_IG(n), Hpar(n), Hpur(n), Vpar(n), dVdT, G

  do j=1,n                 
    xC(j) = xL(j)*0.0177017546+xR(j)*0.9822982454
  end do

  ! setup the initial input
  !sumA=0d0                   
  !do j=1,n-1                 
  !  sumA = sumA + xL(j)      
  !end do                     
                             
  !do j=1,n-1                 
  !  x1(j) = xL(j) / sumA * (1d0-3d-1)
  !  x2(j) = xL(j) / sumA * 1d-4
  !end do                     
  !x1(n) = 3d-1
  !x2(n) = 1d0-1d-4

  do i=1,maxLoop
  !{
    ! save the old values
    x1_old = x1
    x2_old = x2

    ! find max x1: x1m=x1(jm)=max(x1)
    x1m = x1(1)
    jm = 1
    do j=1,n-1
      if (x1m<x1(j)) then
        x1m = x1(j)
        jm = j
      end if
    end do

    !print *, 'jm', jm

    ! calculate initial fugacity
    call thermo_properties(P,T,n,Pc,Tc,w,MW,x1,Tb,SG,H_8,type_k, &
                rho,CP,CV,CP_IG,Hpar,Hpur,Vpar,dVdT,G,lnphi1,am,bm,1)
    call thermo_properties(P,T,n,Pc,Tc,w,MW,x2,Tb,SG,H_8,type_k, &
                rho,CP,CV,CP_IG,Hpar,Hpur,Vpar,dVdT,G,lnphi2,am,bm,1)
    !print *, 'ln(phi)', lnphi1, lnphi2

    phi1 = dexp(lnphi1)
    phi2 = dexp(lnphi2)

    do j=1,n
      K(j) = phi1(j)/phi2(j)
    end do

    do j=1,n
      A(j) = (xC(j )-x2(j ))*(1-K(jm))                  &
           /((xC(jm)-x2(jm))*(1-K(j )))
    end do
!K(j )*
!K(jm)*
    sumA = 0.0
    sumB = 0.0
    do j=1,n-1
      sumA = sumA + A(j)
      sumB = sumB + A(j)*K(j)
    end do

    ! x1(jm)
    x1(jm) = (phi1(n)-phi2(n)) / (sumA*phi1(n) - sumB*phi2(n))

!   x2(n) = K(n)*x1(n)
!   x1(jm) = (1-x2(n))/sumB

    ! x1(j<>jm, n)
    do j=1,n-1
      if (j .ne. jm) x1(j) = A(j)*x1(jm) 
    end do

    ! x2(j<>n)
    do j=1,n-1
      x2(j) = x1(j)*K(j)
    end do

    ! x1(n), x2(n)
    x1(n) = 1 
    x2(n) = 1 
    do j=1,n-1
      x1(n) = x1(n) - x1(j)
      x2(n) = x2(n) - x2(j)
    end do

    ! negative corrections
    do j=1,n
      if (x1(j) .lt. 0) x1(j) = -x1(j)
      if (x2(j) .lt. 0) x2(j) = -x2(j)
    end do

    ! normalize every species
    sum_x1M = 0.0d0
    sum_x2M = 0.0d0
    do j=1,n
      sum_x1M = sum_x1M + x1(j)
      sum_x2M = sum_x2M + x2(j)
    end do
    !print *, sum_x1M, sum_x2M
    do j=1,n
      x1(j) = x1(j)/sum_x1M
      x2(j) = x2(j)/sum_x2M
    end do

    ! check the step size
    dx1 = dabs(x1 - x1_old)
    dx2 = dabs(x2 - x2_old)
    max_dx = -1.0d0
    do j=1,n
      if (dx1(j)>max_dx) max_dx = dx1(j)
      if (dx2(j)>max_dx) max_dx = dx2(j)
    end do
    !print *, i, x1(1), x1(n), x2(1), x2(n) 

!    print *, 'dx', i,max_dx
!    print *, 'X1', i,x1
!    print *, 'X2', i,x2
!    !print *, 'C1', i,(C-x2)/(x1-x2)
!    print *, 'C', i, (0.5*(x1+x2)-xR)/(xL-xR)

    if (max_dx < tol) then
      !print *, 'it takes', i, 'loops'
      !print *, 'fuga diff', x1*dexp(lnphi1) - x2*dexp(lnphi2)
      exit
    end if

  !}
  end do
  
  ! check oil-rich and water-rich
  if (x1(n)>x2(n)) then
    do j=1,n
      max_dx = x1(j)
      x1(j) = x2(j)
      x2(j) = max_dx
    end do
  end if

  !print *, x1, x2

  c0 = 10d0
  do j=1,n
    nc0 = (xC(j)-x2(j))/(x1(j)-x2(j))
    if (nc0<c0) c0=nc0
  end do
  print *, 'c0', c0
  print *, (xC-x2)/(x1-x2)

  n_miscible = -1
  if ((dabs(x1(n) - x2(n)) < 1d-3) .OR. isnan(x1(1))) then
    n_miscible = 1
    print *, 'single phase', x1(n), x2(n)
    return
  end if

!}
end subroutine species_LLE2

! service subroutines
! developed according to J. Karl Johnson's course notes
! puccini.che.pitt.edu/~karlj/Classes/CHE2101/Shell/Ch9_Ideal_gases.pdf
! 1) thermal de Broglie wavelength, unit: meter
subroutine thermalWavelength(T, n, MW, tw)
!{
  integer ::i, n
  real(8) ::T, MW(n), tw(n), &
            h = 6.62606957d-34, & ! Planck constant
            pi= 3.1415926535d0, &
            NA= 6.02214129d23,  & ! Avogadro constant
            kB= 1.3806488e-23,  & ! Boltzmann constant
            mi
          
  do i=1,n
    mi = MW(i)*1d-3/NA
    tw(i) = h / dsqrt(2d0*pi*mi*kB*T)
  end do

  !print *, 'thermal de Broglie wavelength (m):', tw
!}
end subroutine thermalWavelength
! 2) ideal gas chemical potential  
subroutine ideal_gas_chem_potential(P, T, n, MW, x, mu)
!{
  integer ::i, n
  real(8) ::P, T, MW(n), x(n), mu(n), tw(n), xi, &
            kB= 1.3806488e-23,  & ! Boltzmann constant
            R_gas = 8.3144621d0,&
            log_kBT

  log_kBT = dlog(kB*T)

  call thermalWavelength(T, n, MW, tw)

  do i=1,n
!    xi = x(i)
!    if (xi .lt. 1e-7) xi = 1e-7
!    mu(i) = R_gas*T*( dlog(xi*P) - log_kBT + 3d0*dlog(tw(i)) )
    mu(i) = R_gas*T*( - log_kBT + 3d0*dlog(tw(i)) )
    !mu(i) = R_gas*T*( dlog(xi) )
  end do
!}
end subroutine ideal_gas_chem_potential
