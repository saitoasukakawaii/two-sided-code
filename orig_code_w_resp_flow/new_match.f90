module new_match
use f90_tools
implicit none

! The number of generations in the structured tree
integer, parameter :: Maxgen = 100

integer :: done(0:Maxgen,0:Maxgen)
complex(lng) :: comp(0:Maxgen,0:Maxgen,1:2,1:2)

integer alphamax, betamax, count, endbranch

contains

recursive function admit (n,m,omega_k,fa1,fa2,fa3,fv1,fv2,fv3,rho,mu,r_root,r_min,Lr,q,g,asym,expo,lrr_A,lrr_V) result (Y)

    integer, intent(in)   :: n,m  ! n is alpha pow and m is beta pow
    real(lng), intent(in) :: omega_k,fa1,fa2,fa3,fv1,fv2,fv3,rho,mu,r_root,r_min,Lr,q,g,asym,expo,lrr_A,lrr_V

    complex(lng) :: Y(1:2,1:2), YA(1:2,1:2), YV(1:2,1:2), Ymid(1:2,1:2), d1(1:2,1:2), d2(1:2,1:2), i
    complex(lng) :: kappa_A, kappa_V, Fk_A, Fk_V, g_omega_A, c_omega_A, g_omega_V, c_omega_V
    real(lng)    :: alpha, beta, nu, r, r_d, l_A, l_V, A, A_d, D_A, D_V, wom, F0_A, F0_V
    integer :: j,k

    count = count + 1

    do j = 1, 2
      do k = 1, 2
        Y(j,k)    = 0.0
        YA(j,k)   = 0.0
        YV(j,k)   = 0.0
        Ymid(j,k) = 0.0
        d1(j,k)   = 0.0
        d2(j,k)   = 0.0
      end do
    end do
	
    ! Physical constants.
    i     = cmplx(0.0_lng,1.0_lng,lng) ! The complex unit.

    alpha = ((asym**(expo/2)+1.0)**(-1/expo))**n  ! Scaling parameter.
    beta  = (sqrt(asym)*(asym**(expo/2)+1.0)**(-1/expo))**m   ! do.

    if(n >= alphamax) then
      alphamax = n
    endif

    if(m >= betamax) then
      betamax = m
    endif
    
    r_d  = alpha*beta*r_root     ! Radius at root, cm.
    A_d  = pi*r_d**2             ! Cross-sectional area, cm^2.
    r    = r_d                   ! Radius at root, dimension-less.
    A    = A_d                   ! Cross-sectional area, dimension-less.

!    if(r<0.005) then
!      l_A  = 1.88*(r**0.47)
!    else
!      l_A  = lrr_A*r**(1.10)     ! Length of arterial vessel segment.
!    end if

    l_A  = lrr_A*r**(1.10) !mpaun added in

    l_V  = lrr_V*r               ! Length of venous vessel segment.
    nu   = mu/rho                ! Kinematic blood viscosity,cm^2/s.
    D_A  = 1/(fa1*exp(fa2*r_d)+fa3)*3*A_d/2 ! Distensibility.
    D_V  = 1/(fv1*exp(fv2*r_d)+fv3)*3*A_d/2
    wom  = r_d*sqrt(omega_k/nu)  ! Womersleys parameter.

    
    ! Temporary functions of r.
     if (wom > 3.0) then 
 g_omega_A = sqrt(D_A*A/rho)* &
 sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom)) 
c_omega_A = sqrt(A/D_A/rho)* &
              sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom))
                  
  g_omega_V = sqrt(D_V*A/rho)* &
              sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom)) 
  c_omega_V = sqrt(A/D_V/rho)* &
              sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom))
      else
        if (wom > 2.0) then
          g_omega_A = sqrt(D_A*A/rho)*((3.0_lng-wom)* &
                      sqrt(i*wom**2.0/8.0+wom**4.0/48.0) + &
                      (wom-2.0_lng)*&
                      sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom)))
          c_omega_A = sqrt(A/D_A/rho)*((3.0_lng-wom)* &
                      sqrt(i*wom**2.0/8.0+wom**4.0/48.0) + &
                      (wom-2.0_lng)*&
                      sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom)))
                      
          g_omega_V = sqrt(D_V*A/rho)*((3.0_lng-wom)* &
                      sqrt(i*wom**2.0/8.0+wom**4.0/48.0) + &
                      (wom-2.0_lng)*&
                      sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom)))
          c_omega_V = sqrt(A/D_V/rho)*((3.0_lng-wom)* &
                      sqrt(i*wom**2.0/8.0+wom**4.0/48.0) + &
                      (wom-2.0_lng)*&
                      sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom)))
        else
        if (wom == 0) then
          g_omega_A = 0.0
          c_omega_A = 0.0
          g_omega_V = 0.0
          c_omega_V = 0.0
        else
          g_omega_A = sqrt(D_A*A/rho)*sqrt(i*wom**2/8+wom**4/48)
          c_omega_A = sqrt(A/D_A/rho)*sqrt(i*wom**2/8+wom**4/48)
          g_omega_V = sqrt(D_V*A/rho)*sqrt(i*wom**2/8+wom**4/48)
          c_omega_V = sqrt(A/D_V/rho)*sqrt(i*wom**2/8+wom**4/48)
         end if
      end if
    end if
 
    ! Temporary function of omega_k. 
    if (omega_k /= 0) then
      kappa_A  = omega_k*l_A/c_omega_A
      kappa_V  = omega_k*l_V/c_omega_V
else
      kappa_A   = 0.0
      kappa_V   = 0.0
    end if
  
    if (kappa_A == 0) then
      F0_A    = (pi*rho*g*Lr*(r**4))/(8*mu*l_A*q)
!     F0_A    = pi*(r**4)/(8*mu*l_A)
      YA(1,1) =  F0_A
      YA(2,1) = -F0_A
      YA(1,2) = -F0_A
      YA(2,2) =  F0_A
    else
      Fk_A    = i*g_omega_A*rho*g*Lr/(q*sin(kappa_A))
!     Fk_A    = i*g_omega/(q*sin(kappa_A))
      YA(1,1) = -Fk_A*cos(kappa_A)
      YA(2,1) =  Fk_A
      YA(1,2) =  Fk_A
      YA(2,2) = -Fk_A*cos(kappa_A)
    end if

    if (kappa_V == 0) then
      F0_V    = (pi*rho*g*Lr*(r**4))/(8*mu*l_V*q)
!     F0_V    = pi*(r**4)/(8*mu*l_V)
      YV(1,1) =  F0_V
      YV(2,1) = -F0_V
      YV(1,2) = -F0_V
      YV(2,2) =  F0_V
    else
      Fk_V    = i*g_omega_V*rho*g*Lr/(q*sin(kappa_V))
!     Fk_V    = i*g_omega/(q*sin(kappa_V))
      YV(1,1) = -Fk_V*cos(kappa_V)
      YV(2,1) =  Fk_V
      YV(1,2) =  Fk_V
      YV(2,2) = -Fk_V*cos(kappa_V)
    end if
    
    if (r < r_min) then
      endbranch = endbranch + 1
      Y = series(YA,YV)
    else
      if (done(n+1,m) == 1) then
        d1 = comp(n+1,m,:,:)
      else
        d1 = admit(n+1,m,omega_k,fa1,fa2,fa3,fv1,fv2,fv3,rho,mu,r_root,r_min,Lr,q,g,asym,expo,lrr_A,lrr_V)
      end if

      if (done(n,m+1) == 1) then
        d2 = comp(n,m+1,:,:)
      else
        d2 = admit(n,m+1,omega_k,fa1,fa2,fa3,fv1,fv2,fv3,rho,mu,r_root,r_min,Lr,q,g,asym,expo,lrr_A,lrr_V)
      end if
      Ymid = d1 + d2
      Y = series(series(YA,Ymid),YV)
    end if

    done(n,m) = 1
    comp(n,m,:,:) = Y   
end function admit


function series (YA, YB) result(Y)

    complex(lng), intent(in) ::  YA(1:2,1:2), YB(1:2,1:2)
    complex(lng) :: Y(1:2,1:2), DA, DB
	
    DA = YA(1,1)*YA(2,2) - YA(1,2)*YA(2,1)
    DB = YB(1,1)*YB(2,2) - YB(1,2)*YB(2,1)

    Y(1,1) = (DA+YA(1,1)*YB(1,1))/(YA(2,2)+YB(1,1))
    Y(1,2) = (  -YA(1,2)*YB(1,2))/(YA(2,2)+YB(1,1))
    Y(2,1) = (  -YA(2,1)*YB(2,1))/(YA(2,2)+YB(1,1))
    Y(2,2) = (DB+YA(2,2)*YB(2,2))/(YA(2,2)+YB(1,1))

end function series


subroutine impedance (tmstps,Period,rho,mu,r_root,r_min,y11,y12,y21,y22,Lr,q,g,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrr_A,lrr_V)
implicit none

  integer,   intent(in)      :: tmstps
  real(lng), intent(in)      :: Period,rho,mu,Lr,q,g,r_root,r_min,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrr_A,lrr_v
  real(lng)                  :: y11(tmstps),y12(tmstps),y21(tmstps),y22(tmstps)!,a11(tmstps),a12(tmstps),a21(tmstps),a22(tmstps)

  integer          :: j,k
  real(lng)        :: df
  real(lng)        :: Freq(tmstps+1),Omega(tmstps+1)
  complex(lng)     :: Y(1:tmstps,1:2,1:2),tempY(1:2,1:2),YHAT11(tmstps),YHAT12(tmstps),YHAT21(tmstps),YHAT22(tmstps)

!write(*,*)'impedance is called' !mpaun commented out

! Physical parameters
df     = 1/Period                            ! Frequency interval. 
Freq   = (/ (j*df, j=-tmstps/2, tmstps/2) /) ! Frequency-vector (abscissae). 
Omega  = 2*pi*Freq                           ! Freq.-vector scaled by a factor 2pi.

alphamax = 0
betamax  = 0
count = 0
endbranch = 0

do k = tmstps/2+1, tmstps+1
  done = 0
  comp = 0.0_lng
  Y(k-1,:,:) = admit(0,0,Omega(k),fa1,fa2,fa3,fv1,fv2,fv3,rho,mu,r_root,r_min,Lr,q,g,asym,expo,lrr_A,lrr_V)
end do
!write(*,*) 'r_root - ', r_root
!write(*,*) 'Generations (alpha branch) - ', alphamax !mpaun commented out
!write(*,*) 'Generations (beta branch)  - ', betamax
!write(*,*) 'Total vessels  - ', count
!write(*,*) 'End branches   - ', endbranch


tempY(1,1) = Y(tmstps/2,1,1)
tempY(1,2) = Y(tmstps/2,1,2)
tempY(2,1) = Y(tmstps/2,2,1)
tempY(2,2) = Y(tmstps/2,2,2)

Y(1:tmstps/2,1,1) = conjg(flipud(Y(tmstps/2+1:tmstps,1,1)))
Y(1:tmstps/2,1,2) = conjg(flipud(Y(tmstps/2+1:tmstps,1,2)))
Y(1:tmstps/2,2,1) = conjg(flipud(Y(tmstps/2+1:tmstps,2,1)))
Y(1:tmstps/2,2,2) = conjg(flipud(Y(tmstps/2+1:tmstps,2,2)))

Y(tmstps/2+1:tmstps,1,1) = eoshift(Y(tmstps/2+1:tmstps,1,1),-1)
Y(tmstps/2+1:tmstps,1,2) = eoshift(Y(tmstps/2+1:tmstps,1,2),-1)
Y(tmstps/2+1:tmstps,2,1) = eoshift(Y(tmstps/2+1:tmstps,2,1),-1)
Y(tmstps/2+1:tmstps,2,2) = eoshift(Y(tmstps/2+1:tmstps,2,2),-1)

Y(tmstps/2+1,1,1) = tempY(1,1)
Y(tmstps/2+1,1,2) = tempY(1,2)
Y(tmstps/2+1,2,1) = tempY(2,1)
Y(tmstps/2+1,2,2) = tempY(2,2)

YHAT11 = Y(:,1,1)
YHAT12 = Y(:,1,2)
YHAT21 = Y(:,2,1)
YHAT22 = Y(:,2,2)
y11 = real(IFFT(bitreverse(FFTshift(YHAT11/Period))),lng)
y12 = real(IFFT(bitreverse(FFTshift(YHAT12/Period))),lng)
y21 = real(IFFT(bitreverse(FFTshift(YHAT21/Period))),lng)
y22 = real(IFFT(bitreverse(FFTshift(YHAT22/Period))),lng)

!a11 = aimag(YHAT11(k)) ! components of grand admittance matrix
!a12 = aimag(YHAT12(k))
!a21 = aimag(YHAT21(k))
!a22 = aimag(YHAT22(k))


open(unit=1, file='Admit',status='replace')
do k = 1, tmstps
write(1,12) Omega(k)/Lr**3*q, y11(k)*rho*g*Lr/q, y12(k)*rho*g*Lr/q, y21(k)*rho*g*Lr/q, y22(k)*rho*g*Lr/q
12 format(F26.16,2x,F26.16,2x,F26.16,2x,F26.16,2x,F26.16)
end do
close(unit=1)

!open(unit=1, file='Admit.dat',status='replace')
!do k = 1, tmstps
!write(1,13) Omega(k)/Lr**3*q, a11(11), a11(12), a11(21), a11(22)
!13 format(F26.16,2x,F26.16,2x,F26.16,2x,F26.16,2x,F26.16)
!end do
!close(unit=1)


end subroutine impedance

end module new_match
