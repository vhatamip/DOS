!  MeshGen.f90 
!
!  FUNCTIONS:
!  MeshGen - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: MeshGen
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    
MODULE mesh

IMPLICIT NONE


CONTAINS
!*********************************************************************
SUBROUTINE cluster_mesh_left(min_v,max_v,nodes,meshed)

IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: min_v, max_v
INTEGER, INTENT(IN) :: nodes
REAL(KIND=8), DIMENSION(nodes), INTENT(OUT) :: meshed

!local
REAL(KIND=8) :: L ! is the length of the interval
REAL(KIND=8), PARAMETER :: b=1.0001 !closer to the value "one" the more clustered the mesh gets
REAL(KIND=8), DIMENSION(nodes) :: k,a, temp
INTEGER :: i

L=max_v-min_v

DO i=1,nodes
    k(i)=0+(i-1)*(1/(REAL(nodes)-1))
	a(i)=((b+1)/(b-1))**(1-k(i))
    temp(i)=L*((b+1)-(b-1)*a(i))/(1+a(i))
	meshed(i)=temp(i) + min_v
END DO

END SUBROUTINE cluster_mesh_left

!*********************************************************************

SUBROUTINE cluster_mesh_right(min_v,max_v,nodes,meshed)

IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: min_v, max_v
INTEGER, INTENT(IN) :: nodes
REAL(KIND=8), DIMENSION(nodes), INTENT(OUT) :: meshed

!local
REAL(KIND=8) :: L
REAL(KIND=8), PARAMETER :: b=1.0001, alph=0.0
REAL(KIND=8), DIMENSION(nodes) :: k, a, temp
INTEGER :: i

L=max_v-min_v

DO i=1,nodes
    k(i)=((i-1)/(REAL(nodes)-1)) ! 0 to 1 interval is uniformly meshed!
	a(i)=((b+1)/(b-1))**((k(i)-alph)/(1-alph));
    temp(i)=L*((b+2*alph)*a(i)-(b-2*alph))/((2*alph+1)*(1+a(i)))
	meshed(i)=temp(i) + min_v
END DO

END SUBROUTINE cluster_mesh_right

!*********************************************************************

SUBROUTINE cluster_mesh_sides(min_v,max_v,nodes,meshed)

IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: min_v, max_v
INTEGER, INTENT(IN) :: nodes
REAL(KIND=8), DIMENSION(nodes), INTENT(OUT) :: meshed

!local
REAL(KIND=8) :: L
REAL(KIND=8), PARAMETER :: b=1.0001, alpha=0.5
REAL(KIND=8), DIMENSION(nodes) :: k, a, temp
INTEGER :: i

L=max_v-min_v

DO i=1,nodes
    k(i)=0+(i-1)*(1/(REAL(nodes)-1))
	a(i)=((b+1)/(b-1))**((k(i)-alpha)/(1-alpha))
    temp(i)=L*((b+2*alpha)*a(i)-(b-2*alpha))/((2*alpha+1)*(1+a(i)))
	meshed(i)=temp(i) + min_v
END DO

END SUBROUTINE cluster_mesh_sides

!*********************************************************************

SUBROUTINE cluster_mesh_mid(min_v,max_v,nodes,meshed)

IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: min_v, max_v
INTEGER, INTENT(IN) :: nodes
REAL(KIND=8), DIMENSION(nodes), INTENT(OUT) :: meshed

!local
REAL(KIND=8) :: L, D
REAL(KIND=8), PARAMETER :: b=15.0
REAL(KIND=8), DIMENSION(nodes) :: k, a, temp
INTEGER :: i

L=max_v-min_v
D=0.5*L ! clustered mesh in the middle of interval

DO i=1,nodes
    k(i)=0+(i-1)*(1/(REAL(nodes)-1))
	a(i)=(0.5/b)*DLOG((1+(EXP(b)-1)*(D/L))/(1+(EXP(-b)-1)*(D/L)));
    temp(i)=D*(1+(DSINH(b*(k(i)-a(i)))/DSINH(b*a(i))));
    meshed(i)=temp(i) + min_v
END DO

END SUBROUTINE cluster_mesh_mid

!*********************************************************************
SUBROUTINE sort_mesh(meshed1, meshed2, meshed3, sorted) !right, sides, left

IMPLICIT NONE
!nodes_total=nodes_p1+nodes_p2+
REAL(KIND=8), DIMENSION(:), INTENT(IN) :: meshed1, meshed2, meshed3
REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: sorted

!local
!INTEGER, DIMENSION(1) :: n1, n2, n3, n
INTEGER :: n1, n2, n3, n, i
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: temp

n1=UBOUND(meshed1,1) ! number of nodes in the first part 
n2=UBOUND(meshed2,1) ! number of nodes in the middle part
n3=UBOUND(meshed3,1) ! number of nodes in the last part

n=n1+n2+n3-2
ALLOCATE(temp(n)) 

DO i=1,n1
	temp(i)=meshed1(i)
END DO
DO i=1,n2-1
	temp(n1+i)=meshed2(i+1)
END DO
DO i=1,n3-1
	temp(n1+n2-1+i)=meshed3(i+1)
END DO

!ALLOCATE(sorted(n))
sorted=temp

END SUBROUTINE sort_mesh

!****************************************************************************************
!****************************************************************************************

SUBROUTINE Fresnel(j, mu, eps, k_z, r_TM, r_TE)

IMPLICIT NONE

INTEGER, INTENT(IN) :: j
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: k_z
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(OUT) :: r_TM, r_TE
COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: mu, eps

INTEGER :: l, n

n=UBOUND(k_z,2)
                
DO l=1, n
            r_TM(1,l)=-(eps(1) * k_z(2,l) - eps(2)*k_z(1,l))/(eps(1) * k_z(2,l) + eps(2) * k_z(1,l))
	        r_TE(1,l)=-(mu(1)  * k_z(2,l) - mu(2) *k_z(1,l))/(mu(1)  * k_z(2,l) + mu(2)  * k_z(1,l))
			r_TM(2,l)=-(eps(1) * k_z(3,l) - eps(3)*k_z(1,l))/(eps(1) * k_z(3,l) + eps(3) * k_z(1,l))
	        r_TE(2,l)=-(mu(1)  * k_z(3,l) - mu(3) *k_z(1,l))/(mu(1)  * k_z(3,l) + mu(3)  * k_z(1,l))
END DO

END SUBROUTINE Fresnel

SUBROUTINE t_12(r_TM, r_TE, delta, distance, k_z, k_rho, t_TE, t_TM) !t_1(1)-->TE, t_1(2)-->TM

REAL(KIND=8), INTENT(IN) :: delta, distance
REAL(KIND=8), DIMENSION(:), INTENT(IN) :: k_rho
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: k_z
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: r_TM, r_TE
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(OUT) :: t_TE, t_TM !1-->t1, 2-->t2

INTEGER :: l, n

n=UBOUND(k_rho, 1)

DO l=1, n
	      t_TE(1,l)=(((1.0-r_TE(1,l))*EXP(CMPLX(0.0,1.0)*k_z(1,l)*delta))/(1-r_TE(1,l)*r_TE(2,l)*EXP(CMPLX(0.0,1.0)*2*k_z(1,l)*distance)))&
		  -(((1.0-r_TE(1,l))*r_TE(2,l)*EXP(CMPLX(0.0,1.0)*2*k_z(1,l)*(distance-delta)))/(1-r_TE(1,l)*r_TE(2,l)*EXP(CMPLX(0.0,1.0)*2*k_z(1,l)*distance)))
		  t_TE(2,l)=(((1.0-r_TE(1,l))*EXP(CMPLX(0.0,1.0)*k_z(1,l)*delta))/(1-r_TE(1,l)*r_TE(2,l)*EXP(CMPLX(0.0,1.0)*2*k_z(1,l)*distance)))&
		  +(((1.0-r_TE(1,l))*r_TE(2,l)*EXP(CMPLX(0.0,1.0)*2*k_z(1,l)*(distance-delta)))/(1-r_TE(1,l)*r_TE(2,l)*EXP(CMPLX(0.0,1.0)*2*k_z(1,l)*distance)))
		  t_TM(1,l)=(((1.0-r_TM(1,l))*EXP(CMPLX(0.0,1.0)*k_z(1,l)*delta))/(1-r_TM(1,l)*r_TM(2,l)*EXP(CMPLX(0.0,1.0)*2*k_z(1,l)*distance)))&
		  -(((1.0-r_TM(1,l))*r_TM(2,l)*EXP(CMPLX(0.0,1.0)*2*k_z(1,l)*(distance-delta)))/(1-r_TM(1,l)*r_TM(2,l)*EXP(CMPLX(0.0,1.0)*2*k_z(1,l)*distance)))
		  t_TM(2,l)=(((1.0-r_TM(1,l))*EXP(CMPLX(0.0,1.0)*k_z(1,l)*delta))/(1-r_TM(1,l)*r_TM(2,l)*EXP(CMPLX(0.0,1.0)*2*k_z(1,l)*distance)))&
		  +(((1.0-r_TM(1,l))*r_TM(2,l)*EXP(CMPLX(0.0,1.0)*2*k_z(1,l)*(distance-delta)))/(1-r_TM(1,l)*r_TM(2,l)*EXP(CMPLX(0.0,1.0)*2*k_z(1,l)*distance)))
END DO
	  
END SUBROUTINE t_12


SUBROUTINE wave_vector_z(j, mu, eps, k, k_rho, k_z) ! jth omega

IMPLICIT NONE

INTEGER, INTENT(IN) :: j
REAL(KIND=8), DIMENSION(:), INTENT(IN) :: k_rho
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(OUT) :: k_z
COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: k
COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: mu, eps

INTEGER :: l, q, n, m
!COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: k_z0, k_z1, r_TM_10, r_TE_10

n=UBOUND(k_rho,1) !number of nodes
m=UBOUND(eps,1) !number of media
!ALLOCATE(k_z(m,n))

DO l=1, n
	DO q=1,m
			k_z(q,l)=SQRT(mu(q)*eps(q)*k(1)**2-k_rho(l)**2)
	END DO		
END DO

END SUBROUTINE wave_vector_z

SUBROUTINE k_mag(j, omega, mu, eps, k)

IMPLICIT NONE

INTEGER, INTENT(IN) :: j
REAL(KIND=8), INTENT(IN) :: omega
COMPLEX(KIND=8), DIMENSION(:), INTENT(OUT) :: k
COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: mu, eps

INTEGER :: l,n
REAL(KIND=8), PARAMETER :: c0=299792458
!COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: k_z0, k_z1, r_TM_10, r_TE_10
n=UBOUND(eps,1)

DO l=1, n
        k(l)=SQRT(mu(l)*eps(l))*omega/c0    
END DO


END SUBROUTINE k_mag	


SUBROUTINE eps_Basu(j, omega, eps) 

IMPLICIT NONE

INTEGER :: j
REAL(KIND=8), INTENT(IN) :: omega
COMPLEX(KIND=8), DIMENSION(:), INTENT(OUT) :: eps

!INTEGER :: l, k
REAL(KIND=8) :: eps_bl, n_e, e, m_0, m_star, eps_0, mu_1, mu_max, mu_2, c_r, c_s, alpha
REAL(KIND=8) :: betha, T_1, T_2, mu_mob_300, omega_p, mu_mob_T_1, tau_1, mu_mob_T_2, tau_2     

eps_bl=11.7;
n_e=1E20*1E6; !(1/m**3)
e=1.60217662E-19; ! Basu's paper "infrared..."
m_0=9.109E-31;! Basu's paper "infrared..."
m_star=0.27*m_0;
eps_0=8.854187817E-12; ! Basu's paper "infrared..."
mu_1=68.5;
mu_max=1414;
mu_2=56.1;
c_r=9.2E16*1E6; !(1/m**3)
c_s=3.41E20*1E6; !(1/m**3)
alpha=0.711;
betha=1.98;
T_1=400; !Kelvin & emitter
T_2=0.0; !Kelvin & reciever

mu_mob_300=1E-4*(mu_1+((mu_max-mu_1)/(1+(n_e/c_r)**alpha))-(mu_2/(1+(c_s/n_e)**betha))); !m**2/(V.s)
omega_p=SQRT((n_e*e**2)/(m_star*eps_0));

mu_mob_T_1=mu_mob_300*(T_1/300)**1.5; !mob --> mobility
tau_1=(m_star*mu_mob_T_1)/e;


mu_mob_T_2=mu_mob_300*(T_2/300)**1.5; 
tau_2=(m_star*mu_mob_T_2)/e;

IF (tau_1==0.0) THEN 
	eps(1)=DCMPLX(1.0,0.0)
    eps(2)=eps_bl
    eps(3)=eps_bl-(omega_p**2)/(omega*(omega+(DCMPLX(0.0,1.0)/tau_2)));
ELSE IF (tau_2==0.0) THEN
	eps(1)=DCMPLX(1.0,0.0)
    eps(2)=eps_bl-((omega_p**2)/(omega*(omega+(DCMPLX(0.0,1.0)/tau_1))));
    eps(3)=eps_bl
END IF

!eps(3)=eps_inf*(DCMPLX(omega(j)**2-omega_lo**2 , gamma*omega(j))/DCMPLX(omega(j)**2-omega_to**2 , gamma*omega(j)))

END SUBROUTINE eps_Basu

SUBROUTINE eps_Lorentz(j, omega, eps)

IMPLICIT NONE

INTEGER :: j
REAL(KIND=8), INTENT(IN) :: omega
COMPLEX(KIND=8), DIMENSION(:), INTENT(OUT) :: eps

REAL(KIND=8), PARAMETER :: eps_inf=6.7, gamma=8.966E11	!Equation 4, 2009, ()(s^-1)
REAL(KIND=8), PARAMETER :: omega_lo=1.825E14, omega_to=1.494E14  !Equation 4, 2009, (rad s^-1)(rad s^-1)

eps(1)=DCMPLX(1.0,0.0)
eps(2)=eps_inf*(DCMPLX(omega**2-omega_lo**2 , gamma*omega)/DCMPLX(omega**2-omega_to**2 , gamma*omega))
eps(3)=eps(2)

END SUBROUTINE eps_Lorentz

SUBROUTINE mu_mag(omega, mu) ! Magnetic medium

IMPLICIT NONE

REAL(KIND=8), INTENT(IN) :: omega
COMPLEX(KIND=8), DIMENSION(:), INTENT(OUT) :: mu

mu(1)=DCMPLX(1.0,0.0)
mu(2)=DCMPLX(1.0,0.0)
mu(3)=DCMPLX(1.0,0.0)

END SUBROUTINE mu_mag

SUBROUTINE mu_non_mag(j, omega, mu) !Non-Magnetic medium

IMPLICIT NONE

INTEGER :: j
REAL(KIND=8), INTENT(IN) :: omega
COMPLEX(KIND=8), DIMENSION(:), INTENT(OUT) :: mu

mu(1)=DCMPLX(1.0,0.0)
mu(2)=DCMPLX(1.0,0.0)
mu(3)=DCMPLX(1.0,0.0)

END SUBROUTINE mu_non_mag

!****************************************************************************************
!****************************************************************************************

SUBROUTINE integrand_k_rho_e(omega, eps_2, k, k_rho, k_z, t_TE, t_TM, f_k_rho_e)

REAL(KIND=8), INTENT(IN) :: omega
REAL(KIND=8), DIMENSION(:), INTENT(IN) :: k_rho
REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f_k_rho_e
COMPLEX(KIND=8), INTENT(IN) :: eps_2
COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: k
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: k_z
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: t_TE, t_TM

INTEGER :: l, n
REAL(KIND=8), PARAMETER :: c0=299792458
REAL(KIND=8), PARAMETER :: pi=4*ATAN(1.0) !pi=3.14...

n=UBOUND(k_rho,1)

DO l=1, n
	      f_k_rho_e(l)=(((ABS(k_z(1,l)*t_TM(1,l)))**2+(ABS(k_rho(l)*t_TM(2,l)))**2)*(((k_rho(l))**2+k_z(2,l)*CONJG(k_z(2,l)))&
		  /(k(2)*CONJG(k(2))))+(ABS(k(1)*t_TE(2,l)))**2)*((k_rho(l))/(((ABS(k_z(2,l)))**2)*DIMAG(k_z(2,l))))
END DO

f_k_rho_e=((omega*DIMAG(eps_2))/(16*(pi**2)*c0**2))*f_k_rho_e

END SUBROUTINE integrand_k_rho_e

SUBROUTINE integrand_k_rho_m(omega, eps_2, k, k_rho, k_z, t_TE, t_TM, f_k_rho_m)

REAL(KIND=8), INTENT(IN) :: omega
REAL(KIND=8), DIMENSION(:), INTENT(IN) :: k_rho
REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f_k_rho_m
COMPLEX(KIND=8), INTENT(IN) :: eps_2
COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: k
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: k_z
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: t_TE, t_TM

INTEGER :: l, n
REAL(KIND=8), PARAMETER :: c0=299792458
REAL(KIND=8), PARAMETER :: pi=4*ATAN(1.0) !pi=3.14...

n=UBOUND(k_rho,1)

DO l=1, n
	      f_k_rho_m(l)=(((ABS(k(1)*t_TM(2,l)))**2)*(((k_rho(l))**2+k_z(2,l)*CONJG(k_z(2,l)))/(k(2)*CONJG(k(2))))+(ABS(k_z(3,l)&
		  *t_TE(1,l)))**2+(ABS(k_rho(l)*t_TE(2,l)))**2)*((k_rho(l))/(((ABS(k_z(2,l)))**2)*DIMAG(k_z(2,l))))
END DO

f_k_rho_m=((omega*DIMAG(eps_2))/(16*(pi**2)*c0**2))*f_k_rho_m

END SUBROUTINE integrand_k_rho_m

SUBROUTINE integrand_k_rho_book_e(omega, eps_2, k, k_rho, k_z, t_TE, t_TM, f_k_rho_e)

REAL(KIND=8), INTENT(IN) :: omega
REAL(KIND=8), DIMENSION(:), INTENT(IN) :: k_rho
REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f_k_rho_e
COMPLEX(KIND=8), INTENT(IN) :: eps_2
COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: k
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: k_z
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: t_TE, t_TM

INTEGER :: l, n
REAL(KIND=8), PARAMETER :: c0=299792458
REAL(KIND=8), PARAMETER :: pi=4*ATAN(1.0) !pi=3.14...

n=UBOUND(k_rho,1)

DO l=1, n
	      f_k_rho_e(l)=(((ABS(k_z(1,l)*t_TM(1,l)))**2+(ABS(k_rho(l)*t_TM(2,l)))**2)*(((k_rho(l))**2+k_z(2,l)*CONJG(k_z(2,l)))&
		  /(k(2)*CONJG(k(2))))+(ABS(k(1)*t_TE(2,l)))**2)*((k_rho(l))/(((ABS(k_z(2,l)))**2)*DIMAG(k_z(2,l))))
END DO

f_k_rho_e=((omega*DIMAG(eps_2))/(16*(pi**2)*c0**2))*f_k_rho_e

END SUBROUTINE integrand_k_rho_book_e

SUBROUTINE integrand_k_rho_book_m(omega, eps_2, k, k_rho, k_z, t_TE, t_TM, f_k_rho_m)

REAL(KIND=8), INTENT(IN) :: omega
REAL(KIND=8), DIMENSION(:), INTENT(IN) :: k_rho
REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f_k_rho_m
COMPLEX(KIND=8), INTENT(IN) :: eps_2
COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: k
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: k_z
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: t_TE, t_TM

INTEGER :: l, n
REAL(KIND=8), PARAMETER :: c0=299792458
REAL(KIND=8), PARAMETER :: pi=4*ATAN(1.0) !pi=3.14...

n=UBOUND(k_rho,1)

DO l=1, n
	      f_k_rho_m(l)=(((ABS(k(1)*t_TM(2,l)))**2)*(((k_rho(l))**2+k_z(2,l)*CONJG(k_z(2,l)))/(k(2)*CONJG(k(2))))+(ABS(k_z(3,l)&
		  *t_TE(1,l)))**2+(ABS(k_rho(l)*t_TE(2,l)))**2)*((k_rho(l))/(((ABS(k_z(2,l)))**2)*DIMAG(k_z(2,l))))
END DO

f_k_rho_m=((omega*DIMAG(eps_2))/(16*(pi**2)*c0**2))*f_k_rho_m

END SUBROUTINE integrand_k_rho_book_m


SUBROUTINE integrand_k_rho_prop(omega, k, k_z, k_rho, r_TM, r_TE, f_k_rho)

REAL(KIND=8), INTENT(IN) :: omega
REAL(KIND=8), DIMENSION(:), INTENT(IN) :: k_rho
REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f_k_rho
COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: k
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: k_z
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: r_TM, r_TE


INTEGER :: l, n
REAL(KIND=8), PARAMETER :: c0=299792458
REAL(KIND=8), PARAMETER :: pi=4*ATAN(1.0) !pi=3.14...

n=UBOUND(k_rho,1)

DO l=1, n
f_k_rho(l) = (k_rho(l)/ABS(k_z(1,l)))*((1-ABS(r_TE(1,l))**2)+(1-ABS(r_TM(1,l))**2))
END DO

f_k_rho=((k(1)**2)/(omega*4*pi**2))*f_k_rho

END SUBROUTINE integrand_k_rho_prop

SUBROUTINE integrand_k_rho_evan(omega, delta, k, k_z, k_rho, r_TM, r_TE, f_k_rho)

REAL(KIND=8), INTENT(IN) :: omega, delta
REAL(KIND=8), DIMENSION(:), INTENT(IN) :: k_rho
REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f_k_rho
COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: k
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: k_z
COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: r_TM, r_TE


INTEGER :: l, n
REAL(KIND=8), PARAMETER :: c0=299792458
REAL(KIND=8), PARAMETER :: pi=4*ATAN(1.0) !pi=3.14...

n=UBOUND(k_rho,1)

DO l=1, n
f_k_rho(l) = (k_rho(l)**3/ABS(k_z(1,l)))*(IMAG(r_TM(1,l))+IMAG(r_TE(1,l)))*DEXP(-2*ABS(k_z(1,l))*delta) 
END DO

f_k_rho=((2*k(1)**2)/(omega*4*pi**2))*f_k_rho

END SUBROUTINE integrand_k_rho_evan




SUBROUTINE intg_simp (n, h, f, intg)
!
! Subroutine for integration over f(x) with the Simpson rule. f:
! integrand f(x); h: interval(dx); intg: integral.
IMPLICIT NONE
INTEGER, INTENT (IN) :: n
INTEGER :: i
REAL(KIND=8), INTENT (IN) :: h
REAL(KIND=8) :: s0,s1,s2
REAL(KIND=8), INTENT (OUT) :: intg
REAL(KIND=8), INTENT (IN), DIMENSION (:) :: f
!

intg = 0.0
s0 = 0.0
s1 = 0.0
s2 = 0.0
DO i = 2, n-1, 2
s1 = s1+f(i-1)
s0 = s0+f(i)
s2 = s2+f(i+1)
END DO
intg = h*(s1+4.0*s0+s2)/3.0
!
! If n is even, add the last slice separately
!
IF (MOD(n,2).EQ.0) intg = intg &
+h*(5.0*f(n)+8.0*f(n-1)-f(n-2))/12.0
END SUBROUTINE intg_simp
	
!****************************************************************************************
!****************************************************************************************


END MODULE mesh


