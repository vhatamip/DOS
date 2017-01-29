!  LDOS.f90 
!
!  FUNCTIONS:
!  LDOS - Local Density of State.
!

!****************************************************************************
!
!  PROGRAM: LDOS
!
!  PURPOSE:  To calculate energy density U(ω,T) and LDOS(ω,T) at a distance
!	"Δ" from medium zero (Sic) into medium 1 (vacuum).
!
!****************************************************************************
    PROGRAM LDOS
    
	use mesh
	!use properties
	!use integral
	
	IMPLICIT NONE
!****************************************************************************
!
!  Important numerical notes
!  convg_crit=1E10-7 can not handel the resonance omega, convg_crit=1E10-8 fixed the problem
!  
!****************************************************************************    
	INTEGER :: status, number_nodes
	INTEGER :: i, j, l, m
	INTEGER, PARAMETER :: omega_n_1=100, omega_n_2=100, omega_n_3=100
    INTEGER, PARAMETER :: omega_n=omega_n_1+omega_n_2+omega_n_3-2
	INTEGER, PARAMETER :: delta_n=3     !define number of deltas 
	INTEGER, PARAMETER :: nodes=3 !N points of (k_ρ) from 0 to k_vac,
										!simpson 1/3, ""N must be even"", k_vac is excluded from the first integral as k_z1=0.0.
	INTEGER, PARAMETER :: n_k_rho_mesh=50
	INTEGER, PARAMETER :: media_n=3
	
	REAL(KIND=8) :: error
	REAL(KIND=8) :: d_k_rho, d_k_rho_temp, intg_e_t_old, intg_m_t_old
	REAL(KIND=8), DIMENSION(delta_n) :: delta
	REAL(KIND=8), DIMENSION(4) :: omega_range
	REAL(KIND=8) :: int_1_min, int_2_min, int_3_min, int_1_max, int_2_max, int_3_max !omega values (intervals)
	REAL(KIND=8), DIMENSION(omega_n_1) :: omega_1
    REAL(KIND=8), DIMENSION(omega_n_2) :: omega_2
    REAL(KIND=8), DIMENSION(omega_n_3) :: omega_3
	REAL(KIND=8), DIMENSION(omega_n) :: omega, k_rho_final
	REAL(KIND=8), DIMENSION(nodes) :: k_rho
	REAL(KIND=8), DIMENSION(nodes) :: f_k_rho_e, f_k_rho_m
	REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: intg_e, intg_e_t, intg_m, intg_m_t, iter3, LDOSem
	REAL(KIND=8), PARAMETER :: c0=299792458 !Speed of light in vacuum (m/s)
	REAL(KIND=8), PARAMETER :: pi=4*ATAN(1.0) !pi=3.14...
	!REAL(KIND=8), PARAMETER :: distance = 1.0E-8 ! 10nm
	REAL(KIND=8), PARAMETER :: distance = 1000.0 ! 1m
	REAL(KIND=8), PARAMETER :: convg_crit=1.0E-9 ! convergence criterion 
	                                             
	COMPLEX(KIND=8), DIMENSION(media_n, omega_n) ::  eps, mu
	COMPLEX(KIND=8), DIMENSION(media_n) ::  k
	COMPLEX(KIND=8), DIMENSION(media_n,nodes) ::  k_z													
	COMPLEX(KIND=8), DIMENSION(2,nodes) ::  r_TE, r_TM ! 1-->01 & 2-->02 (1:2,l) r01, r02 
	COMPLEX(KIND=8), DIMENSION(2,nodes) ::  t_TE, t_TM ! 1 for t1 & 2 for t2																
	
    CHARACTER(LEN=20) :: filename
!*****************************Initializing***********************************************
	
	omega=0.0; 
	k_rho=0.0;
	f_k_rho_e=0.0;
	f_k_rho_m=0.0;
	r_TE=0.0;
	r_TM=0.0;
	k_z=0.0	
	
!									Body of LDOS
!************************************Data Import*****************************	   	   
       OPEN ( UNIT=9, FILE='data.dat', STATUS='OLD', ACTION='READ', & 
              &IOSTAT=status )
	   ! Was the OPEN successful? 
       fileopen: IF ( status == 0 ) THEN 
       READ(9,90)
       READ(9,90)
       READ (9, 100, IOSTAT=status) (delta(i), i=1,delta_n)
	   !PRINT*, 'delta', delta
	   READ(9,90)
       READ(9,90)
	   !READ(9,90)
       !READ(9,90)
	   READ (9, 100, IOSTAT=status) (omega_range(i), i=1,4)
	   !PRINT*, 'omega range', omega_range
	   90  FORMAT(A80)
       100 FORMAT(30X, E20.15)
       ELSE fileopen 
       ! Else file open failed. Tell user. 
       WRITE (*,1030) status  
       1030 FORMAT (1X,'data.dat open failed--status =', I6)
	   END IF fileopen 
!*******************************Omega discretization*********************************************
int_1_min=omega_range(1)
int_1_max=omega_range(2)
	
int_2_min=omega_range(2)
int_2_max=omega_range(3)
	
int_3_min=omega_range(3)
int_3_max=omega_range(4)
	
CALL cluster_mesh_right(int_1_min, int_1_max, omega_n_1, omega_1)
!print*, omega_1
CALL cluster_mesh_sides(int_2_min, int_2_max, omega_n_2, omega_2)
!print*, omega_2
CALL cluster_mesh_left(int_3_min, int_3_max, omega_n_3, omega_3)
!print*, omega_3
CALL sort_mesh(omega_1, omega_2, omega_3, omega)
	  

!delta=(/0.01*distance, 0.6*distance, distance/)
delta=(/1.0E-8, 5.0E-8, 5.0E-7/)

PRINT*, 'delta', delta
print*, distance
print*, omega

ALLOCATE(LDOSem(delta_n,omega_n), intg_e(delta_n,omega_n), intg_e_t(delta_n,omega_n), intg_m(delta_n,omega_n), intg_m_t(delta_n,omega_n), iter3(delta_n,omega_n))
!***************************************before loop************************************************
            LDOSem=0.0 
            intg_e_t = 0.0
			intg_m_t = 0.0
			intg_e=0.0
			intg_m=0.0
			iter3=0.0

       !DO i=1, delta_n !delta loop
	   !DO i=2,2 !delta loop
	    DO j=1, omega_n !omega loop
			CALL mu_non_mag(j, omega(j), mu(:,j)) !ok
			!CALL eps_Basu(j, omega(j), eps(:,j)) !ok 
			CALL eps_Lorentz(j, omega(j), eps(:,j))
		END DO	
		OPEN(UNIT=26, FILE='eps_mu.dat', STATUS='unknown', ACTION='write', IOSTAT=status)
		DO j=1,omega_n
			WRITE(26,190)	omega(j), (eps(i,j), i=1,media_n), (mu(i,j), i=1, media_n)
	    END DO
	  !WRITE(27,190)	   
!150	FORMAT (14X, '*****  RESULTS FILE FOR PROGRAM LDOS  *****',/,'******************************************************************')	   
!160	FORMAT (T41, 'del1', T69, 'del2', T97, 'del3', T128, 'delta 1', T159, 'delta 2', T190, 'delta 3', T219,'delt1', T248,'delt2', T277,'delt3')
!165	FORMAT (T13,'omega', T41, 'LDOS_e', T69, 'LDOS_e', T97, 'LDOS_e', T128, 'LDOS_m', T159, 'LDOS_m', T190, 'LDOS_m', T219,'iter', T248,'iter', T277,'iter')
!170	FORMAT (15X, 9(E29.2))
!175	FORMAT ('******************************************************************')
!180	FORMAT (E25.15, 6(E30.20), 3(F10.1), 9(E30.20))	   
190	FORMAT (E25.15, 6(2E30.20))
		
		
		
		    !d_k_rho_temp=((omega(1)/c0)-0.0)/(REAL(n_k_rho_mesh)-1) ! dk_ρ in the integral
			
			!d_k_rho=d_k_rho_temp*(1+(0.5)/(REAL(n_k_rho_mesh)-2)) !for eliminating k_0
		
			
	   !DO i=1, delta_n !delta loop
	   DO i=2,2 !delta loop
	    DO j=1, omega_n !omega loop
		   	
			!CALL mu_non_mag(j, omega(j), mu(:,j)) !ok
			!CALL eps_Basu(j, omega(j), eps(:,j)) !ok
			CALL k_mag(j, omega(j), mu(:,j), eps(:,j), k) !ok
			
			!d_k_rho_temp=(k(1)-0.0)/(REAL(n_k_rho_mesh)-1) ! dk_ρ in the integral
			number_nodes=INT((k(1)/((omega(1)/c0)/(n_k_rho_mesh-1)))+1)
			d_k_rho_temp=(k(1)-0.0)/(REAL(number_nodes)-1) ! dk_ρ in the integral
			
			d_k_rho=d_k_rho_temp*(1+(0.5)/(REAL(number_nodes)-2)) !for eliminating k_0
			
			k_rho=0.0
			DO
			
			DO l=1,nodes !k_ρ in the first integral
	        k_rho(l)=k_rho(nodes)+(l-1)*d_k_rho
			END DO
			
			IF 	((k_rho(nodes) - REAL(k(1))) >= 0.0) EXIT
			
			
			CALL wave_vector_z(j, mu(:,j), eps(:,j), k, k_rho, k_z)
			CALL Fresnel(j, mu(:,j), eps(:,j), k_z, r_TM, r_TE)
			CALL t_12(r_TM, r_TE, delta(i), distance, k_z, k_rho, t_TE, t_TM) !t_1(1)-->TE, t_1(2)-->TM
            !!CALL integrand_k_rho_e(omega(j), eps(2,j), k, k_rho, k_z, t_TE, t_TM, f_k_rho_e)
			!!CALL integrand_k_rho_m(omega(j), eps(2,j), k, k_rho, k_z, t_TE, t_TM, f_k_rho_m)
			CALL integrand_k_rho_prop(omega(j), k, k_z, k_rho, r_TM, r_TE, f_k_rho_e)
			CALL intg_simp (nodes, d_k_rho, f_k_rho_e, intg_e(i,j))	
			intg_e_t(i,j) = intg_e_t(i,j) + intg_e(i,j)
			intg_e_t_old = intg_e_t(i,j)
			
			END DO
			print*, REAL(k(1))-k_rho(nodes)
!***********************************************intg1***********************************************
			!PRINT*, 'n_k_rho_1', n_k_rho_1
			!PRINT*, 'n_k_rho_2', n_k_rho_2
			
	   
	        
			
			!coeff_int=1./(4*(pi**2)*omega(j)) !all the terms which are not a function of k_rho
						
			!CALL integrand_prop_k_rho(k_rho_1, mu_0, mu_1, eps_r0, eps_r1, k_vac, f_k_rho_1)	
			
			
			
			
!***********************************************intg2***********************************************
            !d_k_rho_2=(6*k_vac-k_vac)/(n_k_rho_2-1) ! dk_ρ in the integral
	        !d_k_rho_2=d_k_rho_1  ! dk_ρ in the integral
	        !DO k=1,n_k_rho_2 !k_ρ in the second integral from k_ρ to 12k_ρ
	        !k_rho_2(k)=k_vac+(k-1)*d_k_rho_2
			!END DO
			
			!coeff_int2=2.*coeff_int !all the terms which are not a function of k_rho
						
			!CALL integrand_evan_k_rho(k_rho_2,delta(i), mu_0, mu_1, eps_r0, eps_r1, k_vac, f_k_rho_2)	
			!f_k_rho_2_temp(1:n_k_rho_2-1) = f_k_rho_2(2:n_k_rho_2)
			!CALL intg_simp (n_k_rho_2-1, d_k_rho_2, f_k_rho_2_temp, intg2(i,j))
				
			     !intg2(i,j) = coeff_int2*intg2(i,j)					
						
!***********************************************int 3***********************************************		
			
			
			
			!d_k_rho=1*d_k_rho_temp	!to make finer k_rho mesh
	        !k_rho_3=0.0; k_z0_3=0.0; k_z1_3=0.0; r_TM_10_3=0.0; f_k_rho_3=0.0; r_TE_10_3=0.0
			
	        !DO k=1,n_k_rho_3 !k_ρ in the third integral from 6k_ρ to inf
	        !k_rho_3(k)=k_rho_intg2*k_vac+(k-1)*d_k_rho_3 !
			!END DO
			
			!coeff_int3=coeff_int2 !all the terms which are not a function of k_rho
						
			!CALL integrand_evan_k_rho(k_rho_3,delta(i), mu_0, mu_1, eps_r0, eps_r1, k_vac, f_k_rho_3)
							    
			    !CALL intg_simp (n_k_rho_3, d_k_rho_3, f_k_rho_3, intg3(i,j))
			    !intg3(i,j) = coeff_int3*intg3(i,j)
				
				!error=DABS((intg3(i,j)-intg3_old)/(intg1(i,j)+intg2(i,j)+intg3(i,j)))
				!error=DABS(intg3(i,j))
				!intg3_old=intg3(i,j)
			
			!DO l=1,nodes !k_ρ in the first integral
	        k_rho(1)=k_rho(nodes)+d_k_rho
			k_rho(2)=k_rho(nodes)+2*d_k_rho
			k_rho(3)=k_rho(nodes)+3*d_k_rho
			!END DO
			
			IF 	(k_rho(1) <= REAL(k(1))) EXIT
			
			
			CALL wave_vector_z(j, mu(:,j), eps(:,j), k, k_rho, k_z)
			CALL Fresnel(j, mu(:,j), eps(:,j), k_z, r_TM, r_TE)
			CALL t_12(r_TM, r_TE, delta(i), distance, k_z, k_rho, t_TE, t_TM)
				
			CALL integrand_k_rho_evan(omega(j), delta(i), k, k_z, k_rho, r_TM, r_TE, f_k_rho_m)
			CALL intg_simp (nodes, d_k_rho, f_k_rho_m, intg_m(i,j))	
			   		
			
			intg_m_t(i,j) = intg_m_t(i,j) + intg_m(i,j)
			intg_m_t_old = intg_m_t(i,j)
			
		  iter3(i,j)=1								
		  DO !iterative integral
						  
				iter3(i,j)=iter3(i,j)+1		!Counter
				
				DO l=1,nodes !k_ρ in the first integral
	                 k_rho(l)=k_rho(nodes)+(l-1)*d_k_rho
			    END DO
			    !k_rho(1)=k_rho(3)
				!k_rho(2)=k_rho(1)+d_k_rho
				!k_rho(3)=k_rho(2)+d_k_rho
				
				CALL wave_vector_z(j, mu(:,j), eps(:,j), k, k_rho, k_z)
			    CALL Fresnel(j, mu(:,j), eps(:,j), k_z, r_TM, r_TE)
			    CALL t_12(r_TM, r_TE, delta(i), distance, k_z, k_rho, t_TE, t_TM) !t_1(1)-->TE, t_1(2)-->TM
                CALL integrand_k_rho_evan(omega(j), delta(i), k, k_z, k_rho, r_TM, r_TE, f_k_rho_m)
			    CALL intg_simp (nodes, d_k_rho, f_k_rho_m, intg_m(i,j))	
				!CALL integrand_k_rho_e(omega(j), eps(2,j), k, k_rho, k_z, t_TE, t_TM, f_k_rho_e)
			    !CALL integrand_k_rho_m(omega(j), eps(2,j), k, k_rho, k_z, t_TE, t_TM, f_k_rho_m)
				
				!CALL wave_vector_z(j, mu(:,j), eps(:,j), k, k_rho, k_z)
			    !CALL Fresnel(j, mu(:,j), eps(:,j), k_z, r_TM, r_TE)
			    !CALL t_12(r_TM, r_TE, delta(i), distance, k_z, k_rho, t_TE, t_TM) !t_1(1)-->TE, t_1(2)-->TM
                !CALL integrand_k_rho_e(omega(j), eps(2,j), k, k_rho, k_z, t_TE, t_TM, f_k_rho_e)
			    !CALL integrand_k_rho_m(omega(j), eps(2,j), k, k_rho, k_z, t_TE, t_TM, f_k_rho_m)
			    
				!!CALL intg_simp (nodes, d_k_rho, f_k_rho_e, intg_e(i,j))	
			    !!CALL intg_simp (nodes, d_k_rho, f_k_rho_m, intg_m(i,j))	
			   		
			    !!intg_e_t(i,j) = intg_e_t(i,j) + intg_e(i,j)
			    intg_m_t(i,j) = intg_m_t(i,j) + intg_m(i,j)
				!intg3(i,j) = intg3(i,j) + intg3_iter(i,j)	
				
				!error=ABS((intg3(i,j)-intg3_old)/(intg3(i,j)))
				error=MAX(DABS((intg_m_t(i,j)-intg_m_t_old)/(intg_m_t(i,j))), DABS((intg_e_t(i,j)-intg_e_t_old)/(intg_e_t(i,j))))
				error=DABS((intg_m_t(i,j)-intg_m_t_old)/(intg_m_t(i,j)))
				intg_m_t_old = intg_m_t(i,j)
			    !!intg_e_t_old = intg_e_t(i,j)
				
				IF (MOD(iter3(i,j),100000.0)==0) THEN
				print*, 'error', error, 'LDOS', intg_m_t(i,j)+intg_e_t(i,j), 'i', i, 'j', j, 'iter', iter3(i,j), 'k_rho', k_rho(nodes), 'd', distance, 'z', delta(2)
				END IF
				
				IF (error < convg_crit) EXIT
				 				
				
		  END DO ! End  of iterative integral loop
				
			!intg(i,j)=intg1(i,j) +intg2(i,j) + intg3(i,j)	
		    LDOSem(i,j)=intg_m_t(i,j)+intg_e_t(i,j)
			k_rho_final(j)=	k_rho(nodes)	
			
			print*, 'error_finallllll', error, 'LDOS', LDOSem(i,j), 'intg_m', intg_m_t(i,j), 'i', i, 'j', j, 'iter', iter3(i,j), 'k_rho', k_rho(nodes)
			
		END DO ! omega loop
	   END DO  ! delta loop
	   
	   
	    OPEN(UNIT=77, FILE='k_rho.dat', STATUS='unknown', ACTION='write', IOSTAT=status)
		DO j=1,omega_n
			WRITE(77,300)	omega(j), k_rho_final(j)
	    END DO
	  
300	FORMAT (2(E30.20))
	   
	   !DO i= 1,delta_n	!delta
       ! Do j= 1,omega_n !omega
           
          !theta(j)=(h_plank*omega(j))/(DEXP((h_plank*omega(j))/(k_B*T))-1)  ! [J]
          !U_omega(i,j)=theta(j)*intg(i,j)

        !END DO
       !END DO
	   
	   
	   !*****************************************
	   !intg_T=TRANSPOSE(intg)
	   OPEN(UNIT=27, FILE='results.dat', STATUS='unknown', ACTION='write', IOSTAT=status)
	   WRITE(27,150) 
	   WRITE(27,160)
       WRITE(27,170) (delta(i), i=1,delta_n), (delta(i), i=1,delta_n), (delta(i), i=1,delta_n)
	   WRITE(27,165)
	   !WRITE(27,175)
	   
	   
	   
		   DO j=1,omega_n
			   
			WRITE(27,180)	omega(j), (LDOSem(i,j), i=1,delta_n), (intg_e_t(i,j), i=1,delta_n), (intg_m_t(i,j), i=1,delta_n)
			!WRITE(27,180)	omega(j), intg(2,j), intg1(2,j), intg2(2,j), intg3(2,j), iter3(2,j), U_omega(2,j) , (iter3(i,j), i=1,delta_n)   
			   
		   END DO
	  !WRITE(27,190)	   
150	FORMAT (14X, '*****  RESULTS FILE FOR PROGRAM LDOS  *****',/,'******************************************************************')	   
160	FORMAT (T41, 'del1', T69, 'del2', T97, 'del3', T128, 'delta 1', T159, 'delta 2', T190, 'delta 3', T219,'delt1', T248,'delt2', T277,'delt3')
165	FORMAT (T13,'omega', T41, 'LDOS_e', T69, 'LDOS_e', T97, 'LDOS_e', T128, 'LDOS_m', T159, 'LDOS_m', T190, 'LDOS_m', T219,'iter', T248,'iter', T277,'iter')
170	FORMAT (15X, 9(E29.2))
!175	FORMAT ('******************************************************************')
!180	FORMAT (E25.15, 6(E30.20), 3(F10.1), 9(E30.20))	   
180	FORMAT (E25.15, 9(E30.20))	
!190 FORMAT ('******************************************************************',/, 6X, '*****  END OF THE RESULTS FILE FOR PROGRAM LDOS  *****')	  
	
	
	
	END PROGRAM LDOS   
	   