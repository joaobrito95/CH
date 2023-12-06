	SUBROUTINE CH_Eshelby(D_voigt,f,mat_param,sph_tdesign,sph_tdesign_crude,&
		d5_ast,macroSRP,SIGMA_norm)

	USE M_COMP_HOMOG, ONLY : struct_mat_param,struct_sphtdesign

	IMPLICIT NONE
	REAL*8, INTENT(IN) :: D_voigt(6), f
	TYPE (struct_mat_param),  INTENT(IN) :: mat_param
	TYPE (struct_sphtdesign), INTENT(IN) :: sph_tdesign,sph_tdesign_crude
	!LOGICAL, INTENT(IN) :: flag_keep_x0
	REAL*8, INTENT(OUT) :: macroSRP,SIGMA_norm(6)
	REAL*8, INTENT(INOUT) :: d5_ast(5)
	INTEGER maxiter,niter1,niter2,ntypfun,nfun(1),jacmethod
	REAL*8 epsTol(2),epsTol_init(2), d5_ast_init(5),rel_delta(1),rel_delta_init(1)
	LOGICAL updJacobi, J_init_flag
	REAL*8 J_init(5,5)
	real*8 d5_ast_x0(5)
	!INTEGER nf1, nf2
	REAL*8 macroSRP_vect(1),nan
	REAL*8 Dm, D_dev(6)
	INTEGER	soltype
	REAL*8 TolEsh
	LOGICAL 	YesEshelby

	
	IF (SUM(DABS(D_voigt)).EQ.0.d0) THEN
		  macroSRP=0.d0
		  SIGMA_norm=0.d0
		  d5_ast=0.d0
	ENDIF
	
	
	!Eshelby strain-rate field in the cartesian frame
	Dm = (1.d0/3.d0)*SUM(D_voigt(1:3))
	D_dev = D_voigt-Dm*(/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)


	YesEshelby = .TRUE.

	!Solve system non-linear equations
	TolEsh = 1D-6 ! tollx: No need for over-tolerece, due to small sensibity of macroSRP to d_ast
	soltype = 2
	epsTol_init = (/1D-3,-TolEsh/) ! minus flag: force at least one iteration
	epsTol = (/1D-2,-TolEsh/) ! 0 iter possible
	maxiter = 15
	niter1=0
	niter2=0
	jacmethod=3
	ntypfun=1
	nfun=5

	IF (sph_tdesign_crude%N_dsgn.NE.sph_tdesign%N_dsgn) THEN
		! 2 STEP solution
		! Initial solution (crude design)
		IF ((d5_ast(1).EQ.d5_ast(1)).AND.(SUM(d5_ast).NE.0.d0).AND.(DSQRT(SUM(d5_ast*d5_ast)).GT.1.d-14)) THEN
			! test for NaN
			d5_ast_x0 = d5_ast
		ELSE
			d5_ast_x0=0.d0
			!d5_ast_x0(1)=0.01d0*(Dm/f)	
			d5_ast_x0(1)=0.1d0*DSQRT(SUM(D_voigt*D_voigt))*(-1.d0/(2.d0*DLOG(f)))/f 
		ENDIF

		rel_delta_init = 1.D-3
		updJacobi = .FALSE.
		J_init_flag = .FALSE.


		CALL broyden_solver2(fun_min_SRP_init,d5_ast_x0,5,5,ntypfun,nfun,epsTol_init,maxiter,&
			updJacobi,rel_delta_init,jacmethod,J_init_flag,J_init,soltype,1,6,&
			d5_ast_init,macroSRP_vect,SIGMA_norm,niter1)
		IF (d5_ast_init(1).NE.d5_ast_init(1)) THEN
			!Retry with updJacobi
			updJacobi = .TRUE.
			CALL broyden_solver2(fun_min_SRP_init,d5_ast_x0,5,5,ntypfun,nfun,epsTol_init,maxiter,&
				updJacobi,rel_delta_init,jacmethod,J_init_flag,J_init,soltype,1,6,&
				d5_ast_init,macroSRP_vect,SIGMA_norm,niter1)
		ENDIF

		! Actual solution (refined design)
		!IF (flag_keep_x0) THEN
		! Use d5_ast as initial solution (i.e. discart d5_ast_init)
		! use with caution(!) it can speedup computation if d5_ast(in) is near
		!   d5_ast(out), increasing the number of iterations otherwise
		!d5_ast_x0 = d5_ast
		!ELSE
		! More conservative (and robust, general purpuse)
		!d5_ast_x0 = d5_ast_init
		!ENDIF
		IF (niter1.GT.0) J_init_flag = .TRUE.
		updJacobi = .FALSE.
		rel_delta = 1.D-3
		IF (d5_ast_init(1).EQ.d5_ast_init(1).AND.DSQRT(SUM(d5_ast_init*d5_ast_init)).LT.1.d6) THEN
			! d5_ast_init exists and is not NaN
			d5_ast_x0 = d5_ast_init
			CALL broyden_solver2(fun_min_SRP,d5_ast_x0,5,5,ntypfun,nfun,epsTol,maxiter,&
				updJacobi,rel_delta,jacmethod,J_init_flag,J_init,soltype,1,6,&
				d5_ast,macroSRP_vect,SIGMA_norm,niter2)
		ELSE
			updJacobi = .TRUE.
			J_init_flag = .FALSE.
			rel_delta = 1.D-3
			jacmethod=1
			soltype = 3
			CALL broyden_solver2(fun_min_SRP,d5_ast_x0,5,5,ntypfun,nfun,epsTol,maxiter,&
				updJacobi,rel_delta,jacmethod,J_init_flag,J_init,soltype,1,6,&
				d5_ast,macroSRP_vect,SIGMA_norm,niter2)
		ENDIF

	ELSE
		! 1 STEP solution
		niter1=0
		updJacobi = .TRUE.
		J_init_flag = .FALSE.
		rel_delta = 1.d-3
		IF (d5_ast(1).EQ.d5_ast(1)) THEN
			! test for NaN
			d5_ast_x0 = d5_ast
		ELSE
			d5_ast_x0=0.d0
			!d5_ast_x0(1)=0.01d0*(Dm/f)	 1.d0/(2.d0*DLOG(f))
			d5_ast_x0(1)=0.1d0*DSQRT(SUM(D_voigt*D_voigt))*(-1.d0/(2.d0*DLOG(f)))/f 
		ENDIF
		epsTol = (/1D-2,-1D-6/)
		CALL broyden_solver2(fun_min_SRP,d5_ast_x0,5,5,ntypfun,nfun,epsTol,maxiter,&
			updJacobi,rel_delta,jacmethod,J_init_flag,J_init,soltype,1,6,&
			d5_ast,macroSRP_vect,SIGMA_norm,niter2)
	ENDIF

	macroSRP=macroSRP_vect(1)

	! Verify if macroSRP is .LE. than RT.... must be, otherwise fail.. NaN

	IF (d5_ast(1).NE.d5_ast(1)) THEN
		nan=DSQRT(-1.d0)
		macroSRP=nan*macroSRP
		SIGMA_norm=nan*SIGMA_norm
		RETURN
	ENDIF


	RETURN

	CONTAINS

	SUBROUTINE fun_min_SRP_init(x_Esh,F_error,macroSRP,SIGMA_norm)
	REAL*8, INTENT(IN) ::  x_Esh(5)
	REAL*8, INTENT(OUT) :: F_error(5),macroSRP,SIGMA_norm(6)
	REAL*8 rel_delta_jac

	rel_delta_jac = 1.D-3
	CALL eqsystem_Eshelby(D_voigt,x_Esh,f,mat_param,sph_tdesign_crude,rel_delta_jac,&
		macroSRP,SIGMA_norm,F_error)

	ENDSUBROUTINE fun_min_SRP_init

	SUBROUTINE fun_min_SRP(x_Esh,F_error,macroSRP,SIGMA_norm)
	REAL*8, INTENT(IN) ::  x_Esh(5)
	REAL*8, INTENT(OUT) :: F_error(5),macroSRP,SIGMA_norm(6)
	REAL*8 rel_delta_jac

	rel_delta_jac = 1.D-6
	CALL eqsystem_Eshelby(D_voigt,x_Esh,f,mat_param,sph_tdesign,rel_delta_jac,&
		macroSRP,SIGMA_norm,F_error)

	ENDSUBROUTINE fun_min_SRP


	ENDSUBROUTINE CH_Eshelby
	! =============================================================

	SUBROUTINE eqsystem_Eshelby(D_voigt,d5_ast,f,mat_param,sphdesign,rel_delta_jac,&
		macroSRP_i,SIGMA_norm_i,F_error)

	USE M_COMP_HOMOG,  ONLY: b_radius, &
		struct_mat_param,struct_sphtdesign

	IMPLICIT NONE
	REAL*8, INTENT(IN) :: D_voigt(6), d5_ast(5), f, rel_delta_jac
	TYPE (struct_mat_param),  INTENT(IN) :: mat_param
	TYPE (struct_sphtdesign), INTENT(IN) :: sphdesign
	REAL*8, INTENT(OUT) :: macroSRP_i,SIGMA_norm_i(6),F_error(5)
	INTEGER n,m, ntypfun,nfun(1)
	REAL*8 num_gradient_SRP(1,5),y0(1),rel_delta(1)
	LOGICAL YesEshelby
	n=1
	m=5
	ntypfun=1
	nfun=5

	! Get y0 (of iteration i)
	CALL fun_Esh_d5_ast(d5_ast,macroSRP_i,SIGMA_norm_i)

	! Define error funtion:
	y0(1)=macroSRP_i
	rel_delta=rel_delta_jac
	CALL jacobi(fun_Esh_d5_ast,d5_ast,n,m,3,ntypfun,nfun,rel_delta,y0,num_gradient_SRP)
	F_error = [num_gradient_SRP]
	!F_error = F_error/macroSRP_i

	CONTAINS

	SUBROUTINE fun_Esh_d5_ast(x_penta,macroSRP,SIGMA_norm)
	IMPLICIT NONE
	REAL*8  x_penta(5)
	REAL*8  macroSRP,SIGMA_norm(6)

	YesEshelby = .TRUE.
	CALL comphomog_commun(D_voigt,YesEshelby,x_penta,f,mat_param,sphdesign,macroSRP,SIGMA_norm)

	ENDSUBROUTINE fun_Esh_d5_ast

	ENDSUBROUTINE eqsystem_Eshelby

	! =============================================================
	SUBROUTINE local_d_Eshelby(D6,d5_ast,f,r,theta,phi,d6_Esh)
	USE M_COMP_HOMOG,  ONLY: b_radius
	USE M_jit_dgemm,   ONLY: mkl_jit_matmul, &
		mkl_jit_matmul_trans!, &
	!mkl_jit_matmul_sum

	IMPLICIT NONE
	REAL*8, INTENT(IN) :: D6(6), d5_ast(5), f, r, theta, phi
	REAL*8, INTENT(OUT) :: d6_Esh(6)

	REAL*8 Dm, D_dev(6)
	REAL*8 d6_ast_dev(6), d_ast_m, d6_ast_cart(6), A6(6)
	REAL*8 Rotm(3,3), d_ast_cart(3,3), aux1(3,3)
	REAL*8 d_ast_sph(3,3), d_E_sph(3,3), d_E_cart(3,3),d6_E_cart(6)
	REAL*8 d_E_rr,d_E_tt,d_E_pp,d_E_tp,d_E_rp,d_E_rt,a_radius,a_div_r
	REAL*8 sintheta, costheta, sinphi, cosphi

	!Eshelby strain-rate field in the cartesian frame
	Dm = (1.d0/3.d0)*SUM(D6(1:3))
	D_dev = D6-Dm*(/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)

	! Get d_ast and A tensor (deviatoric)
	CALL VE5TE6(d5_ast,d6_ast_dev)
	d_ast_m = Dm/f
	d6_ast_cart = d6_ast_dev + d_ast_m*(/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)

	A6 = D_dev-f*((2.d0/5.d0)*d6_ast_dev)

	! Compute Eshelby strain rate tensor in the spherical frame
	sintheta = DSIN(theta)
	costheta = DCOS(theta)
	sinphi = DSIN(phi)
	cosphi = DCOS(phi)

	! Rotation matrix
	Rotm(1,1) = sintheta*cosphi
	Rotm(1,2) = costheta*cosphi
	Rotm(1,3) = -sinphi
	Rotm(2,1) = sintheta*sinphi
	Rotm(2,2) = costheta*sinphi
	Rotm(2,3) = cosphi
	Rotm(3,1) = costheta
	Rotm(3,2) = -sintheta
	Rotm(3,3) = 0.d0

	CALL VE6TEN(d6_ast_cart,d_ast_cart)

	!d_ast_sph = MATMUL(MATMUL(TRANSPOSE(R)),d_ast_cart),R)
	aux1=0.d0
	CALL mkl_jit_matmul(3,3,3,d_ast_cart,Rotm,aux1)
	CALL mkl_jit_matmul_trans(3,3,3,'T','N',Rotm,aux1,d_ast_sph)

	a_radius = b_radius*f**(1.d0/3.d0)
	a_div_r = a_radius/r

	d_E_rr = -2.d0*(a_div_r)**3.d0*d_ast_sph(1,1) + (4.d0/5.d0)*(a_div_r)**5*&
		( 2.d0*d_ast_sph(1,1) -      d_ast_sph(2,2) -      d_ast_sph(3,3))
	d_E_tt =       (a_div_r)**3.d0*d_ast_sph(1,1) +(1.d0/5.d0)*(a_div_r)**5*&
		(-4.d0*d_ast_sph(1,1) + 3.d0*d_ast_sph(2,2) +      d_ast_sph(3,3))
	d_E_pp =       (a_div_r)**3.d0*d_ast_sph(1,1) +(1.d0/5.d0)*(a_div_r)**5*&
		(-4.d0*d_ast_sph(1,1) +      d_ast_sph(2,2) + 3.d0*d_ast_sph(3,3))
	d_E_tp = (2.d0/5.d0)*(a_div_r)**5*d_ast_sph(2,3)
	d_E_rp = (a_div_r)**3.d0*(1.d0-(8.d0/5.d0)*(a_div_r)**2)*d_ast_sph(1,3)
	d_E_rt = (a_div_r)**3.d0*(1.d0-(8.d0/5.d0)*(a_div_r)**2)*d_ast_sph(1,2)

	CALL VE6TEN((/d_E_rr,d_E_tt,d_E_pp,d_E_tp,d_E_rp,d_E_rt/),d_E_sph)

	! Get the Eshelby strain rate tensor in the original (cartesian) frame
	!d_E_cart = MATMUL(R,MATMUL(d_E_sph,TRANSPOSE(R)))
	aux1=0.d0
	CALL mkl_jit_matmul_trans(3,3,3,'N','T',d_E_sph,Rotm,aux1)
	CALL mkl_jit_matmul(3,3,3,Rotm,aux1,d_E_cart)

	! Total local strain rate field
	CALL TENVE6(d_E_cart,d6_E_cart)
	d6_Esh = A6 + d6_E_cart

	ENDSUBROUTINE local_d_Eshelby