  SUBROUTINE analyt_SRP_RT(D6,f0,mat_param,macroSRP,SIGMA_norm)
	
	! Rice and Tracey pseudo-analytical strain-rate potential and derivative.

	USE M_COMP_HOMOG, ONLY :  hCPB0, &
		struct_mat_param
	IMPLICIT NONE
	REAL*8, INTENT(IN) :: D6(6), f0
	TYPE (struct_mat_param),  INTENT(IN) :: mat_param
	REAL*8, INTENT(OUT) :: macroSRP,SIGMA_norm(6)

	REAL*8 Dm, D6_dev(6)
	REAL*8 h(2), local_srp, s_norm(6), De
	REAL*8 h_case,gamma
	REAL*8 aux1g,auxfg,auxfactor1,auxfactor2,auxnum,auxden
	REAL*8 signDm, dSRPdDm, dSRPdDeq
	REAL*8 testunity, nan,normD6
	REAL*8 macroSDEV(6),SIGMA_m,pwork,dev_pwork,vol_pwork
	REAL*8 eye6(6), I2(6)

    eye6(1:3)=1.d0
    eye6(4:6)=0.d0
    I2(1:3)=1.d0
    I2(4:6)=2.d0

	macroSRP = 0.d0
	SIGMA_norm = 0.d0

	IF (SUM(DABS(D6)).EQ.0.d0) RETURN
	normD6 = DSQRT(SUM(D6*D6))
	IF (D6(1).NE.D6(1).OR.normD6.GT.1.d2) THEN
		nan = DSQRT(-1.d0)
		SIGMA_norm=nan*I2
		macroSRP=nan
		RETURN
	ENDIF

	! Decomposition of D6
	Dm = (1.d0/3.d0)*SUM(D6(1:3))
	D6_dev = D6-Dm*eye6

	!Solve for dir_sdev and compute SRP related with D6_dev
	CALL analyt_SRP_dense(D6_dev,mat_param,local_srp,s_norm)
	IF (local_srp.NE.local_srp) THEN
		!WRITE(*,*) '   Error in analyt_SRP_RT routine. Could not compute spr(D6_dev)'
		nan = DSQRT(-1.d0)
		SIGMA_norm=nan*I2
		macroSRP=nan
		RETURN
	ENDIF

	IF (f0.EQ.0.d0) THEN
		macroSRP=local_srp
		SIGMA_norm=s_norm
		RETURN
	ENDIF

	! local_srp to macroSRP:
	h=hCPB0
	IF (Dm.LE.0.d0) THEN
		h_case=h(1)
	ELSE
		h_case=h(2)
	ENDIF
	De = local_srp

	! analytical strain-rate potential for the porous material
	IF (Dm.EQ.0.d0) THEN
		macroSRP = De*(1.d0-f0)
	ELSEIF (De.EQ.0.d0) THEN
		macroSRP = -h_case*DABS(Dm)*DLOG(f0)
	ELSE
		gamma = h_case*DABS(Dm)/De
		aux1g = DSQRT(1.d0+gamma**2)
		auxfg = DSQRT(f0**2+gamma**2)
		auxfactor1 = (aux1g-auxfg)/gamma
		auxfactor2 = DLOG((1.d0/f0)*((gamma+auxfg)/(gamma+aux1g)))
		macroSRP = h_case*DABS(Dm)*(auxfactor1+auxfactor2)
	ENDIF


	! scalar derivatives of the analytical SRP
	signDm = 0.d0
	IF (Dm.NE.0.d0) signDm = SIGN(1.d0,Dm)
	auxfactor1 = DSQRT(De**2+(Dm*h_case)**2)
	auxfactor2 = DSQRT((De*f0)**2+(Dm*h_case)**2)
	auxnum = Dm*h_case+signDm*auxfactor2
	auxden = Dm*f0*h_case+f0*signDm*auxfactor1

	dSRPdDm = 0.d0 !limit case auxnum/auxden -> 1
	IF (auxden.NE.0.d0) dSRPdDm = h_case*signDm*DLOG(auxnum/auxden)
	dSRPdDeq = 0.d0
	IF (De.NE.0.d0) dSRPdDeq =(auxfactor1-auxfactor2)/De

	! Construct SIGMA_norm fully-analytically
	SIGMA_norm = dSRPdDm*(eye6/3.d0) + dSRPdDeq*s_norm

	testunity= DOT_PRODUCT(SIGMA_norm,D6*I2)/macroSRP
	IF (SIGMA_norm(1).NE.SIGMA_norm(1).OR.DABS((testunity-1.d0)).GT.1.0d-14) THEN
		!Compute hydrostatic macroSIGMA via power considerations
		pwork = macroSRP !Y0_matrix == 1.d0 (unitary)   ! total pwork
		macroSDEV=dSRPdDeq*s_norm
		dev_pwork = DOT_PRODUCT(macroSDEV,D6_dev*I2)    ! 'deviatoric' dissipation
		vol_pwork = (pwork-dev_pwork)                   ! 'volumetric' dissipation
		IF (Dm.NE.0.d0) THEN
			SIGMA_m = vol_pwork/(3.d0*Dm)
		ELSE
			SIGMA_m = 0.d0
		ENDIF
		SIGMA_norm(1:3) = macroSDEV(1:3)+SIGMA_m
		SIGMA_norm(4:6) = macroSDEV(4:6)
		!testunity= DOT_PRODUCT(SIGMA_norm,D6*I2)/macroSRP
		!IF (SIGMA_norm(1).NE.SIGMA_norm(1).OR.DABS((testunity-1.d0)).GT.5.0d-15) THEN
		!    WRITE(*,*) '   Error in analyt_SRP_RT routine. unitary pwork error= ', (1.d0-testunity)
		!    READ(*,*)
		!ENDIF
	ENDIF

	ENDSUBROUTINE analyt_SRP_RT
	! ========================================

	SUBROUTINE analyt_SRP_dense(D6,mat_param,SRP,sigma_norm)
	USE M_COMP_HOMOG, ONLY :  num_ex, &
		struct_mat_param
	!USE M_CONSTANTS,  ONLY: I2,eye6
	IMPLICIT NONE
	REAL*8, intent(IN) :: D6(6)
	TYPE (struct_mat_param),  INTENT(IN) :: mat_param
	REAL*8, INTENT(OUT) :: SRP,sigma_norm(6)

	REAL*8 a, k(num_ex), L(6,6,num_ex), m
	REAL*8 Dm, D6_dev(6)
	REAL*8  sigma_x0(6),dir_sdev(6),local_srp,eqs_dir_sdev
	INTEGER niter
	REAL*8 nan!,testunity
	REAL*8 eye6(6), I2(6)
   eye6(1:3)=1.d0
    eye6(4:6)=0.d0
    I2(1:3)=1.d0
    I2(4:6)=2.d0

	SRP = 0.d0
	sigma_norm = 0.d0

	IF (SUM(DABS(D6)).EQ.0.d0) RETURN
	IF (D6(1).NE.D6(1)) THEN
		nan = DSQRT(-1.d0)
		sigma_norm=nan*I2
		SRP=nan
		RETURN
	ENDIF

	! Unpack material parameters
	a = mat_param%a;
	k = mat_param%k;
	L = mat_param%L;
	m = mat_param%m;

	! Decomposition of D6
	Dm = (1.d0/3.d0)*SUM(D6(1:3))
	D6_dev = D6-Dm*eye6


	! Solve for dir_sdev and compute SRP related with D6_dev
	!sigma_x0 = 0.d0
	sigma_x0=D6_dev/DSQRT(SUM(D6_dev*D6_dev))
	CALL get_local_srp(D6_dev,sigma_x0,a,k,L,m,dir_sdev,eqs_dir_sdev,local_srp,niter)

	IF (local_srp.NE.local_srp) THEN
		!WRITE(*,*) '   Error in analyt_SRP_dense routine. Could not compute spr(D6_dev)'
		!WRITE(*,*) D6_dev, sigma_x0 !, niter
		nan = DSQRT(-1.d0)
		sigma_norm=nan*I2
		SRP=nan
		RETURN
	ENDIF

	sigma_norm = 0.d0
	IF (eqs_dir_sdev.NE.0.d0) sigma_norm = dir_sdev/eqs_dir_sdev

	SRP=local_srp

	!testunity= DOT_PRODUCT(sigma_norm,D6*I2)/SRP
	!IF (sigma_norm(1).NE.sigma_norm(1).OR.DABS((testunity-1.d0)).GT.1.0d-15) THEN
	!	WRITE(*,*) '   Error in analyt_SRP_dense routine. unitary pwork = ', testunity
	!	!READ(*,*)
	!ENDIF

	ENDSUBROUTINE analyt_SRP_dense

	
	! ====================================================

	SUBROUTINE analyt_SP_RT(SIG6,f0,mat_param,flag_EQS,macroEQSx0,D6_normx0)
	
	! Rice and Tracey pseudo-analytical stress potential and derivative.
	
	USE M_COMP_HOMOG,        ONLY : num_ex, &
		struct_mat_param,struct_sphtdesign, &
		hCPB0!, SRP_hydro0,Dm_hydro0,SIGMAm_hydro0
	!USE M_CONSTANTS,  ONLY : eye6, I2

	IMPLICIT NONE
	REAL*8, INTENT(IN)                   :: SIG6(6), f0
	INTEGER, INTENT(IN)						 :: flag_EQS
	TYPE (struct_mat_param),  INTENT(IN) :: mat_param
	REAL*8, INTENT(OUT) :: macroEQSx0,D6_normx0(6)
	REAL*8  SRP_hydro(2), Dm_hydro(2),aproxSIGMAm_hydro(2)
	REAL*8 h(2), signs(2)
	REAL*8 a, k(num_ex), L(6,6,num_ex), m
	REAL*8 SIGeq, dummy, SIGm
	REAL*8 x0_Y0(1), epsTol(2)
	INTEGER maxiter, niter,soltype,jacmethod
	LOGICAL updJacobi
	REAL*8 macroEQSvect(1)
	REAL*8 SIG_penta(5),dummy5(5),dummy6(6), T5(5,5), dir_d5(5),dir_d6_unit(6),h_case
	REAL*8 SIG6dev(6),dir_d6(6),D_dev_x0(6)
	REAL*8 SRPlim,Dm,x0(2),rel_delta(2), xout(2) !,Dmvect(2)
	REAL*8 alpha,beta,auxfrac(2),SRP_D_dev_x0,dummy2(2,2)!,dummy1(1),dummy3
	INTEGER cnt_trial
	REAL*8, PARAMETER :: eps = 2.2204460492503131D-16
	REAL*8 dummyf1(1),dummyf2(1),lb(2),ub(2)
	LOGICAL J_init_flag,bounds
	REAL*8 D6_aux(6),SRPcurrent,Nsig(6)
		REAL*8 eye6(6), I2(6)

        
       eye6(1:3)=1.d0
    eye6(4:6)=0.d0
    I2(1:3)=1.d0
    I2(4:6)=2.d0
    
	macroEQSx0=0.d0
	D6_normx0=0.d0

	IF (SUM(SIG6).EQ.0.d0) RETURN
	IF (SIG6(1).NE.SIG6(1)) RETURN

	! 1) Get initial guess for macroEQS:
	! using gurson-like criteria (stress space)
	! note: only a rought approximation, (quadratic shape-preserving criterion)
	a = mat_param%a;
	k = mat_param%k;
	L = mat_param%L;
	m = mat_param%m;

	CALL local_eqs(SIG6,a,k,L,m,dummy,SIGeq)
	IF (f0.EQ.0.d0) THEN
		! Dense material: no porous phase
		macroEQSx0=SIGeq
	ELSE
		! Porous material
		! index label:   1 - compressive, 2 - tensile
		h = hCPB0
		SIGm = SUM(SIG6(1:3))/3.d0
		IF (DABS(SIGm).LT.5.d0*eps) THEN
			macroEQSx0 = SIGeq*(1.d0-f0)
		ELSE
			x0_Y0 = SIGeq/(1.d0-f0)
			epsTol=(/1.D-3,1.D-6/)
			maxiter = 50
			updJacobi = .TRUE.
			CALL broyden_solver(fun_gurson_Y0,x0_Y0,1,1,epsTol,maxiter,updJacobi,macroEQSvect,niter)
			IF (macroEQSvect(1).NE.macroEQSvect(1)) THEN
				cnt_trial=0
				x0_Y0 = DABS((3.d0*SIGm)/(MINVAL(h)*DLOG(f0)))
				DO WHILE (.true.)
					cnt_trial=cnt_trial+1
					CALL broyden_solver(fun_gurson_Y0,x0_Y0,1,1,epsTol,maxiter,updJacobi,macroEQSvect,niter)
					IF (macroEQSvect(1).NE.macroEQSvect(1)) THEN
						x0_Y0 = x0_Y0/1.5d0
					ELSE
						EXIT
					ENDIF
					IF (cnt_trial.GT.10) THEN
						WRITE(*,*) 'Error in analyt_SP_RT '
						!WRITE(*,*) h
						!READ(*,*)
						macroEQSvect=SIGeq*(1.d0-f0)
						EXIT
					ENDIF

				ENDDO
				!IF (macroEQSvect(1).NE.macroEQSvect(1)) THEN
				! CALL broyden_solver(fun_dgurson_dY0,x0_Y0,1,1,epsTol,maxiter,updJacobi,macroEQSvect,niter)
				! WRITE(*,*) 'check here'
				! READ(*,*)
				!ENDIF
			ENDIF
			macroEQSx0 = macroEQSvect(1)
		ENDIF
	ENDIF

	IF (flag_EQS.EQ.1) RETURN

	! 2) Get initial guess for D6_norm:
	SIG6dev = SIG6
	SIG6dev(1:3) = SIG6(1:3)-SIGm
	CALL TE6VE5(SIG6dev,SIG_penta)
	dummy6 = 0.d0
	CALL penta_Eqsystem(SIG_penta,a,k,L,dummy6,dir_d5,T5,dummy5)
	CALL VE5TE6(dir_d5,dir_d6)
	SRPlim=1.d0
	dir_d6_unit = SRPlim*(m**a)*(SIGeq**(1.d0-a))*dir_d6
	IF (SIGeq.EQ.0.d0.OR.dir_d6_unit(1).NE.dir_d6_unit(1)) dir_d6_unit=0.d0

	IF (f0.EQ.0.d0) THEN
		! Dense material: no porous phase
		D_dev_x0=dir_d6_unit
	ELSE
		! Porous material
		D_dev_x0 = dir_d6_unit/(1.d0-f0)
		SRP_D_dev_x0 = SRPlim/(1.d0-f0)

		IF (SIGm.LT.0.d0) THEN
			h_case=hCPB0(1)
		ELSE
			h_case=hCPB0(2)
		ENDIF

		D6_aux = 2.d0*(SIGeq/macroEQSx0)*(1.d0/macroEQSx0)*dir_d6_unit+ &	 !deviatoric part
			(2.d0/h_case)*(1.d0/macroEQSx0)*f0*dsinh((3.d0/h_case)*(SIGm/macroEQSx0))*eye6	!hydrostatic part
		D6_aux=D6_aux/SQRT(SUM(D6_aux*D6_aux*I2))
		CALL analyt_SRP_RT(D6_aux,f0,mat_param,SRPcurrent,Nsig)
		D6_normx0 = D6_aux/SRPcurrent
	ENDIF

	RETURN

	signs = (/1.d0,-1.d0/)
	Dm_hydro = signs*(1.d0/(h*DLOG(f0)))
	aproxSIGMAm_hydro = signs*h*DLOG(f0)/3.d0
	SRP_hydro = DABS(1.d0/Dm_hydro)

	IF (DABS(SIGm).LT.5.d0*eps) THEN
		alpha = 1.d0
		Dm = 0.d0
	ELSE
		!solve for alpha and Dm
		auxfrac = ((SIGm/macroEQSx0)/aproxSIGMAm_hydro)
		x0(1)=DSQRT(DABS(1.d0-MINVAL(auxfrac)**2)) ! x0 alpha trignometric unit sphere
		alpha = x0(1)
		IF (SIGm.LT.0.d0) THEN
			IF (auxfrac(1).LT.-1.d0) auxfrac=-1.d0
			Dm = auxfrac(1)*Dm_hydro(1) !equal fraction
			x0(2) = auxfrac(1)
		ELSE
			IF (auxfrac(2).GT.1.d0) auxfrac=1.d0
			Dm = auxfrac(2)*Dm_hydro(2) !equal fraction
			x0(2) = auxfrac(2)
		ENDIF

		!Non-symmetric system: broyden_solver2
		epsTol=(/10.d0,1.D-6/)
		maxiter = -25
		updJacobi = .TRUE.
		rel_delta=1.d-3
		jacmethod=3
		soltype=3
		J_init_flag=.false.
		bounds=.true.
		lb=0.d0
		ub=1.d0
		!CALL broyden_solver2(fun_alpha_beta,x0,2,2,2,(/1,1/),epsTol,maxiter,updJacobi,rel_delta,jacmethod,&
		!	.FALSE.,dummy2,soltype,1,1,xout,dummyf1,dummyf2,niter)

		CALL broyden_solver3(fun_alpha_beta,x0,2,2,2,(/1,1/),epsTol,maxiter,updJacobi,rel_delta,jacmethod,&
			J_init_flag,dummy2,soltype,1,1,bounds,lb,ub,&
			xout,dummyf1,dummyf2,niter)

		IF (xout(1).EQ.xout(1)) THEN
			!CALL fun_alpha_beta(xout,dummy2,dummy3,dummy3)
			alpha=xout(1)
			beta=xout(2)
			IF (SIGm.LT.0.d0) THEN
				Dm = beta*Dm_hydro(1)
			ELSE
				Dm = beta*Dm_hydro(2)
			ENDIF
			alpha=DABS(alpha)
			beta=DABS(beta)

			alpha = MIN(alpha,1.d0)
			alpha = MAX(alpha,0.d0)
			beta = MIN(beta,1.d0)
			beta = MAX(beta,0.d0)
		ENDIF

	ENDIF


	D6_normx0(1:3) = alpha*D_dev_x0(1:3)+Dm
	D6_normx0(4:6) = alpha*D_dev_x0(4:6)
	WRITE(*,*) D6_normx0
	CONTAINS

	! =========================================================
	SUBROUTINE fun_alpha_beta(x,F,T,b)
	!USE M_CONSTANTS,        ONLY: I2
	IMPLICIT NONE
	REAL*8 x(2),alpha,beta,F(2),T,b
	!REAL*8  D_norm(6)
	REAL*8 D6(6),Dm_hydro_case,SRPcurrent,Nsig(6), Dm
	REAL*8 eye6(6), I2(6)

    I2(1:3)=1.d0
    I2(4:6)=2.d0

	alpha=(x(1))
	beta=(x(2))

	!alpha = MIN(alpha,1.d0)
	!alpha = MAX(alpha,0.d0)
	!beta = MIN(beta,1.d0)
	!beta = MAX(beta,0.d0)

	IF (SIGm.GT.5.d0*eps) THEN
		Dm_hydro_case=Dm_hydro(2)
	ELSEIF (SIGm.LT.-5.d0*eps) THEN
		Dm_hydro_case=Dm_hydro(1)
	ELSE
		beta=0.d0
	ENDIF

	Dm = 	beta*Dm_hydro_case
	D6(1:3) = alpha*D_dev_x0(1:3)+Dm
	D6(4:6) = alpha*D_dev_x0(4:6)


	F(1) = DOT_PRODUCT(D6*I2,SIG6) - SRPlim*macroEQSx0

	CALL analyt_SRP_RT(D6,f0,mat_param,SRPcurrent,Nsig)
	F(2) = SRPcurrent-SRPlim

	!F = DABS(F)

	T = 1.d0
	b = macroEQSx0

	ENDSUBROUTINE fun_alpha_beta


	! =========================================================
	SUBROUTINE fun_gurson_work(x,F,T,b)
	!USE M_CONSTANTS,        ONLY: I2
	IMPLICIT NONE
	REAL*8 x(2),alpha,Dm,F(2),T,b
	REAL*8  D_norm(6)
	REAL*8 pwork,dev_pwork,vol_pwork,Dmgoal


	alpha=DABS(x(1))
	Dm=DABS(x(2))*SIGN(1.d0,SIGm)

	IF (alpha.GT.1.d0) alpha=1.d0
	IF (Dm.LT.Dm_hydro(1)) Dm = Dm_hydro(1)
	IF (Dm.GT.Dm_hydro(2)) Dm = Dm_hydro(2)

	D_norm(1:3) = alpha*D_dev_x0(1:3)+Dm
	D_norm(4:6) = alpha*D_dev_x0(4:6)
	pwork = SRPlim*macroEQSx0                  ! total pwork
	dev_pwork =DOT_PRODUCT(D_norm*I2,SIG6dev)     ! 'deviatoric' dissipation
	vol_pwork = DABS(pwork-dev_pwork)              ! 'volumetric' dissipation
	IF (SIGm.NE.0.d0.AND.DABS(SIGm).GT.EPS) THEN
		Dmgoal = (vol_pwork/SIGm)/3.d0
	ELSE
		Dmgoal = 0.d0
	ENDIF
	F(1) = DOT_PRODUCT(D_norm*I2,SIG6) - SRPlim*macroEQSx0
	F(2) = Dmgoal-Dm

	!F = DABS(F)

	T = 1.d0
	b = 1.d0

	ENDSUBROUTINE fun_gurson_work

	! =========================================================

	SUBROUTINE fun_gurson_Y0(Y0,F,T,b)
	IMPLICIT NONE
	REAL*8 Y0, F,T,b, h_case

	IF (SIGm.LE.0.d0) THEN
		h_case=h(1)
	ELSE
		h_case=h(2)
	ENDIF
	Y0 = DABS(Y0)
	F = (SIGeq/Y0)**2+2.d0*f0*DCOSH((3.d0/h_case)*(SIGm/Y0))-f0**2-1.d0
	!F = DABS(F)

	T = 1.d0
	b = 1.d0

	ENDSUBROUTINE fun_gurson_Y0
	! =========================================================

	SUBROUTINE fun_dgurson_dY0(Y0,F,T,b)
	IMPLICIT NONE
	REAL*8 Y0, F,T,b, h_case

	IF (SIGm.LE.0.d0) THEN
		h_case=h(1)
	ELSE
		h_case=h(2)
	ENDIF

	F = h_case*SIGeq**2+3.d0*SIGm*Y0*f0*DSINH((3.d0/h_case)*(SIGm/Y0))
	!F = DABS(F)

	T = 1.d0
	b = 1.d0

	ENDSUBROUTINE fun_dgurson_dY0

	ENDSUBROUTINE analyt_SP_RT
	! =========================================================