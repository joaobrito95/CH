    SUBROUTINE CH_INVERSE(SIG6,f,mat_param,sph_tdesign,sph_tdesign_crude,&
        YesEshelby,d5_ast,J_init_flag,J_init,macroEQS,D6_norm)

    USE M_COMP_HOMOG, ONLY : Tol, &
        struct_mat_param,struct_sphtdesign,hCPB0

    IMPLICIT NONE
    REAL*8, INTENT(IN) :: SIG6(6), f
    LOGICAL, INTENT(IN) :: YesEshelby
    TYPE (struct_mat_param),  INTENT(IN) :: mat_param
    TYPE (struct_sphtdesign), INTENT(IN) :: sph_tdesign,sph_tdesign_crude
    REAL*8, INTENT(INOUT) :: macroEQS,D6_norm(6), J_init(7,7)
    REAL*8, INTENT(INOUT) :: d5_ast(5)
    LOGICAL, INTENT(INOUT) :: J_init_flag

    INTEGER maxiter,niter1,niter2,jacmethod
    REAL*8 epsTol(2),epsTol_init(2)!,rel_delta_init(2),rel_delta(2)
    REAL*8 D6_norm_x0(6),Y0_x0,x0(7)!,x0_init(7)
    LOGICAL updJacobi
    INTEGER soltype
    REAL*8 macroEQS_vect_init(1),D6_norm_init(6), x_init(7)
    REAL*8 macroEQS_vect(1), x_out(7), D6_norm_x0_dev(6)

    REAL*8 dummy5(5),nan
    TYPE (struct_sphtdesign) dummydsg
    LOGICAL bounds
    REAL*8 lb(7),ub(7),J_init_51(7,7),x0_51(7),xout_51(7),rel_delta_51(3)
    REAL*8 SIG6_dev(6),SIGm,I2(6)
    REAL*8, PARAMETER :: EPS = 2.2204460492503131D-16
    integer niter2_all, cnt_trial
    
    IF (SUM(DABS(SIG6)).EQ.0.d0) THEN
        macroEQS=0.d0
        D6_norm=0.d0
        RETURN
    ENDIF

    I2(1:3)=1.d0
    I2(4:6)=2.d0
    SIGm = SUM(SIG6(1:3))/3.d0

    !Solve system non-linear equations
    soltype = 3 !use general matrix solver
    epsTol = (/1.D-2,Tol/)
    epsTol_init = (/1.D-2,-1.D-3/)
    niter1=0
    niter2=0
    jacmethod=3

    ! Initial solution handling, test for NaN
    IF ((D6_norm(1).EQ.D6_norm(1)).AND.(macroEQS.EQ.macroEQS).AND.(macroEQS.NE.0.d0)) THEN
        ! initial solution is enforced
        D6_norm_x0 = D6_norm
        Y0_x0 = macroEQS
    ELSE
        ! initial solution is unknown and must be estimated
        CALL get_x0_eqsystem_stress(SIG6,f,mat_param,Y0_x0,D6_norm_x0)
    ENDIF

    ! Solver in the 6x6 strain space works well for high porosities (f~1d-2),
    !  but fails near incompressibility, (maybe due to numerical jacobian error)
    !  therefore the system is solved in the (dev+hydro)/5_1 space.

    ! Initial solution of 5_1_YM space
    D6_norm_x0_dev=D6_norm_x0
    x0_51(6)=SUM(D6_norm_x0(1:3))/3.d0
    CALL TE6VE5(D6_norm_x0_dev,x0_51(1:5))
    x0_51(7)=Y0_x0


    ! (Penta + Hydro) solver: bounded
    bounds=.true.
    lb(1:5)=-5.d0
    ub(1:5)=5.d0
    IF (SIGm.GT.5.d0*EPS) THEN
        lb(6)=0.d0
        ub(6)=(-2.d0/(hCPB0(2)*DLOG(f)))
    ELSEIF (SIGm.LT.5.d0*-EPS) THEN
        lb(6)=(2.d0/(hCPB0(1)*DLOG(f)))
        ub(6)=0.d0
    ELSE
        lb(6)=0.d0
        ub(6)=0.d0
    ENDIF

    !lb(6)=(2.d0/(hCPB0(1)*DLOG(f)))
    !ub(6)=(-2.d0/(hCPB0(2)*DLOG(f)))
    lb(7)=0.50d0*Y0_x0
    ub(7)=2.0d0*Y0_x0

    IF (.NOT.J_init_flag) THEN
        ! Get solution from 'scratch'
        IF (sph_tdesign_crude%N_dsgn.NE.sph_tdesign%N_dsgn) THEN
            ! 2 STEP solution
            ! Initial solution (for crude design)
            updJacobi = .FALSE.
            J_init_flag=.FALSE.
            maxiter=-15 !('crude' error tolerance)
            rel_delta_51 = (/1.D-4,1.D-3,1.D-3/)
            CALL broyden_solver3(fun_macroSIG_51_init,x0_51,7,7,3,(/5,1,1/),epsTol_init,maxiter,&
                updJacobi,rel_delta_51,jacmethod,J_init_flag,J_init_51,soltype,1,6,bounds,lb,ub,&
                x_init,macroEQS_vect_init,D6_norm_init,niter1)

            IF (x_init(1).NE.x_init(1)) THEN !now with 5_1 config.
				J_init_flag=.false.
                ! NaN solution: Retry with updJacobi
                updJacobi = .TRUE.
                CALL broyden_solver3(fun_macroSIG_51_init,x0_51,7,7,3,(/5,1,1/),epsTol_init,maxiter,&
                    updJacobi,rel_delta_51,jacmethod,J_init_flag,J_init_51,soltype,1,6,bounds,lb,ub,&
                    x_init,macroEQS_vect_init,D6_norm_init,niter1)
            ENDIF


            ! Actual solution (refined design)
            J_init_flag = .TRUE.
            updJacobi = .FALSE.
            x0 = x_init
            maxiter=-25
            niter2_all=0
            IF (x_init(1).EQ.x_init(1)) THEN
                CALL broyden_solver3(fun_macroSIG_51,x0,7,7,3,(/5,1,1/),epsTol,maxiter,&
                    updJacobi,rel_delta_51,jacmethod,J_init_flag,J_init_51,soltype,1,6,bounds,lb,ub,&
                    xout_51,macroEQS_vect,D6_norm,niter2)
                niter2_all=niter2_all+niter2
                IF (xout_51(1).NE.xout_51(1)) THEN
				J_init_flag=.false.


            ELSE
                ! NaN solution: Retry with updJacobi
                updJacobi = .TRUE.
                J_init_flag = .FALSE.
                rel_delta_51 = (/1.D-4,1.D-3,1.D-3/)
                x0=x0_51
                x0(6)=SUM(D6_norm_x0(1:3))/3.d0
                CALL broyden_solver3(fun_macroSIG_51,x0,7,7,3,(/5,1,1/),epsTol,maxiter,&
                    updJacobi,rel_delta_51,jacmethod,J_init_flag,J_init_51,soltype,1,6,bounds,lb,ub,&
                    xout_51,macroEQS_vect,D6_norm,niter2)
					  IF (xout_51(1).NE.xout_51(1)) THEN 
						  nan=DSQRT(-1.d0)
				          cnt_trial=0
                    DO WHILE (.true.)
                        maxiter=-20
                        rel_delta_51 = -(/1.D-5,1.D-4,1.D-4/)	 ! rng Jacobi and dx
                        updJacobi=.true.
                        J_init_flag=.false.
                        cnt_trial=cnt_trial+1
                        CALL broyden_solver3(fun_macroSIG_51,x0,7,7,3,(/5,1,1/),epsTol,maxiter,&
                            updJacobi,rel_delta_51,jacmethod,J_init_flag,J_init_51,soltype,1,6,bounds,lb,ub,&
                            xout_51,macroEQS_vect,D6_norm,niter2)
                        niter2_all=niter2_all+niter2
                        IF (.not.ISNAN(xout_51(1)).OR.cnt_trial.GE.3) THEN
                            !WRITE(*,*) iter(1)
                            EXIT
                        ENDIF
                    ENDDO
                ENDIF
                niter2=niter2_all	 
				ENDIF
				

            ENDIF
        ELSE
            ! 1 Step solution
            niter2=0
            updJacobi = .TRUE.
            J_init_flag = .FALSE.
            maxiter=-25
            rel_delta_51 = (/1.D-4,1.D-3,1.D-3/)
            CALL broyden_solver3(fun_macroSIG_51,x0_51,7,7,3,(/5,1,1/),epsTol,maxiter,&
                updJacobi,rel_delta_51,jacmethod,J_init_flag,J_init_51,soltype,1,6,bounds,lb,ub,&
                xout_51,macroEQS_vect,D6_norm,niter2)
            IF (xout_51(1).NE.xout_51(1)) THEN
                CALL broyden_solver3(fun_macroSIG_51,x0_51,7,7,3,(/5,1,1/),epsTol,maxiter,&
                    updJacobi,rel_delta_51,jacmethod,J_init_flag,J_init_51,soltype,1,6,bounds,lb,ub,&
                    xout_51,macroEQS_vect,D6_norm,niter2)
            ENDIF

        ENDIF
    ELSE ! J_init exists

        niter1=0
        updJacobi = .FALSE.
        maxiter=-25

        rel_delta_51 = (/1.D-4,1.D-3,1.D-3/)
        J_init_51=J_init
        CALL broyden_solver3(fun_macroSIG_51,x0_51,7,7,3,(/5,1,1/),epsTol,maxiter,&
            updJacobi,rel_delta_51,jacmethod,J_init_flag,J_init_51,soltype,1,6,bounds,lb,ub,&
            xout_51,macroEQS_vect,D6_norm,niter2)
    ENDIF

    J_init = J_init_51
    macroEQS=macroEQS_vect(1)
    !WRITE(*,*) x_out
    IF (xout_51(1).NE.xout_51(1)) THEN
        J_init_flag=.false.
        nan=DSQRT(-1.d0)
        macroEQS=nan*macroEQS
        D6_norm=nan*D6_norm
        RETURN
    ENDIF

    RETURN

    CONTAINS
    ! ===========================================
    SUBROUTINE fun_macroSIG_init(x,F_error,macroEQSi,D6_normi)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x(7)
    REAL*8, INTENT(OUT) :: F_error(7),macroEQSi,D6_normi(6)
    REAL*8  SRPlim
    ! Note ordering of output arguments is such that in Broyden solver routine,
    !  f1 == macroEQS (and f2 == D6_norm), since f1 is used to 'scale' the error vector F
    !  which as the dimention of the macro stress.


    SRPlim = 1.d0
    !IF (.NOT.YesEshelby) THEN
    !	! Rice and Tracey
    !	CALL eqsystem_stress(x,SIG6,SRPlim,YesEshelby,dummy5,f,mat_param,sph_tdesign_crude,F_error)
    !ELSE
    !	! Eshelby
    !	CALL eqsystem_stress(x,SIG6,SRPlim,YesEshelby,d5_ast,f,mat_param,sph_tdesign_crude,F_error)
    !ENDIF
    CALL eqsystem_stress(x,SIG6,SRPlim,YesEshelby,dummy5,f,mat_param,sph_tdesign_crude,F_error)

    D6_normi = x(1:6)
    macroEQSi = x(7)

    ENDSUBROUTINE fun_macroSIG_init
    ! ==========================================

    SUBROUTINE fun_macroSIG(x,F_error,macroEQS,D6_norm)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x(7)
    REAL*8, INTENT(OUT) :: F_error(7),macroEQS,D6_norm(6)
    REAL*8  SRPlim

    ! Note ordering of output arguments is such that in Broyden solver routine,
    !  f1 == macroEQS (and f2 == D6_norm), since f1 is used to 'scale' the error vector F
    !  which as the dimention of the macro stress.

    SRPlim = 1.d0
    CALL eqsystem_stress(x,SIG6,SRPlim,YesEshelby,dummy5,f,mat_param,sph_tdesign,F_error)

    D6_norm = x(1:6)
    macroEQS = x(7)

    ENDSUBROUTINE fun_macroSIG

    ! ==========================================

    SUBROUTINE fun_macroSIG_51(x,F_error_51,macroEQS,D6_norm)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x(7)
    REAL*8, INTENT(OUT) :: F_error_51(7),macroEQS,D6_norm(6)
    REAL*8  SRPlim,Ddev5(5),D6dev(6),Dm,D6(6),xD6Ym(7),F_error(7)

    ! Note ordering of output arguments is such that in Broyden solver routine,
    !  f1 == macroEQS (and f2 == D6_norm), since f1 is used to 'scale' the error vector F
    !  which as the dimention of the macro stress.

    Ddev5 = x(1:5)
    CALL VE5TE6(Ddev5,D6dev)
    Dm = x(6)
    D6=D6dev
    D6(1:3)=D6(1:3)+Dm

    xD6Ym=(/D6,x(7)/)

    SRPlim = 1.d0
    CALL eqsystem_stress(xD6Ym,SIG6,SRPlim,YesEshelby,d5_ast,f,mat_param,sph_tdesign,F_error)

    CALL TE6VE5(F_error(1:6),F_error_51(1:5))
    F_error_51(6)=SUM(F_error(1:3))/3.d0
    F_error_51(7)=F_error(7)

    D6_norm = xD6Ym(1:6)
    macroEQS = xD6Ym(7)

    ENDSUBROUTINE fun_macroSIG_51

    ! ==========================================

    SUBROUTINE fun_macroSIG_51_init(x,F_error_51,macroEQS,D6_norm)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x(7)
    REAL*8, INTENT(OUT) :: F_error_51(7),macroEQS,D6_norm(6)
    REAL*8  SRPlim,Ddev5(5),D6dev(6),Dm,D6(6),xD6Ym(7),F_error(7)

    ! Note ordering of output arguments is such that in Broyden solver routine,
    !  f1 == macroEQS (and f2 == D6_norm), since f1 is used to 'scale' the error vector F
    !  which as the dimention of the macro stress.

    Ddev5 = x(1:5)
    CALL VE5TE6(Ddev5,D6dev)
    Dm = x(6)
    D6=D6dev
    D6(1:3)=D6(1:3)+Dm

    xD6Ym=(/D6,x(7)/)

    SRPlim = 1.d0
    CALL eqsystem_stress(xD6Ym,SIG6,SRPlim,YesEshelby,d5_ast,f,mat_param,sph_tdesign_crude,F_error)


    CALL TE6VE5(F_error(1:6),F_error_51(1:5))
    F_error_51(6)=SUM(F_error(1:3))/3.d0
    F_error_51(7)=F_error(7)

    D6_norm = xD6Ym(1:6)
    macroEQS = xD6Ym(7)

    ENDSUBROUTINE fun_macroSIG_51_init

    ENDSUBROUTINE CH_INVERSE


    ! ===========================================
    SUBROUTINE eqsystem_stress(x,SIGgoal,SRPlim,YesEshelby,d5_ast,f,mat_param,&
        sph_tdesign,F_error)

    USE M_COMP_HOMOG,  ONLY: struct_mat_param,struct_sphtdesign
    IMPLICIT NONE
    REAL*8, INTENT(IN)  :: x(7), SIGgoal(6), SRPlim, f
    LOGICAL, INTENT(IN) :: YesEshelby
    REAL*8, INTENT(INOUT) :: d5_ast(5)
    TYPE (struct_mat_param),  INTENT(IN) :: mat_param
    TYPE (struct_sphtdesign), INTENT(IN) :: sph_tdesign!,sph_tdesign_crude
    REAL*8, INTENT(OUT) :: F_error(7)
    REAL*8 macroSRP,SIGMA(6)
    REAL*8 D6(6), SIGMA_norm(6), YSmatrix, maxSIGgoal

    D6 = x(1:6)
    YSmatrix = x(7)
    maxSIGgoal = MAXVAL(DABS(SIGgoal))

    IF (SUM(DABS(D6)).EQ.0.d0) THEN
        macroSRP = 0.d0
        SIGMA = 0.d0
        F_error(1:6) = (SIGMA-SIGgoal)/maxSIGgoal
        F_error(7) = (macroSRP-SRPlim)
        RETURN
    ENDIF

    CALL CH_compute(D6,f,mat_param,sph_tdesign,YesEshelby,d5_ast,macroSRP,SIGMA_norm)

    SIGMA = YSmatrix*SIGMA_norm
    F_error(1:6) = (SIGMA-SIGgoal)!/maxSIGgoal
    F_error(7) = (macroSRP-SRPlim)

    ENDSUBROUTINE eqsystem_stress

    ! ======================================================

    SUBROUTINE get_x0_eqsystem_stress(SIG6,f0,mat_param,macroEQSx0,D6_normx0)
    USE M_COMP_HOMOG,        ONLY : num_ex, &
        struct_mat_param,struct_sphtdesign, &
        hCPB0!, SRP_hydro0,Dm_hydro0,SIGMAm_hydro0

    IMPLICIT NONE
    REAL*8, INTENT(IN)                   :: SIG6(6), f0
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
    REAL*8 alpha,beta,auxfrac(2),SRP_D_dev_x0,dummy1(1),dummy2(2,2),dummy3
    INTEGER cnt_trial
    REAL*8, PARAMETER :: eps = 2.2204460492503131D-16
    REAL*8 dummyf1(1),dummyf2(1),lb(2),ub(2)
    LOGICAL J_init_flag,bounds
    REAL*8 I2(6),eye6(6),D6_aux(6),SRPcurrent,Nsig(6)

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

    ! index label:   1 - compressive, 2 - tensile
    h = hCPB0

    signs = (/1.d0,-1.d0/)
    Dm_hydro = signs*(1.d0/(h*DLOG(f0)))
    aproxSIGMAm_hydro = signs*h*DLOG(f0)/3.d0
    SRP_hydro = DABS(1.d0/Dm_hydro)


    SIGm = SUM(SIG6(1:3))/3.d0
    CALL local_eqs(SIG6,a,k,L,m,dummy,SIGeq)

    IF (DABS(SIGm).LT.5.d0*eps) THEN
        macroEQSx0 = SIGeq/(1.d0-f0)
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
                    WRITE(*,*) 'Error in get_x0_eqsystem_stress '
                    !WRITE(*,*) h
                    !READ(*,*)
                    macroEQSvect=SIGeq*(1.d0-f0)
                    EXIT
                ENDIF

            ENDDO
        ENDIF
        macroEQSx0 = macroEQSvect(1)
    ENDIF

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

    D_dev_x0 = dir_d6_unit*(1.d0-f0)
    SRP_D_dev_x0 = SRPlim*(1.d0-f0)

    IF (SIGm.LT.0.d0) THEN
        h_case=hCPB0(1)
    ELSE
        h_case=hCPB0(2)
    ENDIF

    eye6(1:3)=1.d0
    eye6(4:6)=0.d0
    I2(1:3)=1.d0
    I2(4:6)=2.d0
    D6_aux = 2.d0*(SIGeq/macroEQSx0)*(1.d0/macroEQSx0)*dir_d6_unit+ &	 !deviatoric part
        (2.d0/h_case)*(1.d0/macroEQSx0)*f0*dsinh((3.d0/h_case)*(SIGm/macroEQSx0))*eye6	!hydrostatic part
    D6_aux=D6_aux/SQRT(SUM(D6_aux*D6_aux*I2))
    CALL analyt_SRP_RT(D6_aux,f0,mat_param,SRPcurrent,Nsig)
    D6_normx0 = D6_aux/SRPcurrent

    RETURN

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
    REAL*8  D_norm(6)
    REAL*8 D6(6),I2(6),Dm_hydro_case,SRPcurrent,Nsig(6), Dm


    I2(1:3)=1.d0
    I2(4:6)=2.d0

    alpha=(x(1))
    beta=(x(2))

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

    ! =========================================================
    SUBROUTINE fun_gurson_work(x,F,T,b)
    !USE M_CONSTANTS,        ONLY: I2
    IMPLICIT NONE
    REAL*8 x(2),alpha,Dm,F(2),T,b
    REAL*8  D_norm(6)
    REAL*8 pwork,dev_pwork,vol_pwork,Dmgoal,I2(6)

    I2(1:3)=1.d0
    I2(4:6)=2.d0

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

    ENDSUBROUTINE get_x0_eqsystem_stress
    ! =========================================================

    SUBROUTINE get_aprox_hydro_limits(f0,mat_param,sphtdesign,SRP_hydro,Dm_norm_lim,aproxSIGMAm_lim)
    USE M_COMP_HOMOG, ONLY : struct_sphtdesign,struct_mat_param

    IMPLICIT NONE
    REAL*8, INTENT(IN)                   :: f0
    TYPE (struct_sphtdesign), INTENT(IN) :: sphtdesign
    TYPE (struct_mat_param),  INTENT(IN) :: mat_param
    REAL*8, INTENT(OUT) :: SRP_hydro(2),Dm_norm_lim(2),aproxSIGMAm_lim(2)
    LOGICAL  YesEshelby
    REAL*8 dummy(5), D_hydro(6),SRP_hydro_minus, SRP_hydro_plus
    REAL*8  SIGMA_hydro_minus(6), SIGMA_hydro_plus(6)


    D_hydro=0.d0
    D_hydro(1:3)=1.d0

    ! NOTE: SIGMAm_norm_lim are not the actual hydrostatic limits in the stress space
    ! since, generally, a purely hydrostatic macro strain rate does not yield a purelly
    ! hydrostatic macro stress.

    YesEshelby=.FALSE. !hydrostatic limits are the same in both v.f.
    ! Compressive hydrostatic limits
    CALL compHomog_commun(-D_hydro,YesEshelby,dummy,f0,mat_param,sphtdesign,&
        SRP_hydro_minus,SIGMA_hydro_minus)
    ! Tensile hydrostatic limits
    CALL compHomog_commun(D_hydro,YesEshelby,dummy,f0,mat_param,sphtdesign,&
        SRP_hydro_plus,SIGMA_hydro_plus)
    !YesEshelby=YesEshelbyCopy

    SRP_hydro = (/SRP_hydro_minus,SRP_hydro_plus/)
    Dm_norm_lim(1) = (SUM(-D_hydro(1:3))/3.d0)/SRP_hydro_minus
    Dm_norm_lim(2) = (SUM(D_hydro(1:3))/3.d0)/SRP_hydro_plus
    aproxSIGMAm_lim(1) = (SUM(SIGMA_hydro_minus(1:3))/3.d0)
    aproxSIGMAm_lim(2) = (SUM(SIGMA_hydro_plus(1:3))/3.d0)

    ENDSUBROUTINE  get_aprox_hydro_limits
    ! =========================================================