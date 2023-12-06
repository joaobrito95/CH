	SUBROUTINE get_local_srp(d_local,sigma_x0,a,k,L,m,dir_sdev,eqs_dir_sdev,local_srp,niter)

	USE M_COMP_HOMOG,   ONLY: num_ex

	IMPLICIT NONE
	REAL*8, INTENT(IN) :: d_local(6)
	REAL*8, INTENT(IN) :: a,k(num_ex),L(6,6,num_ex),m
	REAL*8, INTENT(OUT) ::  dir_sdev(6), eqs_dir_sdev, local_srp
	REAL*8  norm_d_local, sigma_x0(6),sigma_x0_penta(5)  !F(5),
	REAL*8 epsTol(2), sdir_penta(5),nan
	LOGICAL updJacobi
	INTEGER niter,niter1,niter2, maxiter1,maxiter2
	REAL*8 d_local_pwork(6)!,G
	REAL*8 F
	INTEGER  jacmethod,soltype	!flag_norm_x0
	REAL*8 rel_delta(1),J0(5,5),f1_norm(25),f2_norm(5)!,f1(5*5),f2(5)
	LOGICAL J0_flag
	REAL*8 tol_f
	REAL*8 Tol_local
	norm_d_local = DSQRT(SUM(d_local*d_local))

	IF (norm_d_local.EQ.0.d0) THEN
		dir_sdev=0.d0
		eqs_dir_sdev=0.d0
		local_srp=0.d0
		RETURN
	ENDIF

	!d_local = d_local/norm_d_local

	sigma_x0(3) =-(sigma_x0(1)+sigma_x0(2))
	CALL TE6VE5(sigma_x0,sigma_x0_penta)


	!Solve system non-linear equations d=Ts

    Tol_local=1.d-8
	IF (a.EQ.2.d0) THEN
		tol_f = MAX(Tol_local/1.d3,1.0d-9)
	ELSE
		tol_f = MAX(Tol_local/1.d3,1.0d-8)
    ENDIF
    maxiter1=50
	maxiter2=100
	tol_f = MIN(tol_f,1.0d-3)
	epsTol = (/1.0d-2,tol_f/)
	!epsTol=-1.d0*DABS(epsTol)

	niter=0
	niter1=0
	niter2=0
	!flag_norm_x0=0

	IF (a.EQ.2.d0) THEN
		! System is linear: normalisation with d:d (not helpfull to add non-linearity by normalizing with (m^a)s:d).
		updJacobi = .FALSE.
		CALL broyden_solver(fun_penta,sigma_x0_penta,5,5,epsTol,maxiter1,updJacobi,sdir_penta,niter1)
		IF (sdir_penta(1).NE.sdir_penta(1)) THEN
			! Retry, updating jacobian at each iteration
			updJacobi = .TRUE.
			CALL broyden_solver(fun_penta,sigma_x0_penta,5,5,epsTol,maxiter2,updJacobi,sdir_penta,niter2)
		ENDIF
	ELSE
		! System is always non-linear: solve by normilizing with (m^a)s:d (by definition, EQS(sdir_penta)==1), hence x0 is more stable
		updJacobi = .TRUE.
		CALL broyden_solver(fun_penta,sigma_x0_penta,5,5,epsTol,-maxiter2,updJacobi,sdir_penta,niter1)
		IF (sdir_penta(1).NE.sdir_penta(1)) THEN
			CALL broyden_solver(fun_penta_norm,sigma_x0_penta,5,5,epsTol,maxiter2,updJacobi,sdir_penta,niter2)
			niter =niter+niter2
			IF (sdir_penta(1).NE.sdir_penta(1)) THEN
				! Retry with different solver
				rel_delta(1)=1.d-4
				jacmethod=3 !foward
				J0_flag=.false.
				soltype=3 !general matrix solver
				CALL broyden_solver2(fun_penta_norm,sigma_x0_penta,5,5,1,(/5/),epsTol,maxiter2,updJacobi,rel_delta,&
					jacmethod,J0_flag,J0,soltype,5,5,sdir_penta,f1_norm,f2_norm,niter2)
				niter=niter++niter2
				IF (sdir_penta(1).NE.sdir_penta(1)) THEN
					! Retry, with alternative system, different x0
					updJacobi = .TRUE.
					CALL TE6VE5(d_local/norm_d_local,sigma_x0_penta) ! new sigma_x0_penta
					CALL broyden_solver(fun_penta,sigma_x0_penta,5,5,epsTol,maxiter2,updJacobi,sdir_penta,niter2)
					!IF (sdir_penta(1).EQ.sdir_penta(1)) flag_norm_x0=1
					niter=niter++niter2
				ENDIF
			ENDIF
		ENDIF
	ENDIF
	niter =(niter1+niter2)

	!Return if no convergence:
	IF (sdir_penta(1).NE.sdir_penta(1)) THEN
		write(*,*) 'nan local_srp.f90'
		nan = DSQRT(-1.d0)
		dir_sdev=nan*dir_sdev
		eqs_dir_sdev=nan
		local_srp=nan
		RETURN
	ENDIF

	!Penta to Voigt notation
	CALL VE5TE6(sdir_penta,dir_sdev)

	!Equivalent stress of dir_sdev
	CALL local_eqs(dir_sdev,a,k,L,m,F,eqs_dir_sdev)


	!Compute strain rate potential
	d_local_pwork = d_local
	d_local_pwork(4:6)=2.d0*d_local(4:6)
	local_srp=DOT_PRODUCT(d_local_pwork,dir_sdev)/eqs_dir_sdev


	!IF (flag_norm_x0.EQ.1)  dir_sdev=dir_sdev/eqs_dir_sdev  ! a~=2 system convention for dir_sdev

	!! OLD:
	!!Compute strain rate potential
	!CALL compute_G(dir_sdev,a,k,L,G)
	!local_srp =(1.d0/m)*G**((a-1.d0)/a)
	!
	!!Equivalent stress of dir_sdev
	!d_local_pwork = d_local
	!d_local_pwork(4:6)=2.d0*d_local(4:6)
	!eqs_dir_sdev = DOT_PRODUCT(d_local_pwork,dir_sdev)/local_srp



	RETURN

	CONTAINS

	SUBROUTINE fun_penta_norm(x_penta,F,d_local_penta,d_local_penta_copy)
	IMPLICIT NONE
	REAL*8 x_penta(5), F(5)
	REAL*8 T_penta(5,5), d_local_penta(5),d_local_penta_copy(5)

	CALL penta_norm_Eqsystem(x_penta,a,k,L,m,d_local,F,T_penta,d_local_penta)
	d_local_penta_copy=d_local_penta

	ENDSUBROUTINE fun_penta_norm
	! -------------------------------------------------------------
	SUBROUTINE fun_penta(x_penta,F,T_penta,d_local_penta)
	IMPLICIT NONE
	REAL*8 x_penta(5), F(5)
	REAL*8 T_penta(5,5), d_local_penta(5)

	CALL penta_Eqsystem(x_penta,a,k,L,d_local,F,T_penta,d_local_penta)

	ENDSUBROUTINE fun_penta
	ENDSUBROUTINE get_local_srp
	! =============================================================

	SUBROUTINE penta_norm_Eqsystem(s_penta,a,k,L,m,d_local,F,T_penta,d_local_penta_norm)
	USE M_COMP_HOMOG, ONLY: num_ex
	USE M_CONSTANTS,  ONLY: Trans56, Trans65, I2
	USE M_jit_dgemm,  ONLY: mkl_jit_matmul, &
		mkl_jit_matmul_trans, &
		mkl_jit_matmul_sum

	IMPLICIT NONE
	REAL*8, INTENT(IN) :: d_local(6),s_penta(5),a,k(num_ex),L(6,6,num_ex),m
	REAL*8, INTENT(OUT) ::  F(5),T_penta(5,5), d_local_penta_norm(5)

	INTEGER i!,jj!,j
	REAL*8 s_voigt(6),s_trans_voigt(6), s_trans(3,3),k_i
	REAL*8 eigvect(3,3), eig_s_trans(3)
	REAL*8 Tn_aux_sp(3,3),sign_vect(3)
	REAL*8 Tn33_eff(3,3), Tn66(6,6), T(6,6)

	REAL*8 L_i(6,6)
	REAL*8 Tn_aux_sp_eigvT(3,3)
	REAL*8 Ts(5),T_Trans65(6,5), Tn66L_i(6,6) !T_i(6,6),
	REAL*8 aux_norm


	CALL VE5TE6(s_penta,s_voigt)


	T=0.d0
	DO i=1,num_ex
		L_i = L(:,:,i)
		k_i = k(i)

		!Compute stress transformation, sn
		!s_trans_voigt = MATMUL(L_i,s_voigt)
		CALL mkl_jit_matmul(6,6,1,L_i,s_voigt,s_trans_voigt)

		! Compute eigen decomposition
		CALL VE6TEN(s_trans_voigt,s_trans)
		CALL DSYEVH3(s_trans, eigvect, eig_s_trans)

		!Compute diagonal matrices: Tn, |sp_hat|^(a-2) and their prodct, Tn_aux_sp
		CALL get_aux_Tn_diag(eig_s_trans,k_i,a,sign_vect,Tn_aux_sp)

		!Tn33_eff = MATMUL(eigvect,MATMUL(Tn_aux_sp,TRANSPOSE(eigvect)))
		CALL mkl_jit_matmul_trans(3,3,3,'N','T',Tn_aux_sp,eigvect,Tn_aux_sp_eigvT)
		CALL mkl_jit_matmul(3,3,3,eigvect,Tn_aux_sp_eigvT,Tn33_eff)

		! Equivalent voigt notation:
		CALL dsymm_voigt66(Tn33_eff,Tn66)

		!Construct T
		!T = T + MATMUL(L_i_T,MATMUL(Tn66,L_i))
		CALL mkl_jit_matmul(6,6,6,Tn66,L_i,Tn66L_i)
		CALL mkl_jit_matmul_sum(6,6,6,'T','N',L_i,Tn66L_i,T)
		! Note: T is not symmetric since d_local is not written in *[1 1 1 2 2 2] form

	ENDDO


	!T_penta=MATMUL(Trans56,MATMUL(T,Trans65))
	CALL mkl_jit_matmul(6,6,5,T,Trans65,T_Trans65)
	CALL mkl_jit_matmul(5,6,5,Trans56,T_Trans65,T_penta)
	! Note: T_penta is symmetric (even if T is not)



	!F = MATMUL(T_penta,s_penta)-d_local_penta
	CALL mkl_jit_matmul(5,5,1,T_penta,s_penta,Ts)

	aux_norm = (m**NINT(a))*(SUM((d_local*I2)*s_voigt)) !DABS guarantees positive work (sound solution)
	CALL TE6VE5(d_local/aux_norm,d_local_penta_norm)

	F = Ts-d_local_penta_norm	+ dabs(aux_norm-dabs(aux_norm))

	! TE6VE5(T*s - d_local) == TE6VE5(T)*TE6VE5(s) - TE6VE5(d_local))????? test in matlab...

	ENDSUBROUTINE penta_norm_Eqsystem
	! =============================================================
	SUBROUTINE penta_Eqsystem(s_penta,a,k,L,d_local,F,T_penta,d_local_penta_norm)
	USE M_COMP_HOMOG, ONLY: num_ex
	USE M_CONSTANTS,  ONLY: Trans56, Trans65, I2
	USE M_jit_dgemm,  ONLY: mkl_jit_matmul, &
		mkl_jit_matmul_trans, &
		mkl_jit_matmul_sum

	IMPLICIT NONE
	REAL*8, INTENT(IN) :: d_local(6),s_penta(5),a,k(num_ex),L(6,6,num_ex)
	REAL*8, INTENT(OUT) ::  F(5),T_penta(5,5), d_local_penta_norm(5)

	INTEGER i!,jj!,j
	REAL*8 s_voigt(6),s_trans_voigt(6), s_trans(3,3),k_i
	REAL*8 eigvect(3,3), eig_s_trans(3)
	REAL*8 Tn_aux_sp(3,3),sign_vect(3)
	REAL*8 Tn33_eff(3,3), Tn66(6,6), T(6,6)

	REAL*8 L_i(6,6)
	REAL*8 Tn_aux_sp_eigvT(3,3)
	REAL*8 Ts(5),T_Trans65(6,5), Tn66L_i(6,6) !T_i(6,6),
	REAL*8 aux_norm,d_local_penta(5)

	CALL VE5TE6(s_penta,s_voigt)

	T=0.d0
	DO i=1,num_ex
		L_i = L(:,:,i)
		k_i = k(i)

		!Compute stress transformation, sn
		!s_trans_voigt = MATMUL(L_i,s_voigt)
		CALL mkl_jit_matmul(6,6,1,L_i,s_voigt,s_trans_voigt)

		! Compute eigen decomposition
		CALL VE6TEN(s_trans_voigt,s_trans)
		CALL DSYEVH3(s_trans, eigvect, eig_s_trans)

		!Compute diagonal matrices: Tn, |sp_hat|^(a-2) and their prodct, Tn_aux_sp
		CALL get_aux_Tn_diag(eig_s_trans,k_i,a,sign_vect,Tn_aux_sp)

		!Tn33_eff = MATMUL(eigvect,MATMUL(Tn_aux_sp,TRANSPOSE(eigvect)))
		CALL mkl_jit_matmul_trans(3,3,3,'N','T',Tn_aux_sp,eigvect,Tn_aux_sp_eigvT)
		CALL mkl_jit_matmul(3,3,3,eigvect,Tn_aux_sp_eigvT,Tn33_eff)

		! Equivalent voigt notation:
		CALL dsymm_voigt66(Tn33_eff,Tn66)

		!Construct T
		!T = T + MATMUL(L_i_T,MATMUL(Tn66,L_i))
		CALL mkl_jit_matmul(6,6,6,Tn66,L_i,Tn66L_i)
		CALL mkl_jit_matmul_sum(6,6,6,'T','N',L_i,Tn66L_i,T)
		! Note: T is not symmetric since d_local is not written in *[1 1 1 2 2 2] form

	ENDDO

	!    IF (SUM(s_penta).NE.0.d0) THEN
	!    WRITE(*,*) ' '
	!        DO jj=1,6
	!            WRITE(*,37) T(jj,1),T(jj,2),T(jj,3), &
	!                T(jj,4),T(jj,5),T(jj,6)
	!        ENDDO
	!37      FORMAT(5X,6F15.4)
	!    WRITE(*,*) ' '
	!    ENDIF

	!T_penta=MATMUL(Trans56,MATMUL(T,Trans65))
	CALL mkl_jit_matmul(6,6,5,T,Trans65,T_Trans65)
	CALL mkl_jit_matmul(5,6,5,Trans56,T_Trans65,T_penta)
	! Note: T_penta is symmetric (even if T is not)

	!    IF (SUM(s_penta).NE.0.d0) THEN
	!    WRITE(*,*) ' '
	!        DO jj=1,5
	!            WRITE(*,38) T_penta(jj,1),T_penta(jj,2),T_penta(jj,3), &
	!                T_penta(jj,4),T_penta(jj,5)
	!        ENDDO
	!38      FORMAT(5X,5F15.4)
	!    WRITE(*,*) ' '
	!    ENDIF


	!F = MATMUL(T_penta,s_penta)-d_local_penta
	CALL mkl_jit_matmul(5,5,1,T_penta,s_penta,Ts)


	aux_norm = DSQRT(SUM((d_local*I2)*d_local))
	IF (aux_norm.EQ.0.d0)  aux_norm=1.d0 !(CH_Inverse F==d)
	CALL TE6VE5(d_local,d_local_penta)
	d_local_penta_norm=d_local_penta/aux_norm

	F = Ts-d_local_penta_norm

	! TE6VE5(T*s - d_local) == TE6VE5(T)*TE6VE5(s) - TE6VE5(d_local))????? test in matlab...

	ENDSUBROUTINE penta_Eqsystem
	! =============================================================

	SUBROUTINE get_aux_Tn_diag(eig_s_trans,k,a,sign_vect,Tn_aux_sp)
	USE M_CONSTANTS, ONLY : eye3

	IMPLICIT NONE
	REAL*8, INTENT(IN)  :: eig_s_trans(3),k,a
	REAL*8, INTENT(OUT) :: sign_vect(3),Tn_aux_sp(3,3)
	REAL*8 sign_vect_j,Tn_jj,aux_sp_jj,eig_s_trans_j
	INTEGER j

	Tn_aux_sp=eye3
	IF (a.NE.2.d0) THEN
		DO j=1,3
			eig_s_trans_j = eig_s_trans(j)
			IF (eig_s_trans_j.NE.0.d0) THEN
				sign_vect_j = SIGN(1.d0,eig_s_trans_j)
				aux_sp_jj = DABS(eig_s_trans_j)**NINT(a-2.d0)
			ELSE
				sign_vect_j = 0.d0
				aux_sp_jj=0.d0
			ENDIF

			IF (k.eq.0.d0) THEN
				Tn_jj=1.d0
			ELSE
				sign_vect(j) = sign_vect_j
				Tn_jj = (1.d0-sign_vect_j*k)**NINT(a)
			ENDIF
			Tn_aux_sp(j,j) = Tn_jj*aux_sp_jj
		ENDDO
	ELSE

		IF (k.eq.0.d0) THEN
			!Tn_aux_sp=eye3
			RETURN
		ENDIF

		DO j=1,3
			eig_s_trans_j = eig_s_trans(j)
			IF (eig_s_trans_j.NE.0.d0) THEN
				sign_vect_j = SIGN(1.d0,eig_s_trans_j)
			ELSE
				sign_vect_j = 0.d0

			ENDIF
			sign_vect(j)=sign_vect_j
			Tn_jj = (1.d0-sign_vect_j*k)**2
			Tn_aux_sp(j,j) = Tn_jj
		ENDDO
	ENDIF



	ENDSUBROUTINE

	! =============================================================
	SUBROUTINE compute_G(s_voigt,a,k,L,G)
	USE M_COMP_HOMOG, ONLY: num_ex
	USE M_jit_dgemm,  ONLY: mkl_jit_matmul

	IMPLICIT NONE
	REAL*8, INTENT(IN)  :: s_voigt(6),a,k(num_ex),L(6,6,num_ex)
	REAL*8, INTENT(OUT) :: G

	INTEGER i,j
	REAL*8 s_trans_voigt(6), s_trans(3,3),L_i(6,6),k_i
	REAL*8 eigvect(3,3), eig_s_trans(3)
	REAL*8 sign_vect(3),Tn_aux_sp(3,3)
	REAL*8 beta(3)  !eignv_b(3),
	REAL*8 eignv_b_j, beta_j

	G=0.d0
	DO i=1,num_ex
		L_i = L(:,:,i)
		k_i = k(i)

		!Compute stress transformation, sn
		!s_trans_voigt = MATMUL(L(:,:,i),s_voigt)
		CALL mkl_jit_matmul(6,6,1,L_i,s_voigt,s_trans_voigt)

		! Compute eigen decomposition
		CALL VE6TEN(s_trans_voigt,s_trans)
		CALL DSYEVH3(s_trans, eigvect, eig_s_trans)

		!Compute diagonal matrices: Tn, |sp_hat|^(a-2) and their prodct, Tn_aux_sp
		CALL get_aux_Tn_diag(eig_s_trans,k_i,a,sign_vect,Tn_aux_sp)

		beta=0.d0
		DO j=1,3
			!Compute b_n eigenvalues
			eignv_b_j = Tn_aux_sp(j,j)*eig_s_trans(j)
			IF (DABS(eignv_b_j).EQ.0.d0) CYCLE
			beta_j = (1.d0/(1.d0-sign_vect(j)*k_i))**(a/(a-1.d0))
			G = G + beta_j*DABS(eignv_b_j)**(a/(a-1.d0))
		ENDDO

	ENDDO

	ENDSUBROUTINE compute_G
	! =============================================================

	SUBROUTINE local_eqs(s_voigt,a,k,L,m,F,localEQS)
	USE M_COMP_HOMOG, ONLY: num_ex
	USE M_jit_dgemm,  ONLY: mkl_jit_matmul

	IMPLICIT NONE
	REAL*8, INTENT(IN)  :: s_voigt(6),a,k(num_ex),L(6,6,num_ex),m
	REAL*8, INTENT(OUT) :: F, localEQS

	INTEGER i,j
	REAL*8 s_trans_voigt(6), s_trans(3,3),L_i(6,6),k_i
	REAL*8 eigvect(3,3), eig_s_trans(3),eignv_s_j


	F=0.d0
	DO i=1,num_ex
		L_i = L(:,:,i)
		k_i = k(i)

		!Compute stress transformation, sn
		!s_trans_voigt = MATMUL(L(:,:,i),s_voigt)
		CALL mkl_jit_matmul(6,6,1,L_i,s_voigt,s_trans_voigt)

		! Compute eigen decomposition
		CALL VE6TEN(s_trans_voigt,s_trans)
		CALL DSYEVH3(s_trans, eigvect, eig_s_trans)

		DO j=1,3
			!Compute sum
			eignv_s_j = eig_s_trans(j)
			F = F + (DABS(eignv_s_j)-k_i*eignv_s_j)**NINT(a)
		ENDDO

	ENDDO

	localEQS = m*F**(1.d0/a)

	ENDSUBROUTINE local_eqs

	! =============================================================
	SUBROUTINE dsymm_voigt66(T,T66)
	IMPLICIT NONE
	REAL*8, INTENT(IN) :: T(3,3)
	REAL*8, INTENT(OUT) :: T66(6,6)

	T66=0.d0
	T66(1,1)=T(1,1);
	T66(1,5)=T(1,3);
	T66(1,6)=T(1,2);

	T66(2,2)=T(2,2);
	T66(2,4)=T(2,3);
	T66(2,6)=T(1,2);

	T66(3,3)=T(3,3);
	T66(3,4)=T(2,3);
	T66(3,5)=T(1,3);

	T66(4,2)=T(2,3)/2.d0;
	T66(4,3)=T(2,3)/2.d0;
	T66(4,4)=(T(2,2)+T(3,3))/2.d0;
	T66(4,5)=T(1,2)/2.d0;
	T66(4,6)=T(1,3)/2.d0;

	T66(5,1)=T(1,3)/2.d0;
	T66(5,3)=T(1,3)/2.d0;
	T66(5,4)=T(1,2)/2.d0;
	T66(5,5)=(T(1,1)+T(3,3))/2.d0;
	T66(5,6)=T(2,3)/2.d0;

	T66(6,1)=T(1,2)/2.d0;
	T66(6,2)=T(1,2)/2.d0;
	T66(6,4)=T(1,3)/2.d0;
	T66(6,5)=T(2,3)/2.d0;
	T66(6,6)=(T(1,1)+T(2,2))/2.d0;
	ENDSUBROUTINE dsymm_voigt66

	! =============================================================
	SUBROUTINE VE5TE6(B5,A6)
	USE M_CONSTANTS,  ONLY: DSQRT2,DSQRT3
	IMPLICIT NONE
	REAL*8 A6(6),B5(5)
	
	A6(1)=(B5(1)+B5(2)/DSQRT3)/DSQRT2
	A6(2)=(B5(2)/DSQRT3-B5(1))/DSQRT2
	A6(3)=-DSQRT2/DSQRT3*B5(2)
	A6(4)=B5(3)/DSQRT2
	A6(5)=B5(4)/DSQRT2
	A6(6)=B5(5)/DSQRT2
	
	END
	
	SUBROUTINE TE6VE5(A6,B5)
	USE M_CONSTANTS,  ONLY: DSQRT2,DSQRT3
	IMPLICIT NONE
	REAL*8  A6(6),B5(5)
	
      REAL*8 A6_dev(6)
      
      A6_dev = A6

	 !guarantee dev
      A6_dev(1:3) = A6(1:3)-SUM(A6(1:3))/3.d0
      
      B5(1)=(A6_dev(1)-A6_dev(2))/DSQRT2
      B5(2)=-DSQRT3/DSQRT2*A6_dev(3)
      B5(3)=DSQRT2*A6_dev(4)
      B5(4)=DSQRT2*A6_dev(5)
      B5(5)=DSQRT2*A6_dev(6)
	
	END
	
	SUBROUTINE VE6TEN(B,A)
	IMPLICIT NONE
	REAL*8 A(3,3),B(6)
	
	A(1,1)=B(1)
	A(1,2)=B(6)
	A(1,3)=B(5)
	A(2,1)=B(6)
	A(2,2)=B(2)
	A(2,3)=B(4)
	A(3,1)=B(5)
	A(3,2)=B(4)
	A(3,3)=B(3)
	
	END
	
	SUBROUTINE TENVE6(A,B)
	IMPLICIT NONE
	REAL*8 A(3,3),B(6)
	B(1)=A(1,1)
	B(2)=A(2,2)
	B(3)=A(3,3)
	B(4)=A(2,3)
	B(5)=A(1,3)
	B(6)=A(1,2)
	
	END