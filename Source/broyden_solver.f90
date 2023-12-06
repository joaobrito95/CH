	!INCLUDE 'lapack.f90'
	SUBROUTINE broyden_solver(fun,x0,n,m,eps,maxiter,update_jac,x,i)
	! Solve system of nonlinear equations F(X) = 0 by Broyden's method:
	!  a quasi-Newton method which updates an approximate
	!  Jacobian at each new Newton step.
	! Usage:
	! m: Length of x
	! n: Length of F(x)
	! maxiter: maximum number of iterations of trial step calculation
	! update_jac: logical var, if .TRUE. the jacobian matrix is computed explicitly
	! eps(1): tolx, minimum length of last itertion step, i.e. delta(x)<eps(1)
	! eps(1): tolfun, inimum value for norm(f), i.e. ||f(x)||<eps(2)
	! x0(n): initial solution
	! x(n): converged solution of the problem

	!USE myMKL
	USE M_jit_dgemm,  ONLY: mkl_jit_matmul_sum
	IMPLICIT NONE
	EXTERNAL fun
	INTEGER, INTENT(IN)  :: n,m,maxiter
	REAL*8,  INTENT(IN)  :: x0(m), eps(2)
	LOGICAL, INTENT(IN)  :: update_jac
	REAL*8,  INTENT(OUT) :: x(m)
	REAL*8 f(n), normf,tolx,tolfun,J(n,m),rel_delta(1)
	REAL*8 T(n,m),b(n) !F(x) = 0, T(x)*x - b = 0
	INTEGER  i
	LOGICAL isNaN!, isInf
	INTEGER lda, info
	REAL*8 cholU(n,n)
	INTEGER nrhs, ldb
	REAL*8  dx(n)
	REAL*8  sumdx2,normdx, normx, normb, f_vect(n,1), dx_vect(n,1)
	INTEGER ldx!, iter
	INTEGER ipiv(n),lwork
	REAL*8, ALLOCATABLE :: work(:)
	INTEGER ntypfun,nfun(1)

	IF (m.NE.n) THEN
		WRITE(*,*) 'Not implemented!'
		READ(*,*)
	ENDIF
	tolx = DABS(eps(1))
	tolfun = DABS(eps(2))

	x=x0
	CALL fun(x,f,T,b)
	normf = DSQRT(SUM(f*f))
	normb = DSQRT(SUM(b*b))
	i=0
	IF (normf/normb.LT.tolfun) RETURN

	! Compute initial Jacobian matrix
	ntypfun=1
	nfun=m
	rel_delta = 1.D-6
	IF (update_jac) THEN
		IF (normf/normb.LT.1.D-2) THEN
			! forward differences
			CALL jacobi(fun,x,n,m,3,ntypfun,nfun,rel_delta,f,J)
		ELSE
			! central differences
			CALL jacobi(fun,x,n,m,1,ntypfun,nfun,rel_delta,f,J)
		ENDIF
	ELSE
		J=T
	ENDIF


	lda = n
	ldb=lda
	ldx = n
	nrhs=1
	lwork = n
	ALLOCATE(work(lwork))
	DO i=1,ABS(maxiter)
		isNaN = (J(1,1).NE.J(1,1))
		!isInf = ANY(DABS(J)>HUGE(J))
		IF (isNaN) THEN
			!WRITE(*,'(A,I3)') 'Could not compute jacobian at iteration', i
			x = x*DSQRT(-1.d0) !NaNs
			RETURN
		ENDIF

		!J = 0.5d0*(J+TRANSPOSE(J)) !guarantee symmetric jacobian

		!Solve system of equations (MKL Driver Routine)
		cholU = J
		dx=-f
		CALL dposv('U',n,nrhs,cholU,lda,dx,ldb,info) !

		IF (info.NE.0) THEN

			! try symmetric indefinite solver routine
			cholU = J
			dx=-f
			CALL dsysv('U',n,nrhs,cholU,lda,ipiv,dx,ldb,work,lwork,info)
			IF (info.NE.0) THEN
				!WRITE(*,'(A,I3,A,I3)') 'Error solving system of linear equations at iteration', i,&
				!    '. Error code info = ', info
				x = x*DSQRT(-1.d0) !NaNs
				RETURN
			ENDIF
		ENDIF

		x  = x+dx
		CALL fun(x,f,T,b)
		sumdx2 = SUM(dx*dx)
		normf = DSQRT(SUM(f*f))
		normdx = DSQRT(sumdx2)
		normx = DSQRT(SUM(x*x))

		!IF ((normf/normb.LT.tolfun).AND.(normdx/normx.LT.tolx)) RETURN
		IF (normf/normb.LT.tolfun) RETURN

		IF (update_jac) THEN
			IF (normdx/normx.LE.5.D-2.AND.normf/normb.LE.5.D-3.AND.maxiter.LT.0) THEN
				IF ((m.EQ.5).OR.(m.EQ.6).OR.(m.EQ.7)) THEN
					!J = J + MATMUL(f_vect,TRANSPOSE(dx_vect))/SUM(dx*dx)
					CALL mkl_jit_matmul_sum(m,1,m,'N','N',f,dx/sumdx2,J)
				ELSE
					f_vect(:,1) = f
					dx_vect(:,1) = dx
					J = J + MATMUL(f_vect,TRANSPOSE(dx_vect))/SUM(dx*dx)
				ENDIF
			ELSE
				!forward differences
				CALL jacobi(fun,x,n,m,3,ntypfun,nfun,rel_delta,f,J)
			ENDIF
		ELSE
			IF (normdx/normx.LE.1.D-2.OR.normf/normb.LE.1.D-3) THEN
				IF ((m.EQ.5).OR.(m.EQ.6).OR.(m.EQ.7)) THEN
					!J = J + MATMUL(f_vect,TRANSPOSE(dx_vect))/SUM(dx*dx)
					CALL mkl_jit_matmul_sum(m,1,m,'N','N',f,dx/sumdx2,J)
				ELSE
					f_vect(:,1) = f
					dx_vect(:,1) = dx
					J = J + MATMUL(f_vect,TRANSPOSE(dx_vect))/SUM(dx*dx)
				ENDIF
			ELSE
				J=T
			ENDIF
		ENDIF

	ENDDO
	IF (normf/normb.LT.5d0*tolfun) RETURN
	x = x*DSQRT(-1.d0) !NaNs
	!WRITE(*,'(A,I3)') 'No convergence in Broyden solver, niter > ',maxiter

	ENDSUBROUTINE broyden_solver
	! ========================================================
	SUBROUTINE broyden_solver2(fun,x0,n,m,ntypfun,nfun,eps,maxiter,update_jac,rel_delta0,&
		jacmethod,J_init_flag,J_init,soltype,nf1,nf2,x,f1,f2,i)
	! Solve system of nonlinear equations F(X) = 0 by Broyden's method:
	!  a quasi-Newton method which updates an approximate
	!  Jacobian at each new Newton step.
	! Usage:
	! m: Length of x
	! n: Length of F(x)
	! ntryfun: number of dissimilar equations (e.g. w/ diferent physical dimentions)
	! nfum(ntryfun): number of equations for each type of x(:)
	! maxiter: maximum number of iterations of trial step calculation
	! rel_delta: used in the computation of the jacobian matrix
	! update_jac: logical var, if .TRUE. the jacobian matrix is computed explicitly
	! eps(1): tolx, minimum length of last itertion step, i.e. delta(x)<eps(1)
	! eps(1): tolfun, inimum value for norm(f), i.e. ||f(x)||<eps(2)
	! x0(n): initial solution
	! x(n): converged solution of the problem
	! J_init_flag: boolean flag indicatin if initial jacobian matrix is given
	! J_init(n,m): initial jacobian matrix (if J_init_flag)
	! soltype: flag to indicate which linear solve to employ
	! nf1, nf2: dimention of the auxiliary output vectors of fun
	! f1,f2: auxiliary output vectors of fun
	! i: number of iterations at exit

	USE M_jit_dgemm,  ONLY: mkl_jit_matmul_sum
	IMPLICIT NONE
	EXTERNAL fun
	INTEGER, INTENT(IN)  :: n,m,maxiter,nf1,nf2,ntypfun,nfun(ntypfun),jacmethod
	REAL*8,  INTENT(IN)  :: x0(m), eps(2), rel_delta0(ntypfun)
	LOGICAL, INTENT(INOUT)  :: update_jac,J_init_flag
	REAL*8,  INTENT(OUT) :: x(m),f1(nf1),f2(nf2)
	REAL*8,  INTENT(INOUT) :: J_init(n,m)
	INTEGER, INTENT(INOUT) :: soltype
	REAL*8 f(n),tolx,tolfun,J(n,m)
	INTEGER  i, max_iter, maxnotconvitr, notconvitr
	LOGICAL isNaN, testconv!, isInf
	INTEGER lda, info
	REAL*8 cholU(n,n),LU(n,n)
	INTEGER nrhs, ldb
	REAL*8  dx(n)
	REAL*8  sumdx2,normdx(ntypfun), normx(ntypfun), f_vect(n,1), dx_vect(n,1)
	INTEGER ldx
	REAL*8 best_sol_x(m),f1_best(nf1),f2_best(nf2)
	LOGICAL forceIter
	REAL*8 rng_delta!,rel_delta_org(ntypfun)
	REAL*8  normf, normf1, normf2, min_normf
	INTEGER lwork, ipiv(n)
	REAL*8, ALLOCATABLE :: work2(:)

	INTEGER  jj,mmax,max_prev,idx_aux(ntypfun,2)
	REAL*8  rtnormdxxall,rel_delta(ntypfun)
	REAL*8 jac_eps
	REAL*8 lastdx2(n),lastdx1(n),toldx_diff,norm_dx_diff1(ntypfun),norm_dx_diff2(ntypfun), frac_dx(ntypfun)
	REAL*8 normlast2,normlast1
	LOGICAL bounce_test
	INTEGER cnt_bounce(ntypfun), convitr

	IF (m.NE.n) THEN
		WRITE(*,*) 'Not implemented!'
		READ(*,*)
	ENDIF
	tolx = eps(1)
	tolfun = eps(2)
	IF (tolfun.LT.0.d0.OR.tolx.LT.0.d0) THEN
		forceIter=.true.
		tolfun=DABS(tolfun)
		tolx=DABS(tolx)
	ELSE
		forceIter=.false.
	ENDIF

	IF (soltype.LT.1.OR.soltype.GT.3) soltype=3

	bounce_test=.false.
	cnt_bounce=0
	frac_dx=1.d0
	convitr=0

	testconv = .false.

	IF (maxiter.LT.0) THEN
		max_iter = ABS(maxiter)
		testconv = .true.
		maxnotconvitr = (CEILING(DFLOAT(max_iter)/2.d0))
		notconvitr=0
	ELSEIF (maxiter.EQ.0) THEN
		max_iter=25
	ELSE
		max_iter = maxiter
	ENDIF
	min_normf = HUGE(min_normf)
	best_sol_x = DSQRT(-1.d0)	!nan

	!rel_delta_org=rel_delta0

	mmax = MAXVAL(nfun)
	IF (ntypfun.NE.1) THEN
		! system contains more than one physical quatity
		idx_aux(1,:) = (/1,nfun(1)/)
		DO i=2,ntypfun
			max_prev = (idx_aux(i-1,2))
			idx_aux(i,1) = max_prev+1
			idx_aux(i,2) = max_prev+nfun(i)
		ENDDO
	ELSE
		idx_aux(1,:)= (/1,nfun(1)/)
	ENDIF


	! start computations
	x=x0
	CALL fun(x,f,f1,f2)
	normf = DSQRT(SUM(f*f))
	normf1 = DSQRT(SUM(f1*f1))
	normf2 = DSQRT(SUM(f2*f2))
	IF (normf.NE.normf) THEN
		CALL fun(x,f,f1,f2)
		x = x*DSQRT(-1.d0)
		RETURN
	ENDIF
	i=0
	IF (normf1.EQ.0.d0) normf1=1.d0
	IF (.not.forceIter) THEN
		IF (normf/normf1.LT.tolfun) RETURN
	ENDIF

	IF (J_init_flag) THEN
		J = J_init
	ELSE
		! Compute initial Jacobian matrix
		IF (SUM(rel_delta0).EQ.0.d0) THEN
			! MKL jacobi
			jac_eps = 1.d-6
			CALL mkl_djacobi(n,m,fun,x,J,jac_eps)
		ELSE
			! User made (forward differences)
			CALL jacobi(fun,x,n,m,jacmethod,ntypfun,nfun,rel_delta0,f,J)
		ENDIF
	ENDIF

	lda = n
	ldb=lda
	ldx = n
	nrhs=1
	lwork = n
	ALLOCATE(work2(lwork))

	IF (ANY(rel_delta0.LT.0.d0)) THEN
		max_iter=2*max_iter
	ENDIF

	i=0
	DO WHILE (.true.)
		i=i+1
		!DO i=1,max_iter
		isNaN = (J(1,1).NE.J(1,1))
		!isInf = ANY(DABS(J)>HUGE(J))
		IF (isNaN) THEN
			!WRITE(*,'(A,I3)') 'Could not compute jacobian at iteration', i
			x = x*DSQRT(-1.d0) !NaNs
			RETURN
		ENDIF

		IF (soltype.NE.3.AND.ntypfun.EQ.1) J = 0.5d0*(J+TRANSPOSE(J))


		IF (i.eq.2) THEN
			lastdx1=dx
			toldx_diff=5.d-2
			!bounce_test=.false.
			!cnt_bounce=0
			normlast1=0.d0
			normlast2=0.d0
			IF (bounce_test) normlast1=normf
		ELSEIF (i.gt.2) THEN
			lastdx2=lastdx1
			lastdx1=dx
			IF (bounce_test) normlast2=normlast1
			IF (bounce_test) normlast1=normf
		ENDIF


		!Solve system of equations (MKL Driver Routine)
		IF (soltype.EQ.1) THEN
			cholU = J
			dx=-f
			! symmetric positive-definite solver
			CALL dposv('U',n,nrhs,cholU,lda,dx,ldb,info) !

			IF (info.NE.0) THEN
				! try symmetric indefinite solver routine
				cholU = J
				dx=-f
				call dsysv('U',n,nrhs,cholU,lda,ipiv,dx,ldb,work2,lwork,info)
			ENDIF
		ELSEIF (soltype.EQ.2) THEN
			! symmetric indefinite solver routine
			cholU = J
			dx=-f
			call dsysv('U',n,nrhs,cholU,lda,ipiv,dx,ldb,work2,lwork,info)

		ELSE !soltype.EQ.3
			!general matrix solver
			LU = J
			dx=-f
			CALL dgesv(n,nrhs,LU,lda,ipiv,dx,ldb,info)
		ENDIF
		IF (info.NE.0) THEN
			!WRITE(*,'(A,I3,A,I3)') 'Error solving system of linear equations at iteration', i,&
			!   '. Error code info = ', info
			x = x*DSQRT(-1.d0) !NaNs
			f1 = DSQRT(-1.d0)*f1
			f2 = DSQRT(-1.d0)*f2
			RETURN
		ENDIF

		!        WRITE(*,*) ' '
		!        IF (SIZE(J,1).EQ.6.OR.SIZE(J,1).EQ.7) THEN
		!            DO jj=1,6
		!                WRITE(*,37) J(jj,1),J(jj,2),J(jj,3), &
		!                    J(jj,4),J(jj,5),J(jj,6)
		!            ENDDO
		!37          FORMAT(5X,7F15.4)
		!        ENDIF
		!        WRITE(*,*) ' '
		IF (ANY(rel_delta0.LT.0.d0).AND.i.GE.2) THEN
			IF (normf/normf1.GT.1.d-2) THEN
				CALL RANDOM_NUMBER(rng_delta)
				rng_delta=1.d0+(3.d0*rng_delta)
				!rng_delta=MAX(rng_delta,1.d0)
				dx=dx*(1.d0/rng_delta)
			ELSE
				!update_jac=.false.
			ENDIF

		ENDIF

		IF (bounce_test)	THEN
			DO jj=1,ntypfun
				x(idx_aux(jj,1):idx_aux(jj,2))= x(idx_aux(jj,1):idx_aux(jj,2))+&
					frac_dx(jj)*dx(idx_aux(jj,1):idx_aux(jj,2))
			ENDDO
			!bounce_test=.false.
		ELSE
			x  = x+dx
		ENDIF

		CALL fun(x,f,f1,f2)
		normf1 = DSQRT(SUM(f1*f1))
		IF (normf1.EQ.0.d0) normf1=1.d0

		sumdx2 = SUM(dx*dx)
		normf = DSQRT(SUM(f*f))
		IF (normf.NE.normf) THEN
			x = x*DSQRT(-1.d0)
			RETURN
		ENDIF

		IF (ntypfun.EQ.1) THEN
			normx = DSQRT(SUM(x*x))
			normdx = DSQRT(sumdx2)
		ELSE
			DO jj=1,ntypfun
				normx(jj) = DSQRT(SUM(x(idx_aux(jj,1):idx_aux(jj,2))*x(idx_aux(jj,1):idx_aux(jj,2))))
				normdx(jj) = DSQRT(SUM(dx(idx_aux(jj,1):idx_aux(jj,2))*dx(idx_aux(jj,1):idx_aux(jj,2))))
			ENDDO
		ENDIF
		rtnormdxxall = SUM(normdx/normx)

		J_init = J

		IF ((normf/normf1.LT.tolfun).AND.(rtnormdxxall.LT.tolx)) RETURN

		IF (normf/normf1.LT.tolfun) THEN
			IF (normf.LT.min_normf) THEN
				best_sol_x=x
				f1_best=f1
				f2_best=f2
				!min_normf = normf
			ENDIF
		ENDIF

		IF (update_jac) THEN
			IF (rtnormdxxall.LE.1.D-2.AND.normf/normf1.LE.1.d-4.AND.ALL(rel_delta0.GT.0.d0)) THEN !.AND.ALL(rel_delta0.GT.0.d0)
				IF ((m.EQ.5).OR.(m.EQ.6).OR.(m.EQ.7)) THEN !OR m.EQ.6.OR.m.EQ.7
					!J = J + MATMUL(f_vect,TRANSPOSE(dx_vect))/SUM(dx*dx)
					CALL mkl_jit_matmul_sum(m,1,m,'N','N',f,dx/sumdx2,J)
				ELSE
					f_vect(:,1) = f
					dx_vect(:,1) = dx
					J = J + MATMUL(f_vect,TRANSPOSE(dx_vect))/SUM(dx*dx)
				ENDIF
			ELSE
				IF (SUM(rel_delta0).EQ.0.d0) THEN
					! MKL jacobi
					J=0.d0
					jac_eps=1.d-6
					CALL mkl_djacobi(n,m,fun,x,J,jac_eps)
				ELSE
					! User made (forward differences)
					IF (i.eq.1) rel_delta=rel_delta0
					IF (ANY(rel_delta0.LT.0.d0)) THEN !(normf/normf1.GT.0.10d0.AND.
						!rel_delta=DABS(rel_delta0)
						CALL RANDOM_NUMBER(rng_delta)
						DO jj=1,ntypfun
							!rel_delta(jj)=-MINVAL(DABS(x(idx_aux(jj,1):idx_aux(jj,2))))*1.d-2
							!rel_delta(jj)=normdx(jj)
							rel_delta(jj)=DABS(rel_delta0(jj))*(1000.d0*(rng_delta+0.05d0))
							rel_delta(jj)=MIN(rel_delta(jj),1.d-3)
						ENDDO
					ENDIF
					!IF (normf/normf1.LT.0.10d0.AND.ANY(rel_delta0.LT.0.d0)) THEN
					!	rel_delta=DABS(rel_delta0)
					!ENDIF

					CALL jacobi(fun,x,n,m,jacmethod,ntypfun,nfun,rel_delta,f,J)
				ENDIF
			ENDIF
		ELSE
			IF ((m.EQ.5).OR.(m.EQ.6).OR.(m.EQ.7)) THEN
				!J = J + MATMUL(f_vect,TRANSPOSE(dx_vect))/SUM(dx*dx)
				CALL mkl_jit_matmul_sum(m,1,m,'N','N',f,dx/sumdx2,J)
			ELSE
				f_vect(:,1) = f
				dx_vect(:,1) = dx
				J = J + MATMUL(f_vect,TRANSPOSE(dx_vect))/SUM(dx*dx)
			ENDIF
		ENDIF


		IF (testconv.AND.i.GT.1) THEN
			IF (normf.GE.min_normf) THEN
				notconvitr = notconvitr +1
			ENDIF
			IF (notconvitr.GE.maxnotconvitr) THEN
				EXIT
			ENDIF

		ENDIF
		min_normf = MIN(normf,min_normf)

		! Scale dx
		IF (i.gt.2) THEN
			DO jj=1,ntypfun
				norm_dx_diff1(jj)=DSQRT(SUM((dx(idx_aux(jj,1):idx_aux(jj,2))-lastdx1(idx_aux(jj,1):idx_aux(jj,2)))**2))
				norm_dx_diff2(jj)=DSQRT(SUM((dx(idx_aux(jj,1):idx_aux(jj,2))-lastdx2(idx_aux(jj,1):idx_aux(jj,2)))**2))
				IF ((norm_dx_diff2(jj)/normdx(jj).lt.toldx_diff).AND.(norm_dx_diff1(jj)/normdx(jj).gt.toldx_diff))	THEN
					cnt_bounce(jj)=cnt_bounce(jj)+1
				ELSE
					cnt_bounce(jj)=0
				ENDIF

				IF (cnt_bounce(jj).ge.2) THEN
					bounce_test=.true.
					testconv=.false.
					frac_dx(jj)=frac_dx(jj)*0.20d0

					cnt_bounce(jj)=0
					notconvitr=0
				ENDIF
			ENDDO
			IF (bounce_test) THEN
				IF (normlast1.lt.normlast2.and.normf.lt.normlast1) THEN
					convitr=convitr+1
					IF  (convitr.gt.2.AND.ABS(maxiter).eq.max_iter) THEN
						max_iter=2*max_iter
					ENDIF
					IF (normf/normf1.lt.1.d0.AND.SUM(frac_dx-1.d0).NE.0.d0) THEN
						frac_dx=1.d0
						max_iter=max_iter+NINT(0.5d0*DABS(DFLOAT(maxiter)))
					ENDIF
				ELSE
					IF ((normlast2.ne.0.d0).AND.ALL(frac_dx.ne.1.d0).and.normf.gt.normlast1) THEN
						frac_dx=frac_dx/2.0d0
						frac_dx=MAX(frac_dx,0.05d0)
					ENDIF

				ENDIF
			ENDIF
		ENDIF

		IF (bounce_test) THEN
			max_iter=MIN(max_iter,NINT(2.5d0*(DABS(DFLOAT(maxiter)))))
		ENDIF

		IF (i.ge.max_iter) EXIT

	ENDDO
	IF (best_sol_x(1).EQ.best_sol_x(1)) THEN
		x=best_sol_x
		f1=f1_best
		f2=f2_best
	ELSE
		x = x*DSQRT(-1.d0) !NaNs
		f1 = DSQRT(-1.d0)*f1
		f2 = DSQRT(-1.d0)*f2
		!WRITE(*,'(A,I3)') 'No convergence in Broyden solver, niter > ',maxiter
	ENDIF
	!WRITE(*,'(A,I3)') 'No convergence in Broyden solver, niter > ',maxiter

	ENDSUBROUTINE broyden_solver2

	! =========================================================
	SUBROUTINE jacobi(fun,x,n,m,method,ntypfun,nfun,rel_delta,f0,J)
	! Numerical jacobian matrix of funtion 'fun' in point 'x'
	! Usage:
	! m: Length of x
	! n: Length of F(x)
	! rel_delta: relative step size for computing numerical derivatives. If
	!  negative (flag) then an absolute (constant) value, is used.
	! method: finite differences method: 1: central, 2: backward, 3:forward
	! f0: funtion value in point x, f(x). Used to avoid recompute this value
	!       in the case of backward and forward schemes.
	IMPLICIT NONE
	EXTERNAL fun
	INTEGER, INTENT(IN)  :: n,m,method,ntypfun,nfun(ntypfun)
	REAL*8,  INTENT(IN)  :: x(m), rel_delta(ntypfun),f0(n)
	REAL*8,  INTENT(OUT) :: J(n,m)
	REAL*8   delta(m), dx(m),col(n),normx,f_dx_back(n),f_dx_forw(n)
	REAL*8 T(n,m),b(n) !F(x) = 0, T(x)*x - b = 0
	REAL*8 normxvect(ntypfun)
	INTEGER  i, max_prev, mmax,idx_aux(ntypfun,2),eye(n),jj
	REAL*8 mindelta

	eye=1.d0
	mmax = MAXVAL(nfun)
	IF (mmax.EQ.m) THEN
		normx = DSQRT(SUM(x*x))
		IF (normx.NE.0.D0.AND.rel_delta(1).GT.0.d0) THEN
			delta = rel_delta(1)*normx
		ELSE
			IF (normx.EQ.0.D0) THEN
				delta = DABS(rel_delta(1))
			ELSE
				delta = DABS(rel_delta(1))!*x
			ENDIF
		ENDIF
		IF (method.EQ.4) THEN
			IF (normx.EQ.0.D0) THEN
				delta = DABS(rel_delta(1))
			ELSE
				!delta = 0.d0
				mindelta=MINVAL(DABS(rel_delta(1)*x),1,x.NE.0.d0)
				DO jj=1,m
					IF (x(jj).NE.0.d0) THEN
						delta(jj)=DABS(rel_delta(1)*(x(jj)))
					ELSE
						delta(jj)=mindelta
					ENDIF
				ENDDO
			ENDIF
		ENDIF

	ELSE
		! system contains more than one physical quatity
		idx_aux(1,:) = (/1,nfun(1)/)
		DO i=2,ntypfun
			max_prev = (idx_aux(i-1,2))
			idx_aux(i,1) = max_prev+1
			idx_aux(i,2) = max_prev+nfun(i)
		ENDDO

		DO i=1,ntypfun
			normxvect(i) = DSQRT(SUM(x(idx_aux(i,1):idx_aux(i,2))*x(idx_aux(i,1):idx_aux(i,2))))
			IF (normxvect(i).NE.0.D0.AND.rel_delta(i).GT.0.d0) THEN
				delta(idx_aux(i,1):idx_aux(i,2)) = rel_delta(i)*normxvect(i)
			ELSE
				delta(idx_aux(i,1):idx_aux(i,2)) = DABS(rel_delta(i))
			ENDIF
		ENDDO
	ENDIF


	IF (method.NE.1.AND.SUM(DABS(f0)).EQ.0.d0) THEN
		!f0 is not yet computed
		CALL fun(x,f0,T,b)
	ENDIF


	J=0.d0
	DO i=1,m
		dx = 0.d0
		SELECT CASE (method)
		CASE(1)
			!central
			dx(i) = delta(i)/2.d0
			CALL fun(x-dx,f_dx_back,T,b)
			CALL fun(x+dx,f_dx_forw,T,b)
			col = (f_dx_forw-f_dx_back)/delta(i)
		CASE(2)
			!backward
			dx(i) = delta(i)
			CALL fun(x-dx,f_dx_back,T,b)
			col = (f0-f_dx_back)/delta(i)
			CASE DEFAULT
			!forward
			dx(i) = delta(i)
			CALL fun(x+dx,f_dx_forw,T,b)
			col = (f_dx_forw-f0)/delta(i)
		ENDSELECT
		J(:,i) = col
	ENDDO

	ENDSUBROUTINE jacobi

	! ========================================================
	SUBROUTINE broyden_solver3(fun,x0,n,m,ntypfun,nfun,eps,maxiter,update_jac,rel_delta0,&
		jacmethod,J_init_flag,J_init,soltype,nf1,nf2,bounds,lb,ub,x,f1,f2,i)
	! Solve system of nonlinear equations F(X) = 0 with bounds in X by Broyden's method:
	!  a quasi-Newton method which updates an approximate. lb<=X<=ub
	!  Jacobian at each new Newton step.
	! Usage:
	! m: Length of x
	! n: Length of F(x)
	! ntryfun: number of dissimilar equations (e.g. w/ diferent physical dimentions)
	! nfum(ntryfun): number of equations for each type of x(:)
	! maxiter: maximum number of iterations of trial step calculation
	! rel_delta: used in the computation of the jacobian matrix
	! update_jac: logical var, if .TRUE. the jacobian matrix is computed explicitly
	! eps(1): tolx, minimum length of last itertion step, i.e. delta(x)<eps(1)
	! eps(1): tolfun, inimum value for norm(f), i.e. ||f(x)||<eps(2)
	! x0(n): initial solution
	! x(n): converged solution of the problem
	! J_init_flag: boolean flag indicatin if initial jacobian matrix is given
	! J_init(n,m): initial jacobian matrix (if J_init_flag)
	! soltype: flag to indicate which linear solve to employ
	! nf1, nf2: dimention of the auxiliary output vectors of fun
	! f1,f2: auxiliary output vectors of fun
	! i: number of iterations at exit

	USE M_jit_dgemm,  ONLY: mkl_jit_matmul_sum
	IMPLICIT NONE
	EXTERNAL fun
	INTEGER, INTENT(IN)  :: n,m,maxiter,nf1,nf2,ntypfun,nfun(ntypfun),jacmethod
	REAL*8,  INTENT(IN)  :: x0(m), eps(2), rel_delta0(ntypfun)
	REAL*8, INTENT(IN)   :: lb(m),ub(m)
	LOGICAL, INTENT(IN)  :: update_jac,J_init_flag,bounds
	REAL*8,  INTENT(OUT) :: x(m),f1(nf1),f2(nf2)
	REAL*8,  INTENT(INOUT) :: J_init(n,m)
	INTEGER, INTENT(INOUT) :: soltype
	REAL*8 f(n),tolx,tolfun,J(n,m)
	INTEGER  i, max_iter, maxnotconvitr, notconvitr
	LOGICAL isNaN, testconv!, isInf
	INTEGER lda, info
	REAL*8 cholU(n,n),LU(n,n)
	INTEGER nrhs, ldb
	REAL*8  dx(n)
	REAL*8  sumdx2,normdx(ntypfun), normx(ntypfun), f_vect(n,1), dx_vect(n,1)
	INTEGER ldx
	REAL*8 rng_delta!,rel_delta_org(ntypfun)
	REAL*8  normf, normf1, normf2, min_normf
	INTEGER lwork, ipiv(n)
	REAL*8, ALLOCATABLE :: work2(:)

	INTEGER  jj,mmax,max_prev,idx_aux(ntypfun,2)
	REAL*8  rtnormdxxall,rel_delta(ntypfun)
	REAL*8 jac_eps

	REAL*8 xnew(m)
	LOGICAL lo(m),hi(m)
	INTEGER idx, k
	REAL*8 best_sol_x(m),f1_best(nf1),f2_best(nf2)
	LOGICAL forceIter

	IF (m.NE.n) THEN
		WRITE(*,*) 'Not implemented!'
		READ(*,*)
	ENDIF
	tolx = eps(1)
	tolfun = eps(2)
	IF (tolfun.LT.0.d0.OR.tolx.LT.0.d0) THEN
		forceIter=.true.
		tolfun=DABS(tolfun)
		tolx=DABS(tolx)
	ELSE
		forceIter=.false.
	ENDIF
	IF (soltype.LT.1.OR.soltype.GT.3) soltype=3

	testconv = .false.

	IF (maxiter.LT.0) THEN
		max_iter = ABS(maxiter)
		testconv = .true.
		maxnotconvitr = (CEILING(DFLOAT(max_iter)/2.d0))
		notconvitr=0
	ELSEIF (maxiter.EQ.0) THEN
		max_iter=25
	ELSE
		max_iter = maxiter
	ENDIF
	min_normf = HUGE(min_normf)
	best_sol_x = DSQRT(-1.d0)	!nan

	!rel_delta_org=rel_delta0

	mmax = MAXVAL(nfun)
	IF (ntypfun.NE.1) THEN
		! system contains more than one physical quatity
		idx_aux(1,:) = (/1,nfun(1)/)
		DO i=2,ntypfun
			max_prev = (idx_aux(i-1,2))
			idx_aux(i,1) = max_prev+1
			idx_aux(i,2) = max_prev+nfun(i)
		ENDDO
	ENDIF


	! TEST f(x0) for premature EXIT
	IF (bounds) THEN
		x=x0
		lo = (x < lb);
		hi = (x > ub);

		IF (.not.ANY(lo.OR.hi)) THEN
			xnew=x
		ELSE
			xnew = MERGE(lb,x,lo)
			xnew = MERGE(ub,xnew,hi)
		ENDIF
		x = xnew
	ELSE
		x=x0
	ENDIF
	i=0
	CALL fun(x,f,f1,f2)
	normf = DSQRT(SUM(f*f))
	normf1 = DSQRT(SUM(f1*f1))
	normf2 = DSQRT(SUM(f2*f2))
	IF (normf1.EQ.0.d0) normf1=1.d0
	IF (.not.forceIter) THEN
		IF (normf/normf1.LT.tolfun) RETURN
	ENDIF
	IF (normf.NE.normf) THEN
		CALL fun(x,f,f1,f2)
		x = x*DSQRT(-1.d0)
		RETURN
	ENDIF

	IF (J_init_flag) THEN
		J = J_init
	ELSE
		! Compute initial Jacobian matrix
		IF (SUM(rel_delta0).EQ.0.d0) THEN
			! MKL jacobi
			jac_eps = 1.0d-6
			CALL mkl_djacobi(n,m,fun,x,J,jac_eps)
		ELSE
			! User made (forward differences)
			CALL jacobi(fun,x,n,m,jacmethod,ntypfun,nfun,rel_delta0,f,J)
		ENDIF
	ENDIF

	lda = n
	ldb=lda
	ldx = n
	nrhs=1
	lwork = n
	ALLOCATE(work2(lwork))

	DO i=1,max_iter
		isNaN = (J(1,1).NE.J(1,1))
		!isInf = ANY(DABS(J)>HUGE(J))
		IF (isNaN) THEN
			!WRITE(*,'(A,I3)') 'Could not compute jacobian at iteration', i
			x = x*DSQRT(-1.d0) !NaNs
			RETURN
		ENDIF

		IF (soltype.NE.3.AND.ntypfun.EQ.1) J = 0.5d0*(J+TRANSPOSE(J))


		!Solve system of equations (MKL Driver Routine)
		IF (soltype.EQ.1) THEN
			cholU = J
			dx=-f
			! symmetric positive-definite solver
			CALL dposv('U',n,nrhs,cholU,lda,dx,ldb,info) !

			IF (info.NE.0) THEN
				! try symmetric indefinite solver routine
				cholU = J
				dx=-f
				call dsysv('U',n,nrhs,cholU,lda,ipiv,dx,ldb,work2,lwork,info)
			ENDIF
		ELSEIF (soltype.EQ.2) THEN
			! symmetric indefinite solver routine
			cholU = J
			dx=-f
			call dsysv('U',n,nrhs,cholU,lda,ipiv,dx,ldb,work2,lwork,info)

		ELSE !soltype.EQ.3
			!general matrix solver
			LU = J
			dx=-f
			CALL dgesv(n,nrhs,LU,lda,ipiv,dx,ldb,info)
		ENDIF
		IF (info.NE.0) THEN
			!WRITE(*,'(A,I3,A,I3)') 'Error solving system of linear equations at iteration', i,&
			!   '. Error code info = ', info
			x = x*DSQRT(-1.d0) !NaNs
			f1 = DSQRT(-1.d0)*f1
			f2 = DSQRT(-1.d0)*f2
			RETURN
		ENDIF


		!		 WRITE(*,*) ' '
		!		 IF (SIZE(J,1).EQ.6.OR.SIZE(J,1).EQ.7) THEN
		!		 DO jj=1,6
		!			 WRITE(*,37) J(jj,1),J(jj,2),J(jj,3), &
		!				           J(jj,4),J(jj,5),J(jj,6)
		!		 ENDDO
		!37		 FORMAT(5X,7F15.4)
		!		 ENDIF
		!		 WRITE(*,*) ' '

		! Bounds
		IF (bounds) THEN
			xnew = x + dx;
			lo = (xnew < lb);
			hi = (xnew > ub);
			DO jj=1,20 !max num of truncations
				IF (.not.ANY(lo.OR.hi)) THEN
					EXIT
				ENDIF
				idx = 1
				DO k=1,COUNT(lo)
					idx=FINDLOC(lo(idx:m),.true.,1)
					dx(idx) = (x(idx)-lb(idx))/2.d0;
				ENDDO
				idx = 1
				DO k=1,COUNT(hi)
					idx=FINDLOC(hi(idx:m),.true.,1)
					dx(idx) = (ub(idx)-x(idx))/2.d0;
				ENDDO

				xnew = x+dx;
				lo = (xnew < lb);
				hi = (xnew > ub);
			ENDDO
			DO jj=1,n
				IF (lb(jj).eq.ub(jj)) dx(jj)=0.d0
			ENDDO
		ENDIF

		IF (ANY(rel_delta0.LT.0.d0).AND.i.GE.2) THEN
			IF (normf/normf1.GT.1.d-2) THEN
				CALL RANDOM_NUMBER(rng_delta)
				rng_delta=1.d0+(3.d0*rng_delta)
				!rng_delta=MAX(rng_delta,1.d0)
				dx=dx*(1.d0/rng_delta)
			ELSE
				!update_jac=.false.
			ENDIF
		ENDIF


		x  = x+dx
		CALL fun(x,f,f1,f2)
		normf1 = DSQRT(SUM(f1*f1))
		IF (normf1.EQ.0.d0) normf1=1.d0

		sumdx2 = SUM(dx*dx)
		normf = DSQRT(SUM(f*f))
		IF (normf.NE.normf) THEN
			x = x*DSQRT(-1.d0)
			RETURN
		ENDIF

		IF (ntypfun.EQ.1) THEN
			normx = DSQRT(SUM(x*x))
			normdx = DSQRT(sumdx2)
		ELSE
			DO jj=1,ntypfun
				normx(jj) = DSQRT(SUM(x(idx_aux(jj,1):idx_aux(jj,2))*x(idx_aux(jj,1):idx_aux(jj,2))))
				normdx(jj) = DSQRT(SUM(dx(idx_aux(jj,1):idx_aux(jj,2))*dx(idx_aux(jj,1):idx_aux(jj,2))))
			ENDDO
		ENDIF
		rtnormdxxall = SUM(PACK(normdx/normx,normx.NE.0.d0))

		J_init = J

		IF ((normf/normf1.LT.tolfun).AND.(rtnormdxxall.LT.tolx)) RETURN

		IF (normf/normf1.LT.tolfun) THEN
			IF (normf.LT.min_normf) THEN
				best_sol_x=x
				f1_best=f1
				f2_best=f2
				!min_normf = normf
			ENDIF
		ENDIF

		IF (update_jac) THEN
			IF (rtnormdxxall.LE.1.D-2.AND.normf/normf1.LE.1.D-4.AND.ALL(rel_delta0.GT.0.d0)) THEN
				IF ((m.EQ.5).OR.(m.EQ.6).OR.(m.EQ.7)) THEN
					!J = J + MATMUL(f_vect,TRANSPOSE(dx_vect))/SUM(dx*dx)
					CALL mkl_jit_matmul_sum(m,1,m,'N','N',f,dx/sumdx2,J)
				ELSE
					f_vect(:,1) = f
					dx_vect(:,1) = dx
					J = J + MATMUL(f_vect,TRANSPOSE(dx_vect))/SUM(dx*dx)
					! Create 6x1 1x6 and 7x1 1x7 jit dgemms
				ENDIF
			ELSE
				IF (SUM(rel_delta0).EQ.0.d0) THEN
					! MKL jacobi
					jac_eps = 1.0d-6
					J=0.d0
					CALL mkl_djacobi(n,m,fun,x,J,jac_eps)
				ELSE
					! User made (forward differences)
					IF (i.eq.1) rel_delta=rel_delta0
					IF (ANY(rel_delta0.LT.0.d0)) THEN !(normf/normf1.GT.0.10d0.AND.
						!rel_delta=DABS(rel_delta0)

						DO jj=1,ntypfun
							!rel_delta(jj)=-MINVAL(DABS(x(idx_aux(jj,1):idx_aux(jj,2))))*1.d-2
							!rel_delta(jj)=normdx(jj)
							CALL RANDOM_NUMBER(rng_delta)
							rel_delta(jj)=DABS(rel_delta0(jj))*(1000.d0*(rng_delta+0.05d0))
							rel_delta(jj)=MIN(rel_delta(jj),1.d-3)
						ENDDO
					ENDIF

					CALL jacobi(fun,x,n,m,jacmethod,ntypfun,nfun,rel_delta,f,J)
				ENDIF
			ENDIF

		ELSE
			IF ((m.EQ.5).OR.(m.EQ.6).OR.(m.EQ.7)) THEN
				!J = J + MATMUL(f_vect,TRANSPOSE(dx_vect))/SUM(dx*dx)
				CALL mkl_jit_matmul_sum(m,1,m,'N','N',f,dx/sumdx2,J)
			ELSE
				f_vect(:,1) = f
				dx_vect(:,1) = dx
				J = J + MATMUL(f_vect,TRANSPOSE(dx_vect))/SUM(dx*dx)
			ENDIF
		ENDIF

		IF (testconv.AND.i.GT.1) THEN
			IF (normf.GE.min_normf) THEN
				notconvitr = notconvitr +1
			ENDIF
			IF (notconvitr.GE.maxnotconvitr) THEN
				EXIT
			ENDIF

		ENDIF

		min_normf = MIN(normf,min_normf)

	ENDDO
	IF (best_sol_x(1).EQ.best_sol_x(1)) THEN
		x=best_sol_x
		f1=f1_best
		f2=f2_best
	ELSE
		x = x*DSQRT(-1.d0) !NaNs
		f1 = DSQRT(-1.d0)*f1
		f2 = DSQRT(-1.d0)*f2
		!WRITE(*,'(A,I3)') 'No convergence in Broyden solver, niter > ',maxiter
	ENDIF

	ENDSUBROUTINE broyden_solver3


	! =========================================================
	SUBROUTINE jacobi_partial(fun,x,n,m,method,ntypfun,nfun,col1,colend,rel_delta,f0,J)
	! Numerical jacobian matrix of funtion 'fun' in point 'x', from columm 'col1' to 'colend'
	! Usage:
	! m: Length of x
	! n: Length of F(x)
	! rel_delta: relative step size for computing numerical derivatives. If
	!  negative (flag) then an absolute (constant) value, is used.
	! method: finite differences method: 1: central, 2: backward, 3:forward
	! f0: funtion value in point x, f(x). Used to avoid recompute this value
	!       in the case of backward and forward schemes.
	IMPLICIT NONE
	EXTERNAL fun
	INTEGER, INTENT(IN)  :: n,m,method,ntypfun,nfun(ntypfun),col1,colend
	REAL*8,  INTENT(IN)  :: x(m), rel_delta(ntypfun),f0(n)
	REAL*8,  INTENT(OUT) :: J(n,m)
	REAL*8   delta(m), dx(m),col(n),normx,f_dx_back(n),f_dx_forw(n)
	REAL*8 T(n,m),b(n) !F(x) = 0, T(x)*x - b = 0
	REAL*8 normxvect(ntypfun)
	INTEGER  i, max_prev, mmax,idx_aux(ntypfun,2),eye(n)
	LOGICAL errorcheck
	eye=1.d0
	mmax = MAXVAL(nfun)

	errorcheck=(col1.LT.1.OR.col1.GT.m.OR.colend.LT.1.OR.col1.GT.m.OR.col1.GT.colend)
	IF (errorcheck) THEN
		WRITE(*,*) 'Error in jacobi_partial: check col1 and colend'
		READ(*,*)
	ENDIF

	IF (mmax.EQ.m) THEN
		normx = DSQRT(SUM(x*x))
		IF (normx.NE.0.D0.AND.rel_delta(1).GT.0.d0) THEN
			delta = rel_delta(1)*normx
		ELSE
			IF (normx.EQ.0.D0) THEN
				delta = DABS(rel_delta(1))
			ELSE
				delta = DABS(rel_delta(1))!*x
			ENDIF
		ENDIF
	ELSE
		! system contains more than one physical quatity
		idx_aux(1,:) = (/1,nfun(1)/)
		DO i=2,ntypfun
			max_prev = (idx_aux(i-1,2))
			idx_aux(i,1) = max_prev+1
			idx_aux(i,2) = max_prev+nfun(i)
		ENDDO

		DO i=1,ntypfun
			normxvect(i) = DSQRT(SUM(x(idx_aux(i,1):idx_aux(i,2))*x(idx_aux(i,1):idx_aux(i,2))))
			IF (normxvect(i).NE.0.D0.AND.rel_delta(i).GT.0.d0) THEN
				delta(idx_aux(i,1):idx_aux(i,2)) = rel_delta(i)*normxvect(i)
			ELSE
				delta(idx_aux(i,1):idx_aux(i,2)) = DABS(rel_delta(i))
			ENDIF
		ENDDO
	ENDIF


	IF (method.NE.1.AND.SUM(DABS(f0)).EQ.0.d0) THEN
		!f0 is not yet computed
		CALL fun(x,f0,T,b)
	ENDIF

	J=0.d0
	DO i=col1,colend
		dx = 0.d0
		SELECT CASE (method)
		CASE(1)
			!central
			dx(i) = delta(i)/2.d0
			CALL fun(x-dx,f_dx_back,T,b)
			CALL fun(x+dx,f_dx_forw,T,b)
			col = (f_dx_forw-f_dx_back)/delta(i)
		CASE(2)
			!backward
			dx(i) = delta(i)
			CALL fun(x-dx,f_dx_back,T,b)
			col = (f0-f_dx_back)/delta(i)
			CASE DEFAULT
			!forward
			dx(i) = delta(i)
			CALL fun(x+dx,f_dx_forw,T,b)
			col = (f_dx_forw-f0)/delta(i)
		ENDSELECT
		J(:,i) = col
	ENDDO

	ENDSUBROUTINE jacobi_partial