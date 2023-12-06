    include 'mkl_rci.f90'
    SUBROUTINE mkl_nlsqp(fun,x0,n,m,eps0,maxiter,x,J,st_cr,niterout)
    ! MKL nonlinear least square problem without boundary constraints
    ! using Trust-Region (TR) algorithms
    ! USAGE:
    !   fun: objective function
    !   n: number of function variables
    !   m: dimension of function value
    
	 ! Couln't use this ATM... what is rs and jac_eps? absolute val. or relative val.?
	 !  documentation is lacking... 
	 
    USE MKL_RCI
    USE MKL_RCI_TYPE
    IMPLICIT NONE
    EXTERNAL mkl_free_buffers !routine to track/free memory
    EXTERNAL fun !objective function
    INTEGER, INTENT(IN)  :: n,m,maxiter
    REAL*8,  INTENT(IN)  :: x0(n),eps0(2)
    REAL*8,  INTENT(OUT) :: x(n), J(m,n)
    INTEGER, INTENT(OUT) :: niterout
    INTEGER st_cr
    REAL*8 fvec(m), fjac(m,n)
    REAL*8 eps(6), rs
    !REAL*8 LW(n),UP(n)
    INTEGER iter1, iter2, successful, RCI_request
    INTEGER iter
    REAL*8 r1, r2
    type(handle_tr) :: handle ! TR solver handle
    INTEGER res,info(6)
    REAL*8 tolx,tolfun,jac_eps
    REAL*8 Ix(n), IJ(m,n), nan, epsval
    
    !INTEGER jj
    
    Ix = 1.d0
    IJ = 1.d0
    fvec=0.d0
    fjac=0.d0
    epsval = 1.0D-16
	 
    !set max number of iterations
    iter1 = ABS(maxiter)
    iter2 = CEILING(DFLOAT(iter1)/2.d0)
    
    ! set precisions for stop-criteria 
    ! (eps0 are absolute values!)
    tolx = eps0(1)
    tolfun = eps0(2)
    eps = TINY(tolx)
    !eps(1) = tolx ! minimum trust-region area
    eps(2) = tolfun !norm of the functional, ||F(x)||2
    eps(3) = tolfun !Singularity of the jacobian matrix, ||J(x)||2 
    !eps(4) = tolx !trial step norm, ||s||2
    !eps(5) = tolfun !delta (||F(x)||2 -||F(x)-J(x)s||2)
    jac_eps = 1.D-3
    
    ! set initial values
    x = x0 !initial guess
    !IF (SUM(x).EQ.0.d0) x(1) = x(1)+epsval
    rs = 1.0d1 !initial size of the trust region (boundary of the trial step)
    
    ! Initialize the solver (allocate memory, set initial values)
    res = dtrnlsp_init(handle,n,m,x,eps,iter1,iter2,rs)
    IF (res.ne.TR_SUCCESS) THEN 
        WRITE(*,*) 'Error in dtrnlsp_init MKL solver routine'
        CALL mkl_free_buffers
        READ(*,*)
    ENDIF
    
    !Check the correctness of handle and arrays containing jacobian matrix, 
    !   objective function, lower and upper bounds, and stopping criteria
    res = dtrnlsp_check(handle,n,m,fjac,fvec,eps,info)
    IF (res.ne.TR_SUCCESS) THEN
        WRITE(*,*) 'Error in dtrnlsp_check MKL solver routine'
        CALL mkl_free_buffers
        READ(*,*)
    ELSE
        IF (SUM(info(1:4)).NE.0) THEN
            WRITE(*,*) 'Input parameters for the MKL TR solver are not valid.'
            WRITE(*,'(A,4I2)') '   dtrnlsp_check: info = ', info(1:4)
            CALL mkl_free_buffers
            READ(*,*)
        ENDIF
    ENDIF
    
    ! Initial RCI (reverse communication interface) do-loop variables
    RCI_request = 0
    successful = 0
    st_cr = 1
    DO WHILE (successful.EQ.0)
        ! Call TR solver
        res = dtrnlsp_solve(handle,fvec,fjac,RCI_request)
        
        IF (res.ne.TR_SUCCESS) THEN
            WRITE(*,*) 'Error in dtrnlsp_solve MKL solver routine'
            CALL mkl_free_buffers
            READ(*,*)
        ENDIF
        !RCI_request (in/out): returns a flag converning the next task stage
        SELECT CASE (RCI_request)
            
        CASE (-1)
            ! The algorithm has exceeded the maximum number of iterations
            niterout = iter1
            nan = DSQRT(-1.d0)
            x = nan*Ix
            J = nan*IJ
            RETURN
        CASE (-2,-3,-4,-5,-6)
            ! Exit do-loop
            successful = 1
        CASE (1)
           ! recalculate function value
            CALL fun(x,fvec)
            WRITE(*,*) fvec
        CASE (2)
           ! compute jacobi matrix
            CALL mkl_djacobi(n,m,fun,x,fjac,jac_eps)
!              WRITE(*,*) 'Costum jacobi'
!              CALL jacobi(fun,x,n,m,1,1,(/n/),eps0(1),fvec,fjac)
!
!              IF (SIZE(J,1).EQ.6.OR.SIZE(J,1).EQ.7) THEN
!                  DO jj=1,6
!                      WRITE(*,37) fjac(jj,1),fjac(jj,2),fjac(jj,3), &
!                          fjac(jj,4),fjac(jj,5),fjac(jj,6)
!                  ENDDO
!
!              ENDIF
!              WRITE(*,*) ' '
!
!
!              WRITE(*,*) 'MKL jacobi'
!              CALL mkl_djacobi(n,m,fun,x,fjac,jac_eps)
!              WRITE(*,*) ' '
!              IF (SIZE(J,1).EQ.6.OR.SIZE(J,1).EQ.7) THEN
!                  DO jj=1,6
!                      WRITE(*,37) fjac(jj,1),fjac(jj,2),fjac(jj,3), &
!                          fjac(jj,4),fjac(jj,5),fjac(jj,6)
!                  ENDDO
!37                FORMAT(5X,6F15.4)
!              ENDIF
!              WRITE(*,*) ' '
!
!              !CALL jacobi(fun,x,n,m,1,1,(/n/),eps0(1),fvec,fjac)
!
!              WRITE(*,*) ' '
!              !res = djacobi(fun,n,m,fjac,x,1.d-6)
!              !IF (res.ne.TR_SUCCESS) THEN
!              !    WRITE(*,*) 'Error in djacobi MKL solver routine'
!              !    CALL mkl_free_buffers
!              !    READ(*,*)
!              !ENDIF

        ENDSELECT
    ENDDO
    
    ! Get solution status
    res = dtrnlsp_get(handle,iter,st_cr,r1,r2)
    !r1, r2 :initial and final residuals, resp.
    IF (res.ne.TR_SUCCESS) THEN
        WRITE(*,*) 'Error in dtrnlsp_get MKL solver routine'
        CALL mkl_free_buffers
        READ(*,*)
    ENDIF
    niterout = iter

    
    !Free handle memory
    CALL mkl_free_buffers
    res = dtrnlsp_delete(handle)
    
    ! Test residual of the solution
    IF (r2<tolfun) THEN
        ! success
        J = RESHAPE(fjac,(/m,n/)) ! check for (n,m) or (m,n)
    ELSE
        ! solution fails tolerance requirements
        WRITE(*,*) 'Solution fails tolerance'
        !READ(*,*)
    ENDIF

	ENDSUBROUTINE
	
	! -----------------------------
	
	SUBROUTINE mkl_djacobi(n,m,fun,x,fjac,jac_eps)
	! RCI-based MKL jacobi computation
	USE MKL_RCI
	USE MKL_RCI_TYPE
	IMPLICIT NONE
	EXTERNAL mkl_free_buffers !routine to track/free memory
	EXTERNAL fun !objective function
	INTEGER, INTENT(IN)  :: n,m
	REAL*8,  INTENT(IN)  :: x(n),jac_eps
	REAL*8,  INTENT(OUT) :: fjac(m,n)
	INTEGER res
	INTEGER*8 handle_jac !Jacobi-matrix solver handle
	REAL*8 f1(m), f2(m)
	INTEGER successful_jac, RCI_request_jac
    REAL*8 auxf1(1), auxf2(6)
    
    fjac=0.d0
    
	res=djacobi_init(handle_jac,n,m,x,fjac,jac_eps)
	IF (res.ne.TR_SUCCESS) THEN
		WRITE(*,*) 'Error in djacobi_init MKL solver routine'
		CALL mkl_free_buffers
		READ(*,*)
	ENDIF

	!Set initial RCI_jac
	RCI_request_jac = 0
	successful_jac  = 0
	DO WHILE (successful_jac.EQ.0)
		!Call solver
		res = djacobi_solve(handle_jac,f1,f2,RCI_request_jac)
		IF (res.NE.TR_SUCCESS) THEN
			WRITE(*,*) 'Error in djacobi_solve MKL solver routine'
			CALL mkl_free_buffers
			READ(*,*)
		END IF
		IF (RCI_request_jac.EQ.1) THEN
			!Calculate the function value f1 = f(x+eps)
			CALL fun(x,f1,auxf1,auxf2)
		ELSEIF (RCI_request_jac.EQ.2) THEN
			!Calculate the function value f2 = f(x-eps)
			CALL fun(x,f2,auxf1,auxf2)
		ELSEIF (RCI_request_jac.EQ.0) THEN
			!Exit rci cycle
			successful_jac = 1
		ENDIF
    END DO
    
	! Free handle memory
	res = djacobi_delete(handle_jac)
	IF (res.ne.TR_SUCCESS) THEN
		WRITE(*,*) 'Error in djacobi_delete MKL solver routine'
		CALL mkl_free_buffers
		READ(*,*)
	ENDIF

	ENDSUBROUTINE mkl_djacobi
    