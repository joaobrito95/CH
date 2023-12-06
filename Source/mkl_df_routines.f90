	INCLUDE 'mkl_df.f90'
	SUBROUTINE pp_interpolant(nx,x,ny,y,flag_ppnat,nxq,xq,vq)
	USE MKL_DF_TYPE
	USE MKL_DF

	IMPLICIT NONE
	INTEGER, INTENT(IN) :: nx, ny, nxq
	REAL*8, INTENT(IN) :: x(nx), y(nx*ny)
	REAL*8, INTENT(IN) :: xq(nxq)
	REAL*8, INTENT(OUT) :: vq(nxq,ny)
	INTEGER, INTENT(IN) :: flag_ppnat

	TYPE(DF_TASK) :: task
	INTEGER :: xhint,yhint,status
	INTEGER :: sorder,stype,bctype,ictype,scoeffhint
	INTEGER :: NNSCOEFF
	REAL*8, ALLOCATABLE :: scoeff(:)

	INTEGER :: typee,method,sitehint,ndorder,rhint
	INTEGER :: dorder(1)
	INTEGER :: nsite
	REAL*8 :: site(nxq)
	REAL*8, ALLOCATABLE :: r(:)
	!nx: number of breakpoints in the interval [a,b]
	!ny: umber of functions (dim of vector-valued function y)
	!x: strictly sorted breakpoints from interpolation interval [a,b]
	!y: Vector-valued function y

	xhint = DF_NON_UNIFORM_PARTITION
	yhint = DF_MATRIX_STORAGE_ROWS
	!Create Data Fitting task
	status = dfdnewtask1d(task,nx,x,xhint,ny,y,yhint)
	IF (status.NE.0) CALL CheckDfError(status)

	! Set spline parameters in the Data Fitting task
	IF (flag_ppnat.NE.1) THEN
		IF (nx.LT.5) THEN
			sorder     = DF_PP_LINEAR ! spline order
			stype      = DF_PP_DEFAULT ! spline type
			bctype     = DF_NO_BC  ! boundary conditions type
		ELSE
			sorder     = DF_PP_CUBIC ! spline order
			stype      = DF_PP_AKIMA! !or DF_PP_BESSEL ! spline type
			bctype     = DF_BC_NOT_A_KNOT  ! boundary conditions type
		ENDIF
	ELSE
		sorder     = DF_PP_CUBIC ! spline order
		stype      = DF_PP_NATURAL ! spline type
		bctype     = DF_BC_NOT_A_KNOT  ! boundary conditions type
	ENDIF
	ictype     = DF_NO_IC ! parameters describing internal conditions type
	scoeffhint = DF_NO_HINT! spline coefficients storage format
	NNSCOEFF   = ny*sorder*(nx-1)
	ALLOCATE(scoeff(NNSCOEFF))
	status = dfdeditppspline1d(task,sorder,stype,bctype, &
		ic_type=ictype,scoeff=scoeff,scoeffhint=scoeffhint)
	IF (status.NE.0) CALL CheckDfError(status)

	!Construct spline
	typee = DF_PP_SPLINE
	method = DF_METHOD_STD
	status = dfdconstruct1d(task,typee,method)
	IF (status.NE.0) CALL CheckDfError(status)

	! Initialize interpolation parameters and set site values
	typee = DF_INTERP
	method = DF_METHOD_PP
	sitehint = DF_SORTED_DATA
	ndorder = 1
	dorder = (/0/)
	rhint =	   DF_MATRIX_STORAGE_FUNCS_SITES_DERS
	nsite = nxq
	site = xq
	ALLOCATE (r(ny*nsite))

	!Compute the spline values at the points site
	status = dfdInterpolate1D (task,typee,method,nsite,site,&
		sitehint,ndorder,dorder,r=r,rhint=rhint)
	IF (status.NE.0) CALL CheckDfError(status)

	!De-allocate Data Fitting task resources
	status = dfdeletetask(task);
	IF (status.NE.0) CALL CheckDfError(status)

	vq = RESHAPE(r,(/nxq,ny/))

	ENDSUBROUTINE pp_interpolant

	! ================================================================
	SUBROUTINE pp_integrate1D(nx,x,ny,y,a,b,I)
	USE MKL_DF_TYPE
	USE MKL_DF

	IMPLICIT NONE
	INTEGER, INTENT(IN) :: nx, ny
	REAL*8, INTENT(IN)  :: x(nx), y(nx*ny),a,b
	REAL*8, INTENT(OUT) :: I(1,ny)

	TYPE(DF_TASK) :: task
	INTEGER :: xhint,yhint,status
	INTEGER :: sorder,stype,bctype,ictype,scoeffhint
	INTEGER :: NNSCOEFF
	REAL*8, ALLOCATABLE :: scoeff(:)

	INTEGER :: method,typee,rhint
	!INTEGER :: nsite

	INTEGER :: nlim, llimhint,rlimhint
	REAL*8  :: llim(1), rlim(1)
	REAL*8, ALLOCATABLE :: r(:)

	!IF ((MAXVAL(x).LT.b).OR.(MINVAL(x).GT.a))	THEN
	!	  WRITE(*,*) 'Error in pp_integrate1D, check limits'
	!	  READ(*,*)
	!ENDIF
	
	! Usage:
	!nx: number of breakpoints in the interval [a,b]
	!ny: umber of functions (dim of vector-valued function y)
	!x: strictly sorted breakpoints from interpolation interval [a,b]
	!y: Vector-valued function y

	xhint = DF_NON_UNIFORM_PARTITION
	yhint = DF_MATRIX_STORAGE_ROWS
	!Create Data Fitting task
	status = dfdnewtask1d(task,nx,x,xhint,ny,y,yhint)
	!CALL CheckDfError(status)

	! Set spline parameters in the Data Fitting task
	sorder     = DF_PP_CUBIC ! spline order
	stype      = DF_PP_NATURAL ! spline type, DF_PP_AKIMA, DF_PP_NATURAL
	bctype     = DF_BC_NOT_A_KNOT  ! boundary conditions type
	ictype     = DF_NO_IC ! parameters describing internal conditions type
	scoeffhint = DF_NO_HINT! spline coefficients storage format
	NNSCOEFF   = ny*sorder*(nx-1)
	ALLOCATE(scoeff(NNSCOEFF))
	status = dfdeditppspline1d(task,sorder,stype,bctype, &
		ic_type=ictype,scoeff=scoeff,scoeffhint=scoeffhint)
	!CALL CheckDfError(status)

	!Construct spline
	typee = DF_PP_SPLINE
	method = DF_METHOD_STD
	status = dfdconstruct1d(task,typee,method)
	!CALL CheckDfError(status)

	! Initialize interpolation parameters and set site values
	method = DF_METHOD_PP
	rhint =	   DF_MATRIX_STORAGE_FUNCS_SITES_DERS
	nlim = 1
	llim(1) = a
	rlim(1) = b
	ALLOCATE (r(nlim*ny))

	! Compute integrals ***
	status = dfdintegrate1d(task,method,nlim,llim,          &
		&    llimhint,rlim,rlimhint,r=r,rhint=rhint)
	!CALL CheckDfError(status)

	!De-allocate Data Fitting task resources
	status = dfdeletetask(task);

	I=RESHAPE(r,(/nlim,ny/))

	ENDSUBROUTINE pp_integrate1D

	! ================================================================
	SUBROUTINE CheckDfError(num)
	USE MKL_DF_TYPE

	INTEGER(KIND=4) :: num
	IF ( num == DF_ERROR_NULL_TASK_DESCRIPTOR ) THEN
		PRINT 33, "Error: null task descriptor (code ", num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_MEM_FAILURE ) THEN
		PRINT 33, "Error: memory allocation failure in DF functionality  &
			&             (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_NX ) THEN
		PRINT 33, "Error: the number of breakpoints is invalid (code ",  &
			&             num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_X ) THEN
		PRINT 33, "Error: array of breakpoints is invalid (code ",       &
			&             num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_X_HINT ) THEN
		PRINT 33, "Error: invalid flag describing structure of partition &
			&             (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_NY ) THEN
		PRINT 33, "Error: invalid dimension of vector-valued function y  &
			&             (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_Y ) THEN
		PRINT 33, "Error: array of function values is invalid (code ",   &
			&             num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_Y_HINT ) THEN
		PRINT 33, "Error: invalid flag describing structure of function  &
			&             y (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_SPLINE_ORDER ) THEN
		PRINT 33, "Error: invalid spline order (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_SPLINE_TYPE ) THEN
		PRINT 33, "Error: invalid spline type (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_IC_TYPE ) THEN
		PRINT 33, "Error: invalid type of internal conditions used       &
			&             in the spline construction (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_IC ) THEN
		PRINT 33, "Error: array of internal conditions for spline        &
			&             construction is not defined (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_BC_TYPE ) THEN
		PRINT 33, "Error: invalid type of boundary conditions used       &
			&             in the spline construction (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_BC ) THEN
		PRINT 33, "Error: array of boundary conditions for spline        &
			&             construction is not defined (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_PP_COEFF ) THEN
		PRINT 33, "Error: array of piece-wise polynomial spline          &
			&             coefficients is not defined (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_PP_COEFF_HINT ) THEN
		PRINT 33, "Error: invalid flag describing structure of the       &
			&             piece-wise polynomial spline coefficients (code ",    &
			&             num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_PERIODIC_VAL ) THEN
		PRINT 33, "Error: function values at the end points of the       &
			&             interpolation interval are not equal as required      &
			&             in periodic boundary conditions (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_DATA_ATTR ) THEN
		PRINT 33, "Error: invalid attribute of the pointer to be set or  &
			&             modIFied in Data Fitting task descriptor with         &
			&             EditIdxPtr editor (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_DATA_IDX ) THEN
		PRINT 33, "Error: index of pointer to be set or modIFied in Data &
			&             Fitting task descriptor with EditIdxPtr editor is     &
			&             out of range (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_NSITE ) THEN
		PRINT 33, "Error: invalid number of interpolation sites (code ", &
			&             num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_SITE ) THEN
		PRINT 33, "Error: array of interpolation sites is not defined    &
			&             (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_SITE_HINT ) THEN
		PRINT 33, "Error: invalid flag describing structure of           &
			&             interpolation sites (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_NDORDER ) THEN
		PRINT 33, "Error: invalid size of array that defines order       &
			&             of the derivatives to be computed at the              &
			&             interpolation sites (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_DORDER ) THEN
		PRINT 33, "Error: array defining derivative orders to be         &
			&             computed at interpolation sites is not defined        &
			&             (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_DATA_HINT ) THEN
		PRINT 33, "Error: invalid flag providing a-priori information    &
			&             about partition and/or interpolation sites (code ",   &
			&             num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_INTERP ) THEN
		PRINT 33, "Error: array of spline based interpolation results    &
			&             is not defined (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_INTERP_HINT ) THEN
		PRINT 33, "Error: invalid flag defining structure of spline      &
			&             based interpolation results (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_CELL_IDX ) THEN
		PRINT 33, "Error: array of indices of partition cells containing &
			&             interpolation sites is not defined (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_NLIM ) THEN
		PRINT 33, "Error: invalid size of arrays containing integration  &
			&             limits (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_LLIM ) THEN
		PRINT 33, "Error: array of left integration limits               &
			&             is not defined (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_RLIM ) THEN
		PRINT 33, "Error: array of right integration limits              &
			&             is not defined (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_INTEGR ) THEN
		PRINT 33, "Error: array of spline based integration results      &
			&             is not defined (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_INTEGR_HINT ) THEN
		PRINT 33, "Error: invalid flag defining structure of spline      &
			&             based inetgration results (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_BAD_LOOKUP_INTERP_SITE ) THEN
		PRINT 33, "Error: bad site provided for interpolation with       &
			&             look-up interpolator (code ",num,")."
		READ(*,*)
	end IF
	IF ( num == DF_ERROR_NULL_PTR ) THEN
		PRINT 33, "Error: bad pointer provided in DF function (code ",   &
			&             num,")."
	end IF

33	FORMAT(A,I5,A)
	ENDSUBROUTINE CheckDfError