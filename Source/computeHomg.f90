	SUBROUTINE compHomog_commun(D6,YesEshelby,d5_ast,f0,propMAT,sphtdesign,macroSRP,SIGMA_norm)

	USE M_COMP_HOMOG, ONLY : struct_mat_param,struct_sphtdesign,&
		b_radius, Y0_matrix, Vol_RVE, fmin_adm, fmax_adm
	USE M_CONSTANTS, ONLY : PI

	IMPLICIT NONE

	REAL*8, INTENT(IN)                   :: D6(6), f0, d5_ast(5)
	LOGICAL, INTENT(IN)                  :: YesEshelby
	TYPE (struct_mat_param),  INTENT(IN) :: propMAT
	TYPE (struct_sphtdesign), INTENT(IN) :: sphtdesign
	REAL*8, INTENT(OUT)                  :: SIGMA_norm(6), macroSRP

	INTEGER N_design, N_design_eff , nRad_bar
	INTEGER numel_r, num_r_ref, size_idx_layer(2), size_idx_layer_init(2)
	REAL*8, ALLOCATABLE :: r(:)!, theta(:), phi(:)
	INTEGER, ALLOCATABLE :: idx_layer(:,:), idx_layer_ref(:,:)
	REAL*8  a_radius
	REAL*8, ALLOCATABLE :: srp(:,:), sigma_dev(:,:,:), dir_sigma_dev(:,:,:)
	REAL*8, ALLOCATABLE :: srp_init(:,:), sigma_dev_init(:,:,:), dir_sigma_dev_init(:,:,:)
	LOGICAL, ALLOCATABLE :: mask(:,:)!, masknotNaN(:),mask_vect(:)
	REAL*8, ALLOCATABLE :: surf_aver_srp(:), integrand_srp(:)
	REAL*8, ALLOCATABLE :: surf_aver_srp_init(:), integrand_srp_init(:)
	REAL*8, ALLOCATABLE :: surf_aver_sdev(:,:), integrand_sdev(:,:)
	REAL*8  volint_SRP(1,1), volint_sdev(1,6), macroSDEV(6)
	REAL*8 pwork,dev_pwork,vol_pwork,SIGMA_m,trace_D,D6_dev(6),D6_dev_pwork(6)
	REAL*8 EPS
	REAL*8, ALLOCATABLE :: r_init(:)
	INTEGER numel_r_init, flag_radd, numel_ref
	INTEGER, ALLOCATABLE :: first_layer(:)
	INTEGER n_add_r, ngap_init,ngap_end, n_add_r_sig, i
	REAL*8, ALLOCATABLE :: add_r(:), add_r_sig(:)
	INTEGER, ALLOCATABLE :: idx_r_init(:)
	INTEGER niterall_init, niterall
	REAL*8 niter_avr_init, niter_avr
	LOGICAL test_aymp,test_aymp_sig
	REAL*8, ALLOCATABLE :: r_dummy(:), integrand_srp_dummy(:)
	REAL*8, ALLOCATABLE :: r_all_sig(:), surf_aver_sdev_all(:,:), integrand_sdev_all(:,:)
	REAL*8, ALLOCATABLE :: r_all(:), integrand_srp_all(:), r_all_presort(:)
	INTEGER ndummy, ndummy_sig, numel_r_all, numel_r_all_sig
	INTEGER idx_r_last, idx_last_nz
	INTEGER, ALLOCATABLE :: idx_sort(:)
	REAL*8 ub_factor

	macroSRP = 0.d0
	SIGMA_norm = 0.d0
	IF (SUM(DABS(D6)).EQ.0.d0) RETURN

	! Unpack spherical t-design
	N_design=    sphtdesign%N_dsgn
	N_design_eff=sphtdesign%N_dsgn_eff

	! Get radial discretization
	IF (f0.LE.0.d0.OR.f0.GE.1.d0) THEN
		WRITE(*,*) 'f value not admissible in CH'
		READ(*,*)
	ELSE
		!flag_radd=1
		flag_radd=0
		nRad_bar = sphtdesign%bar_n_rad
		CALL radial_disc_lin(f0,b_radius,a_radius,flag_radd,nRad_bar,numel_r_init,r_init)
		IF (flag_radd.EQ.1) THEN
			IF (FLOOR(dlog10(f0))<=-5) THEN
				numel_ref = (numel_r_init-2)
			ELSE
				numel_ref = (numel_r_init-1)
			ENDIF
		ELSE
			numel_ref = numel_r_init
		ENDIF
	ENDIF

	! Get sets of radial coord (relative ordering of interpolation)
	IF (numel_ref.GT.5) THEN
		num_r_ref = 5
	ELSE
		num_r_ref = 3
	ENDIF
	first_layer=(/0/)
	CALL get_layers(numel_r_init,r_init,num_r_ref,first_layer,idx_layer_ref)
	size_idx_layer_init = SHAPE(idx_layer_ref)

	ALLOCATE(srp_init(numel_r_init,N_design)); srp_init=0.d0
	ALLOCATE(sigma_dev_init(numel_r_init,N_design,6)); sigma_dev_init=0.d0
	ALLOCATE(dir_sigma_dev_init(numel_r_init,N_design,6)); dir_sigma_dev_init=0.d0

	! Compute local fields in the RVE (reference disc)
	CALL local_fields(D6,d5_ast,f0,propMAT,N_design,sphtdesign,numel_r_init,r_init,YesEshelby,&
		size_idx_layer_init,idx_layer_ref,1,size_idx_layer_init(1),&
		srp_init,sigma_dev_init,dir_sigma_dev_init,niterall_init)
	niter_avr_init=DFLOAT(niterall_init)/DFLOAT(N_design_eff*numel_r_init)
	!WRITE(*,*) niter_avr_init

	! Uniform arc-length based reparametrization
	surf_aver_srp_init =(SUM(srp_init,DIM=2))/N_design_eff
	integrand_srp_init = surf_aver_srp_init*(4.d0*PI*r_init*r_init)
	IF (flag_radd.EQ.1) THEN
		IF (FLOOR(dlog10(f0))<=-5) THEN
			ngap_init =3
		ELSE
			ngap_init = 2
		ENDIF
	ELSE
		ngap_init = 1
	ENDIF
	ngap_end=(numel_r_init)

	! Detect jumps of the srp integrand (and/or srp)
	CALL arc_len_resamp(numel_r_init,ngap_init,ngap_end,r_init,surf_aver_srp_init,integrand_srp_init,numel_ref,&
		n_add_r,add_r,test_aymp)

	numel_r=numel_r_init+n_add_r
	ALLOCATE(r(numel_r)); r=0.d0
	ALLOCATE(srp(numel_r,N_design)); srp=0.d0
	ALLOCATE(sigma_dev(numel_r,N_design,6)); sigma_dev=0.d0
	dir_sigma_dev = sigma_dev
	ALLOCATE(idx_r_init(numel_r_init))

	IF (n_add_r.GT.0) THEN
		! Compute local fields in the new radial surfaces
		r(1:numel_r_init) = r_init
		r(numel_r_init+1:numel_r)=add_r
		CALL DSORT(r,1,numel_r)
		DO i=1,numel_r_init
			idx_r_init(i)=MINLOC(DABS(r-r_init(i)),1)
		ENDDO
		srp(idx_r_init,:)=srp_init
		sigma_dev(idx_r_init,:,:) = sigma_dev_init
		dir_sigma_dev(idx_r_init,:,:) = dir_sigma_dev_init

		IF (ALLOCATED(first_layer)) DEALLOCATE(first_layer)
		first_layer=idx_r_init
		CALL get_layers(numel_r,r,numel_r_init,first_layer,idx_layer)
		size_idx_layer = SHAPE(idx_layer)

		CALL local_fields(D6,d5_ast,f0,propMAT,N_design,sphtdesign,numel_r,r,YesEshelby,&
			size_idx_layer,idx_layer,2,size_idx_layer(1),&
			srp,sigma_dev,dir_sigma_dev,niterall)
		niter_avr=DFLOAT(niterall)/DFLOAT(N_design_eff*n_add_r)
		!WRITE(*,*) niter_avr
	ELSE
		r=r_init
		srp=srp_init
		sigma_dev=sigma_dev_init
		dir_sigma_dev=dir_sigma_dev_init
	ENDIF


	numel_r=numel_r_init+n_add_r
	ngap_end=numel_r

	! Detect jumps of sigma_dev
	ALLOCATE(surf_aver_sdev(numel_r,6))
	DO i=1,6
		surf_aver_sdev(:,i) = (SUM(sigma_dev(:,:,i),DIM=2))/N_design_eff
	ENDDO
	CALL arc_len_resamp_sig(numel_r,ngap_init,ngap_end,r,surf_aver_sdev,numel_ref,&
		n_add_r_sig,add_r_sig,test_aymp_sig,idx_last_nz)

	IF (test_aymp_sig.AND.n_add_r_sig.GT.0) THEN

		DEALLOCATE(r_init,srp_init,sigma_dev_init,dir_sigma_dev_init)
		numel_r_init=numel_r
		r_init=r
		srp_init=srp
		sigma_dev_init=sigma_dev
		dir_sigma_dev_init=dir_sigma_dev
		DEALLOCATE(r,srp,sigma_dev,dir_sigma_dev,idx_r_init,surf_aver_sdev)

		! Compute local fields in the new radial surfaces
		numel_r=numel_r+n_add_r_sig
		ALLOCATE(r(numel_r)); r=0.d0
		ALLOCATE(srp(numel_r,N_design)); srp=0.d0
		ALLOCATE(sigma_dev(numel_r,N_design,6)); sigma_dev=0.d0
		dir_sigma_dev = sigma_dev
		ALLOCATE(idx_r_init(numel_r_init)); idx_r_init=0

		r(1:numel_r_init) = r_init
		r(numel_r_init+1:numel_r)=add_r_sig
		CALL DSORT(r,1,numel_r)
		DO i=1,numel_r_init
			idx_r_init(i)=MINLOC(DABS(r-r_init(i)),1)
		ENDDO
		srp(idx_r_init,:)=srp_init
		sigma_dev(idx_r_init,:,:) = sigma_dev_init
		dir_sigma_dev(idx_r_init,:,:) = dir_sigma_dev_init

		IF (ALLOCATED(first_layer)) DEALLOCATE(first_layer)
		first_layer=idx_r_init
		CALL get_layers(numel_r,r,numel_r_init,first_layer,idx_layer)
		size_idx_layer = SHAPE(idx_layer)

		CALL local_fields(D6,d5_ast,f0,propMAT,N_design,sphtdesign,numel_r,r,YesEshelby,&
			size_idx_layer,idx_layer,2,size_idx_layer(1),&
			srp,sigma_dev,dir_sigma_dev,niterall)
		niter_avr=DFLOAT(niterall)/DFLOAT(N_design_eff*n_add_r_sig)
		!WRITE(*,*) niter_avr
	ELSE


	ENDIF

	ngap_end=numel_r

	! Get surface averages and integrand of srp and stress components
	mask = (srp.EQ.srp) !non-NaN points (failed computation) in the sphere
	IF (count(mask).ne.SIZE(mask)) WRITE(*,*) 'NaN values in the local/sphere SRP computation'
	surf_aver_srp =	(SUM(srp,DIM=2,MASK=mask))/N_design_eff
	integrand_srp = surf_aver_srp*(4.d0*PI*r*r)
	IF (.not.ALLOCATED(surf_aver_sdev))	THEN
		ALLOCATE(surf_aver_sdev(numel_r,6))
		DO i=1,6
			surf_aver_sdev(:,i) = (SUM(sigma_dev(:,:,i),DIM=2))/N_design_eff
		ENDDO
	ENDIF
	!surf_aver_sdev = (SUM(sigma_dev,DIM=2,MASK=SPREAD(mask,3,6)))/N_design_eff
	integrand_sdev = surf_aver_sdev*SPREAD((4.d0*PI*r*r),2,6)

	! Add dummy points to stabilize the jumps of dev_sig near inner surf
	! (always done, regardeless of (test_aymp_sig))
	!IF (test_aymp_sig) THEN
	!	idx_r_last = MAX(idx_last_nz+n_add_r_sig,4)
	!ELSE
	!	ub_factor=(a_radius+(b_radius-a_radius)*0.25d0)
	!	idx_r_last=FINDLOC(r<=ub_factor,.true.,DIM=1,BACK=.true.)
	!ENDIF
	!IF (test_aymp_sig.AND.n_add_r_sig.GT.0) THEN
	!	CALL add_asymptote_pts_sig(numel_r,idx_r_last,r,surf_aver_sdev,flag_radd,test_aymp_sig,&
	!		numel_r_all_sig,r_all_sig,surf_aver_sdev_all)
	!ELSE
	!	numel_r_all_sig=numel_r
	!	r_all_sig=r
	!	surf_aver_sdev_all=surf_aver_sdev
	!ENDIF
	numel_r_all_sig=numel_r
	r_all_sig=r
	surf_aver_sdev_all=surf_aver_sdev
	integrand_sdev_all = surf_aver_sdev_all*(4.d0*PI)*SPREAD(r_all_sig*r_all_sig,2,6)


	! Add dummy points to stabilize the spline near inner surf y~a(1/r)
	! (avoid runge phenomenon/overshoot)
	IF (test_aymp.AND.n_add_r.GT.1) THEN
		CALL add_asymptote_pts(numel_r,r,integrand_srp,flag_radd,ndummy,r_dummy,integrand_srp_dummy)
		numel_r_all=numel_r+ndummy
		ALLOCATE(r_all(numel_r_all))
		r_all(1:numel_r)=r
		ALLOCATE(integrand_srp_all(numel_r_all))
		integrand_srp_all(1:numel_r)=integrand_srp
		IF (ndummy.GT.0) THEN
			r_all(numel_r+1:numel_r_all)=r_dummy
			r_all_presort=r_all
			CALL DSORT(r_all,1,numel_r_all)
			ALLOCATE(idx_sort(numel_r_all))
			DO i=1,numel_r_all
				idx_sort(i)=MINLOC(DABS(r_all(i)-r_all_presort),1)
			ENDDO
			integrand_srp_all(numel_r+1:numel_r_all)=integrand_srp_dummy
			integrand_srp_all=integrand_srp_all(idx_sort)
		ENDIF
	ELSE
		numel_r_all=numel_r
		r_all=r
		integrand_srp_all=integrand_srp
	ENDIF


	! Integrate (1D) surf average quantities: limits [a,b], yet using spline with added points
	CALL pp_integrate1D(numel_r_all,r_all,1,integrand_srp_all,a_radius,b_radius,volint_SRP)
	macroSRP=volint_SRP(1,1)/Vol_RVE !eqaul to 3.d0
	!CALL pp_integrate1D(numel_r,r,6,integrand_sdev,a_radius,b_radius,volint_sdev)
	CALL pp_integrate1D(numel_r_all_sig,r_all_sig,6,integrand_sdev_all,a_radius,b_radius,volint_sdev)

	macroSDEV = RESHAPE(volint_sdev,(/6/))/Vol_RVE
	macroSDEV(1:3)=macroSDEV(1:3)-SUM(macroSDEV(1:3))/3.d0

	!Compute hydrostatic macroSIGMA
	trace_D = SUM(D6(1:3))
	D6_dev(1:3) = D6(1:3)-(trace_D/3.d0)
	D6_dev(4:6) = D6(4:6)
	D6_dev_pwork = D6_dev
	D6_dev_pwork(4:6) = 2.d0*D6_dev(4:6) !Voigt correction for dot product

	pwork = Y0_matrix*macroSRP                      ! total pwork
	dev_pwork = DOT_PRODUCT(macroSDEV,D6_dev_pwork) ! 'deviatoric' dissipation
	vol_pwork = (pwork-dev_pwork)                   ! 'volumetric' dissipation
	EPS = 2.2204460492503131D-16
	IF (trace_D.NE.0.d0.AND.DABS(trace_D/3.d0).GT.EPS) THEN
		SIGMA_m = vol_pwork/trace_D
	ELSE
		SIGMA_m = 0.d0
	ENDIF
	SIGMA_norm(1:3) = macroSDEV(1:3)+SIGMA_m
	SIGMA_norm(4:6) = macroSDEV(4:6)

	RETURN

	CONTAINS

	! =============================================================
	SUBROUTINE add_asymptote_pts_sig(n,nend,x,y,flag_radd,test_aymp,nnew,xnew,ynew)
	USE M_CONSTANTS, ONLY : PI
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, flag_radd
	INTEGER nend
	LOGICAL, INTENT(IN) :: test_aymp
	REAL*8, INTENT(IN) :: x(n),y(n,6)
	INTEGER, INTENT(OUT) :: nnew
	REAL*8, ALLOCATABLE, INTENT(OUT) :: xnew(:),ynew(:,:)

	INTEGER, PARAMETER :: ndisc=100

	REAL*8 x_red(nend),y_red(nend,6)
	REAL*8 xq(ndisc), x1, x2, x_all(ndisc+nend)
	REAL*8, ALLOCATABLE :: x_all_red(:)
	INTEGER r_ref_idx(ndisc)
	INTEGER n_x_all_red
	LOGICAL idx_del(ndisc+nend)
	REAL*8, ALLOCATABLE :: xadd_vect(:)
	REAL*8, ALLOCATABLE :: xnew_aux(:), y_prod(:), vq_aux(:,:),y_aux(:,:)
	INTEGER flag_ppnat, last_idx

	x_red=x(1:nend)
	y_red=y(1:nend,:)

	!ALLOCATE (r_xq(ndisc))
	r_ref_idx = (/(i,i=0,ndisc-1)/)
	x1	= x(1)
	x2 = x(nend)
	xq = x1+((x2-x1)/DFLOAT(ndisc-1))*DFLOAT(r_ref_idx(:))

	x_all=0.d0
	x_all(1:nend)=x_red
	x_all(nend+1:ndisc+nend)=xq
	CALL DSORT(x_all,1,ndisc+nend)
	idx_del=.false.
	idx_del(1:SIZE(x_all)-1)= x_all(2:ndisc+nend)-x_all(1:ndisc+nend-1).eq.0.d0
	x_all_red = PACK(x_all,.not.idx_del)

	n_x_all_red = SIZE(x_all_red)
	ALLOCATE(xnew(n-(nend)+n_x_all_red))
	xnew(1:n_x_all_red)=x_all_red
	xnew(n_x_all_red+1:n_x_all_red+(n-nend))=x((nend+1):n)

	!xnew=x_all_red
	nnew=SIZE(xnew)

	ALLOCATE (vq_aux(n_x_all_red,6))
	IF (test_aymp) THEN
		flag_ppnat=0 !Makima spline (C1), not C2
	ELSE
		flag_ppnat=1 !Natural spline (C2)
	ENDIF

	CALL pp_interpolant(n,x,6,y,flag_ppnat,n_x_all_red,x_all_red,vq_aux)
	y_aux=vq_aux !*(4.d0*PI)*SPREAD(xnew*xnew,2,6)
	ALLOCATE(ynew(nnew,6))
	ynew(1:n_x_all_red,:)=y_aux
	ynew((n_x_all_red+1):nnew,:)=y((nend+1):n,:) !*(4.d0*PI)*SPREAD(x(nend+1:n)*x(nend+1:n),2,6)

	! Verify uniqueness and sort ascending: implicit in code
	! CALL DSORT(xnew,1,cnt)

	ENDSUBROUTINE add_asymptote_pts_sig

	! =============================================================
	SUBROUTINE add_asymptote_pts(n,x,y,flag_radd,nnew,xnew,ynew)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, flag_radd
	REAL*8, INTENT(IN) :: x(n),y(n)
	INTEGER, INTENT(OUT) :: nnew
	REAL*8, ALLOCATABLE, INTENT(OUT) :: xnew(:),ynew(:)
	REAL*8, ALLOCATABLE :: x_red(:),y_red(:), error(:)
	LOGICAL idxtest(n)
	REAL*8 xmin,xmax,x0,x1,y0,y1,apower,error_max
	INTEGER idx_last
	REAL*8  x_normlizer, y_normlizer, L_aprox, maxLi
	REAL*8, ALLOCATABLE :: x_norm(:), y_norm(:), diffvec(:,:), arc_leng(:)
	INTEGER, ALLOCATABLE :: testmaxLi(:),numptsxq(:)
	INTEGER n_discf, cnt, igap, npoints_add
	REAL*8 f0
	INTEGER, PARAMETER :: ndisc=50
	REAL*8, ALLOCATABLE :: xadd_vect(:)
	REAL*8, ALLOCATABLE :: xnew_aux(:), y_prod(:), vq_aux(:)
	INTEGER flag_ppnat

	xmin=MINVAL(x)
	xmax=MAXVAL(x)

	idxtest= x < xmin+(xmax-xmin)/4.d0
	!idxtest= x < 0.25d0*xmax
	IF (COUNT(idxtest).LT.3) THEN
		!idxtest=.true.
		nnew=0
		RETURN
	ENDIF

	x_red=PACK(x,idxtest)
	y_red=PACK(y,idxtest)
	IF (flag_radd.EQ.1)THEN
		IF (FLOOR(dlog10((x(3)/x(n))**3))<=-5) THEN
			x0=x_red(3)
			x1=x_red(4)
			y0=y_red(3)
			y1=y_red(4)
			f0=(x(3)/x(n))**3
		ELSE
			x0=x_red(2)
			x1=x_red(3)
			y0=y_red(2)
			y1=y_red(3)
			f0=(x(2)/x(n))**3
		ENDIF

	ELSE
		x0=x_red(1)
		x1=x_red(2)
		y0=y_red(1)
		y1=y_red(2)
		f0=(x(1)/x(n))**3
	ENDIF
	apower=y0*(y1/y0)**(DLOG(x0)/(DLOG(x0/x1)))
	error=DABS(apower/x_red-y_red)/y_red
	error_max=5.0d-2   !error_max=5.d-2
	idx_last = FINDLOC(error.LT.error_max,.true.,DIM=1,BACK=.true.)

	!idx_last=0

	IF (idx_last.EQ.0) THEN
		nnew=0
		RETURN
	ENDIF

	IF (flag_radd.EQ.1) THEN
		idx_last=MAX(idx_last,5)
	ELSE
		idx_last=MAX(idx_last,4)
	ENDIF

	x_red=x_red(1:idx_last)
	y_red=y_red(1:idx_last)

	x_normlizer = MAXVAL(x_red)
	y_normlizer = MAXVAL(y_red)
	x_norm=x_red/x_normlizer
	y_norm=y_red/y_normlizer
	ALLOCATE (diffvec(SIZE(x_norm)-1,2))
	diffvec(:,1)=	x_norm(1:idx_last-1)-x_norm(2:idx_last)
	diffvec(:,2)=	y_norm(1:idx_last-1)-y_norm(2:idx_last)
	arc_leng=DSQRT(SUM(diffvec*diffvec,DIM=2))
	L_aprox = SUM(arc_leng)

	n_discf=(-NINT(DLOG10(f0))-2)*(2*idx_last) !*10 !less important at high prosity val
	n_discf=MAX(n_discf,10)
	n_discf=MIN(n_discf,2*SIZE(x))

	maxLi=(L_aprox/n_discf)
	testmaxLi = NINT(arc_leng/maxLi)
	numptsxq = (testmaxLi-1)
	cnt=0

	ALLOCATE(xnew_aux(2*n_discf))
	xnew_aux=0.d0
	DO igap=1,SIZE(numptsxq)
		npoints_add=numptsxq(igap)
		IF (npoints_add.GT.0) THEN

			IF (ALLOCATED(xadd_vect)) DEALLOCATE(xadd_vect)
			xadd_vect = x_red(igap)+((x_red(igap+1)-x_red(igap))/DFLOAT(npoints_add+1))*DFLOAT((/(i,i=0,npoints_add+1)/))
			xnew_aux(cnt+1:cnt+npoints_add)=xadd_vect(2:npoints_add+1)
			cnt=cnt+npoints_add
		ENDIF
	ENDDO

	xnew=xnew_aux(1:cnt)
	nnew=cnt

	IF (nnew.GT.0) THEN
		ALLOCATE (vq_aux(nnew))
		! Get interpolation
		y_prod = x_red*y_red
		!IF  (SIZE(x_red).GE.5) THEN
		!	flag_ppnat=0 !Makima spline
		!ELSE
		!	flag_ppnat=1 !Natural spline (C2)
		!ENDIF

		IF (SIZE(x_red).LE.3) THEN
			WRITE(*,*) 'check add_asymptote_pts'
			READ(*,*)
		ENDIF


		flag_ppnat=1 !Natural spline (C2)
		CALL pp_interpolant(SIZE(x_red),x_red,1,y_prod,flag_ppnat,&
			nnew,xnew,vq_aux)
		ynew=vq_aux/xnew
	ELSE
		RETURN
	ENDIF

	! Verify uniqueness and sort ascending: implicit in code
	! CALL DSORT(xnew,1,cnt)

	ENDSUBROUTINE add_asymptote_pts

	! =============================================================

	SUBROUTINE get_layers(num_points,r,num_r_ref,first_layer,idx_layer)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: num_points
	REAL*8,  INTENT(IN) :: r(num_points)
	INTEGER, ALLOCATABLE, INTENT(IN) :: first_layer(:)
	INTEGER             :: num_r_ref
	INTEGER, ALLOCATABLE, INTENT(OUT) :: idx_layer(:,:)
	INTEGER, ALLOCATABLE :: idx_layer_aux(:,:), vect_idx_layer(:), &
		row_idx_layer(:), nnz_idx_layer(:), aux_vect_idx_layer_aux(:)

	INTEGER n_layer, i, idx_i, target_vect(num_points), cnt, max_j
	LOGICAL, ALLOCATABLE :: mask(:,:)
	REAL*8  a_aux, b_aux
	REAL*8, ALLOCATABLE :: aux_sort(:)
	INTEGER auxnum

	IF (num_r_ref.GT.num_points) num_r_ref = num_points

	ALLOCATE (idx_layer_aux(num_points,num_points))
	ALLOCATE (aux_vect_idx_layer_aux(num_points*num_points))
	idx_layer_aux = 0

	n_layer = 0
	target_vect = (/(i,i=1,num_points)/)
	DO WHILE (.TRUE.)
		n_layer = n_layer+1
		IF (n_layer.EQ.1) THEN
			IF (ALL(first_layer.EQ.0)) THEN
				DO i=1,(num_r_ref-1)
					a_aux = r(1) ! r already sorted
					b_aux = r(num_points)
					idx_i = MINLOC(DABS(r-(a_aux + DFLOAT(i)*(b_aux-a_aux)/(DFLOAT(num_r_ref)-1.d0))),1)
					idx_layer_aux(1,i+1) = idx_i
				ENDDO
				idx_layer_aux(1,1) = 1
				vect_idx_layer = PACK(idx_layer_aux(1,:),idx_layer_aux(1,:).NE.0)
			ELSE
				vect_idx_layer = first_layer;
				idx_layer_aux(1,1:SIZE(first_layer)) = first_layer
			ENDIF
		ELSE
			DO i=1,(SIZE(vect_idx_layer)-1)
				IF (vect_idx_layer(i+1)-vect_idx_layer(i).GT.1) THEN
					a_aux = r(vect_idx_layer(i))
					b_aux = r(vect_idx_layer(i+1))
					idx_i = MINLOC(DABS(r-(a_aux + (b_aux-a_aux)/2.d0)),1)
					IF (ALL(vect_idx_layer.NE.idx_i)) idx_layer_aux(n_layer,i) = idx_i
				ENDIF
			ENDDO
			aux_vect_idx_layer_aux=RESHAPE(idx_layer_aux,(/num_points*num_points/))
			vect_idx_layer = PACK(aux_vect_idx_layer_aux,aux_vect_idx_layer_aux.NE.0)
			!vect_idx_layer = ISORT(SIZE(vect_idx_layer),vect_idx_layer)
			CALL ISORT(SIZE(vect_idx_layer),vect_idx_layer,vect_idx_layer)

			IF (SIZE(vect_idx_layer).EQ.SIZE(target_vect)) THEN
				IF (ALL(vect_idx_layer.EQ.target_vect)) EXIT
			ENDIF
		ENDIF
		IF (n_layer.GT.100) THEN
			WRITE(*,*) 'ERROR in get_layers, infinite loop'
			READ(*,*)
		ENDIF

	ENDDO


	mask = (idx_layer_aux.NE.0)
	nnz_idx_layer = COUNT(mask,DIM=2)
	max_j = MAXVAL(nnz_idx_layer);
	ALLOCATE (idx_layer(COUNT(nnz_idx_layer.NE.0),max_j))
	idx_layer = 0

	cnt=0
	DO i=1,SIZE(idx_layer_aux,1)
		IF (nnz_idx_layer(i).EQ.0) CYCLE
		cnt = cnt +1

		aux_sort = DFLOAT(idx_layer_aux(i,:))
		CALL DSORT (aux_sort,1,SIZE(aux_sort))

		row_idx_layer = NINT(aux_sort)
		!idx_layer(i,1:nnz_idx_layer(i)) = row_idx_layer((max_j+1)-nnz_idx_layer(i):max_j)
		row_idx_layer = PACK(row_idx_layer,row_idx_layer.NE.0)
		auxnum = SIZE(row_idx_layer)
		idx_layer(i,1:auxnum) = row_idx_layer
	ENDDO

	! final test and return
	row_idx_layer = PACK(idx_layer,idx_layer.NE.0)
	!row_idx_layer = ISORT(SIZE(row_idx_layer),row_idx_layer)
	CALL ISORT(SIZE(row_idx_layer),row_idx_layer,row_idx_layer)
	IF (ALL(vect_idx_layer.EQ.target_vect)) THEN
		RETURN
	ELSE
		WRITE(*,*) 'ERROR in get_layers, not all radial points where sorted.'
		READ(*,*)
	ENDIF

	ENDSUBROUTINE get_layers

	! =============================================================
	SUBROUTINE radial_disc_lin(f0,b_RVE,a_RVE,flag_radd,nRad_bar,numel_r,r)
	IMPLICIT NONE
	REAL*8, INTENT(IN)               :: f0, b_RVE
	REAL*8, INTENT(OUT)              :: a_RVE
	INTEGER, INTENT(IN)              :: flag_radd,nRad_bar
	INTEGER, INTENT(OUT)             :: numel_r
	REAL*8, ALLOCATABLE, INTENT(OUT) :: r(:)
	REAL*8 TypRadStp
	INTEGER, ALLOCATABLE :: r_ref_idx(:)

	INTEGER  i, numel_ref

	!IF (.not.YesHOMOG) RETURN

	a_RVE = b_RVE*f0**(1.d0/3.d0); !inner radius
	!PI=4.d0*DATAN(1.d0)
	!Vol_RVE = (4.d0*PI/3.d0)*b_RVE**3.d0;


	! Set radial discretization (TO DO ONCE!!, not adaptative...)
	TypRadStp=b_RVE/DFLOAT(nRad_bar)
	numel_ref = MAX(CEILING((b_RVE-a_RVE)/TypRadStp),3)

	IF (flag_radd.EQ.1) THEN
		IF (FLOOR(dlog10(f0))<=-5) THEN
			numel_r=numel_ref+2
		ELSE
			numel_r=numel_ref+1
		ENDIF

		!numel_r=numel_ref+2
		ALLOCATE (r(numel_r))

		! Point interior to a_radius (to better describe sline BC at r=a_radius)
		r_ref_idx = (/(i,i=0,numel_ref-1)/)

		IF (FLOOR(dlog10(f0))<=-5) THEN
			r(1)=0.975d0*a_RVE
			r(2)=0.990d0*a_RVE
			r(3:numel_r) = a_RVE+((b_RVE-a_RVE)/DFLOAT(numel_ref-1))*DFLOAT(r_ref_idx(:))
		ELSE
			r(1)=0.95d0*a_RVE
			r(2:numel_r) = a_RVE+((b_RVE-a_RVE)/DFLOAT(numel_ref-1))*DFLOAT(r_ref_idx(:))
		ENDIF

	ELSE
		numel_r=numel_ref
		ALLOCATE (r(numel_r))
		r_ref_idx = (/(i,i=0,numel_ref-1)/)
		r = a_RVE+((b_RVE-a_RVE)/DFLOAT(numel_ref-1))*DFLOAT(r_ref_idx(:))
	ENDIF

	! Sort r ascending
	!CALL DSORT(r,1,numel_r) ! No need, r is necessarily sorted


	ENDSUBROUTINE radial_disc_lin

	! =============================================================
	SUBROUTINE arc_len_resamp_sig(n,n0,nend,x,y,ngap_uni,nxnew,xnew,test_aymp,idx_last_nz)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n,n0,nend,ngap_uni
	REAL*8, INTENT(IN) :: x(n),y(n)
	INTEGER, INTENT(OUT) :: nxnew
	REAL*8, ALLOCATABLE, INTENT(OUT):: xnew(:)
	LOGICAL, INTENT(OUT):: test_aymp

	REAL*8  x_normlizer, y_normlizer, L_aprox, maxLi_ref
	REAL*8, ALLOCATABLE :: x_norm(:), y_norm(:), diffvec(:,:), arc_leng(:)
	INTEGER, ALLOCATABLE :: testmaxLi(:)
	LOGICAL test_aymp1,test_aymp2
	REAL*8, ALLOCATABLE :: xnew_aux(:)
	INTEGER cnt, igap, ngap
	INTEGER i, ndisc
	REAL*8 goal_arc_norm,xtest
	INTEGER isum, idxmin, idxmin_old, npoints_add
	REAL*8, ALLOCATABLE :: radd_vect(:), maxLi(:)
	LOGICAL, ALLOCATABLE :: idx_muti(:)
	INTEGER firsttrue, lasttrue
	REAL*8 faux, factor
	REAL*8  a_radius,b_radius,ub_factor,lb_factor,tol_diff
	INTEGER, ALLOCATABLE :: maxLi_val_sig(:,:), maxLi_diff_sig(:,:), dummy_Li(:,:)
	LOGICAL test_val_sig,test_diff_sig
	LOGICAL test_asym_sig
	LOGICAL dummy_tf
	LOGICAL, ALLOCATABLE ::	dummy_alloc(:),nz_grad(:),test(:)
	INTEGER flag_job,ny,last_idx
	INTEGER, ALLOCATABLE :: maxLi_devsig(:),both_maxLi(:,:)
	INTEGER idx_last_nz, idx_last_nz_aux, idx_fst_plat, j, mingap

	nxnew=0
	test_aymp=.false.

	a_radius=x(n0)
	b_radius=x(nend)
	faux = (a_radius/b_radius)**3.d0
	ub_factor=(a_radius+(b_radius-a_radius)*0.25d0)
	lb_factor=a_radius
	tol_diff=1.d-6
	ny=6
	last_idx = FINDLOC(x<=ub_factor.and.x>=lb_factor,.true.,DIM=1,BACK=.true.)

	! sigma_dev (y1) funtion based jump detection:
	! - value based:
	flag_job=1
	factor=1.75d0
	CALL arc_based_detect_step(flag_job,n,n0,x,ny,y,faux,ngap_uni,lb_factor,ub_factor,factor,tol_diff,&
		test_val_sig,maxLi_val_sig,dummy_tf,dummy_Li,dummy_alloc)
	! - diff (dy/dx) based:
	flag_job=2
	factor=1.0d0
	CALL arc_based_detect_step(flag_job,n,n0,x,ny,y,faux,ngap_uni,lb_factor,ub_factor,factor,tol_diff,&
		dummy_tf,dummy_Li,test_diff_sig,maxLi_diff_sig,nz_grad)


	test_asym_sig=.false.
	ALLOCATE(maxLi_devsig(n-2))
	ALLOCATE(both_maxLi(n-2,12))
	maxLi_devsig=0
	IF (test_val_sig.OR.test_diff_sig) THEN
		test_asym_sig=.true.
		maxLi_diff_sig=SPREAD(NINT(SUM(maxLi_diff_sig,DIM=2)/DFLOAT(6)),2,6)
		maxLi_diff_sig = MERGE(2,maxLi_diff_sig,maxLi_diff_sig>2)
		maxLi_diff_sig(last_idx:n-2,:)=0
		both_maxLi(:,1:6)=maxLi_val_sig(1:n-2,:)
		both_maxLi(:,7:12)=maxLi_diff_sig
		maxLi_devsig = MAXVAL(both_maxLi,DIM=2)
	ENDIF

	test_aymp=test_asym_sig
	IF (.not.test_aymp) THEN
		nxnew=0
		xnew=0.d0
		idx_last_nz=0
		RETURN
	ENDIF


	IF (test_asym_sig.OR.COUNT(nz_grad(1:last_idx)).GT.0) THEN
		IF (test_val_sig) THEN
			mingap=3
		ELSE
			mingap=2
		ENDIF
		idx_last_nz=0
		DO j=1,1 !6 (equal cols)
			idx_last_nz_aux=FINDLOC(maxLi_diff_sig(1:last_idx,j).ne.0,.true.,DIM=1,BACK=.true.)
			IF (COUNT(maxLi_diff_sig(1:idx_last_nz_aux,j).ne.0).EQ.idx_last_nz_aux) THEN
				idx_last_nz=MAX(idx_last_nz,FINDLOC(maxLi_diff_sig(1:last_idx,j).gt.0,.true.,DIM=1,BACK=.true.))
			ENDIF
		ENDDO
		idx_fst_plat = FINDLOC(nz_grad(1:last_idx),.true.,DIM=1,BACK=.true.)+1
		ALLOCATE(test(n-2))
		test=.false.
		IF (MAX(idx_fst_plat,idx_last_nz).LE.n-2) THEN
			test(1:MAX(idx_fst_plat,idx_last_nz))=NINT(DFLOAT(SUM(maxLi_val_sig(1:MAX(idx_fst_plat,idx_last_nz),:),DIM=2))/DFLOAT(6))<=1
		ENDIF
		maxLi_devsig  = MERGE(1,maxLi_devsig,test)

		IF (idx_fst_plat.NE.1) 	idx_fst_plat=MAX(idx_fst_plat,3)
		IF (idx_last_nz.NE.0) 	idx_last_nz=MAX(idx_last_nz,2)

		! add points in the jump and near this gap
		IF (idx_last_nz.NE.0) THEN
			IF (idx_last_nz+1>=1.and.idx_last_nz+1<=n-2) maxLi_devsig(idx_last_nz+1)=MAX(maxLi_devsig(idx_last_nz+1)+1,2)
			IF (idx_last_nz>=1.and.idx_last_nz<=n-2) maxLi_devsig(idx_last_nz)=MAX(maxLi_devsig(idx_last_nz)+1,mingap)
			IF (idx_last_nz>1.and.idx_last_nz<=n-2) maxLi_devsig(idx_last_nz-1)=MAX(maxLi_devsig(idx_last_nz-1)+1,2)
		ELSE
			IF (idx_fst_plat>1) maxLi_devsig(idx_fst_plat-1)=MAX(maxLi_devsig(idx_fst_plat-1),2)
			IF (idx_fst_plat>2) maxLi_devsig(idx_fst_plat-2)=MAX(maxLi_devsig(idx_fst_plat-2),mingap)
			IF (idx_fst_plat<=n-2) maxLi_devsig(idx_fst_plat)=MAX(maxLi_devsig(idx_fst_plat),2)
		ENDIF
		idx_last_nz=MAX(1,MAX(idx_last_nz+1,idx_fst_plat))
	ENDIF


	! Add points based on a quasi-uniform arc lenght partition
	ALLOCATE(xnew_aux(ngap_uni*n))
	xnew_aux=0.d0

	cnt=0
	DO igap=n0,last_idx!nend-1
		ngap=maxLi_devsig(igap)
		IF (ngap<=1) CYCLE
		ngap=MIN(ngap,MAX(2,NINT(DFLOAT(ngap_uni)/4.0d0)))

		npoints_add=ngap-1
		IF (ALLOCATED(radd_vect)) DEALLOCATE(radd_vect)
		radd_vect = x(igap)+((x(igap+1)-x(igap))/DFLOAT(ngap))*DFLOAT((/(i,i=0,ngap)/))
		xnew_aux(cnt+1:cnt+npoints_add)=radd_vect(2:ngap)
		cnt=cnt+npoints_add
	ENDDO

	xnew=xnew_aux(1:cnt)
	nxnew=cnt

	! Verify uniqueness and sort ascending: implicit in code
	! CALL DSORT(xnew,1,cnt)


	ENDSUBROUTINE arc_len_resamp_sig

	! =============================================================
	SUBROUTINE arc_len_resamp(n,n0,nend,x,y,yr2,ngap_uni,nxnew,xnew,test_aymp)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n,n0,nend,ngap_uni
	REAL*8, INTENT(IN) :: x(n),y(n),yr2(n)
	INTEGER, INTENT(OUT) :: nxnew
	REAL*8, ALLOCATABLE, INTENT(OUT):: xnew(:)
	LOGICAL, INTENT(OUT):: test_aymp

	REAL*8  x_normlizer, y_normlizer, L_aprox, maxLi_ref
	REAL*8, ALLOCATABLE :: x_norm(:), y_norm(:), diffvec(:,:), arc_leng(:)
	INTEGER, ALLOCATABLE :: testmaxLi(:)
	LOGICAL test_aymp1,test_aymp2
	REAL*8, ALLOCATABLE :: xnew_aux(:)
	INTEGER cnt, igap, ngap
	REAL*8 x0,x1,x2, y0,y1, apower,ppower
	!INTEGER, PARAMETER :: ndisc=1000
	INTEGER i !, xq_vect(ndisc)
	!REAL*8 xq(ndisc),vq(ndisc),xq_norm(ndisc),vq_norm(ndisc)
	!REAL*8 diffvec_xqvq(ndisc-1,2),L_gap_norm(ndisc-1)
	REAL*8 goal_arc_norm,xtest !,cumsumL_norm(ndisc-1)

	INTEGER ndisc
	INTEGER, ALLOCATABLE :: xq_vect(:)
	REAL*8, ALLOCATABLE ::	xq(:),vq(:),xq_norm(:),vq_norm(:)
	REAL*8, ALLOCATABLE :: diffvec_xqvq(:,:),L_gap_norm(:),cumsumL_norm(:)

	INTEGER isum, idxmin, idxmin_old, npoints_add
	REAL*8, ALLOCATABLE :: radd_vect(:), maxLi(:)
	LOGICAL, ALLOCATABLE :: idx_muti(:)
	INTEGER firsttrue, lasttrue

	REAL*8 faux, factor
	REAL*8  a_radius,b_radius,ub_factor,lb_factor,tol_diff
	INTEGER, ALLOCATABLE :: maxLi_val_srp(:,:), maxLi_diff_srp(:,:)
	INTEGER, ALLOCATABLE :: maxLi_val_srpr2(:,:), maxLi_diff_srpr2(:,:)
	LOGICAL test_val_srp,test_diff_srp
	LOGICAL test_val_srpr2,test_diff_srpr2
	LOGICAL test_asym_srp,test_asym_srpr2
	INTEGER max_gap_srp,max_gap_srpr2
	INTEGER mingap,flag_job,ny,last_idx
	LOGICAL, ALLOCATABLE ::	dummy_alloc(:)

	nxnew=0
	test_aymp=.false.

	a_radius=x(n0)
	b_radius=x(nend)
	faux = (a_radius/b_radius)**3.d0
	ub_factor=(a_radius+(b_radius-a_radius)*0.25d0)
	lb_factor=a_radius
	tol_diff=1.d-8
	flag_job=0	! both
	factor = 1.d0
	ny=1
	! SRP (y1) funtion based jump detection
	CALL arc_based_detect_step(flag_job,n,n0,x,ny,y,faux,ngap_uni,lb_factor,ub_factor,factor,tol_diff,&
		test_val_srp,maxLi_val_srp,test_diff_srp,maxLi_diff_srp,dummy_alloc)
	test_asym_srp=.false.
	max_gap_srp=0
	IF (test_val_srp.or.test_diff_srp) THEN
		test_asym_srp=.true.
		max_gap_srp=MAX(MAXVAL(maxLi_val_srp),NINT(MAXVAL(0.5d0*maxLi_diff_srp)))
	ENDIF

	! Integrand (~y1*r^2) funtion based jump detection
	CALL arc_based_detect_step(flag_job,n,n0,x,ny,yr2,faux,ngap_uni,lb_factor,ub_factor,factor,tol_diff,&
		test_val_srpr2,maxLi_val_srpr2,test_diff_srpr2,maxLi_diff_srpr2,dummy_alloc)
	test_asym_srpr2=.false.
	max_gap_srpr2=0
	IF (test_val_srpr2.or.test_diff_srpr2) THEN
		test_asym_srpr2=.true.
		max_gap_srpr2=MAX(MAXVAL(maxLi_val_srpr2),NINT(MAXVAL(0.5d0*maxLi_diff_srpr2)))
	ENDIF


	! Add points based on a quasi-uniform arc lenght partition
	ALLOCATE(xnew_aux(ngap_uni*n))
	xnew_aux=0.d0
	test_aymp1=.false.
	test_aymp2=.false.

	! power regression of points near r=a, i.e.
	! is y2~srp*r^2 it an asymptote y~a*(1/x)?
	x0=x(1)
	x1=x(2)
	y0=yr2(1)
	y1=yr2(2)
	apower=y0*(y1/y0)**(DLOG(x0)/(DLOG(x0/x1)))
	ppower=-DLOG(y1/y0)/(DLOG(x0/x1))
	test_aymp1 = (DABS(ppower-(-1.d0)).LT.1.d-1)
	IF (.not.test_aymp1) test_aymp2 = (DABS(apower/x(3)-y(3))/y(3).LT.5.d-2)

	test_aymp=(test_aymp1.OR.test_aymp2.OR.test_asym_srp.or.test_asym_srpr2)
	IF (.not.test_aymp) THEN
		nxnew=0
		xnew=0.d0
		RETURN
	ENDIF

	x1=x(n0); !a_radius
	x2=x(n0+1)
	x_normlizer=MAXVAL(x)
	y_normlizer=MAXVAL(yr2)
	IF (FLOOR(dlog10(faux))<=-5) THEN
		mingap=5
	ELSEIF ((FLOOR(log10(faux))>-5).AND.(FLOOR(log10(faux))<=-4)) THEN
		mingap=4
	ELSE
		mingap=2
	ENDIF
	cnt=0
	last_idx = FINDLOC(x<=ub_factor.and.x>=lb_factor,.true.,DIM=1,BACK=.true.)

	DO igap=n0,last_idx!nend-1
		ngap=maxLi_val_srpr2(igap,1)
		IF (igap==n0) THEN
			IF (test_aymp) THEN
				! it is an asymptote
				IF (test_asym_srp.AND.((.not.test_asym_srpr2).or.(.not.test_aymp1))) THEN
					ndisc=1000
				ELSE
					ndisc=100
				ENDIF
				ALLOCATE(xq_vect(ndisc))
				ALLOCATE(xq(ndisc),vq(ndisc),xq_norm(ndisc),vq_norm(ndisc))
				ALLOCATE(diffvec_xqvq(ndisc-1,2),L_gap_norm(ndisc-1),cumsumL_norm(ndisc-1))
				xq_vect=(/(i,i=0,ndisc-1)/)
				xq=x1+((x2-x1)/DFLOAT(ndisc-1))*DFLOAT(xq_vect(:))
				vq=apower/xq
				xq_norm=xq/x_normlizer
				vq_norm=vq/y_normlizer
				diffvec_xqvq(:,1)=	xq_norm(1:ndisc-1)-xq_norm(2:ndisc)
				diffvec_xqvq(:,2)=vq_norm(1:ndisc-1)-vq_norm(2:ndisc)
				L_gap_norm=DSQRT(SUM(diffvec_xqvq*diffvec_xqvq,DIM=2))
				IF (test_asym_srp.AND.((.not.test_asym_srpr2).or.(.not.test_aymp1))) THEN
					ngap=MAX(ngap,MIN(max_gap_srp,NINT(DFLOAT(ngap_uni)/2.5d0)))
				ELSE
					ngap=MIN(ngap,NINT(DFLOAT(ngap_uni)/2.0d0))
				ENDIF
				ngap=MAX(ngap,mingap)
				goal_arc_norm=MAX(SUM(L_gap_norm)/DFLOAT(ngap),MAXVAL(L_gap_norm))
				CALL cumsum(L_gap_norm,cumsumL_norm)
				idxmin_old=0
				DO isum=1,(ngap-1)
					idxmin = MINLOC(DABS(cumsumL_norm-DFLOAT(isum)*goal_arc_norm),1)
					xtest=xq(idxmin)
					IF ((idxmin.NE.idxmin_old).AND.(xtest.GT.x1).AND.(xtest.LT.x2)& !THEN !&
						.AND.((DABS(xtest-x2)/DABS(x1-x2)).GT.0.01d0)) THEN
						cnt=cnt+1
						xnew_aux(cnt)=xtest
						idxmin_old=idxmin
					ENDIF
				ENDDO

			ENDIF
		ELSE ! not the first gap (near r=a_radius)
			IF ((igap.EQ.(n0+1)).AND.(cnt.GT.0).AND.(DABS(xnew_aux(cnt)-x2)/DABS(x1-x2).LT.0.25d0)) THEN
				ngap=MAX(ngap,2)
			ENDIF
			IF ((ngap.GT.1).AND.(x(igap+1)<ub_factor).AND.(test_aymp)) THEN
				npoints_add=ngap-1
				IF (ALLOCATED(radd_vect)) DEALLOCATE(radd_vect)
				radd_vect = x(igap)+((x(igap+1)-x(igap))/DFLOAT(ngap))*DFLOAT((/(i,i=0,ngap)/))
				xnew_aux(cnt+1:cnt+npoints_add)=radd_vect(2:ngap)
				cnt=cnt+npoints_add
			ENDIF
		ENDIF
	ENDDO

	xnew=xnew_aux(1:cnt)
	nxnew=cnt

	! Verify uniqueness and sort ascending: implicit in code
	! CALL DSORT(xnew,1,cnt)


	ENDSUBROUTINE arc_len_resamp
	! =============================================================
	!FUNCTION ISORT(n,a)
	!INTEGER, INTENT(IN)  :: n
	!INTEGER, INTENT(IN)  :: a(n)
	!INTEGER  ISORT(n)
	!
	!INTEGER :: acopy(n),i,j,tmp
	!
	!acopy=a
	!do i=1,n
	!	do j=i+1,n
	!		if (acopy(j) < acopy(i)) then
	!			tmp = acopy(j)
	!			acopy(j) = acopy(i)
	!			acopy(i) = tmp
	!		end if
	!	end do
	!end do
	!
	!ISORT = acopy
	!ENDFUNCTION ISORT
	!
	!! =============================================================
	SUBROUTINE cumsum(a,b)
	REAL*8, INTENT(IN) :: a(:)
	REAL*8, INTENT(OUT):: b(SIZE(a))
	INTEGER :: i
	b(1)=a(1)
	DO i=2,SIZE(a)
		b(i) = b(i-1) + a(i)
	ENDDO
	ENDSUBROUTINE cumsum
	! =============================================================

	SUBROUTINE arc_based_detect_step(flag_job,n,n0,x,ny,y,f,nref,lb_factor,ub_factor,factor,tol_diff,&
		jump_val,eq_gaps_val,jump_diff,eq_gaps_diff,nz_grad)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n,n0,nref,flag_job,ny
	REAL*8, INTENT(IN) :: x(n),y(n,ny),f,lb_factor,ub_factor,tol_diff,factor
	INTEGER, ALLOCATABLE, INTENT(OUT) :: eq_gaps_val(:,:), eq_gaps_diff(:,:)
	LOGICAL, INTENT(OUT) :: jump_val,jump_diff
	LOGICAL, ALLOCATABLE, INTENT(OUT) :: nz_grad(:)

	REAL*8  x_normlizer, y_normlizer, L_aprox, maxLi_ref
	REAL*8 x_norm(n), y_norm(n)
	REAL*8 diffvec(n-1,2), arc_leng(n-1), maxLi(n-1)
	LOGICAL idx_muti(n)
	REAL*8 faux
	INTEGER firsttrue, lasttrue,j
	REAL*8 x_norm_diff(n-1), y_norm_diff(n-1)
	REAL*8 diffvec2(n-2,2),arc_leng2(n-2),maxLi2(n-1),dy_dx(n-1)

	IF (COUNT(flag_job.NE.(/0,1,2/)).EQ.3)	THEN
		WRITE(*,*) 'Wrong flag_job={0,1,2} in arc_based_detect_step'
		READ(*,*)
	ENDIF

	x_normlizer = MAXVAL(x)
	x_norm=x/x_normlizer

	idx_muti = (x<=ub_factor.and.x>=lb_factor)
	!idx_muti(1:n0) = .false.
	faux = (x(n0)/MAXVAL(x))**3.d0
	firsttrue=FINDLOC(idx_muti(n0:n-1),.true.,DIM=1,BACK=.false.)
	lasttrue=FINDLOC(idx_muti(n0:n-1),.true.,DIM=1,BACK=.true.)

	! 1st) Funtion value based
	IF (flag_job.EQ.0.OR.flag_job.EQ.1) THEN
		ALLOCATE(eq_gaps_val(n-1,ny))

		DO j=1,ny
			y_normlizer = MAXVAL(DABS(y(:,j)))
			y_norm=DABS(y(:,j))/y_normlizer
			!ALLOCATE (diffvec(n-1,2))
			diffvec(:,1)=	x_norm(1:n-1)-x_norm(2:n)
			diffvec(:,2)=	y_norm(1:n-1)-y_norm(2:n)
			arc_leng=DSQRT(SUM(diffvec*diffvec,DIM=2))
			L_aprox = SUM(arc_leng)
			maxLi_ref=(L_aprox/numel_ref)
			maxLi = maxLi_ref
			IF (firsttrue.NE.0.AND.lasttrue.NE.0) THEN
				maxLi(firsttrue:lasttrue)=factor*maxLi_ref
			ENDIF
			eq_gaps_val(:,j) = NINT(arc_leng/maxLi)
		ENDDO

		jump_val = .false.
		IF (any(eq_gaps_val(firsttrue:lasttrue,:)>1)) jump_val = .true.

		IF (flag_job.EQ.1) THEN
			jump_diff=.false.
			ALLOCATE(eq_gaps_diff(n-2,ny))
			eq_gaps_diff=0
			ALLOCATE(nz_grad(n-1))
			nz_grad=.false.
			RETURN
		ENDIF
	ENDIF

	! 2nd) diff based
	IF (flag_job.EQ.0.OR.flag_job.EQ.2) THEN
		ALLOCATE(eq_gaps_diff(n-2,ny))

		x_norm_diff=dabs(x_norm(n0+1:n-n0)-x_norm(n0:n-1))

		ALLOCATE(nz_grad(n-1))
		nz_grad=.false.



		!x_norm_diff=x_norm_diff/MAXVAL(x_norm_diff)
		DO j=1,ny

			y_normlizer = MAXVAL(DABS(y(:,j)))
			y_norm=DABS(y(:,j))/y_normlizer

			y_norm_diff=(y_norm(n0+1:n)-y_norm(n0:n-1))

			IF (MAXVAL(DABS(y_norm_diff))<=tol_diff) THEN
				jump_diff=.false.
				eq_gaps_diff(:,j)=0
				CYCLE
			ENDIF

			dy_dx=y_norm_diff/x_norm_diff
			nz_grad = (dabs(dy_dx)/MAXVAL(dabs(dy_dx))>1.0d-3.or.nz_grad)

			y_norm_diff=dy_dx
			y_norm_diff=y_norm_diff/MAXVAL(dabs(y_norm_diff))

			diffvec2(:,1)=	x_norm(1:n-2)-x_norm(2:n-1)
			diffvec2(:,2)=	y_norm_diff(1:n-2)-y_norm_diff(2:n-1)
			arc_leng2=DSQRT(SUM(diffvec2*diffvec2,DIM=2))
			L_aprox = SUM(arc_leng2)
			maxLi_ref=(L_aprox/numel_ref)
			maxLi2 = maxLi_ref
			IF (firsttrue.NE.0.AND.lasttrue.NE.0) THEN
				maxLi2(firsttrue:lasttrue)=factor*maxLi_ref
			ENDIF

			jump_diff = .false.
			eq_gaps_diff(:,j)=NINT(arc_leng2/maxLi2)
		ENDDO
		IF (any(eq_gaps_diff(firsttrue:lasttrue,:)>1)) jump_diff = .true.

		IF (flag_job.EQ.2) THEN
			jump_val=.false.
			ALLOCATE(eq_gaps_val(n-2,ny))
			eq_gaps_val=0
			RETURN
		ENDIF
	ENDIF

	ENDSUBROUTINE arc_based_detect_step
	! =============================================================
	ENDSUBROUTINE compHomog_commun

	! =============================================================


	SUBROUTINE local_fields(D6,d5_ast,f0,propMAT,N_design,sphtdesign,numel_r,r,YesEshelby,&
		size_idx_layer,idx_layer,first_layer,last_layer,&
		srp,sigma_dev,dir_sigma_dev,niter_all)

	USE M_COMP_HOMOG, ONLY : struct_mat_param,struct_sphtdesign,&
		num_ex, Y0_matrix

	IMPLICIT NONE
	INTEGER, INTENT(IN)                  :: numel_r, N_design, size_idx_layer(2)
	INTEGER, INTENT(IN)                  :: idx_layer(size_idx_layer(1),size_idx_layer(2)),&
		first_layer,last_layer
	REAL*8, INTENT(IN)                   :: D6(6), f0, d5_ast(5),r(numel_r)
	TYPE (struct_mat_param),  INTENT(IN) :: propMAT
	TYPE (struct_sphtdesign), INTENT(IN) :: sphtdesign
	LOGICAL, INTENT(IN)                  :: YesEshelby
	REAL*8,  INTENT(INOUT)                 :: srp(numel_r,N_design), &
		sigma_dev(numel_r,N_design,6), &
		dir_sigma_dev(numel_r,N_design,6)
	INTEGER, INTENT(OUT) :: niter_all

	REAL*8  theta(N_design), phi(N_design)
	REAL*8 a, k(num_ex), L(6,6,num_ex), m
	LOGICAL flag_sym
	REAL*8, ALLOCATABLE :: design_sph(:,:)
	INTEGER  N_design_eff


	INTEGER j, i_idx, j_idx
	REAL*8  r_ij, d_local(6)
	REAL*8  sigma_x0(6),dir_sdev(6),local_srp,eqs_dir_sdev

	INTEGER niter
	INTEGER, ALLOCATABLE :: idx_to_interp(:), idx_avail(:)
	LOGICAL, ALLOCATABLE :: mask(:,:), masknotNaN(:),mask_vect(:)
	REAL*8, ALLOCATABLE ::  xinterp(:), yinterp(:)
	INTEGER, ALLOCATABLE :: next_layer(:)
	REAL*8, ALLOCATABLE :: sigma_x0_next(:,:)
	REAL*8, ALLOCATABLE :: x_all(:),y_all(:,:)
	INTEGER npts_delta, nyinterp, idx_previous, aux_sizeidx_avail, aux_countnonnan
	REAL*8 	aux_nantest

	! Unpack material parameters
	a = propMAT%a;
	k = propMAT%k;
	L = propMAT%L;
	m = propMAT%m;

	! Unpack spherical t-design

	N_design_eff=sphtdesign%N_dsgn_eff
	flag_sym=    sphtdesign%flagsym
	design_sph=  sphtdesign%dsgn_sph
	theta = design_sph(:,2)
	phi = design_sph(:,3)



	niter_all=0
	DO j=1,N_design_eff

		! Build intepolation based on previously computed fields
		IF (first_layer.NE.1) THEN
			! Select data to interpolate, ignoring NaNs (failed convergence)
			idx_previous = (first_layer-1)
			mask = (idx_layer(1:idx_previous,:).NE.0)
			!idx_avail = ISORT(COUNT(mask),PACK(idx_layer(1:idx_previous,:),mask))
			ALLOCATE(idx_avail(COUNT(mask)))
			CALL ISORT(COUNT(mask),PACK(idx_layer(1:idx_previous,:),mask),idx_avail)
			masknotNaN = (srp(idx_avail,j).EQ.srp(idx_avail,j))
			idx_avail = PACK(idx_avail,masknotNaN)
			! Define query points
			mask_vect = (idx_layer(idx_previous+1,:).NE.0)
			!next_layer = ISORT(COUNT(mask_vect),PACK(idx_layer(idx_previous+1,:),mask_vect))
			ALLOCATE(next_layer(COUNT(mask_vect)))
			CALL ISORT(COUNT(mask_vect),PACK(idx_layer(idx_previous+1,:),mask_vect),next_layer)
			x_all=r
			y_all=dir_sigma_dev(:,j,(/1,2,4,5,6/))
			npts_delta=2
			nyinterp=5 ! only 5 independent ppvals
			aux_sizeidx_avail=MAX(0,size(idx_avail))
			aux_countnonnan=COUNT(masknotNaN)
			! MKL spline interpolant
			IF (aux_countnonnan.LT.2.OR.aux_countnonnan.NE.COUNT(mask)) THEN	!(aux_sizeidx_avail.LT.2).OR.
				ALLOCATE(sigma_x0_next(SIZE(idx_to_interp),nyinterp))
				WRITE(*,*) 'here'
				READ(*,*)
				sigma_x0_next(:,4:6)=SPREAD((/0.d0,0.d0,0.d0/),1,SIZE(idx_to_interp))
				sigma_x0_next(:,1:2)=SPREAD(SIGN((/1.d0,1.d0/),d_local(1:2)),1,SIZE(idx_to_interp))
			ELSE
				CALL interp_sigma(aux_sizeidx_avail,idx_avail,size(next_layer),next_layer,&
					x_all,nyinterp,y_all,npts_delta,idx_to_interp,sigma_x0_next)
			ENDIF

			DEALLOCATE(idx_avail,next_layer)
		ENDIF

		DO i_idx=first_layer,last_layer
			DO j_idx=1,size_idx_layer(2)
				IF (idx_layer(i_idx,j_idx).EQ.0) EXIT

				r_ij = r(idx_layer(i_idx,j_idx));

				!Get local strain rate field
				IF (YesEshelby) THEN
					CALL local_d_Eshelby(D6,d5_ast,f0,r_ij,theta(j),phi(j),d_local)
				ELSE
					! Rice and Tracey v.f.
					CALL local_d_RT(D6,r_ij,theta(j),phi(j),d_local)
				ENDIF
				IF (i_idx.EQ.1) THEN
					! Set x0 based on sign(d_local) or previous solutions
					IF (j.GT.1.AND.flag_sym) THEN
						sigma_x0 = dir_sigma_dev(idx_layer(i_idx,j_idx),j-1,:)
					ELSE
						sigma_x0 = 0.d0
						sigma_x0(1:3) = SIGN((/1.d0,1.d0,1.d0/),d_local(1:3))
					ENDIF
				ELSE
					! Set x0 based on ppval
					sigma_x0((/1,2,4,5,6/)) = sigma_x0_next(j_idx,:)
					sigma_x0(3) = -(sigma_x0(1)+sigma_x0(2))
					aux_nantest=sum(sigma_x0)
					IF (aux_nantest.ne.aux_nantest) THEN
						sigma_x0 = 0.d0
						sigma_x0(1:3) = SIGN((/1.d0,1.d0,1.d0/),d_local(1:3))
					ENDIF
				ENDIF

				! Solve for dir_sdev and compute local srp
				CALL get_local_srp(d_local,sigma_x0,a,k,L,m,dir_sdev,eqs_dir_sdev,local_srp,niter)
				niter_all=niter_all+niter
				!WRITE(*,*) niter

				! Data to integrate
				srp(idx_layer(i_idx,j_idx),j) = local_srp
				sigma_dev(idx_layer(i_idx,j_idx),j,:) = (Y0_matrix/eqs_dir_sdev)*dir_sdev

				! Useful to set x0:
				dir_sigma_dev(idx_layer(i_idx,j_idx),j,:) = dir_sdev
			ENDDO

			! Perform interpolation to improve x0
			IF (i_idx.LT.SIZE(idx_layer,1)) THEN
				! Select data to interpolate, ignoring NaNs (failed convergence)
				mask = (idx_layer(1:i_idx,:).NE.0)
				!idx_avail = ISORT(COUNT(mask),PACK(idx_layer(1:i_idx,:),mask))
				ALLOCATE(idx_avail(COUNT(mask)))
				CALL ISORT(COUNT(mask),PACK(idx_layer(1:i_idx,:),mask),idx_avail)
				masknotNaN = (srp(idx_avail,j).EQ.srp(idx_avail,j))
				idx_avail = PACK(idx_avail,masknotNaN)
				! Define query points
				mask_vect = (idx_layer(i_idx+1,:).NE.0)
				!next_layer = ISORT(COUNT(mask_vect),PACK(idx_layer(i_idx+1,:),mask_vect))
				ALLOCATE(next_layer(COUNT(mask_vect)))
				CALL ISORT(COUNT(mask_vect),PACK(idx_layer(i_idx+1,:),mask_vect),next_layer)
				x_all=r
				y_all=dir_sigma_dev(:,j,(/1,2,4,5,6/))
				npts_delta=2
				nyinterp=5 ! only 5 independent ppvals
				! MKL spline interpolant
				CALL interp_sigma(size(idx_avail),idx_avail,size(next_layer),next_layer,&
					x_all,nyinterp,y_all,npts_delta,idx_to_interp,sigma_x0_next)
				DEALLOCATE(idx_avail,next_layer)
			ENDIF
		ENDDO
	ENDDO
	CONTAINS

	! =============================================================
	!
	!FUNCTION ISORT(n,a)
	!INTEGER, INTENT(IN)  :: n
	!INTEGER, INTENT(IN)  :: a(n)
	!INTEGER  ISORT(n)
	!
	!INTEGER :: acopy(n),i,j,tmp
	!
	!acopy=a
	!do i=1,n
	!	do j=i+1,n
	!		if (acopy(j) < acopy(i)) then
	!			tmp = acopy(j)
	!			acopy(j) = acopy(i)
	!			acopy(i) = tmp
	!		end if
	!	end do
	!end do
	!
	!ISORT = acopy
	!ENDFUNCTION ISORT

	! =============================================================

	SUBROUTINE interp_sigma(nx,idx_x,nxq,idx_xq,&
		x,nytyp,y,npts_delta,idx2interp,vq)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: nx, nxq, nytyp,npts_delta
	INTEGER, INTENT(INOUT) :: idx_x(nx), idx_xq(nxq)
	REAL*8, ALLOCATABLE, INTENT(IN) :: x(:), y(:,:)
	INTEGER, ALLOCATABLE, INTENT(OUT) :: idx2interp(:)
	REAL*8, ALLOCATABLE, INTENT(OUT) :: vq(:,:)

	INTEGER num_inf, num_sup
	LOGICAL, ALLOCATABLE :: test_max(:),test_min(:)
	INTEGER idx_test_min,idx_test_max,idx_lim_sup,idx_lim_inf
	REAL*8, ALLOCATABLE ::  xinterp(:), yinterp(:)
	INTEGER nxinterp, nyinterp, flag_ppnat, i

	IF (nx.GE.3) THEN

		! Reduce sample points (x,y) based on query points span
		test_min=idx_x<MINVAL(idx_xq)
		test_max=idx_x>MAXVAL(idx_xq)
		num_inf=COUNT(test_min)
		num_sup=COUNT(test_max)
		idx_lim_inf=1
		idx_lim_sup=SIZE(idx_x)
		IF (num_inf.gt.npts_delta) THEN
			idx_test_min = FINDLOC(test_min,.true.,DIM=1,BACK=.true.)
			idx_lim_inf = (idx_test_min-npts_delta)
		ENDIF
		IF (num_sup.gt.npts_delta) THEN
			idx_test_max = FINDLOC(test_max,.true.,DIM=1)
			idx_lim_sup = (idx_test_max+npts_delta-1)
		ENDIF

		idx2interp=idx_x(idx_lim_inf:idx_lim_sup)
		! Define sample points (x,y)
		nxinterp = SIZE(idx2interp)
		xinterp=x(idx2interp)
		yinterp=[y(idx2interp,:)]

	ELSE
		xinterp=x
		yinterp=[y(:,:)]
		nxinterp=SIZE(x)
	ENDIF

	!IF (nxinterp.LE.3) THEN
	!	WRITE(*,*) nx,idx_x,nxq,idx_xq
	!	WRITE(*,*) x
	!	WRITE(*,*) y
	!	!nxinterp=nx
	!	!xinterp=x
	!	!yinterp=[y(:,:)]
	!ENDIF

	IF (ALLOCATED(vq)) DEALLOCATE (vq)
	ALLOCATE (vq(SIZE(idx_xq),nytyp))
	! Get interpolation
	flag_ppnat=0
	CALL pp_interpolant(nxinterp,xinterp,nytyp,yinterp,flag_ppnat,&
		SIZE(idx_xq),x(idx_xq),vq)


	ENDSUBROUTINE interp_sigma
	ENDSUBROUTINE local_fields
	! =============================================================


	RECURSIVE SUBROUTINE DSORT(a,first,last)
	IMPLICIT NONE
	REAL*8  a(*), x, t
	INTEGER first, last
	INTEGER i, j

	x = a((first+last)/2)
	i = first
	j = last
	DO
		DO WHILE (a(i)<x)
			i=i+1
		ENDDO
		DO WHILE (x<a(j))
			j=j-1
		ENDDO
		IF (i >= j) EXIT
		t = a(i);  a(i)=a(j);a(j) = t
		i=i+1
		j=j-1
	end do
	IF (first<i-1) CALL DSORT(a,first,i-1)
	IF (j+1<last) CALL DSORT(a,j+1,last)

	END SUBROUTINE DSORT
	! =============================================================

	SUBROUTINE ISORT(n,a,asort)
	INTEGER, INTENT(IN)  :: n
	INTEGER, INTENT(IN)  :: a(n)
	INTEGER, INTENT(OUT) :: asort(n)

	INTEGER :: acopy(n),i,j,tmp

	acopy=a
	do i=1,n
		do j=i+1,n
			if (acopy(j) < acopy(i)) then
				tmp = acopy(j)
				acopy(j) = acopy(i)
				acopy(i) = tmp
			end if
		end do
	end do

	asort = acopy
	ENDSUBROUTINE ISORT

	! =============================================================