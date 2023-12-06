
    SUBROUTINE CH_INIT

    USE M_COMP_HOMOG
    USE M_CONSTANTS
    USE M_jit_dgemm, ONLY: create_jit_dgemms

    IMPLICIT NONE

    ! Refined curbature:
    INTEGER                N_design, t_design, N_design_eff
    LOGICAL                flag_sym
    REAL*8, ALLOCATABLE :: design_cart(:,:), design_sph(:,:)
    ! Crude curbature (used in Eshelby):
    INTEGER                N_design_coarse, t_design_coarse, N_design_eff_coarse
    LOGICAL                flag_sym_coarse
    REAL*8, ALLOCATABLE :: design_cart_coarse(:,:), design_sph_coarse(:,:)
    REAL*8 signs(2)

    CALL checkCOMPHOMOG
    IF (.not.YesHOMOG) RETURN

    WRITE(*,*) 'Importing material data and spherical t-design'

    ! Import data from 'DD3oPorous.dat' file
    CALL getCOMPHOMOGdata

    ! Assign material data to structure mat_param
    mat_param%a = a_CPB
    mat_param%k = k_CPB
    mat_param%L = L_CPB
    mat_param%m = m_CPB

    ! Import spherical t-design (curbature rule)
    CALL import_design(authorID,nPtsDesign,N_design,t_design,N_design_eff,&
        flag_sym,design_sph,design_cart)
    sph_tdesign%N_dsgn= N_design
    sph_tdesign%t_dsgn = t_design
    sph_tdesign%N_dsgn_eff = N_design_eff
    sph_tdesign%flagsym = flag_sym
    sph_tdesign%dsgn_sph = design_sph
    sph_tdesign%bar_n_rad = nRad_bar


    !IF (YesEshelby) THEN
    nPtsDesign_coarse = 120
    bar_nrad_crude = 10

    bar_nrad_crude=MIN(nRad_bar,bar_nrad_crude)
    nPtsDesign_coarse=MIN(nPtsDesign,nPtsDesign_coarse)
    CALL import_design(1,nPtsDesign_coarse,N_design_coarse,t_design_coarse,N_design_eff_coarse,&
        flag_sym_coarse,design_sph_coarse,design_cart_coarse)
    !ENDIF
    sph_tdesign_crude%N_dsgn = N_design_coarse
    sph_tdesign_crude%t_dsgn = t_design_coarse
    sph_tdesign_crude%N_dsgn_eff = N_design_eff_coarse
    sph_tdesign_crude%flagsym = flag_sym_coarse
    sph_tdesign_crude%dsgn_sph = design_sph_coarse
    sph_tdesign_crude%bar_n_rad = bar_nrad_crude

    ! Create tailored GEMM kernels for small matrix operations
    CALL create_jit_dgemms


    Vol_RVE = (4.d0*PI/3.d0)*b_radius**3.d0;
    Y0_matrix = 1; !(unitary and constant in the matrix)

    CALL set_penta(Trans56,Trans65)

    ! Initial hydro limits (and compute h)
    ! Constant CPB_mat_param =>  constant h(2) parameter
    ! index label:   1 - compressive, 2 - tensile
    hCPB0(1)=2.d0
    hCPB0(2)=2.d0
    CALL get_aprox_hydro_limits(f_read,mat_param,sph_tdesign,&
        SRP_hydro0,Dm_hydro0,SIGMAm_hydro0)
    signs = (/1.d0,-1.d0/)
    hCPB0=signs*(3.d0/DLOG(f_read))*(SIGMAm_hydro0)

    fmin_adm=1.d-7
    fmax_adm=0.95d0
    


    CONTAINS

    SUBROUTINE import_design(designID,num_points,N_design,t_design,N_design_eff,flag_sym,&
        design_sph,design_cart)

    USE M_COMP_HOMOG, ONLY: IavailDESIGN,IDESIGN
    IMPLICIT NONE

    INTEGER, INTENT(IN)   :: num_points, designID
    INTEGER, INTENT(OUT)  :: N_design,t_design,N_design_eff
    REAL*8, INTENT(OUT), ALLOCATABLE :: design_sph(:,:), design_cart(:,:)
    LOGICAL, INTENT(OUT)  :: flag_sym

    LOGICAL      checkExist, cd_checkExist
    CHARACTER*80 filename, cd_filename
    INTEGER, ALLOCATABLE :: avail_designs(:,:)
    INTEGER  idx_closest
    !INTEGER N_design,t_design,N_design_eff

    IavailDESIGN = 40
    IDESIGN = 41

    ! Get list of available desings:
    SELECT CASE (designID)
    CASE(1) !'WomersleySym'
        filename = 'ss_avail_designs.dat'
    CASE(2) !'WomersleyNonSym'
        filename = 'sf_avail_designs.dat'
    CASE(3) !'Hardin'
        filename = 'hs_avail_designs.dat'
        CASE DEFAULT
        WRITE(*,*) 'designID not valid. Must be one of the following: (1,2,3).'
        filename =''
        READ(*,*)
    ENDSELECT
    cd_filename = ('sphdesigns/' // filename)

    INQUIRE(FILE=TRIM(filename),EXIST=checkExist)
    IF (.NOT.checkExist) THEN
        INQUIRE(FILE=cd_filename,EXIST=cd_checkExist)
    ENDIF
    IF (checkExist) THEN
        OPEN(IavailDESIGN,FILE=TRIM(filename),STATUS='OLD',FORM='FORMATTED')
        CALL read_avail_designs(avail_designs)
    ELSEIF (cd_checkExist) THEN
        OPEN(IavailDESIGN,FILE=TRIM(cd_filename),STATUS='OLD',FORM='FORMATTED')
        CALL read_avail_designs(avail_designs)
    ELSE
        WRITE(*,*) 'Error: avail_designs file missing.'
        READ(*,*)
    ENDIF

    ! Import design points such that 'num_points' is closest to avail_designs(:,2)
    idx_closest = select_design(SIZE(avail_designs,DIM=1),avail_designs,num_points)
    N_design = avail_designs(idx_closest,2);
    t_design = avail_designs(idx_closest,1);
    SELECT CASE (designID)
    CASE(1) !'WomersleySym'
        filename = (TRIM('ss' // str(t_design,'(I0.3)')) // '.' // str(N_design,'(I0.5)'))
        flag_sym=.TRUE.
        N_design_eff = N_design/2
    CASE(2) !'WomersleyNonSym'
        filename = (TRIM('sf' // str(t_design,'(I0.3)')) // '.' // str(N_design,'(I0.5)'))
        flag_sym = .FALSE.
        N_design_eff = N_design
    CASE(3) !'Hardin'
        filename = (TRIM('hs' // str(t_design,'(I0.3)')) // '.' // str(N_design,'(I0.5)'))
        flag_sym = .FALSE.
        N_design_eff = N_design
    ENDSELECT
    cd_filename = ('sphdesigns/' // filename)

    IF (checkExist) THEN
        OPEN(IDESIGN,FILE=TRIM(filename),STATUS='OLD',FORM='FORMATTED')
    ELSEIF (cd_checkExist) THEN
        OPEN(IDESIGN,FILE=TRIM(cd_filename),STATUS='OLD',FORM='FORMATTED')
        CALL read_design(N_design,design_sph,design_cart)
    ENDIF

    ENDSUBROUTINE import_design
    !---------------------------------------

    SUBROUTINE read_design(N_design,design_sph,design_cart)       ! NEW
    !USE M_COMP_HOMOG,       ONLY: IDESIGN, N_design
    IMPLICIT NONE
    REAL*8, ALLOCATABLE :: design_cart(:,:), design_sph(:,:)
    INTEGER i,N_design
    REWIND(IDESIGN)

    ALLOCATE (design_sph(N_design,3))
    ALLOCATE (design_cart(N_design,3))
    DO i=1,N_design
        READ(IDESIGN,'(3(F25.0))',END=100) design_cart(i,:)
    ENDDO
100 CONTINUE

    CALL cart2sph(N_design,design_cart,design_sph)

    ENDSUBROUTINE read_design
    !---------------------------------------

    SUBROUTINE read_avail_designs(avail_designs)
    USE M_COMP_HOMOG,       ONLY: IavailDESIGN
    IMPLICIT NONE
    INTEGER, ALLOCATABLE :: avail_designs(:,:)
    INTEGER  avail_designs_aux(1000,2), cnt

    REWIND(IavailDESIGN)

    cnt=0
    DO WHILE (.TRUE.)
        cnt = cnt+1
        READ(IavailDESIGN,'(2I10)',END=100) avail_designs_aux(cnt,1), avail_designs_aux(cnt,2)
    ENDDO
100 CONTINUE

    ALLOCATE (avail_designs(cnt-1,2))
    avail_designs = avail_designs_aux(1:cnt-1,:)

    ENDSUBROUTINE read_avail_designs
    !---------------------------------------
    FUNCTION select_design(nAvail,avail_designs,num_points)
    IMPLICIT NONE
    INTEGER nAvail, avail_designs(nAvail,2), num_points
    INTEGER select_design, aux_diff(nAvail)

    aux_diff=(ABS(avail_designs(:,2)-num_points))
    select_design = MINLOC(aux_diff,1)
    ENDFUNCTION
    !---------------------------------------
    FUNCTION str(k,formatspec)
    !   "Convert an integer to string."
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k
    CHARACTER(len=6), INTENT(IN) :: formatspec

    CHARACTER(len=20)      str
    WRITE(str, formatspec) k
    str = adjustl(str)
    ENDFUNCTION str
    !---------------------------------------
    SUBROUTINE cart2sph(nPoints,cart_coord,sph_coord)
    USE M_CONSTANTS, ONLY : PI
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: nPoints
    REAL*8, INTENT(IN)   :: cart_coord(nPoints,3)
    REAL*8, INTENT(OUT)  :: sph_coord(nPoints,3)
    REAL*8  x(nPoints), y(nPoints), z(nPoints)!, PI

    !spherical = @(x,y,z) [(x.^2 + y.^2 + z.^2).^(1/2), atan2((x.^2 + y.^2).^(1/2),z), mod(atan2(y,x),2*pi)]; % (r,theta,phi)
    !% note: mod(atan2(y,x),2*pi) returns [0,2*pi], while atan2(y,x) returns [-pi,pi];

    x = cart_coord(:,1)+1e-16 !... to avoid atan2(0,0) ~= 0
    y = cart_coord(:,2)
    z = cart_coord(:,3)

    !PI=4.D0*DATAN(1.D0)

    sph_coord(:,1) = DSQRT(SUM(cart_coord*cart_coord,2)) !r
    sph_coord(:,2) = DATAN2(DSQRT(x*x+y*y),z) !theta [0,pi]
    sph_coord(:,3) = MODULO(DATAN2(y,x),2.d0*PI) !phi:  [0,2*pi]

    ENDSUBROUTINE cart2sph

    !---------------------------------------
    SUBROUTINE sph2cart(nPoints,sph_coord,cart_coord)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: nPoints
    REAL*8, INTENT(IN)   :: sph_coord(nPoints,3)
    REAL*8, INTENT(OUT)  :: cart_coord(nPoints,3)
    REAL*8  r(nPoints), theta(nPoints), phi(nPoints)

    !cartesian = @(r,theta,phi) [r.*sin(theta).*cos(phi), r.*sin(theta).*sin(phi), r.*cos(theta)];
    r = sph_coord(:,1)
    theta = sph_coord(:,2)
    phi = sph_coord(:,3)

    cart_coord(:,1) = r*DSIN(theta)*DCOS(phi)
    cart_coord(:,2) = r*DSIN(theta)*DSIN(phi)
    cart_coord(:,3) = r*DCOS(theta)
    ENDSUBROUTINE sph2cart

    ENDSUBROUTINE CH_INIT
    !---------------------------------------
    SUBROUTINE checkCOMPHOMOG

    USE M_COMP_HOMOG
    IMPLICIT NONE

    !Verifies if the external file "DD3oPorous.dat" exists
    YesHOMOG=.FALSE.

    CHfilename = 'matmodel.dat'
    INQUIRE(FILE=TRIM(CHfilename),EXIST=YesHOMOG)
    IF (.not.YesHOMOG) RETURN
    ICOMPHOMOG = 42

    ENDSUBROUTINE checkCOMPHOMOG

    SUBROUTINE getCOMPHOMOGdata
    USE M_COMP_HOMOG
    IMPLICIT NONE
    CHARACTER*80 TITL1,TITL2,T80
    REAL*8       C_aux(6,6)
    INTEGER      I,J

    IF (.not.YesHOMOG) RETURN

    OPEN(ICOMPHOMOG,FILE=TRIM(CHfilename),STATUS='OLD',FORM='FORMATTED')

    READ(ICOMPHOMOG,'(A80)') TITL1
    READ(ICOMPHOMOG,'(A80)') TITL2
    READ(ICOMPHOMOG,'(/)')
    READ(ICOMPHOMOG,'(20X,2I10)')   authorID,nPtsDesign
    READ(ICOMPHOMOG,'(/)')
    READ(ICOMPHOMOG,'(20X,F10.0,I10,F10.0)') Tol,nRad_bar!,xMaxDiff ! NEW!
    READ(ICOMPHOMOG,'(/)')
    READ(ICOMPHOMOG,'(20X,F10.0)')  f_read
    READ(ICOMPHOMOG,'(/)')
    READ(ICOMPHOMOG,'(20X,I10)')    num_ex
    READ(ICOMPHOMOG,'(/)')

    !nRad_bar=20
    !f_read=1.d-6
    !WRITE(*,*) 'f_read = ',f_read

    ! MATERIAL:
    CALL ALLOC_CPB_MATER(num_ex)
    IF (num_ex.EQ.1) THEN
        READ(ICOMPHOMOG,'(20X,2F10.0)')  a_CPB, (k_CPB(I),I=1,num_ex)
    ELSEIF (num_ex.EQ.2) THEN
        READ(ICOMPHOMOG,'(20X,3F10.0)')  a_CPB, (k_CPB(I),I=1,num_ex)
    ELSEIF (num_ex.EQ.3) THEN
        READ(ICOMPHOMOG,'(20X,4F10.0)')  a_CPB, (k_CPB(I),I=1,num_ex)
    ELSEIF (num_ex.EQ.4) THEN
        READ(ICOMPHOMOG,'(20X,5F10.0)')  a_CPB, (k_CPB(I),I=1,num_ex)
    ELSE
        WRITE(*,*) 'Error reading material data. Max number of linear transf: 4'
    ENDIF
    DO J=1,num_ex
        READ(ICOMPHOMOG,'(/20X,6F10.0)') (C_aux(I,I),I=1,6)
        READ(ICOMPHOMOG,'(/20X,6F10.0)') C_aux(2,3),C_aux(1,3),C_aux(1,2)
        C_CPB(:,:,J)=CijsymmTensor(6,C_aux)
        C_aux = 0.d0
    ENDDO

    CLOSE(ICOMPHOMOG)


    CALL CPBexn_param(num_ex,k_CPB,C_CPB,a_CPB,m_CPB,L_CPB)

    CONTAINS
    !---------------------------------------

    FUNCTION CijsymmTensor(n,C_aux)
    IMPLICIT NONE
    INTEGER  n, i
    REAL*8   C_aux(n,n), CijsymmTensor(n,n)


    CijsymmTensor =  (C_aux+TRANSPOSE(C_aux))
    DO i=1,n
        CijsymmTensor(i,i) = CijsymmTensor(i,i)/2.d0
    ENDDO
    ENDFUNCTION

    ENDSUBROUTINE getCOMPHOMOGdata
    !---------------------------------------

    SUBROUTINE CPBexn_param(num_ex,k,C,a,m,L)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: num_ex
    REAL*8, INTENT(IN)  :: k(num_ex), C(6,6,num_ex), a
    REAL*8, INTENT(OUT) :: m, L(6,6,num_ex)
    REAL*8  Kdev(6,6), R13, R23
    INTEGER i
    REAL*8  phi(3,num_ex), aux_denom

    L=0.d0
    m=0.d0

    ! 4th order deviatoric unit tensor
    Kdev=0.d0
    R13 = (1.d0/3.d0)
    R23 = (2.d0/3.d0)
    Kdev(1,1:3) = (/R23,-R13,-R13/)
    Kdev(2,1:3) = (/-R13,R23,-R13/)
    Kdev(3,1:3) = (/-R13,-R13,R23/)
    Kdev(4,4) = 1.d0
    Kdev(5,5) = 1.d0
    Kdev(6,6) = 1.d0

    DO i=1,num_ex
        L(:,:,i)=MATMUL(C(:,:,i),Kdev);
        phi(:,i)=L(1:3,1,i);
    ENDDO

    aux_denom = 0.d0
    DO i=1,num_ex
        aux_denom = aux_denom  &
            + ((DABS(phi(1,i))-k(i)*phi(1,i))**a+(DABS(phi(2,i))-k(i)*phi(2,i))**a+(DABS(phi(3,i))-k(i)*phi(3,i))**a);
    ENDDO

    m=(1/aux_denom)**(1/a)

    ENDSUBROUTINE CPBexn_param
    !---------------------------------------


    SUBROUTINE checkLoading

    USE M_COMP_HOMOG
    IMPLICIT NONE

    YesLoad=.false.
    loading_filename = 'loading.dat'
    INQUIRE(FILE=TRIM(loading_filename),EXIST=YesLoad)

    IF (.not.YesLoad) RETURN
    ILoading=55

    ENDSUBROUTINE checkLoading

    SUBROUTINE getLoadingdata(D,SIG)
    USE M_COMP_HOMOG
    IMPLICIT NONE
    CHARACTER*80 TITL1,TITL2,T80
    REAL*8       D(6), SIG(6)
    INTEGER      I,J

    CALL checkLoading
    IF (.not.YesLoad) RETURN

    OPEN(ILoading,FILE=TRIM(loading_filename),STATUS='OLD',FORM='FORMATTED')

    READ(ILoading,'(A80)') TITL1
    READ(ILoading,'(A80)') TITL2
    READ(ILoading,'(/)')
    READ(ILoading,'(/)')
    READ(ILoading,'(20X,6F10.0)')   (D(I),I=1,6)
    READ(ILoading,'(/)')
    READ(ILoading,'(/)')
    READ(ILoading,'(20X,6F10.0)')   (SIG(I),I=1,6)
    !READ(ILoading,'(/)')

    CLOSE(ILoading)

    ENDSUBROUTINE getLoadingdata