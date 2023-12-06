
    MODULE M_COMP_HOMOG

    ! MODULE == Workspace in Matlab

    IMPLICIT NONE

    CHARACTER*80           CHfilename, loading_filename
    LOGICAL                YesHOMOG, YesLoad
    ! Read from input file:
    INTEGER                authorID, nPtsDesign, nPtsDesign_coarse
    REAL*8                 Tol, f_read, fmin_adm, fmax_adm 
    INTEGER                nRad_bar, bar_nrad_crude
    ! Matrix material data
    INTEGER                num_ex
    REAL*8, ALLOCATABLE :: C_CPB(:,:,:), k_CPB(:), L_CPB(:,:,:)
    REAL*8                 a_CPB, m_CPB
    ! I/O
    INTEGER                IavailDESIGN, IDESIGN, ICompHomog, ILoading
    ! Eshelby flag:
    !LOGICAL                 YesEshelby
    ! Spherical t-designs (surface curbature)
    TYPE struct_sphtdesign
        INTEGER              :: N_dsgn, t_dsgn, N_dsgn_eff, bar_n_rad
        LOGICAL              :: flagsym
        REAL*8, ALLOCATABLE ::  dsgn_sph(:,:)!,dsgn_sph_cart(:,:)
    ENDTYPE
    TYPE (struct_sphtdesign) :: sph_tdesign, sph_tdesign_crude
    ! Material data structure
    TYPE struct_mat_param
        REAL*8              :: a
        REAL*8, ALLOCATABLE :: k(:)
        REAL*8, ALLOCATABLE :: L(:,:,:)
        REAL*8              :: m
    ENDTYPE struct_mat_param
    TYPE (struct_mat_param) :: mat_param
    ! RVE geometry
    REAL*8 a_radius, Vol_RVE, Y0_matrix
    REAL*8, PARAMETER :: b_radius = 1.d0 !,PI=4.D0*DATAN(1.D0)

    ! Initial hydrostatic limits (aproximation)
    REAL*8  SRP_hydro0(2), Dm_hydro0(2),SIGMAm_hydro0(2), hCPB0(2)

    CONTAINS

    SUBROUTINE ALLOC_CPB_MATER(num_ex)

    IMPLICIT NONE
    INTEGER  num_ex

	 WRITE(*,*) ALLOCATED(C_CPB)
	 IF (ALLOCATED(C_CPB)) DEALLOCATE(C_CPB)
	 IF (ALLOCATED(L_CPB)) DEALLOCATE(L_CPB)
	 IF (ALLOCATED(k_CPB)) DEALLOCATE(k_CPB)

    ALLOCATE(C_CPB(6,6,num_ex))
	 C_CPB=    0.d0
    ALLOCATE(L_CPB(6,6,num_ex))
	 L_CPB=    0.d0
    ALLOCATE(k_CPB(num_ex))
	 k_CPB=    0.d0

    ENDSUBROUTINE ALLOC_CPB_MATER

    ENDMODULE M_COMP_HOMOG

    MODULE M_CONSTANTS
    IMPLICIT NONE

    REAL*8, PARAMETER :: PI=4.D0*DATAN(1.D0)
    REAL*8, PARAMETER :: DSQRT2=DSQRT(2.d0)
    REAL*8, PARAMETER :: DSQRT3=DSQRT(3.d0)

    REAL*8 Trans56(5,6), Trans65(6,5)
    REAL*8 eye3(3,3),eye6(6)
    REAL*8 I2(6), I6(6)


    CONTAINS

    SUBROUTINE set_penta(Trans56,Trans65)
    REAL*8 Trans56(5,6), Trans65(6,5)

    Trans56 = 0.d0
    Trans56(1,1:2) = (/1/DSQRT2, -1/DSQRT2/)
    !Trans56(2,1:2) = (/DSQRT3/DSQRT2, DSQRT3/DSQRT2/)
    Trans56(2,1:3) = (/1/(DSQRT2*DSQRT3),1/(DSQRT2*DSQRT3),-2/(DSQRT2*DSQRT3)/)
    Trans56(3,4) = DSQRT2
    Trans56(4,5) = DSQRT2
    Trans56(5,6) = DSQRT2

    Trans65=0.d0
    Trans65(1,1:2) = (/1/DSQRT2, 1/(DSQRT2*DSQRT3)/)
    Trans65(2,1:2) = (/-1/DSQRT2, 1/(DSQRT2*DSQRT3)/)
    Trans65(3,2) = -DSQRT2/DSQRT3
    Trans65(4,3) = 1/DSQRT2
    Trans65(5,4) = 1/DSQRT2
    Trans65(6,5) = 1/DSQRT2
    
    I2(1:3)=1.d0
	I2(4:6)=2.d0
	! Identity matrix
	eye3=0.d0
	eye3(1,1)=1.d0
	eye3(2,2)=1.d0
	eye3(3,3)=1.d0
	eye6 = (/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)

    ENDSUBROUTINE

    ENDMODULE M_CONSTANTS

