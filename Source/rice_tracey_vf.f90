    SUBROUTINE CH_RiceTracey(D_voigt,f,mat_param,sph_tdesign,macroSRP,SIGMA_norm)
    
    USE M_COMP_HOMOG, ONLY : struct_mat_param,struct_sphtdesign
    
    IMPLICIT NONE
    REAL*8, intent(IN) :: D_voigt(6), f 
    TYPE (struct_mat_param),  INTENT(IN) :: mat_param
    TYPE (struct_sphtdesign), INTENT(IN) :: sph_tdesign
    REAL*8, INTENT(OUT) :: macroSRP,SIGMA_norm(6)
    LOGICAL YesEshelby
    REAL*8 dummy(5)
    
    YesEshelby = .FALSE.
    CALL comphomog_commun(D_voigt,YesEshelby,dummy,f,mat_param,sph_tdesign,macroSRP,SIGMA_norm)
    
    ENDSUBROUTINE CH_RiceTracey
    
    ! =============================================================
    SUBROUTINE local_d_RT(D6,r,theta,phi,d6_RT)
    USE M_COMP_HOMOG,  ONLY: b_radius
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: D6(6), r, theta, phi
    REAL*8, INTENT(OUT) :: d6_RT(6)
    
    REAL*8 Dm, D_dev(6)
    REAL*8 dv_11,dv_22,dv_33,dv_23,dv_13,dv_12,aux_factor
    REAL*8 sintheta, costheta, sinphi, cosphi
    
    !Rice and Tracey strain-rate field in the cartesian frame
    Dm = (1.d0/3.d0)*SUM(D6(1:3))
    D_dev = D6-Dm*(/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)
    
	sintheta = DSIN(theta)
    costheta = DCOS(theta)
    sinphi = DSIN(phi)
    cosphi = DCOS(phi)
    
    aux_factor = (Dm*(b_radius/r)**(3))
    dv_11 = aux_factor*(1.d0-3.d0*cosphi**2*sintheta**2);
    dv_22 = aux_factor*(1.d0-3.d0*sinphi**2*sintheta**2);
    dv_33 = -(dv_11+dv_22)
    dv_23 = aux_factor*(-3.d0*costheta*sinphi*sintheta);
    dv_13 = aux_factor*(-3.d0*sintheta*costheta*cosphi);
    dv_12 = aux_factor*(-3.d0*sintheta**2*sinphi*cosphi);

    d6_RT = D_dev+(/dv_11,dv_22,dv_33,dv_23,dv_13,dv_12/)
    
    ENDSUBROUTINE local_d_RT

