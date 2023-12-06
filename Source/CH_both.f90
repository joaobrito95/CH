	SUBROUTINE CH_compute(D_voigt,f,mat_param,sph_tdesign,YesEshelby,d5_ast,macroSRP,SIGMA_norm)
	USE M_COMP_HOMOG, ONLY :  struct_mat_param,struct_sphtdesign,&
		sph_tdesign_crude
	IMPLICIT NONE
	REAL*8, intent(IN) :: D_voigt(6)
	REAL*8, intent(IN) :: f
	LOGICAL, INTENT(IN) :: YesEshelby
	TYPE (struct_mat_param),  INTENT(IN) :: mat_param
	TYPE (struct_sphtdesign), INTENT(IN) :: sph_tdesign
	REAL*8, INTENT(INOUT) :: d5_ast(5)
	REAL*8, INTENT(OUT) :: macroSRP,SIGMA_norm(6)
	REAL*8  macroSRP_RT,SIGMA_norm_RT(6)
	REAL*8 norm_D_voigt, D_voigt_norm(6), dast_norm(5)
	
	norm_D_voigt = DSQRT(SUM(D_voigt*D_voigt))
	IF (norm_D_voigt.NE.0.d0) THEN
		!D_voigt = D_voigt/norm_D_voigt
		!d5_ast = d5_ast/norm_D_voigt
		D_voigt_norm= D_voigt/norm_D_voigt
		dast_norm=d5_ast/norm_D_voigt
	ELSE
		D_voigt_norm=0.d0
		dast_norm=0.d0
	ENDIF
	
	
	IF (.not.YesEshelby) THEN
		CALL CH_RiceTracey(D_voigt_norm,f,mat_param,sph_tdesign,macroSRP,SIGMA_norm)
		dast_norm=0.d0
	ELSE
		CALL CH_RiceTracey(D_voigt_norm,f,mat_param,sph_tdesign,macroSRP_RT,SIGMA_norm_RT)
		CALL CH_Eshelby(D_voigt_norm,f,mat_param,sph_tdesign,sph_tdesign_crude,&
			dast_norm,macroSRP,SIGMA_norm)
		IF (macroSRP.EQ.macroSRP) THEN
			  IF (macroSRP.GT.macroSRP_RT) THEN
				  	 SIGMA_norm=SIGMA_norm_RT
					 macroSRP=macroSRP_RT
					 dast_norm=0.d0
			  ENDIF
		ENDIF
		
	ENDIF
	
	!D_voigt = D_voigt*norm_D_voigt
	!d5_ast = d5_ast*norm_D_voigt
	macroSRP = macroSRP*norm_D_voigt
	
	ENDSUBROUTINE

