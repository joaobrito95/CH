    program CH_main

    USE M_COMP_HOMOG
    USE MKL_SERVICE
    USE M_jit_dgemm, ONLY: create_jit_dgemms, del_jit_dgemms

    IMPLICIT NONE
    REAL*8  D_voigt(6), f, macroSRP, SIGMA_norm(6), D_voigt_dev(6),Dm,SIGMA(6)
    INTEGER i
    INTEGER isa, irc
    INTEGER ::  count, rate
    REAL*8  :: timeAtStart, timeAtEnd,time
    real*8 d5_ast(5),dummy(5)
    LOGICAL keep
    REAL*8 macroEQS, D6_norm(6), nan, J0(7,7)
    LOGICAL J0_exists
    LOGICAL YesEshelby
    
    
    NaN = DSQRT(-1.d0)

#ifdef _AVX512
    isa = MKL_ENABLE_AVX512_E4
#else
    isa = MKL_ENABLE_AVX2
#endif

    irc=mkl_enable_instructions(isa)

    ! Import control parameters, material data and integration rule
    CALL CH_INIT
    
    f = f_read
    
    CALL getLoadingdata(D_voigt,SIGMA)     
     
	 D_voigt_dev=D_voigt      
	 Dm = SUM(D_voigt(1:3))/3.d0	        
	 D_voigt_dev(1:3)=D_voigt(1:3)-Dm
	 
   
    WRITE(*,'(A,F10.6)') 'f = ',f
    WRITE(*,'(A)') ' '

    WRITE(*,*) 'DIRECT PROBLEM'
    !------------------------------------
    ! Rice and Tracey v.f. 
    YesEshelby = .FALSE.
    CALL SYSTEM_CLOCK(count = count,count_rate = rate)
    timeAtStart = DFLOAT(count)/DFLOAT(rate)
    
    WRITE(*,*) 'Rice and Tracey v.f'
    CALL CH_RiceTracey(D_voigt,f,mat_param,sph_tdesign,macroSRP,SIGMA_norm)
    
    WRITE(*,'(A,F10.6)') 'macroSRP = ',macroSRP
    WRITE(*,'(A,6F10.6)') 'macroSIGMA = ',SIGMA_norm

    CALL SYSTEM_CLOCK(count = count, count_rate = rate)
    timeAtEnd = DFLOAT(count)/DFLOAT(rate)
    time=timeAtEnd-timeAtStart
    WRITE(*,'(A,F6.3,A)') 'Elapsed time: ',time, ' sec.'
    !READ(*,*) 

	 !------------------------------------
    !! Eshelby v.f. 
    CALL SYSTEM_CLOCK(count = count,count_rate = rate)
    timeAtStart = DFLOAT(count)/DFLOAT(rate)
    
    WRITE(*,'(A)') ' '
    WRITE(*,*) 'Eshelby v.f'
    
    d5_ast = nan
    CALL CH_Eshelby(D_voigt,f,mat_param,sph_tdesign,sph_tdesign_crude,&
            d5_ast,macroSRP,SIGMA_norm)
    
    WRITE(*,'(A,F10.6)') 'macroSRP = ',macroSRP
    WRITE(*,'(A,6F10.6)') 'macroSIGMA = ',SIGMA_norm
        
    CALL SYSTEM_CLOCK(count = count, count_rate = rate)
    timeAtEnd = DFLOAT(count)/DFLOAT(rate)
    time=timeAtEnd-timeAtStart
    WRITE(*,'(A,F6.3,A)') 'Elapsed time: ',time, ' sec.'
    WRITE(*,'(A)') ' '
    
    !------------------------------------
    ! Inverse problem: 
    !------------------------------------
    WRITE(*,*) 'INVERSE PROBLEM'
    ! Rice and Tracey v.f. 

    CALL SYSTEM_CLOCK(count = count,count_rate = rate)
    timeAtStart = DFLOAT(count)/DFLOAT(rate)
    
    WRITE(*,*) 'Rice and Tracey v.f.'
    YesEshelby = .FALSE.

    macroEQS = nan
    D6_norm = nan*D6_norm
    J0_exists = .FALSE.
    J0=0.d0
    CALL CH_INVERSE(SIGMA,f,mat_param,sph_tdesign,sph_tdesign_crude,YesEshelby,dummy,&
                         J0_exists,J0,macroEQS,D6_norm)
    WRITE(*,'(A,F10.6)') 'macroEQS = ',macroEQS
    WRITE(*,'(A,6F10.6)') 'D6_norm = ',D6_norm

    CALL SYSTEM_CLOCK(count = count, count_rate = rate)
    timeAtEnd = DFLOAT(count)/DFLOAT(rate)
    time=timeAtEnd-timeAtStart
    WRITE(*,'(A,F6.3,A)') 'Elapsed time: ',time, ' sec.'
    
    !------------------------------------
    ! Eshelby v.f. 
    CALL SYSTEM_CLOCK(count = count,count_rate = rate)
    timeAtStart = DFLOAT(count)/DFLOAT(rate)
    
    WRITE(*,'(A)') ' '
    WRITE(*,*) 'Eshelby v.f'
    
    !macroEQS = nan*macroEQS
    !D6_norm = nan*D6_norm
    ! Note: if macroEQS and D6_norm are unknown, firtly compute these in the case
    d5_ast=nan*d5_ast
    YesEshelby = .TRUE.
    J0_exists = .FALSE.
    J0=0.d0
    CALL CH_INVERSE(SIGMA,f,mat_param,sph_tdesign,sph_tdesign_crude,YesEshelby,d5_ast,&
                        J0_exists,J0,macroEQS,D6_norm)
    WRITE(*,'(A,F10.6)') 'macroEQS = ',macroEQS
    WRITE(*,'(A,6F10.6)') 'D6_norm = ',D6_norm
    
    CALL SYSTEM_CLOCK(count = count, count_rate = rate)
    timeAtEnd = DFLOAT(count)/DFLOAT(rate)
    time=timeAtEnd-timeAtStart
    WRITE(*,'(A,F6.3,A)') 'Elapsed time: ',time, ' sec.'
    
    CALL del_jit_dgemms
    
    READ(*,*)
    
    
    end program CH_main

