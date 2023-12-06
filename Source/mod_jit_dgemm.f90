    INCLUDE 'mkl_blas.f90'
    MODULE M_jit_dgemm
    ! MKL's JIT optimized GEMM kernels for small matrix MATMUL operations
    USE ISO_C_BINDING,  ONLY:  c_ptr, c_int, c_intptr_t, &
                               c_funptr, c_double, c_f_procpointer
    USE MKL_JIT_BLAS
    IMPLICIT NONE
    
    PROCEDURE(dgemm_jit_kernel_t), POINTER :: mkl_jit_dgemm_6_6,mkl_jit_dgemm_6_6T, mkl_jit_dgemm_6T_6_sum, & !used
                                              !mkl_jit_dgemm_5_6_6, & !not used
                                              mkl_jit_dgemm_5_6_5, & !used         
                                              mkl_jit_dgemm_6_6_5, & !used
                                              mkl_jit_dgemm_5_5_1, & !used
                                              mkl_jit_dgemm_6_6_1, & !used
                                              mkl_jit_dgemm_5_1_5_sum, & !used
                                              mkl_jit_dgemm_6_1_6_sum, & !used
                                              mkl_jit_dgemm_7_1_7_sum, & !used
                                              mkl_jit_dgemm_3_3,mkl_jit_dgemm_3_3T,mkl_jit_dgemm_3T_3,& !used
                                              mkl_jit_dgemm_7_7
                                              
    
    TYPE(C_PTR) :: jitter_6_6,jitter_6_6T, jitter_6T_6_sum,  &
                   !jitter_5_6_6, &   
                   jitter_5_6_5, & 
                   jitter_6_6_5, &    
                   jitter_5_5_1, & 
                   jitter_6_6_1, & 
                   jitter_5_1_5_sum, & 
                   jitter_6_1_6_sum, & 
                   jitter_7_1_7_sum, & 
                   jitter_3_3,jitter_3_3T,jitter_3T_3, &
                   jitter_7_7
    CONTAINS
    
    SUBROUTINE create_jit_dgemms
    IMPLICIT NONE
    
    WRITE(*,*) 'Creating JIT optimized GEMM kernels for small matrix MATMUL operations'
    CALL get_custom_jit_dgemm(6,6,6,6,6,6,1.d0,0.d0,'N','N',jitter_6_6,mkl_jit_dgemm_6_6)       !A_66*B_66
    CALL get_custom_jit_dgemm(6,6,6,6,6,6,1.d0,0.d0,'N','T',jitter_6_6T,mkl_jit_dgemm_6_6T)     !A_66*B_66
    CALL get_custom_jit_dgemm(6,6,6,6,6,6,1.d0,1.d0,'T','N',jitter_6T_6_sum,mkl_jit_dgemm_6T_6_sum)       !A_66T*B_66+C_66
    CALL get_custom_jit_dgemm(5,6,5,5,6,5,1.d0,0.d0,'N','N',jitter_5_6_5,mkl_jit_dgemm_5_6_5)   !A_56*B_65
    CALL get_custom_jit_dgemm(6,6,5,6,6,6,1.d0,0.d0,'N','N',jitter_6_6_5,mkl_jit_dgemm_6_6_5)   !A_66*B_65
    CALL get_custom_jit_dgemm(5,5,1,5,5,5,1.d0,0.d0,'N','N',jitter_5_5_1,mkl_jit_dgemm_5_5_1)   !A_55*B_51
    CALL get_custom_jit_dgemm(6,6,1,6,6,6,1.d0,0.d0,'N','N',jitter_6_6_1,mkl_jit_dgemm_6_6_1)   !A_66*B_61
    CALL get_custom_jit_dgemm(5,1,5,5,1,5,1.d0,1.d0,'N','N',jitter_5_1_5_sum,mkl_jit_dgemm_5_1_5_sum)   !A_51*B_15 + C_55
    CALL get_custom_jit_dgemm(6,1,6,6,1,6,1.d0,1.d0,'N','N',jitter_6_1_6_sum,mkl_jit_dgemm_6_1_6_sum)   !A_61*B_16 + C_66
    CALL get_custom_jit_dgemm(7,1,7,7,1,7,1.d0,1.d0,'N','N',jitter_7_1_7_sum,mkl_jit_dgemm_7_1_7_sum)   !A_71*B_17 + C_77 
    CALL get_custom_jit_dgemm(3,3,3,3,3,3,1.d0,0.d0,'N','N',jitter_3_3,mkl_jit_dgemm_3_3)       !A_33*B_33
    CALL get_custom_jit_dgemm(3,3,3,3,3,3,1.d0,0.d0,'N','T',jitter_3_3T,mkl_jit_dgemm_3_3T)     !A_33*B_33T
    CALL get_custom_jit_dgemm(3,3,3,3,3,3,1.d0,0.d0,'T','N',jitter_3T_3,mkl_jit_dgemm_3T_3)     !A_33*B_33T
    CALL get_custom_jit_dgemm(6,6,6,6,6,6,1.d0,0.d0,'N','N',jitter_6_6,mkl_jit_dgemm_6_6)
    CALL get_custom_jit_dgemm(7,7,7,7,7,7,1.d0,0.d0,'N','N',jitter_7_7,mkl_jit_dgemm_7_7)       !A_77*B_77
    
    
    ENDSUBROUTINE create_jit_dgemms
    
    SUBROUTINE del_jit_dgemms
    IMPLICIT NONE
    INTEGER     :: status

    ! delete the DGEMM kernel and the code generator to free memory
    status = mkl_jit_destroy(jitter_6_6)
    status = mkl_jit_destroy(jitter_6_6T)
    status = mkl_jit_destroy(jitter_3_3)
    status = mkl_jit_destroy(jitter_6_6_5)
    status = mkl_jit_destroy(jitter_5_5_1)
    status = mkl_jit_destroy(jitter_6_6_1)
    status = mkl_jit_destroy(jitter_5_1_5_sum)
    status = mkl_jit_destroy(jitter_6_1_6_sum)
    status = mkl_jit_destroy(jitter_7_1_7_sum)
    status = mkl_jit_destroy(jitter_6T_6_sum)
    status = mkl_jit_destroy(jitter_3_3T)
    status = mkl_jit_destroy(jitter_3T_3)
    status = mkl_jit_destroy(jitter_7_7)
    
    ENDSUBROUTINE del_jit_dgemms
    
    SUBROUTINE mkl_jit_matmul(m,k,n,a,b,c)
    IMPLICIT NONE
    INTEGER :: m,n,k
    REAL*8 a(m,k),b(k,n),c(m,n)
    
    ! Custom operation of the generated mkl_jit_matmul kernel:
    !   C_mn := A_mk*B_kn
    IF ((m.eq.6).AND.(k.eq.6).AND.(n.eq.6)) THEN
        CALL mkl_jit_dgemm_6_6(jitter_6_6,a,b,c)
    ELSEIF ((m.eq.3).AND.(k.eq.3).AND.(n.eq.3)) THEN
        CALL mkl_jit_dgemm_3_3(jitter_3_3,a,b,c)
    ELSEIF ((m.eq.5).AND.(k.eq.6).AND.(n.eq.5)) THEN
        CALL mkl_jit_dgemm_5_6_5(jitter_5_6_5,a,b,c)
    ELSEIF ((m.eq.5).AND.(k.eq.6).AND.(n.eq.6)) THEN
        !CALL mkl_jit_dgemm_5_6_6(jitter_5_6_6,a,b,c)
    ELSEIF ((m.eq.6).AND.(k.eq.6).AND.(n.eq.5)) THEN
        CALL mkl_jit_dgemm_6_6_5(jitter_6_6_5,a,b,c)
    ELSEIF ((m.eq.5).AND.(k.eq.5).AND.(n.eq.1)) THEN
        CALL mkl_jit_dgemm_5_5_1(jitter_5_5_1,a,b,c)
    ELSEIF ((m.eq.6).AND.(k.eq.6).AND.(n.eq.1)) THEN
        CALL mkl_jit_dgemm_6_6_1(jitter_6_6_1,a,b,c)
    ELSEIF ((m.eq.5).AND.(k.eq.1).AND.(n.eq.5)) THEN
        !CALL mkl_jit_dgemm_5_1_5(jitter_5_1_5,a,b,c)
    ELSEIF ((m.eq.7).AND.(k.eq.7).AND.(n.eq.7)) THEN
        CALL mkl_jit_dgemm_7_7(jitter_7_7,a,b,c)
    ELSE
        WRITE(*,*) 'JIT GEMM not implementd!'
        READ(*,*)
    ENDIF
    
    
    ENDSUBROUTINE mkl_jit_matmul
    
    SUBROUTINE mkl_jit_matmul_trans(m,k,n,transa,transb,a,b,c)
    IMPLICIT NONE
    INTEGER :: m,n,k
    CHARACTER*1 :: transa,transb
    REAL*8 a(m,k),b(k,n),c(m,n)
    
     IF ((m.eq.3).AND.(k.eq.3).AND.(n.eq.3) &
        .AND.(transa.eq.'N').AND.(transb.eq.'T')) THEN
        CALL mkl_jit_dgemm_3_3T(jitter_3_3T,a,b,c)
    ELSEIF ((m.eq.3).AND.(k.eq.3).AND.(n.eq.3) &
        .AND.(transa.eq.'T').AND.(transb.eq.'N')) THEN
         CALL mkl_jit_dgemm_3T_3(jitter_3T_3,a,b,c)
    ELSEIF ((m.eq.6).AND.(k.eq.6).AND.(n.eq.6) &
        .AND.(transa.eq.'N').AND.(transb.eq.'T')) THEN
        CALL mkl_jit_dgemm_6_6T(jitter_6_6T,a,b,c)
    ELSE
        WRITE(*,*) 'JIT GEMM not implementd!'
        READ(*,*)
    ENDIF
    
    ENDSUBROUTINE
    
    
    SUBROUTINE mkl_jit_matmul_sum(m,k,n,transa,transb,a,b,c)
    IMPLICIT NONE
    INTEGER :: m,n,k
    CHARACTER*1 :: transa,transb
    REAL*8 a(m,k),b(k,n),c(m,n)
    
    IF ((m.eq.6).AND.(k.eq.6).AND.(n.eq.6) &
        .AND.(transa.eq.'T').AND.(transb.eq.'N')) THEN
        CALL mkl_jit_dgemm_6T_6_sum(jitter_6T_6_sum,a,b,c)
    ELSEIF ((m.eq.5).AND.(k.eq.1).AND.(n.eq.5) &
        .AND.(transa.eq.'N').AND.(transb.eq.'N')) THEN
        CALL mkl_jit_dgemm_5_1_5_sum(jitter_5_1_5_sum,a,b,c)
    ELSEIF ((m.eq.6).AND.(k.eq.1).AND.(n.eq.6) &
        .AND.(transa.eq.'N').AND.(transb.eq.'N')) THEN
        CALL mkl_jit_dgemm_6_1_6_sum(jitter_6_1_6_sum,a,b,c)
    ELSEIF ((m.eq.7).AND.(k.eq.1).AND.(n.eq.7) &
        .AND.(transa.eq.'N').AND.(transb.eq.'N')) THEN
        CALL mkl_jit_dgemm_7_1_7_sum(jitter_7_1_7_sum,a,b,c)
    ELSE
        WRITE(*,*) 'JIT GEMM not implementd!'
        READ(*,*)
    ENDIF
    
    
    ENDSUBROUTINE mkl_jit_matmul_sum
    
    SUBROUTINE get_custom_jit_dgemm(m,k,n,lda,ldb,ldc,alpha,beta,transa,transb,jitter,f_kernel) 
    
    IMPLICIT NONE
    INTEGER :: m,n,k
    CHARACTER*1 :: transa,transb
    INTEGER     :: lda,ldb,ldc
    REAL(KIND=c_double)  :: alpha, beta
    INTEGER     :: status
    TYPE(C_PTR) :: jitter ! C pointer
    TYPE(C_FUNPTR) :: c_kernel
    PROCEDURE(dgemm_jit_kernel_t),POINTER :: f_kernel
    
    
    ! General operation of the generated GEMM kernel:
    ! C_mn := alpha*op(A_mk)*op(B_kn) + beta*C_mn
    
    status = mkl_jit_create_dgemm(jitter,transa,transb,m,n,k,alpha,lda,ldb,beta,ldc)
    
    IF (status.EQ.MKL_JIT_ERROR) THEN
        WRITE(*,*) 'jitter was not created'
        READ(*,*)
    ELSEIF (status.EQ.MKL_NO_JIT) THEN
        WRITE(*,*) 'jitter has been created, but a JIT GEMM kernel was not created &
                   & because JIT is not beneficial for the given input parameters. &
                   & Standard GEMM will be used.'
    ELSE
        !MKL_JIT_SUCCESS
        !WRITE(*,*) 'jitter has been created and the GEMM kernel was successfully created'
    ENDIF
    
    ! C procedure pointer to the generated GEMM kernel
    ! this pointer needs to be converted to a Fortran procedure pointer 
    c_kernel=mkl_jit_get_dgemm_ptr(jitter)
    
    ! convert C function pointer to fortran procedure pointer
    CALL C_F_PROCPOINTER(c_kernel,f_kernel)
    
    ! usage:
    !CALL f_kernel(jitter,a,b,c)


    ENDSUBROUTINE get_custom_jit_dgemm
    
    ENDMODULE
