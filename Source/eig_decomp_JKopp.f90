    ! ----------------------------------------------------------------------------
    SUBROUTINE DSYEVC3(A, W)
    ! ----------------------------------------------------------------------------
    ! Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
    ! analytical algorithm.
    ! Only the diagonal and upper triangular parts of A are accessed. The access
    ! is read-only.
    ! ----------------------------------------------------------------------------
    ! Parameters:
    !   A: The symmetric input matrix
    !   W: Storage buffer for eigenvalues
    ! ----------------------------------------------------------------------------
    DOUBLE PRECISION A(3,3)
    DOUBLE PRECISION W(3)

    DOUBLE PRECISION DSQRT3
    PARAMETER        ( DSQRT3 = 1.73205080756887729352744634151D0 )

    DOUBLE PRECISION M, C1, C0
    DOUBLE PRECISION DE, DD, EE, FF
    DOUBLE PRECISION P, DSQRTP, Q, C, S, PHI


    !     Determine coefficients of characteristic poynomial. We write
    !           | A   D   F  |
    !      A =  | D*  B   E  |
    !           | F*  E*  C  |

    DE    = A(1,2) * A(2,3)
    DD    = A(1,2)**2
    EE    = A(2,3)**2
    FF    = A(1,3)**2
    M     = A(1,1) + A(2,2) + A(3,3)
    C1    = ( A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3) ) &
        - (DD + EE + FF)
    C0    = A(3,3)*DD + A(1,1)*EE + A(2,2)*FF - A(1,1)*A(2,2)*A(3,3) &
        - 2.0D0 * A(1,3)*DE

    P     = M**2 - 3.0D0 * C1
    Q     = M*(P - (3.0D0/2.0D0)*C1) - (27.0D0/2.0D0)*C0
    DSQRTP = DSQRT(DABS(P))

    PHI   = 27.0D0 * ( 0.25D0 * C1**2 * (P - C1) &
        + C0 * (Q + (27.0D0/4.0D0)*C0) )
    PHI   = (1.0D0/3.0D0) * DATAN2(DSQRT(DABS(PHI)), Q)

    C     = DSQRTP * DCOS(PHI)
    S     = (1.0D0/DSQRT3) * DSQRTP * DSIN(PHI)

    W(2) = (1.0D0/3.0D0) * (M - C)
    W(3) = W(2) + S
    W(1) = W(2) + C
    W(2) = W(2) - S

    END SUBROUTINE

    ! ----------------------------------------------------------------------------
    SUBROUTINE DSYEVH3(A, Q, W)
    ! ----------------------------------------------------------------------------
    ! Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
    ! matrix A using Cardano's method for the eigenvalues and an analytical
    ! method based on vector cross products for the eigenvectors. However,
    ! if conditions are such that a large error in the results is to be
    ! expected, the routine falls back to using the slower, but more
    ! accurate QL algorithm. Only the diagonal and upper triangular parts of A need
    ! to contain meaningful values. Access to A is read-only.
    ! ----------------------------------------------------------------------------
    ! Parameters:
    !   A: The symmetric input matrix
    !   Q: Storage buffer for eigenvectors
    !   W: Storage buffer for eigenvalues
    ! ----------------------------------------------------------------------------
    ! Dependencies:
    !   DSYEVC3(), DSYTRD3(), DSYEVQ3()
    ! ----------------------------------------------------------------------------

    DOUBLE PRECISION A(3,3)
    DOUBLE PRECISION Q(3,3)
    DOUBLE PRECISION W(3)

    DOUBLE PRECISION EPS
    PARAMETER        ( EPS = 2.2204460492503131D-16 )

    DOUBLE PRECISION NORM
    DOUBLE PRECISION ERROR
    DOUBLE PRECISION T, U
    INTEGER          J

    EXTERNAL         DSYEVC3, DSYEVQ3

    !     Calculate eigenvalues
    CALL DSYEVC3(A, W)

    !     --- The rest of this subroutine can be omitted if only the eigenvalues are desired ---

    !     Prepare calculation of eigenvectors
    !      N1      = A(1, 1)**2 + A(1, 2)**2 + A(1, 3)**2
    !      N2      = A(1, 2)**2 + A(2, 2)**2 + A(2, 3)**2
    T       = MAX(DABS(W(1)), DABS(W(2)), DABS(W(3)))
    U       = MAX(T, T**2)
    ERROR   = 256.0D0 * EPS * U**2

	 !      ERROR   = 256.0D0 * EPS * (N1 + U) * (N2 + U)
    Q(1, 2) = A(1, 2) * A(2, 3) - A(1, 3) * A(2, 2)
    Q(2, 2) = A(1, 3) * A(1, 2) - A(2, 3) * A(1, 1)
    Q(3, 2) = A(1, 2)**2

    !     Calculate first eigenvector by the formula
    !       v[0] = (A - lambda[0]).e1 x (A - lambda[0]).e2
    Q(1, 1) = Q(1, 2) + A(1, 3) * W(1)
    Q(2, 1) = Q(2, 2) + A(2, 3) * W(1)
    Q(3, 1) = (A(1,1) - W(1)) * (A(2,2) - W(1)) - Q(3,2)
    NORM    = Q(1, 1)**2 + Q(2, 1)**2 + Q(3, 1)**2
	 
    !     If vectors are nearly linearly dependent, or if there might have
    !     been large cancellations in the calculation of A(I,I) - W(1), fall
    !     back to QL algorithm
    !     Note that this simultaneously ensures that multiple eigenvalues do
    !     not cause problems: If W(1) = W(2), then A - W(1) * I has rank 1,
    !     i.e. all columns of A - W(1) * I are linearly dependent.
    IF (NORM .LE. ERROR) THEN
        CALL DSYEVQ3(A, Q, W)
        RETURN
        !*     This is the standard branch
    ELSE
        NORM = DSQRT(1.0D0 / NORM)
        DO 20, J = 1, 3
            Q(J, 1) = Q(J, 1) * NORM
20      CONTINUE
	 END IF
	 
	 !*     Calculate second eigenvector by the formula
    !*       v[1] = (A - lambda[1]).e1 x (A - lambda[1]).e2
    Q(1, 2) = Q(1, 2) + A(1, 3) * W(2)
    Q(2, 2) = Q(2, 2) + A(2, 3) * W(2)
    Q(3, 2) = (A(1,1) - W(2)) * (A(2,2) - W(2)) - Q(3, 2)
    NORM    = Q(1, 2)**2 + Q(2, 2)**2 + Q(3, 2)**2
    IF (NORM .LE. ERROR) THEN
        CALL DSYEVQ3(A, Q, W)
        RETURN
    ELSE
        NORM = DSQRT(1.0D0 / NORM)
        DO 40, J = 1, 3
            Q(J, 2) = Q(J, 2) * NORM
40      CONTINUE
	 END IF

	 !*     Calculate third eigenvector according to
    !*       v[2] = v[0] x v[1]
    Q(1, 3) = Q(2, 1) * Q(3, 2) - Q(3, 1) * Q(2, 2)
    Q(2, 3) = Q(3, 1) * Q(1, 2) - Q(1, 1) * Q(3, 2)
    Q(3, 3) = Q(1, 1) * Q(2, 2) - Q(2, 1) * Q(1, 2)
     
	 IF (DABS(SUM(Q(:, 3)*Q(:, 3))-1.d0).GT.(10.0*EPS)) THEN
		 CALL DSYEVQ3(A, Q, W)
		 RETURN
	 ENDIF
	 
	 
    END SUBROUTINE


    ! ----------------------------------------------------------------------------
    SUBROUTINE DSYEVQ3(A, Q, W)
    ! ----------------------------------------------------------------------------
    ! Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
    ! matrix A using the QL algorithm with implicit shifts, preceded by a
    ! Householder reduction to real tridiagonal form.
    ! The function accesses only the diagonal and upper triangular parts of
    ! A. The access is read-only.
    ! ----------------------------------------------------------------------------
    ! Parameters:
    !   A: The symmetric input matrix
    !   Q: Storage buffer for eigenvectors
    !   W: Storage buffer for eigenvalues
    ! ----------------------------------------------------------------------------
    ! Dependencies:
    !   DSYTRD3()
    ! ----------------------------------------------------------------------------
    DOUBLE PRECISION A(3,3)
    DOUBLE PRECISION Q(3,3)
    DOUBLE PRECISION W(3)

    INTEGER          N
    PARAMETER        ( N = 3 )

    DOUBLE PRECISION E(3)
    DOUBLE PRECISION G, R, P, F, B, S, C, T
    INTEGER          NITER
    INTEGER          L, M, I, J, K

    EXTERNAL         DSYTRD3

    !     Transform A to real tridiagonal form by the Householder method
    CALL DSYTRD3(A, Q, W, E)

    !     Calculate eigensystem of the remaining real symmetric tridiagonal
    !     matrix with the QL method

    !     Loop over all off-diagonal elements
    DO 10 L = 1, N-1
        NITER = 0

        !       Iteration loop
        DO 11 I = 1, 50
            !         Check for convergence and exit iteration loop if off-diagonal
            !         element E(L) is zero
            DO 20 M = L, N-1
                G = DABS(W(M)) + DABS(W(M+1))
                IF (DABS(E(M)) + G .EQ. G) THEN
                    GO TO 30
                END IF
20          CONTINUE
30          IF (M .EQ. L) THEN
                GO TO 10
            END IF

            NITER = NITER + 1
            IF (NITER >= 30) THEN
                PRINT *, 'DSYEVQ3: No convergence.'
                RETURN
            END IF

            !         Calculate G = D(M) - K
            G = (W(L+1) - W(L)) / (2.0D0 * E(L))
            R = DSQRT(1.0D0 + G**2)
            IF (G .GE. 0.0D0) THEN
                G = W(M) - W(L) + E(L)/(G + R)
            ELSE
                G = W(M) - W(L) + E(L)/(G - R)
            END IF

            S = 1.0D0
            C = 1.0D0
            P = 0.0D0
            DO 40 J = M - 1, L, -1
                F = S * E(J)
                B = C * E(J)
                IF (DABS(F) .GT. DABS(G)) THEN
                    C      = G / F
                    R      = DSQRT(1.0D0 + C**2)
                    E(J+1) = F * R
                    S      = 1.0D0 / R
                    C      = C * S
                ELSE
                    S      = F / G
                    R      = DSQRT(1.0D0 + S**2)
                    E(J+1) = G * R
                    C      = 1.0D0 / R
                    S      = S * C
                END IF

                G      = W(J+1) - P
                R      = (W(J) - G) * S + 2.0D0 * C * B
                P      = S * R
                W(J+1) = G + P
                G      = C * R - B

                !           Form eigenvectors
                !           --- This loop can be omitted if only the eigenvalues are desired ---
                DO 50 K = 1, N
                    T         = Q(K, J+1)
                    Q(K, J+1) = S * Q(K, J) + C * T
                    Q(K, J)   = C * Q(K, J) - S * T
50              CONTINUE
40          CONTINUE
            W(L) = W(L) - P
            E(L) = G
            E(M) = 0.0D0
11      CONTINUE
10  CONTINUE

    END SUBROUTINE


    ! ----------------------------------------------------------------------------
    SUBROUTINE DSYTRD3(A, Q, D, E)
    ! ----------------------------------------------------------------------------
    ! Reduces a symmetric 3x3 matrix to real tridiagonal form by applying
    ! (unitary) Householder transformations:
    !            [ D[1]  E[1]       ]
    !    A = Q . [ E[1]  D[2]  E[2] ] . Q^T
    !            [       E[2]  D[3] ]
    ! The function accesses only the diagonal and upper triangular parts of
    ! A. The access is read-only.
    ! ---------------------------------------------------------------------------
    DOUBLE PRECISION A(3,3)
    DOUBLE PRECISION Q(3,3)
    DOUBLE PRECISION D(3)
    DOUBLE PRECISION E(2)

    INTEGER          N
    PARAMETER        ( N = 3 )

    DOUBLE PRECISION U(N), P(N)
    DOUBLE PRECISION OMEGA, F
    DOUBLE PRECISION K, H, G
    INTEGER          I, J

    !     Initialize Q to the identitity matrix
    !     --- This loop can be omitted if only the eigenvalues are desired ---
    DO 10 I = 1, N
        Q(I,I) = 1.0D0
        DO 11, J = 1, I-1
            Q(I, J) = 0.0D0
            Q(J, I) = 0.0D0
11      CONTINUE
10  CONTINUE

    !     Bring first row and column to the desired form
    H = A(1,2)**2 + A(1,3)**2
    IF (A(1,2) .GT. 0.0D0) THEN
        G = -DSQRT(H)
    ELSE
        G = DSQRT(H)
    END IF
    E(1)  = G
    F     = G * A(1,2)
    U(2)  = A(1,2) - G
    U(3)  = A(1,3)

    OMEGA = H - F
    IF (OMEGA > 0.0D0) THEN
        OMEGA = 1.0D0 / OMEGA
        K     = 0.0D0
        DO 20 I = 2, N
            F    = A(2,I)*U(2) + A(I,3)*U(3)
            P(I) = OMEGA * F
            K    = K + U(I) * F
20      CONTINUE
        K = 0.5D0 * K * OMEGA**2

        DO 30 I = 2, N
            P(I) = P(I) - K * U(I)
30      CONTINUE

        D(1) = A(1,1)
        D(2) = A(2,2) - 2.0D0 * P(2) * U(2)
        D(3) = A(3,3) - 2.0D0 * P(3) * U(3)

        !       Store inverse Householder transformation in Q
        !       --- This loop can be omitted if only the eigenvalues are desired ---
        DO 40, J = 2, N
            F = OMEGA * U(J)
            DO 41 I = 2, N
                Q(I,J) = Q(I,J) - F * U(I)
41          CONTINUE
40      CONTINUE

        !       Calculated updated A(2, 3) and store it in E(2)
        E(2) = A(2, 3) - P(2) * U(3) - U(2) * P(3)
    ELSE
        DO 50 I = 1, N
            D(I) = A(I, I)
50      CONTINUE
        E(2) = A(2, 3)
    END IF

    END SUBROUTINE


