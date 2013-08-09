MODULE ROBUST_UTILITIES
  USE DATATYPES
  IMPLICIT NONE
  !This module contains some procedures used to compute the truncated
  !score used in the REMM estimation algorithm.
  !The non truncated score is referenced as "score" while the truncated
  !and scaled score function, i.e. A*score*Weights_Huber, is referenced
  !as "psi".


  
CONTAINS

  SUBROUTINE WEIGHTS_HUBER(psi,m,n,weights)
    USE ROB_COMMON, ONLY: bound
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: m, n
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:m,1:n) :: psi
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:m,1) :: weights

    INTEGER :: i,j
    REAL (KIND=realkind), DIMENSION(1:m) :: L

    !Compute the norm of the m rows of the matrix h
    L = 0._realkind
    DO j=1,n
       DO i=1,m
          L(i)=L(i)+psi(i,j)**2
       END DO
    END DO

    L=SQRT(L)

    DO i=1,m
       IF (L(i) .GT. 0._realkind) THEN
          weights(i,1)=MIN(1._realkind,bound/L(i))
       ELSE
          weights(i,1)=1._realkind
       END IF
    END DO

    RETURN
  END SUBROUTINE WEIGHTS_HUBER


  SUBROUTINE UPDATE_SCORE_N1(beta_restricted)
    USE S1_COMMON,   ONLY: AM_nb_free_par, T1, N1, grad_i, modified_series
    USE ROB_COMMON,  ONLY: score_N1
    USE LIKELIHOOD_UTILITIES, ONLY: SNP_ORT_FUNC, LIKELIHOOD_UPDATE

    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_free_par) :: beta_restricted

    INTEGER :: i

    !compute the ortogonality function of the modified_series at beta
    !and save the grad_i in psi_matrix_tmp for later use
    i = 1
    CALL LIKELIHOOD_UPDATE(beta_restricted,T1,modified_series,1,i)
    CALL SNP_ORT_FUNC(T1)

    score_N1(1:N1,:) = grad_i(1:N1,:)

    RETURN
  END SUBROUTINE UPDATE_SCORE_N1



  SUBROUTINE UPDATE_PSI_N1
    USE S1_COMMON,   ONLY: AM_nb_not_fixed_par, N1
    USE S2_COMMON,   ONLY: SM_nb_free_par, M_rho
    USE ROB_COMMON,  ONLY: rob_A, rob_A_tilde, S, rob_weights_N1, &
         psi_N1, score_N1, psi_N1_tilde, weights_metric
    USE MYFUNCTIONS, ONLY: MM
    IMPLICIT NONE

    INTEGER :: i,j

    IF (weights_metric .EQ. 0) THEN

       !scale the gradient using rob_A
       psi_N1 = MM(score_N1,rob_A,N1,AM_nb_not_fixed_par,AM_nb_not_fixed_par)

       !compute the weights
       CALL WEIGHTS_HUBER(psi_N1, N1, AM_nb_not_fixed_par, rob_weights_N1)

       !now compute the criterion
       DO j=1,AM_nb_not_fixed_par
          DO i=1,N1
             psi_N1(i,j) = psi_N1(i,j)*rob_weights_N1(i,1)
          END DO
       END DO

    ELSE

       psi_N1 = MM(score_N1,S,N1,AM_nb_not_fixed_par,AM_nb_not_fixed_par)
       psi_N1_tilde = MM(psi_N1,M_rho,N1,AM_nb_not_fixed_par,SM_nb_free_par)
       psi_N1_tilde = MM(psi_N1_tilde,rob_A_tilde,N1,SM_nb_free_par,SM_nb_free_par)
       CALL WEIGHTS_HUBER(psi_N1_tilde, N1, SM_nb_free_par, rob_weights_N1)

       !now compute the criterion
       DO j=1,AM_nb_not_fixed_par
          DO i=1,N1
             psi_N1(i,j) = score_N1(i,j)*rob_weights_N1(i,1)
          END DO
       END DO

    END IF

    RETURN
  END SUBROUTINE UPDATE_PSI_N1


  SUBROUTINE UPDATE_PSI_TMP_N2(beta_restricted)
    USE S1_COMMON,   ONLY: AM_nb_free_par, grad_i
    USE S2_COMMON,   ONLY: T2, N2, simulated_series_step_2
    USE ROB_COMMON,  ONLY: score_N2
    USE LIKELIHOOD_UTILITIES, ONLY: SNP_ORT_FUNC, LIKELIHOOD_UPDATE

    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_free_par) :: beta_restricted

    INTEGER :: i

    !compute the ortogonality function of the modified_series at the beta
    !and save the grad_i in psi_matrix_tmp for later use
    i=1
    CALL LIKELIHOOD_UPDATE(beta_restricted,T2,simulated_series_step_2,2,i)
    CALL SNP_ORT_FUNC(T2)

    score_N2(1:N2,:) = grad_i(1:N2,:)

    RETURN
  END SUBROUTINE UPDATE_PSI_TMP_N2


  SUBROUTINE UPDATE_PSI_N2()
    USE DATATYPES
    USE S1_COMMON,   ONLY: AM_nb_not_fixed_par
    USE S2_COMMON,   ONLY: T2, N2, SM_nb_par, rho_SM_end, rho_rest_SM_old, SM_nb_free_par, &
         nb_discard_step_2, nb_iid_series, simulated_series_step_2, random_iid_step_2, M_rho
    USE ROB_COMMON,  ONLY: rob_A, rob_A_tilde, S, rob_weights_N2, &
         psi_N2, score_N2, psi_N2_tilde, weights_metric
    USE MYFUNCTIONS,          ONLY: MM

    IMPLICIT NONE

    INTEGER :: i,j

    IF (weights_metric .EQ. 0) THEN
       !scale the gradient using rob_A
       psi_N2 = MM(score_N2,rob_A,N2,AM_nb_not_fixed_par,AM_nb_not_fixed_par)

       !compute the weights
       CALL WEIGHTS_HUBER(psi_N2, N2, AM_nb_not_fixed_par, rob_weights_N2)

       !now compute the criterion
       DO j=1,AM_nb_not_fixed_par
          DO i=1,N2
             psi_N2(i,j) = psi_N2(i,j)*rob_weights_N2(i,1)
          END DO
       END DO

    ELSE

       psi_N2 = MM(score_N2,S,N2,AM_nb_not_fixed_par,AM_nb_not_fixed_par)
       psi_N2_tilde = MM(psi_N2,M_rho,N2,AM_nb_not_fixed_par,SM_nb_free_par)
       psi_N2_tilde = MM(psi_N2_tilde,rob_A_tilde,N2,SM_nb_free_par,SM_nb_free_par)
       CALL WEIGHTS_HUBER(psi_N2_tilde, N2, SM_nb_free_par, rob_weights_N2)

       DO j=1,AM_nb_not_fixed_par
          DO i=1,N2
             psi_N2(i,j) = score_N2(i,j)*rob_weights_N2(i,1)
          END DO
       END DO

    END IF

    RETURN
  END SUBROUTINE UPDATE_PSI_N2


  SUBROUTINE UPDATE_A(beta_restricted,do_update_score_N1)
    USE S1_COMMON,   ONLY: AM_nb_free_par, AM_nb_not_fixed_par, sel_score_cov_mat_s1, N1, &
         lags_score_cov_mat_s1, write_output_3
    USE S2_COMMON,   ONLY: SM_nb_free_par, M_rho 
    USE MYFUNCTIONS, ONLY: INV, CHOL, MM, MAT_CROSS_LAG
    USE ROB_COMMON,  ONLY: rob_weights_N1, psi_N1, score_N1, weights_metric, &
         nb_A_updates, rob_A, rob_A_tilde, B_0, S 

    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_free_par) :: beta_restricted
    INTEGER, INTENT(IN) :: do_update_score_N1

    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par,AM_nb_not_fixed_par) :: A
    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par,AM_nb_not_fixed_par) :: SM, gamma
    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par,SM_nb_free_par) :: A_tmp
    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par) :: mean
    INTEGER :: i, j, k


100 FORMAT(I2,TR1,25(ES10.3,TR1))

    IF (weights_metric .EQ. 0) THEN

       IF (do_update_score_N1 .EQ. 1) THEN
          CALL UPDATE_SCORE_N1(beta_restricted)
       END IF

       DO k=1, nb_A_updates
          !scale the gradient using rob_A
          psi_N1 = MM(score_N1,rob_A,N1,AM_nb_not_fixed_par,AM_nb_not_fixed_par)

          !compute the weights
          CALL WEIGHTS_HUBER(psi_N1, N1, AM_nb_not_fixed_par, rob_weights_N1)

          !now compute the criterion
          DO j=1,AM_nb_not_fixed_par
             DO i=1,N1
                psi_N1(i,j) = score_N1(i,j)*rob_weights_N1(i,1)
             END DO
          END DO
          !introdotto il 20.9.2002
        !  mean = 0._realkind
        !  DO j=1,AM_nb_not_fixed_par
        !     DO i=1,N1
        !        mean(j) = mean(j) + psi_N1(i,j)
        !     END DO
        !  END DO
        !  mean = mean / N1

        !  DO j=1,AM_nb_not_fixed_par
        !     DO i=1,N1
        !        psi_N1(i,j) = psi_N1(i,j) - mean(j)
        !     END DO
        !  END DO
          !introdotto il 20.9.2002

          SM=MM(TRANSPOSE(psi_N1),psi_N1,AM_nb_not_fixed_par,N1,AM_nb_not_fixed_par)

          IF (sel_score_cov_mat_s1 .EQ. 2) THEN
             i=1
             DO
                IF (i .GT. lags_score_cov_mat_s1) EXIT
                gamma = MAT_CROSS_LAG(psi_N1,N1,AM_nb_not_fixed_par,i)
                SM = SM+(1.0_realkind-REAL(i)/REAL(lags_score_cov_mat_s1+1))*(gamma+TRANSPOSE(gamma))
                i=i+1
             END DO
          END IF

          A=SM/REAL(N1)

          rob_A = CHOL(INV(A,AM_nb_not_fixed_par),AM_nb_not_fixed_par)

          IF (write_output_3) THEN
             WRITE(UNIT=37,fmt='(A,I4)') "Update: ", k
             DO j=1, AM_nb_not_fixed_par
                WRITE(UNIT=37,fmt=100) j, rob_A(j,:)
             END DO
             WRITE(UNIT=37,fmt=*) 
          END IF

       END DO

    ELSE

       DO k=1, nb_A_updates
          A=MM(B_0,S,AM_nb_not_fixed_par,AM_nb_not_fixed_par,AM_nb_not_fixed_par)
          A=MM(S,A,AM_nb_not_fixed_par,AM_nb_not_fixed_par,AM_nb_not_fixed_par)
          A_tmp=MM(A,M_rho,AM_nb_not_fixed_par,AM_nb_not_fixed_par,SM_nb_free_par)
          rob_A_tilde=MM(TRANSPOSE(M_rho),A_tmp,SM_nb_free_par,AM_nb_not_fixed_par,SM_nb_free_par)
          rob_A_tilde = CHOL(INV(rob_A_tilde,SM_nb_free_par),SM_nb_free_par)

          IF (write_output_3) THEN
             WRITE(UNIT=37,fmt='(A,I4)') "Update: ", k
             DO j=1, SM_nb_free_par
                WRITE(UNIT=37,fmt=100) j, rob_A_tilde(j,:)
             END DO
             WRITE(UNIT=37,fmt=*)
          END IF

       END DO

    END IF

    RETURN
  END SUBROUTINE UPDATE_A


  SUBROUTINE UPDATE_B(beta_restricted, B, do_update_score_N1)
    USE DATATYPES
    USE S1_COMMON,   ONLY: N1, AM_nb_not_fixed_par, AM_nb_free_par, &
         sel_score_cov_mat_s1, lags_score_cov_mat_s1, write_output_3
    USE S2_COMMON,   ONLY: SM_nb_free_par, M_rho
    USE ROB_COMMON,  ONLY: weights_metric, rob_A, rob_A_tilde, rob_weights_N1, &
         psi_N1, score_N1, S, psi_N1_tilde, nb_B_updates
    USE MYFUNCTIONS, ONLY: MM, MAT_CROSS_LAG, INV
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_free_par) :: beta_restricted
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:AM_nb_not_fixed_par,1:AM_nb_not_fixed_par) :: B
    INTEGER, INTENT(IN) :: do_update_score_N1

    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par,1:AM_nb_not_fixed_par) :: gamm
    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par) :: mean

    INTEGER :: i , j, k


100 FORMAT(I2,TR1,25(ES10.3,TR1))

    IF (do_update_score_N1 .EQ. 1) THEN
       CALL UPDATE_SCORE_N1(beta_restricted)
    END IF

    DO k=1, nb_B_updates

       IF (weights_metric .EQ. 0) THEN
          !scale the gradient using rob_A
          psi_N1 = MM(score_N1,rob_A,N1,AM_nb_not_fixed_par,AM_nb_not_fixed_par)

          !compute the weights
          CALL WEIGHTS_HUBER(psi_N1, N1, AM_nb_not_fixed_par, rob_weights_N1)

          !now compute the criterion
          DO j=1,AM_nb_not_fixed_par
             DO i=1,N1
                psi_N1(i,j) = psi_N1(i,j)*rob_weights_N1(i,1)
             END DO
          END DO
          !introdotto il 20.9.2002
        !  mean = 0._realkind
        !  DO j=1,AM_nb_not_fixed_par
        !     DO i=1,N1
        !        mean(j) = mean(j) + psi_N1(i,j)
        !     END DO
        !  END DO
        !  mean = mean / N1

        !  DO j=1,AM_nb_not_fixed_par
        !     DO i=1,N1
        !        psi_N1(i,j) = psi_N1(i,j) - mean(j)
        !     END DO
        !  END DO
          !introdotto il 20.9.2002

       ELSE

          psi_N1 = MM(score_N1,S,N1,AM_nb_not_fixed_par,AM_nb_not_fixed_par)
          psi_N1_tilde = MM(psi_N1,M_rho,N1,AM_nb_not_fixed_par,SM_nb_free_par)
          psi_N1_tilde = MM(psi_N1_tilde,rob_A_tilde,N1,SM_nb_free_par,SM_nb_free_par)
          CALL WEIGHTS_HUBER(psi_N1_tilde, N1, SM_nb_free_par, rob_weights_N1)

          !now compute the criterion
          DO j=1,AM_nb_not_fixed_par
             DO i=1,N1
                psi_N1(i,j) = score_N1(i,j)*rob_weights_N1(i,1)
             END DO
          END DO
          mean = 0._realkind
          DO j=1,AM_nb_not_fixed_par
             DO i=1,N1
                mean(j) = mean(j) + psi_N1(i,j)
             END DO
          END DO
          mean = mean / N1

          DO j=1,AM_nb_not_fixed_par
             DO i=1,N1
                psi_N1(i,j) = psi_N1(i,j) - mean(j)
             END DO
          END DO
       END IF

       B=MM(TRANSPOSE(psi_N1),psi_N1,AM_nb_not_fixed_par,N1,AM_nb_not_fixed_par)

       IF (sel_score_cov_mat_s1 .EQ. 2) THEN
          i=1
          DO
             IF (i .GT. lags_score_cov_mat_s1) EXIT
             gamm = MAT_CROSS_LAG(psi_N1,N1,AM_nb_not_fixed_par,i)
             B = B + (1.0_realkind-REAL(i)/REAL(lags_score_cov_mat_s1+1))*(gamm+TRANSPOSE(gamm))
             i=i+1
          END DO
       END IF

       B=B/REAL(N1)

       S = INV(B,AM_nb_not_fixed_par)

       IF (write_output_3) THEN
          WRITE(UNIT=39,fmt='(A,I4)') "Update: ", k
          DO j=1, AM_nb_not_fixed_par
             WRITE(UNIT=39,fmt=100) j, B(j,:)
          END DO
          WRITE(UNIT=39,fmt=*) 
       END IF

       !further updates are not necessary if weights_metric .EQ. 0
       IF (weights_metric .EQ. 0) THEN 
          exit 
       END IF

    END DO


    RETURN
  END SUBROUTINE UPDATE_B


  SUBROUTINE COMPUTE_WEIGHTS_N1(beta_restricted, do_update_psi)
    USE DATATYPES
    USE S1_COMMON,   ONLY: N1, AM_nb_not_fixed_par, AM_nb_free_par
    USE S2_COMMON,   ONLY: SM_nb_free_par, M_rho 
    USE ROB_COMMON,  ONLY: weights_metric, rob_A, rob_A_tilde, rob_weights_N1, &
         psi_N1, score_N1, S, psi_N1_tilde
    USE MYFUNCTIONS, ONLY: MM

    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_free_par) :: beta_restricted
    INTEGER, INTENT(IN) :: do_update_psi

    IF (do_update_psi .EQ. 1) THEN
       CALL UPDATE_SCORE_N1(beta_restricted)
    END IF

    IF (weights_metric .EQ. 0) THEN
       !scale the gradient using rob_A
       psi_N1 = MM(score_N1,rob_A,N1,AM_nb_not_fixed_par,AM_nb_not_fixed_par)

       !compute the weights
       CALL WEIGHTS_HUBER(psi_N1, N1, AM_nb_not_fixed_par, rob_weights_N1)

    ELSE

       psi_N1 = MM(score_N1,S,N1,AM_nb_not_fixed_par,AM_nb_not_fixed_par)
       psi_N1_tilde = MM(psi_N1,M_rho,N1,AM_nb_not_fixed_par,SM_nb_free_par)
       psi_N1_tilde = MM(psi_N1_tilde,rob_A_tilde,N1,SM_nb_free_par,SM_nb_free_par)
       CALL WEIGHTS_HUBER(psi_N1_tilde, N1, SM_nb_free_par, rob_weights_N1)

    END IF

    RETURN
  END SUBROUTINE COMPUTE_WEIGHTS_N1

END MODULE ROBUST_UTILITIES


MODULE M_RHO_MATRIX
  USE DATATYPES
  IMPLICIT NONE

CONTAINS

  SUBROUTINE PSI_FOR_M_RHO_ESTIMATION(beta_restricted, rho_restricted, H, is_rob_estimation)
    USE DATATYPES
    USE S1_COMMON,   ONLY: AM_nb_not_fixed_par, AM_nb_free_par, &
         write_output_3, dev_unit_nb
    USE S2_COMMON,   ONLY: T2, N2, SM_nb_par, SM_nb_free_par, &
         nb_discard_step_2, nb_iid_series, simulated_series_step_2, random_iid_step_2, &
         SM_nb_not_fixed_par
    USE ROB_COMMON,  ONLY: psi_N2, score_N2
    USE Z_MATRIX,             ONLY: is_new2
    USE UTILITIES,            ONLY: COMPLETE_RHO_RESTRICTED 
    USE SIMULATION,           ONLY: STRUCTURAL_MODEL_SIMULATION
    USE MYFUNCTIONS,          ONLY: MM
    USE ROBUST_UTILITIES,     ONLY: UPDATE_PSI_TMP_N2, UPDATE_PSI_N2
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_free_par) :: beta_restricted
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:SM_nb_free_par) :: rho_restricted
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:AM_nb_not_fixed_par) :: H
    INTEGER, INTENT(IN) :: is_rob_estimation

    REAL (KIND=realkind), DIMENSION(1:SM_nb_par) :: rho
    INTEGER :: i,j

    !simulate the process
    CALL COMPLETE_RHO_RESTRICTED(rho_restricted,rho)
    CALL STRUCTURAL_MODEL_SIMULATION(SM_nb_par,rho,T2,simulated_series_step_2, &
         nb_discard_step_2,nb_iid_series,random_iid_step_2)

    is_new2 = 1

    CALL UPDATE_PSI_TMP_N2(beta_restricted)

    IF (is_rob_estimation .EQ. 1) THEN
       CALL UPDATE_PSI_N2()
    ELSE
       psi_N2 = score_N2
    END IF

    !sum the columns of psi_N2
    H=0._realkind
    DO j=1, AM_nb_not_fixed_par
       DO i=1, N2
          H(j) = H(j) + psi_N2(i,j)
       END DO
    END DO

    H=H/REAL(N2)


  END SUBROUTINE PSI_FOR_M_RHO_ESTIMATION


  SUBROUTINE APPROX_M_RHO(beta, rho, M_rho, is_rob_estimation)
    USE DATATYPES
    USE S1_COMMON,   ONLY: AM_nb_par, AM_nb_not_fixed_par, AM_nb_free_par, grad_i, &
         write_output_3, dev_unit_nb
    USE S2_COMMON,   ONLY: T2, SM_nb_par,  SM_nb_free_par, delta_M_rho, method_derivative_M_rho
    USE ROB_COMMON,  ONLY: M_rho_delta_selection
    USE UTILITIES,   ONLY: RESTRICT_BETA, RESTRICT_RHO
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:SM_nb_par) :: rho
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_par) :: beta
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:AM_nb_not_fixed_par,1:SM_nb_free_par) :: M_rho
    INTEGER, INTENT(IN) :: is_rob_estimation

    REAL (KIND=realkind), DIMENSION(1:AM_nb_free_par) :: beta_restricted
    REAL (KIND=realkind), DIMENSION(1:SM_nb_free_par) :: rho_restricted
    REAL (KIND=realkind), DIMENSION(1:SM_nb_free_par) :: rho_restricted_1
    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par) :: H_1
    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par) :: H_2
    REAL (KIND=realkind) :: delta

    INTEGER :: i_SM_nb_par

100 FORMAT(25(ES11.4,TR1))

    IF (write_output_3) THEN
       WRITE(UNIT=dev_unit_nb,FMT='(A)') "Computation of M_rho"
       WRITE(UNIT=dev_unit_nb,FMT='(A,EN10.3)') " Delta_M_rho:", delta_M_rho      
    END IF

    CALL RESTRICT_BETA(beta,beta_restricted)
    CALL RESTRICT_RHO(rho, rho_restricted)

    !loop over all SM_parameters (compute the directional derivatives)
    DO i_SM_nb_par=1, SM_nb_free_par

       IF (M_rho_delta_selection .EQ. 2) THEN
          delta = MAX(delta_M_rho,ABS(rho_restricted(i_SM_nb_par))*0.0001_realkind)
       ELSE
          delta = delta_M_rho
       END IF

       IF (write_output_3) THEN
          WRITE(UNIT=dev_unit_nb,FMT='(A,I4)') "Call to PSI_FOR_M_RHO_ESTIMATION, coefficient nb. ", i_SM_nb_par
       END IF

       rho_restricted_1 = rho_restricted

       IF ((method_derivative_M_rho .EQ. 1) .OR. (method_derivative_M_rho .EQ. 2))  THEN 
          rho_restricted_1(i_SM_nb_par) = rho_restricted_1(i_SM_nb_par) + delta
       END IF

       CALL PSI_FOR_M_RHO_ESTIMATION(beta_restricted, rho_restricted_1,H_1, is_rob_estimation)

       rho_restricted_1 = rho_restricted

       IF ((method_derivative_M_rho .EQ. 1) .OR. (method_derivative_M_rho .EQ. 3))  THEN 
          rho_restricted_1(i_SM_nb_par) = rho_restricted_1(i_SM_nb_par) - delta
       END IF

       CALL PSI_FOR_M_RHO_ESTIMATION(beta_restricted, rho_restricted_1,H_2, is_rob_estimation)

       IF (method_derivative_M_rho .EQ. 1)  THEN 
          M_rho(:,i_SM_nb_par) = (H_1(:) - H_2(:) ) / (2._realkind * delta)
       ELSE
          M_rho(:,i_SM_nb_par) = (H_1(:) - H_2(:) ) / delta  
       END IF

       IF (write_output_3) THEN
          WRITE(UNIT=dev_unit_nb,FMT='(A11,TR1,A)') "Delta","Derivative:"
          WRITE(UNIT=dev_unit_nb,FMT=100) delta, M_rho(:,i_SM_nb_par)
          WRITE(UNIT=dev_unit_nb,FMT=*) 
       END IF

    END DO  !i_SM_nb_par

  END SUBROUTINE APPROX_M_RHO

  SUBROUTINE WRITE_FUNCTION_M(beta, rho, is_rob_estimation)
    USE DATATYPES
    USE S1_COMMON,   ONLY: AM_nb_par, AM_nb_not_fixed_par, AM_nb_free_par, grad_i, &
         write_output_3, dev_unit_nb
    USE ROB_COMMON,  ONLY: M_rho_delta_selection
    USE S2_COMMON,   ONLY: T2, SM_nb_par,  SM_nb_free_par, delta_M_rho, method_derivative_M_rho
    USE UTILITIES,   ONLY: RESTRICT_BETA, RESTRICT_RHO
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:SM_nb_par) :: rho
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_par) :: beta
    INTEGER, INTENT(IN) :: is_rob_estimation

    REAL (KIND=realkind), DIMENSION(1:AM_nb_free_par) :: beta_restricted
    REAL (KIND=realkind), DIMENSION(1:SM_nb_free_par) :: rho_restricted
    REAL (KIND=realkind), DIMENSION(1:SM_nb_free_par) :: rho_restricted_1
    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par) :: H
    REAL (KIND=realkind) :: delta

    INTEGER :: i,i_SM_nb_par

100 FORMAT(25(ES11.4,TR1))

    IF (write_output_3) THEN
       WRITE(UNIT=dev_unit_nb,FMT='(A)') "Computation of M"
       WRITE(UNIT=dev_unit_nb,FMT='(A,TR1,E9.3)') "the:", beta
       WRITE(UNIT=dev_unit_nb,FMT='(A,TR1,E9.3)') "rho:", rho
    END IF

    CALL RESTRICT_BETA(beta,beta_restricted)
    CALL RESTRICT_RHO(rho, rho_restricted)

    !loop over all SM_parameters (compute the directional derivatives)
    DO i_SM_nb_par=1, SM_nb_free_par

       IF (M_rho_delta_selection .EQ. 2) THEN
          delta = MAX(delta_M_rho,ABS(rho_restricted(i_SM_nb_par))*0.0001_realkind)
       ELSE
          delta = delta_M_rho
       END IF

       IF (write_output_3) THEN
          WRITE(UNIT=dev_unit_nb,FMT='(A,I4)') "Call to COMPUTE_M, coefficient nb. ", i_SM_nb_par
       END IF
       DO i=-20, 20
          rho_restricted_1 = rho_restricted
          rho_restricted_1(i_SM_nb_par) = rho_restricted_1(i_SM_nb_par) + i*delta

          CALL PSI_FOR_M_RHO_ESTIMATION(beta_restricted, rho_restricted_1,H, is_rob_estimation)

          IF (write_output_3) THEN
             WRITE(UNIT=dev_unit_nb,FMT='(A11,TR1,A)') "Delta","Derivative:"
             WRITE(UNIT=dev_unit_nb,FMT=*) rho_restricted_1(i_SM_nb_par), H 
          END IF
       END DO
    END DO  !i_SM_nb_par

  END SUBROUTINE WRITE_FUNCTION_M

END MODULE M_RHO_MATRIX



MODULE STEPS
  USE DATATYPES
  IMPLICIT NONE

CONTAINS

  SUBROUTINE CONSTRAINTS_STEP_1(mode, ncnln, n, ldJ, needc, x, c, cJac, nstate)
    IMPLICIT NONE

    INTEGER, INTENT(IN OUT) :: mode
    INTEGER, INTENT(IN) :: ncnln, n, ldJ, nstate
    INTEGER, INTENT(IN), DIMENSION(1:n) :: needc

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: x
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:ncnln) :: c
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:ldJ,1:n) :: cJac

    RETURN
  END SUBROUTINE CONSTRAINTS_STEP_1


  SUBROUTINE CONSTRAINTS_STEP_2(mode, ncnln, n, ldJ, needc, x, c, cJac, nstate)
    IMPLICIT NONE

    INTEGER, INTENT(IN OUT) :: mode
    INTEGER, INTENT(IN) :: ncnln, n, ldJ, nstate
    INTEGER, INTENT(IN), DIMENSION(1:n) :: needc

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: x
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:ncnln) :: c
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:ldJ,1:n) :: cJac

    RETURN
  END SUBROUTINE CONSTRAINTS_STEP_2


  SUBROUTINE SNP_CRITERION_STEP_1(mode, n, beta_restricted, likelihood, g, nstate)
    USE S1_COMMON, ONLY: T1, modified_series
    USE LIKELIHOOD_UTILITIES, ONLY: LIKELIHOOD_UPDATE, CDF_CONSTRUCTION, GRADIENT_STEP_1
    IMPLICIT NONE

    INTEGER, INTENT(IN OUT) :: mode
    INTEGER, INTENT(IN) :: n, nstate
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: beta_restricted
    REAL (KIND=realkind), INTENT(OUT) :: likelihood
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:n) :: g

    CALL LIKELIHOOD_UPDATE(beta_restricted,T1,modified_series,1,mode)
    CALL CDF_CONSTRUCTION(T1,likelihood)
    CALL GRADIENT_STEP_1(mode,n,g,T1)

    RETURN
  END SUBROUTINE SNP_CRITERION_STEP_1


  SUBROUTINE SNP_CRITERION_STEP_1_UPDATES(mode, n, beta_restricted, likelihood, g, nstate)
    USE S2_COMMON, ONLY: T2, simulated_series_step_2
    USE LIKELIHOOD_UTILITIES, ONLY: LIKELIHOOD_UPDATE, CDF_CONSTRUCTION, GRADIENT_STEP_1
    IMPLICIT NONE

    INTEGER, INTENT(IN OUT) :: mode
    INTEGER, INTENT(IN) :: n, nstate
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: beta_restricted
    REAL (KIND=realkind), INTENT(OUT) :: likelihood
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:n) :: g

    CALL LIKELIHOOD_UPDATE(beta_restricted,T2,simulated_series_step_2,1,mode)
    CALL CDF_CONSTRUCTION(T2,likelihood)
    CALL GRADIENT_STEP_1(mode,n,g,T2)

    RETURN
  END SUBROUTINE SNP_CRITERION_STEP_1_UPDATES


  FUNCTION MY_SNP_CRITERION_STEP_1(beta_restricted) RESULT(function_value)
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(:) :: beta_restricted
    REAL (KIND=realkind) :: function_value

    INTEGER :: mode, n
    REAL (KIND=realkind), DIMENSION(1:SIZE(beta_restricted,1)) :: g 
    mode = 0
    n = SIZE(beta_restricted,1)

    CALL SNP_CRITERION_STEP_1(mode, n, beta_restricted, function_value, g, 1)

    RETURN
  END FUNCTION MY_SNP_CRITERION_STEP_1


  SUBROUTINE SNP_CRITERION_STEP_2(mode, n, rho, q, g, nstate)
    USE S1_COMMON,   ONLY: AM_nb_par, AM_nb_free_par, beta_AM_end, t_start, grad_i, &
         inv_score_cov, AM_nb_not_fixed_par, write_output_3, dev_unit_nb
    USE S2_COMMON,   ONLY: T2, SM_nb_par, rho_SM_start, rho_SM_end, rho_rest_SM_old, &
         nb_discard_step_2, random_iid_step_2, nb_iid_series, &
         simulated_series_step_2
    USE Z_MATRIX,    ONLY: is_new2
    USE UTILITIES,   ONLY: RESTRICT_BETA, COMPLETE_RHO_RESTRICTED
    USE SIMULATION, ONLY: STRUCTURAL_MODEL_SIMULATION
    USE LIKELIHOOD_UTILITIES, ONLY: SNP_ORT_FUNC, LIKELIHOOD_UPDATE
    USE DEVEL, ONLY: count
    IMPLICIT NONE

    INTEGER, INTENT(IN OUT) :: mode
    INTEGER, INTENT(IN) :: n, nstate
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: rho
    REAL (KIND=realkind), INTENT(OUT) :: q
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:n) :: g

    REAL (KIND=realkind), DIMENSION(1:AM_nb_free_par) :: beta_restricted
    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par) :: h
    INTEGER :: i, j, s
    character(len=10) :: nome

    !simulate the series but only if rho has changed
    DO s=1, n
       IF (rho(s) .NE. rho_rest_SM_old(s)) THEN
          rho_rest_SM_old = rho
          CALL COMPLETE_RHO_RESTRICTED(rho, rho_SM_end)
          IF (write_output_3) THEN
             WRITE(UNIT=dev_unit_nb,FMT=*) "chiamata a SNP_CRITERION_STEP_2"
             WRITE(UNIT=dev_unit_nb,FMT=*) "   rho:", rho_SM_end      
          END IF

          i = T2+nb_discard_step_2
          CALL STRUCTURAL_MODEL_SIMULATION(SM_nb_par,rho_SM_end,T2,simulated_series_step_2, &
               nb_discard_step_2,nb_iid_series,random_iid_step_2)
          is_new2 = 1
          if (write_output_3) THEN
             WRITE(nome,fmt="(I3)") count
             OPEN (UNIT=22,FILE=nome,FORM="FORMATTED")
             WRITE(unit=22,fmt=*) "   rho:", rho_SM_end 
             DO j=1, T2
                write(unit=22,fmt=*) j, simulated_series_step_2(j)
             END DO
             CLOSE(unit=22)
             count = count + 1
          end if
          EXIT
       END IF
    END DO

    !remove the restricted components from the vector of estimates
    CALL RESTRICT_BETA(beta_AM_end, beta_restricted)

    !IF (write_output_3) THEN
    !   WRITE(UNIT=dev_unit_nb,FMT=*) "  beta:", beta_AM_end      
    !END IF
    j = 1
    CALL LIKELIHOOD_UPDATE(beta_restricted,T2,simulated_series_step_2,2,j)
    CALL SNP_ORT_FUNC(T2)

    !compute the vector h for the quadratic form 
    h = 0._realkind
    DO j=1,AM_nb_not_fixed_par
       DO i=1, T2-t_start+1
          h(j) = h(j) + grad_i(i,j)
       END DO
    END DO

    h = h/REAL(T2-t_start+1)
    q = 0._realkind
    DO j=1, AM_nb_not_fixed_par
       q = q + inv_score_cov(j,j)*h(j)**2
       DO i=j+1, AM_nb_not_fixed_par
          q = q + 2._realkind*inv_score_cov(i,j)*h(i)*h(j)
       END DO
    END DO

    IF (write_output_3) THEN
       WRITE(UNIT=dev_unit_nb,FMT=*) "      h=", h
       WRITE(UNIT=dev_unit_nb,FMT=*) "   h'Mh=", q
       WRITE(UNIT=dev_unit_nb,FMT=*)
    END IF

    RETURN
  END SUBROUTINE SNP_CRITERION_STEP_2


  FUNCTION MY_SNP_CRITERION_STEP_2(rho_restricted) RESULT(function_value)
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(:) :: rho_restricted
    REAL (KIND=realkind) :: function_value

    INTEGER :: mode, n
    REAL (KIND=realkind), DIMENSION(1:SIZE(rho_restricted,1)) :: g 
    mode = 0
    n = SIZE(rho_restricted,1)

    CALL SNP_CRITERION_STEP_2(mode, n, rho_restricted, function_value, g, 1)

    RETURN
  END FUNCTION MY_SNP_CRITERION_STEP_2


  SUBROUTINE ROB_CRITERION_STEP_1(mode, n, beta_restricted, q, g, nstate)
    USE S1_COMMON,   ONLY: N1, AM_nb_not_fixed_par
    USE S2_COMMON,   ONLY: SM_nb_free_par
    USE ROB_COMMON,  ONLY: psi_N1
    USE MYFUNCTIONS, ONLY: MM
    USE ROBUST_UTILITIES, ONLY: WEIGHTS_HUBER, UPDATE_SCORE_N1, UPDATE_PSI_N1
    IMPLICIT NONE

    INTEGER, INTENT(IN OUT) :: mode
    INTEGER, INTENT(IN) :: n, nstate
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: beta_restricted
    REAL (KIND=realkind), INTENT(OUT) :: q
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:n) :: g

    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par) :: H
    REAL (KIND=realkind) :: norm
    INTEGER :: i,j

    !compute the gradient for T1
    CALL UPDATE_SCORE_N1(beta_restricted)

    CALL UPDATE_PSI_N1()

    !sum the columns of psi_N1
    H=0._realkind
    DO j=1, AM_nb_not_fixed_par
       DO i=1, N1
          H(j) = H(j) + psi_N1(i,j)
       END DO
    END DO

    H=H/REAL(N1)

    norm=0._realkind

    DO i=1,AM_nb_not_fixed_par
       norm=norm+h(i)*h(i)
    END DO

    q = norm

    RETURN
  END SUBROUTINE ROB_CRITERION_STEP_1


  FUNCTION MY_ROB_CRITERION_STEP_1(beta_restricted) RESULT(function_value)
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(:) :: beta_restricted
    REAL (KIND=realkind) :: function_value

    INTEGER :: mode, n
    REAL (KIND=realkind), DIMENSION(1:SIZE(beta_restricted,1)) :: g 
    mode = 0
    n = SIZE(beta_restricted,1)

    CALL ROB_CRITERION_STEP_1(mode, n, beta_restricted, function_value, g, 1)

    RETURN
  END FUNCTION MY_ROB_CRITERION_STEP_1


  SUBROUTINE ROB_CRITERION_STEP_2(mode, n, rho, q, g, nstate)
    USE DATATYPES
    USE S1_COMMON,   ONLY: AM_nb_not_fixed_par, AM_nb_free_par, beta_AM_end, &
         write_output_3, dev_unit_nb
    USE S2_COMMON,   ONLY: T2, N2, SM_nb_par, rho_SM_end, rho_rest_SM_old, SM_nb_free_par, &
         nb_discard_step_2, nb_iid_series, simulated_series_step_2, random_iid_step_2
    USE ROB_COMMON,  ONLY: psi_N2, S, weights_metric
    USE Z_MATRIX,             ONLY: is_new2
    USE UTILITIES,            ONLY: RESTRICT_BETA, COMPLETE_RHO_RESTRICTED
    USE SIMULATION,           ONLY: STRUCTURAL_MODEL_SIMULATION
    USE MYFUNCTIONS,          ONLY: MM
    USE ROBUST_UTILITIES,     ONLY: WEIGHTS_HUBER, UPDATE_PSI_TMP_N2, UPDATE_PSI_N2

    IMPLICIT NONE

    INTEGER, INTENT(IN OUT) :: mode
    INTEGER, INTENT(IN) :: n, nstate
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: rho
    REAL (KIND=realkind), INTENT(OUT) :: q
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:n) :: g


    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par) :: H
    REAL (KIND=realkind), DIMENSION(1:AM_nb_free_par) :: beta_restricted
    REAL (KIND=realkind) :: norm
    INTEGER :: i,j,l

    !simulate the series but only if rho has changed
    DO l=1, n
       IF (rho(l) .NE. rho_rest_SM_old(l)) THEN
          rho_rest_SM_old = rho
          !IF (write_output_3) THEN
          !   WRITE(UNIT=dev_unit_nb,FMT=*) "chiamata a ROB_CRITERION_STEP_2"
          !   WRITE(UNIT=dev_unit_nb,FMT=*) "   rho:", rho      
          !END IF

          CALL COMPLETE_RHO_RESTRICTED(rho,rho_SM_end)

          CALL STRUCTURAL_MODEL_SIMULATION(SM_nb_par,rho_SM_end,T2,simulated_series_step_2, &
               nb_discard_step_2,nb_iid_series,random_iid_step_2)

          is_new2 = 1
          EXIT
       END IF
    END DO

    CALL RESTRICT_BETA(beta_AM_end, beta_restricted)

    !compute the gradient for T2
    CALL UPDATE_PSI_TMP_N2(beta_restricted)

    CALL UPDATE_PSI_N2()

    !sum the columns of psi_N2
    H=0._realkind
    DO j=1, AM_nb_not_fixed_par
       DO i=1, N2
          H(j) = H(j) + psi_N2(i,j)
       END DO
    END DO

    H=H/REAL(N2)

    norm=0._realkind

    IF (weights_metric .EQ. 0) THEN
       DO i=1,AM_nb_not_fixed_par
         norm=norm+H(i)**2
       END DO
    ELSE
       DO j=1, AM_nb_not_fixed_par
          norm = norm + S(j,j)*H(j)**2
          DO i=j+1, AM_nb_not_fixed_par
             norm = norm + 2._realkind*S(i,j)*H(i)*H(j)
          END DO
       END DO
    END IF

    q = norm

    !IF (write_output_3) THEN
    !   WRITE(UNIT=dev_unit_nb,FMT=*) "ROB_CRITERION_STEP_2 -> q:", q      
    !END IF

    RETURN
  END SUBROUTINE ROB_CRITERION_STEP_2



  FUNCTION MY_ROB_CRITERION_STEP_2(rho_restricted) RESULT(function_value)
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(:) :: rho_restricted
    REAL (KIND=realkind) :: function_value

    INTEGER :: mode, n
    REAL (KIND=realkind), DIMENSION(1:SIZE(rho_restricted,1)) :: g 
    mode = 0
    n = SIZE(rho_restricted,1)

    CALL ROB_CRITERION_STEP_2(mode, n, rho_restricted, function_value, g, 1)

    RETURN
  END FUNCTION MY_ROB_CRITERION_STEP_2



  SUBROUTINE DO_AM_GRID_SEARCH(beta,which_criterion,grid_search_log_lik)
    USE DATATYPES
    USE S1_COMMON, ONLY: T1, AM_nb_par, AM_nb_free_par, write_output_2, tmp_unit_nb
    USE UTILITIES, ONLY: RESTRICT_BETA, CHECK_CONSTRAINTS_STEP_1
    USE GRID_SEARCH_AM
    IMPLICIT NONE

    CHARACTER (LEN=3), INTENT(IN) :: which_criterion
    REAL (KIND=realkind), INTENT(OUT), DIMENSION(1:AM_nb_par) :: beta
    REAL (KIND=realkind), INTENT(OUT) :: grid_search_log_lik    

    INTEGER :: i,p,s,first_grid_par,last_grid_par,first_completed, last_completed
    INTEGER, DIMENSION(1:AM_nb_par) :: k
    INTEGER :: mode, nstate
    REAL (KIND=realkind) :: minimum, likelihood
    REAL (KIND=realkind), DIMENSION(1:AM_nb_free_par) :: g, beta_restricted
    REAL (KIND=realkind), DIMENSION(1:AM_nb_par) :: beta_grid
    LOGICAL :: fullfilled, first_time

    mode = 0
    nstate = 1
    k=1
    first_time = .TRUE.
    grid_search_AM_successfull = .FALSE.

    IF (write_output_2) THEN
       IF (which_criterion .EQ. "snp") THEN
          OPEN (UNIT=tmp_unit_nb,FILE="grid_search_step_1.dat",FORM="FORMATTED")
       ELSE
          OPEN (UNIT=tmp_unit_nb,FILE="grid_search_step_1_rob.dat",FORM="FORMATTED")
       END IF
    END IF

    DO i=1,AM_nb_par
       IF (is_grid_search_AM(i)) THEN
          first_grid_par = i
          EXIT
       END IF
    END DO

    i=AM_nb_par
    DO
       IF (is_grid_search_AM(i)) THEN
          last_grid_par = i
          EXIT
       END IF
       i=i-1
    END DO

    DO i=1,nb_grid_searches_AM
       DO s=1, AM_nb_par
          p = cum_grid_points_AM(s)+k(s)
          beta_grid(s) = grid_points_values_AM(p)
       END DO
       CALL CHECK_CONSTRAINTS_STEP_1(beta_grid,fullfilled)
       IF (fullfilled) THEN

          grid_search_AM_successfull = .TRUE.

          CALL RESTRICT_BETA(beta_grid, beta_restricted)
          IF (which_criterion .EQ. "snp") THEN
             CALL SNP_CRITERION_STEP_1(mode, AM_nb_free_par, beta_restricted, likelihood, g, nstate)
          ELSE
             CALL ROB_CRITERION_STEP_1(mode, AM_nb_free_par, beta_restricted, likelihood, g, nstate)
          END IF

          grid_search_results_AM(i) = likelihood

          IF (first_time) THEN
             beta = beta_grid
             minimum = likelihood
          ELSE
             IF (minimum .GT. likelihood) THEN
                minimum = likelihood
                beta =  beta_grid
             END IF
          END IF
          first_time = .FALSE.
          IF (write_output_2) THEN
             WRITE(UNIT=tmp_unit_nb,fmt=*) "grid search nb:", i, grid_search_results_AM(i), beta_grid
          END IF
       END IF !IF (fullfilled)

       k(first_grid_par) = k(first_grid_par) + 1
       IF (k(first_grid_par) .GT. nb_grid_points_AM(first_grid_par)) THEN
          k(first_grid_par) = 1

          DO s=1+first_grid_par, AM_nb_par
             IF (is_grid_search_AM(s)) THEN
                k(s) = k(s) + 1
                IF (k(s) .GT. nb_grid_points_AM(s)) THEN
                   k(s) = 1
                ELSE
                   EXIT
                END IF
             END IF
          END DO
       END IF
    END DO

    IF (write_output_2) THEN
       CLOSE (UNIT=tmp_unit_nb)
    END IF

    IF (grid_search_AM_successfull) THEN
       grid_search_log_lik = minimum
    END IF

    RETURN
  END SUBROUTINE DO_AM_GRID_SEARCH


  SUBROUTINE DO_SM_GRID_SEARCH(rho,which_criterion,grid_search_minimum)
    USE DATATYPES
    USE S1_COMMON, ONLY: write_output_2, tmp_unit_nb
    USE S2_COMMON, ONLY: T2, SM_nb_par, SM_nb_free_par
    USE UTILITIES, ONLY: RESTRICT_RHO, CHECK_CONSTRAINTS_STEP_2
    USE GRID_SEARCH_SM
    IMPLICIT NONE

    CHARACTER (LEN=3), INTENT(IN) :: which_criterion
    REAL (KIND=realkind), INTENT(OUT), DIMENSION(1:SM_nb_par) :: rho
    REAL (KIND=realkind), INTENT(OUT) :: grid_search_minimum    

    INTEGER :: i, p, s, first_grid_par, last_grid_par
    INTEGER, DIMENSION(1:SM_nb_par) :: k
    INTEGER :: mode, nstate
    REAL (KIND=realkind) :: minimum, criterion
    REAL (KIND=realkind), DIMENSION(1:SM_nb_free_par) :: g, rho_restricted
    REAL (KIND=realkind), DIMENSION(1:SM_nb_par) :: rho_grid

    LOGICAL :: fullfilled, first_time

    mode = 0
    nstate = 1
    k=1
    first_time = .TRUE.
    grid_search_SM_successfull = .FALSE.


    IF (write_output_2) THEN
       IF (which_criterion .EQ. "snp") THEN
          OPEN (UNIT=tmp_unit_nb,FILE="grid_search_step_2.dat",FORM="FORMATTED")
       ELSE
          OPEN (UNIT=tmp_unit_nb,FILE="grid_search_step_2_rob.dat",FORM="FORMATTED")
       END IF
    END IF

    DO i=1,SM_nb_par
       IF (is_grid_search_SM(i)) THEN
          first_grid_par = i
          EXIT
       END IF
    END DO

    i=SM_nb_par
    DO
       IF (is_grid_search_SM(i)) THEN
          last_grid_par = i
          EXIT
       END IF
       i=i-1
    END DO

    DO i=1,nb_grid_searches_SM
       DO s=1, SM_nb_par
          p = cum_grid_points_SM(s)+k(s)
          rho_grid(s) = grid_points_values_SM(p)
       END DO

       CALL CHECK_CONSTRAINTS_STEP_2(rho_grid,fullfilled)

       IF (fullfilled) THEN
          grid_search_SM_successfull = .TRUE.
          CALL RESTRICT_RHO(rho_grid, rho_restricted)

          IF (which_criterion .EQ. "snp") THEN

             CALL SNP_CRITERION_STEP_2(mode, SM_nb_free_par, rho_restricted, criterion, g, nstate)

          ELSE
             CALL ROB_CRITERION_STEP_2(mode, SM_nb_free_par, rho_restricted, criterion, g, nstate)
          END IF

          grid_search_results_SM(i) = criterion

          IF (first_time) THEN
             rho = rho_grid
             minimum = criterion
          ELSE
             IF (minimum .GT. criterion) THEN
                minimum = criterion
                rho = rho_grid
             END IF
          END IF
          first_time = .FALSE.

          IF (write_output_2) THEN
             WRITE(UNIT=tmp_unit_nb,fmt=*) "grid search nb:", i, grid_search_results_SM(i), rho_grid
          END IF
       END IF !IF (fullfilled)

       k(first_grid_par) = k(first_grid_par) + 1
       IF (k(first_grid_par) .GT. nb_grid_points_SM(first_grid_par)) THEN
          k(first_grid_par) = 1

          DO s=1+first_grid_par, SM_nb_par
             IF (is_grid_search_SM(s)) THEN
                k(s) = k(s) + 1
                IF (k(s) .GT. nb_grid_points_SM(s)) THEN
                   k(s) = 1
                ELSE
                   EXIT
                END IF
             END IF
          END DO
       END IF
    END DO

    IF (write_output_2) THEN
       CLOSE (UNIT=tmp_unit_nb)
    END IF

    IF (grid_search_SM_successfull) THEN
       grid_search_minimum = minimum
    END IF

    RETURN
  END SUBROUTINE DO_SM_GRID_SEARCH


  SUBROUTINE STEP_1()
    USE S1_COMMON,   ONLY: AM_nb_par, AM_nb_free_par, beta_AM_end, beta_AM_start, &
         nclin, ncnln, ldA, ldJ, ldR, A, bl, bu, inform, iter, istate, &
         cJac, clamda, RR, iw, leniw, w, lenw, npsol_inform_step_1, npsol_iter_step_1, &
         tmp_unit_nb, log_unit_nb, write_output_1, write_output_2, space_lenght, use_npsol, &
         minimum_step_1_vector, do_grid_search
    USE UTILITIES,        ONLY: RESTRICT_BETA, COMPLETE_BETA_RESTRICTED
    USE STARTUP,          ONLY: GRID_POINTS_AM
    USE Z_MATRIX,         ONLY: beta_old
    USE MONTE_CARLO,      ONLY: nb_rep_monte_carlo, mc_counter
    USE MYINTERFACE,      ONLY: npsol
    USE GRID_SEARCH_AM,   ONLY: do_grid_search_AM, grid_search_AM_successfull
    USE OUTPUT_UTILITIES, ONLY: OUTPUT_SNP_STEP_1, NEW_RUN
    USE MY_NPSOL
    IMPLICIT NONE

    REAL (KIND=realkind), DIMENSION(1:AM_nb_free_par) :: beta_restricted

    !variables necessary in npsol
    REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: c
    REAL (KIND=realkind), DIMENSION(1:AM_nb_free_par) :: g

    !other variables
    INTEGER :: j
    REAL (KIND=realkind) :: f, start_log_lik, grid_search_log_lik
    REAL (KIND=realkind), DIMENSION(1:AM_nb_par) :: beta_grid_search
    CHARACTER (LEN=space_lenght) :: space

100 FORMAT(A,A,TR1,ES11.4," | ",25(ES11.4,TR1))
101 FORMAT(A5,TR1,ES11.4," | ",25(ES11.4,TR1))

    IF (ncnln .EQ. 0) THEN
       ALLOCATE (c(1:1))
    ELSE
       ALLOCATE (c(1:ncnln))
    END IF

    space = " "
    space_lenght = space_lenght + 3
    grid_search_AM_successfull = .FALSE.
    !compute the criterion at the proposed starting value
    CALL RESTRICT_BETA(beta_AM_start,beta_restricted)
    beta_old = -999999.9_realkind
    j=0
    CALL SNP_CRITERION_STEP_1(j,AM_nb_free_par,beta_restricted,start_log_lik,g,1)
    beta_AM_end = beta_AM_start
    !write the output to the screen
    write(*,fmt=101) "the: ", start_log_lik, beta_AM_start

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=100) space, "log-L & the:", start_log_lik, &
            beta_restricted
    END IF

    !if required perform a first grid search
    IF (do_grid_search_AM .AND. ((do_grid_search .EQ. 1) .OR. (do_grid_search .EQ. 3))) THEN
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "START GRID_SEARCH"
       END IF
       CALL GRID_POINTS_AM(AM_nb_par,beta_AM_start)
       CALL DO_AM_GRID_SEARCH(beta_grid_search,"snp",grid_search_log_lik)
       IF (grid_search_AM_successfull) THEN  !i.e. if at least 1 grid_point is feasible
          IF (grid_search_log_lik .LT. start_log_lik) THEN
             CALL RESTRICT_BETA(beta_grid_search,beta_restricted)
             IF (write_output_2) THEN
                WRITE(UNIT=log_unit_nb,FMT=*) space, "GRID SEARCH SUCCESSFULL: NEW START VALUES"
                WRITE(UNIT=log_unit_nb,FMT=100) space, "log-L & the:", grid_search_log_lik, &
                     beta_restricted
             END IF
             beta_AM_end = beta_grid_search
          ELSE
             IF (write_output_2) THEN
                WRITE(UNIT=log_unit_nb,FMT=*) space, "GRID SEARCH NOT SUCCESSFULL"
             END IF
          END IF
       ELSE
          IF (write_output_2) THEN
             WRITE(UNIT=log_unit_nb,FMT=*) space, "GRID SEARCH NOT SUCCESSFULL"
          END IF
       END IF
    ELSE
       !set the grid_search_log_lik = 99999999.9
       grid_search_log_lik = 99999999._realkind
    END IF

    IF ((nb_rep_monte_carlo .LT. 2) .AND. write_output_1) THEN
       !write the model.dat file in the estimation_path.txt file
       OPEN (UNIT=tmp_unit_nb,FILE="estimation_path.txt",STATUS="UNKNOWN",POSITION="APPEND",FORM="FORMATTED")
       !       CALL NEW_RUN(tmp_unit_nb)     
       WRITE(UNIT=tmp_unit_nb,FMT="(A)") "          *******************************************************" 
       CLOSE (UNIT=tmp_unit_nb)
    END IF

    beta_old = -999999.9_realkind

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "START NPSOL"
    END IF

    IF (use_npsol(1)) THEN
       CALL npsol(AM_nb_free_par,nclin,ncnln,ldA,ldJ,ldR,A,bl,bu,CONSTRAINTS_STEP_1, &
            SNP_CRITERION_STEP_1,inform,iter,istate,c,cJac,clamda,f,g,RR, &
            beta_restricted,iw,leniw,w,lenw )
    ELSE
       CALL minimum(MY_SNP_CRITERION_STEP_1,beta_restricted,f,iter)
       inform=iter
    END IF

    minimum_step_1_vector(mc_counter) = f

    !add the restricted components to the parameter vector
    CALL COMPLETE_BETA_RESTRICTED(beta_restricted, beta_AM_end)
    !write to the screen the estimated values
    write(*,fmt=101) "the: ", f, beta_AM_end

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=100) space, "log-L &^the:", f, beta_restricted
    END IF

    IF ((nb_rep_monte_carlo .GT. 0) .AND. (write_output_1)) THEN
       npsol_inform_step_1(mc_counter)=inform
       npsol_iter_step_1(mc_counter)=iter
    END IF

    IF ((nb_rep_monte_carlo .LT. 2) .AND. (write_output_1)) THEN
       IF (do_grid_search_AM) THEN
          IF (grid_search_log_lik .LT. start_log_lik) THEN
             CALL OUTPUT_SNP_STEP_1(f,grid_search_log_lik)
          ELSE
             CALL OUTPUT_SNP_STEP_1(f,start_log_lik)
          END IF
       ELSE
          CALL OUTPUT_SNP_STEP_1(f,start_log_lik)
       END IF
    END IF

    space_lenght = space_lenght - 3

    DEALLOCATE(c)

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "EXIT STEP_1"
    END IF

    RETURN

  END SUBROUTINE STEP_1


  SUBROUTINE STEP_2()
    USE S1_COMMON,   ONLY: nclin, ncnln, ldA, ldJ, ldR, A, bl, bu, &
         inform, iter, istate, cJac, clamda, RR, iw, leniw, w, &
         lenw, write_output_1, write_output_2, log_unit_nb, T1, space_lenght, &
         use_npsol, do_grid_search
    USE S2_COMMON,   ONLY: rho_SM_end, rho_SM_start, SM_nb_free_par, &
         SM_nb_par, rho_rest_SM_start, rho_rest_SM_end, rho_rest_SM_old, &
         npsol_inform_step_2, npsol_iter_step_2, hansen_test_results, &
         minimum_step_2_vector
    USE STARTUP,     ONLY: GRID_POINTS_SM
    USE MONTE_CARLO, ONLY: nb_rep_monte_carlo, mc_counter
    USE MYINTERFACE, ONLY: npsol
    USE UTILITIES,   ONLY: RESTRICT_RHO, COMPLETE_RHO_RESTRICTED
    USE OUTPUT_UTILITIES, ONLY: OUTPUT_SNP_STEP_2
    USE GRID_SEARCH_SM,   ONLY: do_grid_search_SM, grid_search_SM_successfull
    USE MY_NPSOL
    IMPLICIT NONE

    !variables necessary in npsol
    REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: c
    REAL (KIND=realkind) :: f, criterion
    REAL (KIND=realkind), DIMENSION(1:SM_nb_free_par) :: g

    !other variables
    REAL (KIND=realkind), DIMENSION(1:SM_nb_par) :: rho_grid_search
    INTEGER :: j
    CHARACTER (LEN=space_lenght) :: space

100 FORMAT(A,A,TR1,ES11.4," | ",25(ES11.4,TR1))
101 FORMAT(A5,TR1,ES11.4," | ",25(ES11.4,TR1))

    IF (ncnln .EQ. 0) THEN
       ALLOCATE (c(1:1))
    ELSE
       ALLOCATE (c(1:ncnln))
    END IF

    space = " "
    space_lenght = space_lenght + 3
    grid_search_SM_successfull = .FALSE.

    !optimize the objective function with respect to rho_restricted
    CALL RESTRICT_RHO(rho_SM_start,rho_rest_SM_start)
    rho_rest_SM_old = -999999.9_realkind
    j=0

    CALL SNP_CRITERION_STEP_2(j, SM_nb_free_par, rho_rest_SM_start, f, g, 1)
    !write to the screen the starting values
    write(*,fmt=101) "rho: ", f, rho_SM_start

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=100) space, "q'Mq |  rho:", f, rho_rest_SM_start
    END IF

    !if required perform a first grid search
    IF (do_grid_search_SM .AND. ((do_grid_search .EQ. 1) .OR. (do_grid_search.EQ. 3))) THEN
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "START GRID_SEARCH"
       END IF
       CALL GRID_POINTS_SM(SM_nb_par,rho_SM_start)
       CALL DO_SM_GRID_SEARCH(rho_grid_search,"snp",criterion)
       IF (grid_search_SM_successfull) THEN
          IF (criterion .LT. f) THEN
             CALL RESTRICT_RHO(rho_grid_search,rho_rest_SM_start)
             IF (write_output_2) THEN
                WRITE(UNIT=log_unit_nb,FMT=*) space,  "GRID SEARCH SUCCESSFULL: NEW START VALUES"
                WRITE(UNIT=log_unit_nb,FMT=100) space,  "q'Mq |  rho:", criterion, &
                     rho_rest_SM_start
             END IF
          ELSE
             IF (write_output_2) THEN
                WRITE(UNIT=log_unit_nb,FMT=*) space, "GRID SEARCH NOT SUCCESSFULL"
             END IF
          END IF
       ELSE
          IF (write_output_2) THEN
             WRITE(UNIT=log_unit_nb,FMT=*) space, "GRID SEARCH NOT SUCCESSFULL"
          END IF
       END IF
    END IF

    rho_rest_SM_old = -999999.9_realkind
    rho_rest_SM_end = rho_rest_SM_start

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "START NPSOL"
    END IF

    IF (use_npsol(2)) THEN
       CALL npsol( SM_nb_free_par,nclin,ncnln,ldA,ldJ,ldR,A,bl,bu,CONSTRAINTS_STEP_2, &
            SNP_CRITERION_STEP_2,inform,iter,istate,c,cJac,clamda,f,g,RR, &
            rho_rest_SM_end,iw,leniw,w,lenw )
    ELSE
       CALL minimum(MY_SNP_CRITERION_STEP_2,rho_rest_SM_end,f,iter)
       inform=iter
    END IF

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=100) space, "q'Mq | ^rho:", f, rho_rest_SM_end
    END IF

    minimum_step_2_vector(mc_counter) = f
    CALL COMPLETE_RHO_RESTRICTED(rho_rest_SM_end,rho_SM_end)
    !write to the screen the estimated values
    write(*,fmt=101) "rho: ", f, rho_SM_end

    IF ((nb_rep_monte_carlo .GT. 0) .AND. (write_output_1)) THEN
       hansen_test_results(mc_counter,1) = T1*f
       npsol_inform_step_2(mc_counter,1) = inform
       npsol_iter_step_2(mc_counter,1) = iter
    END IF


    IF ((nb_rep_monte_carlo .LT. 2) .AND. (write_output_1)) THEN
       CALL OUTPUT_SNP_STEP_2(f)
    END IF

    space_lenght = space_lenght - 3

    DEALLOCATE(c)

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "EXIT STEP_2"
    END IF

    RETURN
  END SUBROUTINE STEP_2


  SUBROUTINE EXPECTATION_MC()
    USE S1_COMMON, ONLY:  T1, nclin, ncnln, ldA, ldJ, ldR, A, bl, bu, &
         inform, iter, istate, cJac, clamda, RR, iw, leniw, w, &
         lenw, write_output_1, write_output_2, log_unit_nb, write_output_3, dev_unit_nb, &
         use_npsol, beta_AM_end, AM_nb_par, AM_nb_not_fixed_par, inv_score_cov, &
         space_lenght
    USE S2_COMMON, ONLY: SM_nb_par, SM_nb_free_par, rho_rest_SM_old, & 
         nb_discard_step_2, nb_iid_series, random_iid_step_2, T2, max_T2, &
         npsol_inform_step_2, npsol_iter_step_2, hansen_test_results
    USE MONTE_CARLO, ONLY: mc_results_step_1, mc_results_step_2, & 
         nb_rep_monte_carlo, inv_score_cov_mc, nb_rep_exp_analysis
    USE STARTUP,     ONLY: NPSOL_SETUP
    USE UTILITIES,   ONLY: RESTRICT_RHO, COMPLETE_RHO_RESTRICTED
    USE MYFUNCTIONS, ONLY: IID_SIM
    USE MY_NPSOL
    USE OUTPUT_UTILITIES, ONLY: OUTPUT_EXPECTATION_MC_STEP_2
    IMPLICIT NONE

    INTEGER :: mc_counter, exp_counter, i, j
    REAL (KIND=realkind), DIMENSION(1:SM_nb_free_par) :: rho_restricted
    REAL (KIND=realkind), DIMENSION(1:(max_T2+nb_discard_step_2),1:nb_iid_series) :: iid_tmp 
    CHARACTER (LEN=space_lenght) :: space

    !variables necessary in npsol
    REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: c
    REAL (KIND=realkind) :: f
    REAL (KIND=realkind), DIMENSION(1:SM_nb_free_par) :: g

    IF (ncnln .EQ. 0) THEN
       ALLOCATE (c(1:1))
    ELSE
       ALLOCATE (c(1:ncnln))
    END IF

    space = " "
    space_lenght = space_lenght + 3

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "NPSOL_SETUP(2)"
    END IF
    CALL NPSOL_SETUP(2)

    DO exp_counter=1, nb_rep_exp_analysis
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "SIMULATE THE random_iid_step_2 SERIES:", exp_counter
       END IF

       !simulate again the iid series for the approximation of the expectation  
       iid_tmp = IID_SIM(max_T2+nb_discard_step_2,nb_iid_series,"NORMAL_")
       DO j=1, nb_iid_series
          DO i=1, T2+nb_discard_step_2
             random_iid_step_2(i,j) = iid_tmp(i,j)
          END DO
       END DO

       DO mc_counter=1, nb_rep_monte_carlo
          rho_rest_SM_old = -999999.9_realkind
          CALL RESTRICT_RHO(mc_results_step_2(1:SM_nb_par,mc_counter,1),rho_restricted)
          DO i=1, AM_nb_par
             beta_AM_end(i) = mc_results_step_1(i,mc_counter)
          END DO
          IF (write_output_3) THEN
             WRITE(UNIT=dev_unit_nb,FMT=*) "beta_am_end for exp:", mc_counter, exp_counter, beta_am_end(:)
             WRITE(UNIT=dev_unit_nb,FMT=*) "rho_rest.   for exp:", rho_restricted(:)
          END IF
          inv_score_cov(:,:) = inv_score_cov_mc(:,:,nb_rep_monte_carlo)

          IF (use_npsol(2)) THEN
             CALL npsol( SM_nb_free_par,nclin,ncnln,ldA,ldJ,ldR,A,bl,bu,CONSTRAINTS_STEP_2, &
                  SNP_CRITERION_STEP_2,inform,iter,istate,c,cJac,clamda,f,g,RR, &
                  rho_restricted,iw,leniw,w,lenw )
          ELSE
             CALL minimum(MY_SNP_CRITERION_STEP_2,rho_restricted,f,iter)
             inform=iter
          END IF
          !save the results for every mc_repetition and exp_repetitions

          CALL COMPLETE_RHO_RESTRICTED(rho_restricted,mc_results_step_2(:,mc_counter,1+exp_counter))

          IF (write_output_1) THEN
             hansen_test_results(mc_counter,1+exp_counter) = T1*f
             npsol_inform_step_2(mc_counter,1+exp_counter) = inform
             npsol_iter_step_2(mc_counter,1+exp_counter) = iter
          END IF
       END DO
    END DO

    !write the results
    IF ((write_output_2) .AND. (write_output_1)) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "WRITE THE OUTPUT"
    END IF
    IF (write_output_1) THEN
       CALL OUTPUT_EXPECTATION_MC_STEP_2()
    END IF

    space_lenght = space_lenght - 3

    RETURN
  END SUBROUTINE EXPECTATION_MC



  SUBROUTINE STEP_1_ROB(beta_restricted, AM_nb_free_par)
    USE S1_COMMON,   ONLY: AM_nb_par, AM_nb_not_fixed_par, nclin, ncnln, ldA, ldJ, &
         ldR, A, bl, bu, inform, iter, istate, cJac, clamda, RR, iw, leniw, w, lenw, &
         log_unit_nb, write_output_2, space_lenght, use_npsol, step, derivative_level, &
         do_grid_search
    USE Z_MATRIX,       ONLY: beta_old
    USE UTILITIES,      ONLY: RESTRICT_BETA
    USE ROB_COMMON,     ONLY: npsol_inform_step_1_rob, npsol_iter_step_1_rob, &
         minimum_step_1_rob, minimum_step_1_rob_vector
    USE MYINTERFACE,    ONLY: npsol
    USE MONTE_CARLO,    ONLY: nb_rep_monte_carlo, mc_counter
    USE GRID_SEARCH_AM, ONLY: do_grid_search_AM, grid_search_AM_successfull
    USE MY_NPSOL
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: AM_nb_free_par
    REAL (KIND = realkind), INTENT(IN OUT), DIMENSION(1:AM_nb_free_par) :: beta_restricted

    !variables necessary in npsol
    REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: c
    REAL (KIND=realkind), DIMENSION(1:AM_nb_free_par) :: g

    !other variables
    INTEGER :: j
    REAL (KIND=realkind) :: f, start_log_lik, grid_search_log_lik
    REAL (KIND=realkind), DIMENSION(1:AM_nb_par) :: beta_grid_search
    CHARACTER (LEN=space_lenght) :: space

100 FORMAT(A,A,TR1,ES11.4," | ",25(ES11.4,TR1))
    
    IF (ncnln .EQ. 0) THEN
       ALLOCATE (c(1:1))
    ELSE
       ALLOCATE (c(1:ncnln))
    END IF

    space = " "
    space_lenght = space_lenght + 3
    grid_search_AM_successfull = .FALSE.

    !compute the criterion at the proposed starting value
    beta_old = -999999.9_realkind
    j=0
    CALL ROB_CRITERION_STEP_1(j,AM_nb_free_par,beta_restricted,start_log_lik,g,1)
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=100) space, "q'Mq |  the:", start_log_lik, &
            beta_restricted
    END IF

    !if required perform a first grid search
    IF (do_grid_search_AM  .AND. ((do_grid_search .EQ. 2) .OR. (do_grid_search.EQ. 3))) THEN
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "START GRID_SEARCH"
       END IF
       CALL DO_AM_GRID_SEARCH(beta_grid_search,"rob",grid_search_log_lik)
       IF (grid_search_AM_successfull) THEN  !i.e. if at least 1 grid_point is feasible
          IF (grid_search_log_lik .LT. start_log_lik) THEN
             CALL RESTRICT_BETA(beta_grid_search,beta_restricted)
             IF (write_output_2) THEN
                WRITE(UNIT=log_unit_nb,FMT=*) space, "GRID SEARCH SUCCESSFULL: NEW START VALUES"
                WRITE(UNIT=log_unit_nb,FMT=100) space, "q'Mq |  the:", grid_search_log_lik, &
                     beta_restricted
             END IF
          ELSE
             IF (write_output_2) THEN
                WRITE(UNIT=log_unit_nb,FMT=*) space, "GRID SEARCH NOT SUCCESSFULL"
             END IF
          END IF
       ELSE
          IF (write_output_2) THEN
             WRITE(UNIT=log_unit_nb,FMT=*) space, "GRID SEARCH NOT SUCCESSFULL"
          END IF
       END IF
    END IF

    beta_old = -999999.9_realkind

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "START NPSOL"
    END IF

    !here the subroutine STEP_1_ROB() starts
    IF (use_npsol(3)) THEN
       CALL npsol(AM_nb_free_par,nclin,ncnln,ldA,ldJ,ldR,A,bl,bu,CONSTRAINTS_STEP_1, &
            ROB_CRITERION_STEP_1,inform,iter,istate,c,cJac,clamda,f,g,RR, &
            beta_restricted,iw,leniw,w,lenw )
       IF (nb_rep_monte_carlo .GT. 0) THEN
          !save the results
          npsol_inform_step_1_rob(mc_counter) = inform
          npsol_iter_step_1_rob(mc_counter) = iter
       END IF
    ELSE
       CALL minimum(MY_ROB_CRITERION_STEP_1,beta_restricted,f,iter)
       inform=iter
    END IF

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=100) space, "q'Mq | ^the:", f, beta_restricted
    END IF

    IF (nb_rep_monte_carlo .GE. 1) THEN
       npsol_inform_step_1_rob(mc_counter) = inform
       npsol_iter_step_1_rob(mc_counter) = iter
    END IF

    minimum_step_1_rob_vector(mc_counter) = f
    minimum_step_1_rob = f

    space_lenght = space_lenght - 3

    DEALLOCATE(c)

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "EXIT STEP_1_ROB"
    END IF

    RETURN

  END SUBROUTINE STEP_1_ROB


  SUBROUTINE STEP_2_ROB(rho_restricted,SM_nb_free_par)
    USE S1_COMMON,   ONLY: nclin, ncnln, ldA, ldJ, ldR, A, bl, bu, inform, iter, istate, &
         cJac, clamda, RR, iw, leniw, w, lenw, write_output_2, log_unit_nb, space_lenght, &
         use_npsol, step, derivative_level,T1, do_grid_search
    USE S2_COMMON,      ONLY: SM_nb_par, rho_rest_SM_end, rho_rest_SM_start, rho_rest_SM_old
    USE UTILITIES,      ONLY: RESTRICT_RHO
    USE ROB_COMMON,     ONLY: npsol_inform_step_2_rob, npsol_iter_step_2_rob, &
         minimum_step_2_rob, minimum_step_2_rob_vector, hansen_test_results_rob 
    USE MYINTERFACE,    ONLY: npsol
    USE MONTE_CARLO,    ONLY: nb_rep_monte_carlo, mc_counter
    USE GRID_SEARCH_SM, ONLY: do_grid_search_SM, grid_search_SM_successfull
    USE MY_NPSOL
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: SM_nb_free_par
    REAL (KIND = realkind), INTENT(IN OUT), DIMENSION(1:SM_nb_free_par) :: rho_restricted

    !variables necessary in npsol
    REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: c
    REAL (KIND=realkind) :: f, criterion
    REAL (KIND=realkind), DIMENSION(1:SM_nb_free_par) :: g

    !other variables
    REAL (KIND=realkind), DIMENSION(1:SM_nb_par) :: rho_grid_search
    INTEGER :: j
    CHARACTER (LEN=space_lenght) :: space

100 FORMAT(A,A,TR1,ES11.4," | ",25(ES11.4,TR1))

    IF (ncnln .EQ. 0) THEN
       ALLOCATE (c(1:1))
    ELSE
       ALLOCATE (c(1:ncnln))
    END IF

    space = " "
    space_lenght = space_lenght + 3
    grid_search_SM_successfull = .FALSE.

    rho_rest_SM_old = -999999.9_realkind
    j=0
    CALL ROB_CRITERION_STEP_2(j, SM_nb_free_par, rho_restricted, f, g, 1)

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=100) space, "q'Mq |  rho:", f, rho_restricted
    END IF

    !if required perform a first grid search
    IF (do_grid_search_SM  .AND. ((do_grid_search .EQ. 2) .OR. (do_grid_search.EQ. 3))) THEN
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "START GRID_SEARCH"
       END IF
       CALL DO_SM_GRID_SEARCH(rho_grid_search,"rob",criterion)
       IF (grid_search_SM_successfull) THEN
          IF (criterion .LT. f) THEN
             CALL RESTRICT_RHO(rho_grid_search,rho_restricted)
             IF (write_output_2) THEN
                WRITE(UNIT=log_unit_nb,FMT=*) space,  "GRID SEARCH SUCCESSFULL: NEW START VALUES"
                WRITE(UNIT=log_unit_nb,FMT=100) space,  "q'Mq &  rho:", criterion, &
                     rho_restricted
             END IF
          ELSE
             IF (write_output_2) THEN
                WRITE(UNIT=log_unit_nb,FMT=*) space, "GRID SEARCH NOT SUCCESSFULL"
             END IF
          END IF
       ELSE
          IF (write_output_2) THEN
             WRITE(UNIT=log_unit_nb,FMT=*) space, "GRID SEARCH NOT SUCCESSFULL"
          END IF
       END IF
    END IF

    rho_rest_SM_start = rho_restricted
    rho_rest_SM_end = rho_restricted
    rho_rest_SM_old = -999999.0_realkind

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "START NPSOL"
    END IF

    IF (use_npsol(4)) THEN
       CALL npsol(SM_nb_free_par,nclin,ncnln,ldA,ldJ,ldR,A,bl,bu,CONSTRAINTS_STEP_2, &
            ROB_CRITERION_STEP_2,inform,iter,istate,c,cJac,clamda,f,g,RR, &
            rho_rest_SM_end,iw,leniw,w,lenw )
       IF (nb_rep_monte_carlo .GT. 0) THEN
          !save the results
          IF ((step .EQ. 3) .AND. (derivative_level .EQ. 1)) THEN
             npsol_inform_step_2_rob(mc_counter) = inform
             npsol_iter_step_2_rob(mc_counter) = iter
          END IF
       END IF
    ELSE
       CALL minimum(MY_ROB_CRITERION_STEP_2,rho_rest_SM_end,f,iter)
       inform=iter
    END IF

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=100) space, "q'Mq | ^rho:", f, rho_rest_SM_end
    END IF
    
    minimum_step_2_rob_vector(mc_counter) = f
    minimum_step_2_rob = f

    IF (nb_rep_monte_carlo .GE. 1) THEN
       hansen_test_results_rob(mc_counter,1) = T1*f
       npsol_inform_step_2_rob(mc_counter) = inform
       npsol_iter_step_2_rob(mc_counter) = iter
    END IF

    rho_restricted = rho_rest_SM_end

    space_lenght = space_lenght - 3

    DEALLOCATE(c)

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "EXIT STEP_2_ROB"
    END IF
    RETURN

  END SUBROUTINE STEP_2_ROB


  SUBROUTINE UPDATES(beta,rho,update_kind,nb_loop)
    USE S1_COMMON,   ONLY: N1, AM_nb_par, AM_nb_not_fixed_par, AM_nb_free_par, &
         write_output_3, sel_score_cov_mat_s1, lags_score_cov_mat_s1
    USE S2_COMMON,   ONLY: T2, SM_nb_par, SM_nb_free_par, M_rho, &
         nb_discard_step_2, nb_iid_series, random_iid_step_2, simulated_series_step_2, &
         compute_accuracy_M_rho
    USE STARTUP,     ONLY: NPSOL_SETUP
    USE Z_MATRIX,    ONLY: is_new3, beta_old
    USE UTILITIES,   ONLY: RESTRICT_BETA, RESTRICT_RHO
    USE ROB_COMMON,  ONLY: rob_A, rob_A_tilde, B_0, S, psi_N1, score_N1, &
         rob_weights_N1, nb_A_updates, nb_B_updates, nb_M_rho_updates, wait_update_M_rho, &
         count_update_M_rho
    USE MYFUNCTIONS, ONLY: MM, CHOL, INV, MAT_CROSS_LAG
    USE SIMULATION,  ONLY: STRUCTURAL_MODEL_SIMULATION
    USE ROBUST_UTILITIES
    USE M_RHO_MATRIX, ONLY: APPROX_M_RHO, WRITE_FUNCTION_M
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_par) :: beta
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:SM_nb_par) :: rho
    INTEGER, INTENT(IN) :: update_kind, nb_loop


    INTEGER :: i, j, k
    REAL (KIND=realkind), DIMENSION(1:AM_nb_free_par) :: beta_restricted
    REAL (KIND=realkind), DIMENSION(1:SM_nb_free_par) :: rho_restricted
    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par,1:AM_nb_not_fixed_par) :: SM
    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par,1:AM_nb_not_fixed_par) :: gamm
    REAL (KIND=realkind), DIMENSION(1:SM_nb_free_par, 1:SM_nb_free_par) :: A_tilde_check

100 FORMAT(I2,TR1,25(ES10.3,TR1))

    !remove the restricted components from the vector of estimates of beta and rho
    CALL RESTRICT_BETA(beta,beta_restricted)
    CALL RESTRICT_RHO(rho,rho_restricted)

    CALL UPDATE_SCORE_N1(beta_restricted)

    IF (update_kind .EQ. 0) THEN

       IF (write_output_3) THEN
          WRITE(UNIT=37,FMT='(A,I4)') "Loop: ", nb_loop
       END IF
       !update the A matrix
       CALL UPDATE_A(beta_restricted,0) !attribute 1 refers to do_update_score_N1

       IF (write_output_3) THEN

          psi_N1 = MM(score_N1,rob_A,N1,AM_nb_not_fixed_par,AM_nb_not_fixed_par)

          CALL WEIGHTS_HUBER(psi_N1, N1, AM_nb_not_fixed_par, rob_weights_N1)
          DO j=1,AM_nb_not_fixed_par
             DO k=1,N1
                psi_N1(k,j)=psi_N1(k,j)*rob_weights_N1(k,1)
             END DO
          END DO

          !compute the covariance matrix of the orth. conditions
          SM=MM(TRANSPOSE(psi_N1),psi_N1,AM_nb_not_fixed_par,N1,AM_nb_not_fixed_par)

          IF (sel_score_cov_mat_s1 .EQ. 2) THEN
             k=1
             DO
                IF (k .GT. lags_score_cov_mat_s1) EXIT
                gamm = MAT_CROSS_LAG(psi_N1,N1,AM_nb_not_fixed_par,k)
                SM = SM+(1.0_realkind-REAL(k)/REAL(lags_score_cov_mat_s1+1))*(gamm+TRANSPOSE(gamm))
                k=k+1
             END DO
          END IF

          SM=SM/REAL(N1)
          WRITE(UNIT=40,FMT='(A,I4)') "Check Var(A*psi) = I at loop: ", nb_loop 
          DO j=1, AM_nb_not_fixed_par
             WRITE(UNIT=40,fmt=100) j, SM(j,:)
          END DO
          WRITE(UNIT=40,fmt=*)

       END IF

    ELSE

       IF (write_output_3) THEN
          WRITE(UNIT=39,fmt='(A,I4)') "Loop: ", nb_loop
       END IF

       !update matrix B:
       CALL UPDATE_B(beta_restricted, B_0, 0) !attribute 1 refers to do_update_psi_matrix_tmp_N

       !update S
       !Qui andrebbe fatto l'update della matrice S. Noi utilizziamo la matrice
       !S ottimale, cioe' S=B^(-1)
       S = INV(B_0,AM_nb_not_fixed_par)

       IF (write_output_3) THEN
          WRITE(UNIT=38,fmt='(A,I4)') "Loop: ", nb_loop
       END IF
       IF (nb_M_rho_updates .NE. -1) THEN
          IF (count_update_M_rho .EQ. wait_update_M_rho) THEN
             DO i=1, nb_M_rho_updates
                !update matrix M_rho:
                CALL APPROX_M_RHO(beta, rho, M_rho, 1)
                IF (write_output_3) THEN
                   DO j=1, AM_nb_not_fixed_par
                      WRITE(UNIT=38,fmt=100) j, M_rho(j,:)
                   END DO
                   WRITE(UNIT=38,fmt=*)
                END IF
             END DO
             count_update_M_rho = 1
          ELSE
             IF (write_output_3) THEN
                WRITE(UNIT=38,fmt='(A)') "No update!"
                count_update_M_rho = count_update_M_rho + 1
             END IF
          END IF
       ELSE
          IF (nb_loop .LE. 1) THEN
             !update matrix M_rho:
             CALL APPROX_M_RHO(beta, rho, M_rho, 1)
             IF (write_output_3) THEN
                DO j=1, AM_nb_not_fixed_par
                   WRITE(UNIT=38,fmt=100) j, M_rho(j,:)
                END DO
                WRITE(UNIT=38,fmt=*)
             END IF
          ELSE
             IF (write_output_3) THEN
                WRITE(UNIT=38,fmt='(A)') "No update!"
             END IF
          END IF
       END IF

       IF (compute_accuracy_M_rho) THEN
          CALL WRITE_FUNCTION_M(beta, rho, 1)
       END IF

       IF (write_output_3) THEN
          WRITE(UNIT=37,fmt='(A,I4)') "Loop: ", nb_loop
       END IF
       !update A_tilde
       CALL UPDATE_A(beta_restricted,0)

       IF (write_output_3) THEN

          !update matrix M_rho:
          SM = MM(S, B_0,  AM_nb_not_fixed_par, AM_nb_not_fixed_par, AM_nb_not_fixed_par)
          SM = MM(SM, S, AM_nb_not_fixed_par, AM_nb_not_fixed_par, AM_nb_not_fixed_par)
          A_tilde_check = MM(TRANSPOSE(M_rho), &
               MM(SM, M_rho, AM_nb_not_fixed_par, AM_nb_not_fixed_par, SM_nb_free_par), &
               SM_nb_free_par, AM_nb_not_fixed_par, SM_nb_free_par)

          WRITE(UNIT=40,fmt='(A,I4,A)') "Loop ", nb_loop, " check M_rho'*S*B*S*M_rho:"
          DO j=1, SM_nb_free_par
             WRITE(UNIT=40,fmt=100) j, A_tilde_check(j,:)
          END DO

          A_tilde_check = INV(MM(rob_A_tilde, TRANSPOSE(rob_A_tilde), &
               SM_nb_free_par, SM_nb_free_par, SM_nb_free_par), SM_nb_free_par)
          WRITE(UNIT=40,fmt='(A)') "Check (A*A')^-1:"
          DO j=1, SM_nb_free_par
             WRITE(UNIT=40,fmt=100) j, A_tilde_check(j,:)
          END DO
          WRITE(UNIT=40,fmt=*)

       END IF

    END IF

    RETURN
  END SUBROUTINE UPDATES


  SUBROUTINE ROB_MAIN(beta_start, rho_start)
    USE DATATYPES
    USE S1_COMMON,   ONLY: T1, AM_nb_par, AM_nb_not_fixed_par, AM_nb_free_par, beta_AM_end, &
         grad_i, modified_series, space_lenght, write_output_2, write_output_3, &
         log_unit_nb, step, derivative_level, do_grid_search_first_n_loops
    USE S2_COMMON,   ONLY: T2, SM_nb_par, SM_nb_free_par, rho_SM_end, M_rho, &
         nb_discard_step_2, nb_iid_series, random_iid_step_2, simulated_series_step_2
    USE STARTUP,     ONLY: NPSOL_SETUP
    USE UTILITIES
    USE ROB_COMMON
    USE SIMULATION,  ONLY: STRUCTURAL_MODEL_SIMULATION
    USE MYFUNCTIONS, ONLY: EYE, INV
    USE MONTE_CARLO, ONLY: nb_rep_monte_carlo, mc_counter, mc_results_step_1_rob, mc_results_step_2_rob
    USE M_RHO_MATRIX, ONLY: APPROX_M_RHO
    USE GRID_SEARCH_AM, ONLY: do_grid_search_AM
    USE GRID_SEARCH_SM, ONLY: do_grid_search_SM
    USE ROBUST_UTILITIES, ONLY: UPDATE_A, UPDATE_B
    USE LIKELIHOOD_UTILITIES, ONLY: SNP_ORT_FUNC, LIKELIHOOD_UPDATE
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_par) :: beta_start
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:SM_nb_par) :: rho_start

    INTEGER :: nb_loops, i
    REAL (KIND=realkind) :: bound_tmp
    REAL (KIND=realkind), DIMENSION(1:AM_nb_free_par) :: beta_restricted_start
    REAL (KIND=realkind), DIMENSION(1:AM_nb_free_par) :: beta_restricted_old
    REAL (KIND=realkind), DIMENSION(1:AM_nb_free_par) :: beta_restricted
    REAL (KIND=realkind), DIMENSION(1:SM_nb_free_par) :: rho_restricted_start
    REAL (KIND=realkind), DIMENSION(1:SM_nb_free_par) :: rho_restricted_old
    REAL (KIND=realkind), DIMENSION(1:SM_nb_free_par) :: rho_restricted
    REAL (KIND=realkind), DIMENSION(1:AM_nb_free_par) :: max1
    REAL (KIND=realkind), DIMENSION(1:SM_nb_free_par) :: max2
    CHARACTER (LEN=space_lenght) :: space
    LOGICAL :: do_grid_search_AM_tmp
    LOGICAL :: do_grid_search_SM_tmp


100 FORMAT(A5,TR1,ES11.4," | ",25(ES11.4,TR1))
101 FORMAT(I2,TR1,25(ES10.3,TR1))

    space = " "
    space_lenght = space_lenght + 3

    IF (write_output_3) THEN
       OPEN (UNIT=33,FILE="out_theta.dat",FORM="FORMATTED")
       OPEN (UNIT=34,FILE="out_rho.dat",FORM="FORMATTED")
       OPEN (UNIT=35,FILE="out_dtheta.dat",FORM="FORMATTED")
       OPEN (UNIT=36,FILE="out_drho.dat",FORM="FORMATTED")
       OPEN (UNIT=37,FILE="out_A_rob.dat",FORM="FORMATTED")
       OPEN (UNIT=38,FILE="out_M_rob.dat",FORM="FORMATTED")
       OPEN (UNIT=39,FILE="out_B_rob.dat",FORM="FORMATTED")
       OPEN (UNIT=40,FILE="out_check_A.dat",FORM="FORMATTED")
       OPEN (UNIT=41,FILE="out_minimum.dat",FORM="FORMATTED")     
    END IF

    !remove the restricted components from the vector of estimates of beta and rho
    CALL RESTRICT_BETA(beta_start,beta_restricted_start)
    CALL RESTRICT_RHO(rho_start,rho_restricted_start)
    nb_loops = -1  !initialisation
    IF (weights_metric .EQ. 0) THEN
       !perform a first update of A, alpha and the weights. We do not use 
       !a bound = infinity because the matrix A is used for both the computation of
       !the weights and the quadratic form psi'*A*psi. It is then important to
       !dowload the effect of outliers in estimation of A. 
       rob_A = EYE(AM_nb_not_fixed_par)
       DO i=1,nb_A_updates_start
          IF (rob_start_value(1)) THEN
             CALL UPDATES(beta_start, rho_start,weights_metric,0)
          END IF
       END DO

    ELSE
       !calcola la derivata di M_rho con un c=inf
       M_rho = 0._realkind
       rob_A_tilde = EYE(SM_nb_free_par)
       B_0 = EYE(AM_nb_not_fixed_par)
       S = EYE(AM_nb_not_fixed_par)
       !In the first and second steps we do not need a weighting matrix for the quadratic
       !form so that we can use a bound of infinity in order to calculate the
       !matrices used as first approximation for the computation of the weights.
       bound_tmp = bound
       bound = 1000000._realkind
       !the following UPDATE_B and APPROX_M_RHO are necessary to have a valid M_rho starting matrix
       !update matrix M_rho:
       IF (write_output_3) THEN
          WRITE(UNIT=38,fmt='(A)') "ROB_MAIN: Initialising M matrix for robust estimation"
       END IF
       CALL APPROX_M_RHO(beta_start, rho_start, M_rho, 1)
       count_update_M_rho = wait_update_M_rho
       IF (write_output_3) THEN
          DO i=1, AM_nb_not_fixed_par
             WRITE(UNIT=38,fmt=101) i, M_rho(i,:)
          END DO
          WRITE(UNIT=38,fmt=*) 
       END IF

       IF (write_output_3) THEN
          WRITE(UNIT=39,fmt='(A)') "ROB_MAIN: Initialising B matrix for robust estimation" 
       END IF
       CALL UPDATE_B(beta_restricted_start, B_0,1)
       S = INV(B_0,AM_nb_not_fixed_par)

       IF (write_output_3) THEN
          WRITE(UNIT=37,fmt='(A)') "ROB_MAIN: Initialising A matrix for robust estimation" 
       END IF
       CALL UPDATE_A(beta_restricted_start,0)

       bound = bound_tmp
        CALL UPDATES(beta_start, rho_start,weights_metric,0)
    END IF

    !set the starting values for the parameters
    beta_restricted = beta_restricted_start
    beta_restricted_old = beta_restricted_start
    rho_restricted = rho_restricted_start
    rho_restricted_old = rho_restricted_start

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "START LOOP OF ROBUST OPTIMIZATIONS"
    END IF
    IF (write_output_3) THEN
       WRITE(UNIT=41, FMT='(A,I9)') "Starting Monte Carlo nb. ", mc_counter
    END IF

    do_grid_search_AM_tmp = do_grid_search_AM
    do_grid_search_SM_tmp = do_grid_search_SM
    
    nb_loops = 1
    DO
       IF (nb_rep_monte_carlo .GT. 0) THEN
          iter_rob(mc_counter) = nb_loops
       END IF
       IF (nb_loops .GT. rob_iteration_limit) THEN
          WRITE (*,FMT='(A)') "The maximum number of iterations in ROB_GMM has been reached."
          IF (nb_rep_monte_carlo .GT. 0) THEN
             inform_rob(mc_counter) = 1
          END IF
          IF (write_output_2) THEN
             WRITE(UNIT=log_unit_nb,FMT=*) space, "The maximum number of iterations in ROB_GMM has been reached."
          END IF
          EXIT
       END IF

       !set the parameters of the NPSOL subroutine
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "NPSOL_SETUP(3)"
       END IF
       CALL NPSOL_SETUP(3)

       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "START STEP_1_ROB"
       END IF
       CALL STEP_1_ROB(beta_restricted,AM_nb_free_par)
       
       IF ((do_grid_search_first_n_loops .NE. 0) .AND. (nb_loops .GE. do_grid_search_first_n_loops)) THEN
         do_grid_search_AM = .FALSE.
       END IF 
         
       !add the restricted components to the parameter vector beta_AM_end
       CALL COMPLETE_BETA_RESTRICTED(beta_restricted,beta_AM_end)
       write(*,FMT=100) "the: ", minimum_step_1_rob, beta_restricted

       !update matrix A but only for the first method
       IF (weights_metric .EQ. 0) THEN
         CALL UPDATES(beta_AM_end,rho_SM_end,weights_metric,nb_loops)
       END IF
       
       IF (weights_metric .EQ. 1) THEN
          !set the parameters of the NPSOL subroutine
          IF (write_output_2) THEN
             WRITE(UNIT=log_unit_nb,FMT=*) space, "NPSOL_SETUP(4)"
          END IF
          CALL NPSOL_SETUP(4)

          !perform the second step
          IF (write_output_2) THEN
             WRITE(UNIT=log_unit_nb,FMT=*) space, "START STEP_2_ROB"
          END IF
          CALL STEP_2_ROB(rho_restricted,SM_nb_free_par)

          IF ((do_grid_search_first_n_loops .NE. 0) .AND. (nb_loops .GE. do_grid_search_first_n_loops)) THEN
             do_grid_search_SM = .FALSE.
          END IF

          !add the restricted components to the parameter vector rho_SM_end
          CALL COMPLETE_RHO_RESTRICTED(rho_restricted,rho_SM_end)
          write(*,FMT=100) "rho: ", minimum_step_2_rob, rho_restricted
          write(*,*) 
       END IF
       
       !check if exit condition 
       max1 = ABS(beta_restricted-beta_restricted_old)
       max2 = ABS(rho_restricted-rho_restricted_old)

       max1 = ABS(beta_restricted-beta_restricted_old)/(1._realkind+ABS(beta_restricted_old))
       max2 = ABS(rho_restricted-rho_restricted_old)/(1._realkind+ABS(rho_restricted_old))
       IF (write_output_3) THEN
          WRITE(UNIT=33,FMT=*) nb_loops, beta_AM_end
          WRITE(UNIT=34,FMT=*) nb_loops, rho_SM_end
          WRITE(UNIT=35,FMT=*) nb_loops, max1
          WRITE(UNIT=36,FMT=*) nb_loops, max2
          WRITE(UNIT=41,FMT='(I2,TR1,E9.3,TR1,E9.3)') nb_loops, minimum_step_1_rob, minimum_step_2_rob
       END IF
       IF ((MAXVAL(max1) .LT. rob_parameter_tolerance) .AND. &
            (MAXVAL(max2) .LT. rob_parameter_tolerance)) THEN
          IF (nb_rep_monte_carlo .GT. 0) THEN
             inform_rob(mc_counter) = 0
          END IF
          EXIT
       END IF

       beta_restricted_old = beta_restricted
       rho_restricted_old = rho_restricted

       IF (weights_metric .EQ. 1) THEN
         CALL UPDATES(beta_AM_end,rho_SM_end,weights_metric,nb_loops)
       END IF

       nb_loops = nb_loops + 1
    END DO


    IF (weights_metric .EQ. 0) THEN
       !set the nb_loops = 1 in order for the grid search to be effective
       nb_loops = 1
       !set the parameters of the NPSOL subroutine
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "NPSOL_SETUP(4)"
       END IF
       CALL NPSOL_SETUP(4)

       !perform the second step
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "START STEP_2_ROB"
       END IF
       CALL STEP_2_ROB(rho_restricted,SM_nb_free_par)

       !add the restricted components to the parameter vector rho_SM_end
       CALL COMPLETE_RHO_RESTRICTED(rho_restricted,rho_SM_end)
       write(*,FMT=100) "rho: ", minimum_step_2_rob, rho_restricted
       write(*,*) 
       rho_restricted_old = rho_restricted
    END IF
    
    do_grid_search_AM = do_grid_search_AM_tmp
    do_grid_search_SM = do_grid_search_SM_tmp
    
    IF (write_output_3) THEN
       CLOSE (UNIT=33) !out_theta.dat
       CLOSE (UNIT=34) !out_rho.dat
       CLOSE (UNIT=35) !out_dtheta.dat
       CLOSE (UNIT=36) !out_drho.dat
       CLOSE (UNIT=37) !out_A_rob.dat
       CLOSE (UNIT=38) !out_M_rob.dat
       CLOSE (UNIT=39) !out_B_rob.dat
       CLOSE (UNIT=40) !out_check_A.dat
       CLOSE (UNIT=41) !out_minimium.dat      
    END IF

    IF (nb_rep_monte_carlo .GT. 0) THEN
       !save the results
       mc_results_step_1_rob(:,mc_counter)=beta_AM_end(:)
       IF ((step .EQ. 3) .AND. (derivative_level .EQ. 1)) THEN
          mc_results_step_2_rob(:,mc_counter) = rho_SM_end(:)
       END IF
    END IF

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "EXIT ROB_MAIN"
    END IF

    space_lenght = space_lenght - 3

    RETURN
  END SUBROUTINE ROB_MAIN

END MODULE STEPS

