
SUBROUTINE SNP()
  USE DATATYPES
  USE S1_COMMON,   ONLY: T1, step, beta_AM_start, beta_AM_end, modified_series, &
       derivative_level, inv_score_cov, score_cov, AM_nb_not_fixed_par, write_output_1, &
       write_output_2, log_unit_nb, space_lenght, expectation_analysis

  USE S2_COMMON,   ONLY: change_random_iid
  USE MONTE_CARLO, ONLY: mc_results_step_1, mc_results_step_2, nb_rep_monte_carlo, &
       mc_counter, inv_score_cov_mc, estimate_all_mc_simulations, &
       selected_simulations, skip_first_n
  USE STEPS,       ONLY: STEP_1, STEP_2, DO_AM_GRID_SEARCH
  USE STARTUP,     ONLY: IMPORT_SERIES, MODIFY_SERIES, NPSOL_SETUP
  USE Z_MATRIX,    ONLY: beta_old, is_new1, is_new2, series_old
  USE UTILITIES,   ONLY: RESTRICT_BETA, UPDATE_RANDOM_IID
  USE S2_COMMON,   ONLY: rho_SM_end, rho_rest_SM_old
  ! tolto perché non utilizzato USE ROB_COMMON,  ONLY: random_iid_rob_update
  USE MYFUNCTIONS, ONLY: INV
  USE OUTPUT_UTILITIES,     ONLY: OUTPUT_SNP_MC
  USE COVARIANCE_ESTIMATES, ONLY: EST_COV_SCORE_AM_STEP1
  IMPLICIT NONE

  INTEGER :: sub_i, j
  REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par) :: beta_restricted
  CHARACTER (LEN=space_lenght) :: space
  LOGICAL :: ok

  space = " "
  space_lenght = space_lenght + 3

  IF ((step .EQ. 1) .OR. (step .EQ. 3)) THEN
     !start a loop over all monte carlo repetitions but at least once
     j=MAX(1,nb_rep_monte_carlo)
     sub_i=1
     DO mc_counter = 1, j
        !first import/generate the "observed" time series
        IF (write_output_2) THEN
           WRITE(UNIT=log_unit_nb,FMT=*) "MC nb.", mc_counter
           WRITE(UNIT=log_unit_nb,FMT=*) space, "IMPORT_SERIES"
        END IF
        CALL IMPORT_SERIES(10)

        IF (write_output_2) THEN
           WRITE(UNIT=log_unit_nb,FMT=*) space, "MODIFY_SERIES"
        END IF
        CALL MODIFY_SERIES()

        IF (estimate_all_mc_simulations .EQ. 1) THEN
           IF (skip_first_n .LT. mc_counter) THEN
              ok=.true.
           ELSE
              ok=.false.
           END IF
        ELSE
           IF (mc_counter .EQ. selected_simulations(sub_i)) THEN
              sub_i = sub_i + 1
              ok=.true.
           ELSE
              ok=.false.
           END IF
        END IF

        IF (ok) THEN
           series_old = 0  
           is_new1 = 1
           is_new2 = 1
           beta_old = -999999.0_realkind

           !Now set the parameters of the NPSOL subroutine
           IF (write_output_2) THEN
              WRITE(UNIT=log_unit_nb,FMT=*) space, "NPSOL_SETUP(1)"
           END IF
           CALL NPSOL_SETUP(1)

           !perform the first step
           IF (write_output_2) THEN
              WRITE(UNIT=log_unit_nb,FMT=*) space, "START STEP_1"
           END IF
           CALL STEP_1()

           !if step=3 perform the second step
           IF (step .EQ. 3) THEN
              IF (derivative_level .EQ. 1) THEN
                 beta_old = -999999.9_realkind

                 !estimate the covariance matrix of the score used in the quadratic form defining the
                 !objective function in step 2
                 IF (write_output_2) THEN
                    WRITE(UNIT=log_unit_nb,FMT=*) space, "EST_COV_SCORE_AM_STEP1"
                 END IF
                 CALL RESTRICT_BETA(beta_AM_end,beta_restricted)
                 CALL EST_COV_SCORE_AM_STEP1(beta_restricted,T1,modified_series,1)
                 !the "1" is which_series
                 inv_score_cov = INV(score_cov,AM_nb_not_fixed_par)
                 IF (expectation_analysis) THEN
                    inv_score_cov_mc(:,:,mc_counter) = inv_score_cov(:,:)
                 END IF
                 rho_rest_SM_old = -999999.9_realkind  !necessary only from the second repetition on

                 IF (write_output_2) THEN
                    WRITE(UNIT=log_unit_nb,FMT=*) space, "NPSOL_SETUP(2)"
                 END IF
                 CALL NPSOL_SETUP(2)

                 IF (write_output_2) THEN
                    WRITE(UNIT=log_unit_nb,FMT=*) space, "START STEP_2"
                 END IF
                 CALL STEP_2()
              ELSE
                 WRITE(*,*) "Derivative level .NE. 1: step 2 skipped"
              END IF
           END IF
           !if nb_rep_monte_carlo > 0 save the results

           IF (nb_rep_monte_carlo .GT. 0) THEN
              !save the results
              mc_results_step_1(:,mc_counter)=beta_AM_end(:)
              IF ((step .EQ. 3) .AND. (derivative_level .EQ. 1)) THEN
                 mc_results_step_2(:,mc_counter,1) = rho_SM_end(:)
              END IF
           END IF
           write(*,*) "Fine step1-montecarlo nr. ", mc_counter
        END IF
        !generate a new series of iid innovations used to approximate the integrals
        IF (change_random_iid) THEN
           CALL UPDATE_RANDOM_IID()
        END IF
     END DO
  ELSE    !IF ((step .EQ. 1) .OR. (step .EQ. 3)) 
     IF (derivative_level .EQ. 1) THEN
        series_old = 0  
        is_new2 = 1
        beta_old = -99999.0_realkind
        beta_AM_end = beta_AM_start

        !estimate the covariance matrix of the score
        IF (write_output_2) THEN
           WRITE(UNIT=log_unit_nb,FMT=*) space, "EST_COV_SCORE_AM_STEP1"
        END IF

        CALL RESTRICT_BETA(beta_AM_end,beta_restricted)
        CALL EST_COV_SCORE_AM_STEP1(beta_restricted,T1,modified_series,1)
        !the "1" is which_series
        inv_score_cov = INV(score_cov,AM_nb_not_fixed_par)

        rho_rest_SM_old = -999999.9_realkind
        !setup the parameters for the NPSOL subroutine
        IF (write_output_2) THEN
           WRITE(UNIT=log_unit_nb,FMT=*) space, "NPSOL_SETUP(2)"
        END IF
        CALL NPSOL_SETUP(2)

        IF (write_output_2) THEN
           WRITE(UNIT=log_unit_nb,FMT=*) space, "START STEP_2"
        END IF
        CALL STEP_2()
     ELSE
        WRITE(*,*) "Derivative level .NE. 1: impossible to perform step 2"
        STOP
     END IF
  END IF

  IF ((nb_rep_monte_carlo .GT. 0) .AND. (write_output_1)) THEN
     CALL OUTPUT_SNP_MC()
  END IF

  space_lenght = space_lenght - 3

  RETURN

END SUBROUTINE SNP

