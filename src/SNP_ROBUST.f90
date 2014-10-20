SUBROUTINE SNP_ROBUST()
  USE DATATYPES
  USE S1_COMMON,   ONLY: T1, step, AM_nb_par, beta_AM_start, beta_AM_end, modified_series, &
       derivative_level, inv_score_cov, score_cov, AM_nb_not_fixed_par, &
       write_output_1, write_output_2, log_unit_nb, space_lenght

  USE S2_COMMON,   ONLY: rho_SM_end, rho_rest_SM_old, rho_SM_start, SM_nb_par, &
       change_random_iid
  USE STEPS,       ONLY: STEP_1, STEP_2, ROB_MAIN
  USE STARTUP,     ONLY: IMPORT_SERIES, MODIFY_SERIES, NPSOL_SETUP, GRID_POINTS_AM, GRID_POINTS_SM
  USE Z_MATRIX,    ONLY: beta_old, series_old, is_new1, is_new2, is_new3
  USE UTILITIES,   ONLY: RESTRICT_BETA, UPDATE_RANDOM_IID
  USE MYFUNCTIONS, ONLY: SNP_DENSITY_FUNC, INV, IID_SIM
  USE OUTPUT_UTILITIES,     ONLY: OUTPUT_SNP_ROB_MC, OUTPUT_SNP_ROB_STEP_1, &
       OUTPUT_SNP_ROB_STEP_2
  USE COVARIANCE_ESTIMATES, ONLY: EST_COV_SCORE_AM_STEP1
  USE GRID_SEARCH_AM,   ONLY: do_grid_search_AM
  USE GRID_SEARCH_SM,   ONLY: do_grid_search_SM
  USE ROB_COMMON
  USE MONTE_CARLO
  IMPLICIT NONE

  INTEGER :: sub_i,j
  REAL (KIND=realkind), DIMENSION(1:AM_nb_par) :: beta
  REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par) :: beta_restricted
  REAL (KIND=realkind), DIMENSION(1:SM_nb_par) :: rho
  REAL (KIND=realkind) :: dummy_real
  CHARACTER (LEN=space_lenght) :: space
  LOGICAL :: ok

  !100 FORMAT(A5,25(ES10.3,TR1))
100 FORMAT(A5,TR1,ES11.4," | ",25(ES11.4,TR1))
101 FORMAT(A,A,TR1,25(ES10.3,TR1))

  IF (step .NE. 3) THEN
     IF (step .EQ. 1) THEN
        WRITE(*,*) "The robust estimation of the first step only is not possible with this methodology!"
        STOP
     ELSE
        WRITE(*,*) "The robust estimation of the second step only has not been implemeted yet!"
        STOP
     END IF
  END IF

  IF (derivative_level .NE. 1) THEN
     WRITE(*,*) "derivative_level must be 1 in order to do a robust estimation"
     STOP
  END IF

  space = " "
  space_lenght = space_lenght + 3
  beta = beta_AM_start
  rho = rho_SM_start

  !start a loop over all monte carlo repetitions but at least once
  j=MAX(1,nb_rep_monte_carlo)
  sub_i=1
  DO mc_counter=1, j
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

     series_old = 0
     is_new1 = 1
     is_new2 = 1
     is_new3 = 1
     beta_old = -999999.0_realkind

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

     IF (ok) THEN ! esegui la stima
        write(*,fmt='(A,I6)') "starting values and classical estimates ", mc_counter
        !if it's required, perform a first non robust estimation and use the results
        !as starting values
        IF (rob_start_value(1)) THEN
           IF (write_output_2) THEN
              WRITE(UNIT=log_unit_nb,FMT=*) space, "CLASSICAL ESTIMATION OF THE AUXILIARY MODEL"
           END IF
           !Now set the parameters of the NPSOL subroutine
           IF (write_output_2) THEN
              WRITE(UNIT=log_unit_nb,FMT=*) space, "NPSOL_SETUP(1)"
           END IF
           CALL NPSOL_SETUP(1)

           !perform the first step
           IF (write_output_2) THEN
              WRITE(UNIT=log_unit_nb,FMT=*) space, "START STEP_1()"
           END IF
           CALL STEP_1()
           !set the results as the new starting values
           beta = beta_AM_end
           !save the results of the non robust estimation
           IF (nb_rep_monte_carlo .GT. 0) THEN
              mc_results_step_1(:,mc_counter) = beta_AM_end(:)
           END IF
           IF (do_grid_search_AM) THEN
              CALL GRID_POINTS_AM(AM_nb_par,beta)
           END IF
        ELSE
           dummy_real=-99999.999_realkind
           write(*,fmt=100) "the: ", dummy_real, beta_AM_start
        END IF

        IF (rob_start_value(2)) THEN
           !now the second step
           IF (write_output_2) THEN
              WRITE(UNIT=log_unit_nb,FMT=*) space, "CLASSICAL ESTIMATION OF THE STRUCTURAL MODEL"
           END IF
           IF (.NOT. rob_start_value(1)) THEN
                beta_AM_end = beta_AM_start
           END IF
           !estimate the covariance matrix of the score used in the quadratic form defining the
           !objective function in step 2
           IF (write_output_2) THEN
              WRITE(UNIT=log_unit_nb,FMT=*) space, "EST_COV_SCORE_AM_STEP1()"
           END IF

           CALL RESTRICT_BETA(beta_AM_end,beta_restricted)

           beta_old = -999999.9_realkind
           CALL EST_COV_SCORE_AM_STEP1(beta_restricted,T1,modified_series,1)

           !the "1" is "which_series"
           inv_score_cov = INV(score_cov,AM_nb_not_fixed_par)

           rho_rest_SM_old = -999999.9_realkind  !necessary only from the second repetition on

           !Now set the parameters of the NPSOL subroutine
           IF (write_output_2) THEN
              WRITE(UNIT=log_unit_nb,FMT=*) space, "NPSOL_SETUP(2)"
           END IF
           CALL NPSOL_SETUP(2)

           IF (write_output_2) THEN
              WRITE(UNIT=log_unit_nb,FMT=*) space, "START STEP_2"
           END IF
           CALL STEP_2()
           !set the results as the new starting values
           rho = rho_SM_end
           !save the results of the non robust estimation
           IF (nb_rep_monte_carlo .GT. 0) THEN
              mc_results_step_2(:,mc_counter,1) = rho_SM_end(:)
           END IF
           IF (do_grid_search_SM) THEN
              CALL GRID_POINTS_SM(SM_nb_par,rho)
           END IF
        ELSE
           dummy_real=-99999.999_realkind
           write(*,fmt=100) "rho: ", dummy_real, rho_SM_start
        END IF

        IF (write_output_2) THEN
           WRITE(UNIT=log_unit_nb,FMT=*) space, "START ROB_MAIN"
           WRITE(UNIT=log_unit_nb,FMT=101) space, "beta start values: ", beta
           WRITE(UNIT=log_unit_nb,FMT=101) space, "rho  start values: ", rho
        END IF

        write(*,*)
        write(*,fmt='(A)') "start the robust estimation subroutine:"

        !beta and rho are "intent(in)". The results are stored in beta_AM_end, rho_SM_end
        CALL ROB_MAIN(beta, rho)

        IF ((nb_rep_monte_carlo .LT. 2) .AND. (write_output_1)) THEN
           CALL OUTPUT_SNP_ROB_STEP_1()
           CALL OUTPUT_SNP_ROB_STEP_2()
        END IF

        IF (nb_rep_monte_carlo .GT. 0) THEN
           mc_results_step_1_rob(:,mc_counter) = beta_AM_end(:)
           mc_results_step_2_rob(:,mc_counter) = rho_SM_end(:)
        END IF

        beta = beta_AM_start
        rho = rho_SM_start
     END IF
     !generate a new series of iid innovations used to approximate the integrals
     IF (change_random_iid) THEN
        CALL UPDATE_RANDOM_IID()
     END IF
  END DO

  IF ((nb_rep_monte_carlo .GT. 0) .AND. (write_output_1)) THEN
     CALL OUTPUT_SNP_ROB_MC()
  END IF

  space_lenght = space_lenght - 3

  IF (write_output_2) THEN
     WRITE(UNIT=log_unit_nb,FMT=*) space, "EXIT SNP_ROBUST"
  END IF

  RETURN

END SUBROUTINE SNP_ROBUST

