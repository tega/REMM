MODULE SIMULATION
  USE DATATYPES
  IMPLICIT NONE
CONTAINS

  SUBROUTINE DGP_SIMULATION(T, time_series, unit_nb)
    USE DGP_DEFINITION
    USE CONTAMINATION_MODEL
    USE S1_COMMON,   ONLY: write_output_2
    USE MONTE_CARLO, ONLY: nb_rep_monte_carlo
    USE MYFUNCTIONS, ONLY: IID_SIM, CONTAMINATE, ARMA_SIMULATION, SVM1_SIMULATION, &
         ARMA_GARCH, AR_GARCH_SNP_LAG, SWITCH_AR
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: T
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:T) :: time_series
    INTEGER, INTENT(IN) :: unit_nb

    INTEGER :: i, is_contaminated, is_simulated
    ! remove this line, start, dim_mu, dim_arch, dim_garch, dim_sigma
    REAL (KIND=realkind), DIMENSION(1:T+DGP_nb_discard,DGP_nb_iid_series) :: iid
    REAL (KIND=realkind), DIMENSION(1:T+DGP_nb_discard,1) :: iid_temp
    ! remove this line REAL (KIND=realkind), DIMENSION(1,1) :: iid_temp_bis
    REAL (KIND=realkind), DIMENSION(1:nb_contam_points(1)) :: additive_contam_points
    REAL (KIND=realkind), DIMENSION(1:nb_contam_points(1)) :: zero_additive_contam_points
    REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: simulated_not_observable
    CHARACTER (LEN=2) :: simulation_kind

    IF ((DGP_name .EQ. "STOC_VOLATILITY") .OR. (DGP_name .EQ. "SWITCH_AR_MODEL")) THEN
       ALLOCATE(simulated_not_observable(1:T,1))
       simulation_kind = "AR"
    END IF


    zero_additive_contam_points=0._realkind
    !verify if we want to contaminate the model
    is_simulated = 0
    is_contaminated = 0
    DO i=1, DGP_nb_iid_series + 1
       IF (contam_model(i) .EQ. 1) THEN
          is_contaminated = 1
       END IF
    END DO

    IF (DGP_pol_dim .LE. 1) THEN
       !first simulate the iid series used to construct the process
       DO i=1, DGP_nb_iid_series
          IF (DGP_innov_dist(i) .EQ. "SNP_DST") THEN
             iid_temp = IID_SIM(T+DGP_nb_discard,1,DGP_innov_dist(i),DGP_deg_freedom(i), &
                  DGP_pol_coef,DGP_nb_pol_coef,DGP_epsilon_0)
          ELSE
             iid_temp = IID_SIM(T+DGP_nb_discard,1,DGP_innov_dist(i),DGP_deg_freedom(i))
          END IF
          iid(1:T+DGP_nb_discard,i) = iid_temp(1:T+DGP_nb_discard,1)
       END DO

       IF ((nb_rep_monte_carlo .GT. 0) .AND. (write_output_2)) THEN 
          !write the generated series in a file
          OPEN (UNIT=unit_nb,FILE="series_iid.dat",FORM="FORMATTED")
          DO i=1, T+DGP_nb_discard   
             write (unit=unit_nb,fmt=*) iid(i,:) 
          END DO
          CLOSE(UNIT=unit_nb)
       END IF

       IF ((nb_rep_monte_carlo .GT. 0) .OR. is_contaminated .EQ. 0) THEN 
          !simulate the uncontaminated process
          is_simulated = 1
          IF (DGP_name .EQ. "ARMA(p,q)_MODEL") THEN
             CALL ARMA_SIMULATION(DGP_nb_par,DGP_par_vec,DGP_order_model_dim(1), DGP_order_model_dim(2), &
                  T, time_series, DGP_nb_discard, DGP_nb_iid_series, iid)
          END IF

          IF (DGP_name .EQ. "STOC_VOLATILITY") THEN
             CALL SVM1_SIMULATION(DGP_nb_par,DGP_par_vec,DGP_order_model_dim(1), DGP_order_model_dim(2), &
                  T, time_series, DGP_nb_discard, DGP_nb_iid_series, iid, simulation_kind, &
                  simulated_not_observable=simulated_not_observable)
          END IF

          IF (DGP_name .EQ. "ARMA_GARCH_MOD_") THEN
             CALL ARMA_GARCH(DGP_nb_par,DGP_par_vec,DGP_order_model_dim(1), DGP_order_model_dim(2), &
                  DGP_order_model_dim(3), DGP_order_model_dim(4), T, time_series, DGP_nb_discard, &
                  DGP_nb_iid_series, iid)
          END IF
          
          IF (DGP_name .EQ. "SWITCH_AR_MODEL") THEN
             CALL SWITCH_AR(DGP_order_model_dim(1),DGP_order_model_dim(2), DGP_nb_par, DGP_par_vec, &
                T, time_series, DGP_nb_discard, DGP_nb_iid_series, iid, &
                states_series=simulated_not_observable)
          END IF
 
          IF ((nb_rep_monte_carlo .GT. 0) .AND. (write_output_2)) THEN 
             !write the generated series in a file
             OPEN (UNIT=unit_nb,FILE="series_uncontaminated.dat",FORM="FORMATTED")
             DO i=1, T   
                write (unit=unit_nb,fmt=*) time_series(i) 
             END DO
             CLOSE(UNIT=unit_nb)
             IF (ALLOCATED(simulated_not_observable)) THEN
                !write the inobservable generated series in a file
                OPEN (UNIT=unit_nb,FILE="series_not_observable_uncontaminated.dat",FORM="FORMATTED")
                DO i=1, T   
                   WRITE (unit=unit_nb,fmt=*) simulated_not_observable(i,:) 
                END DO
                CLOSE(UNIT=unit_nb)
             END IF
          END IF
       END IF

       !Perform the true or a zero contamination. Because of the zero contamination the 
       !iid sequences used along many MC simulations are the same!

       !first for the innovation outliers 
       DO i=1, DGP_nb_iid_series
          IF (contam_model(i+1) .EQ. 1) THEN
             IF (contam_dist(i+1) .EQ. "_DIRAC_") THEN
                additive_contam_points(:)=contam_x0(1+i,:)
                CALL CONTAMINATE(T+DGP_nb_discard,iid(:,i),epsilon_p(1+i),nb_contam_points(i+1),x0=additive_contam_points)
             END IF
             IF (contam_dist(i+1) .EQ. "UNIFORM") THEN
                iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(i+1))
                iid_temp = iid_temp * (contam_upper_bound(i+1) -  contam_lower_bound(i+1)) + &
                     contam_lower_bound(i+1)
                CALL CONTAMINATE(T+DGP_nb_discard,iid(:,i),epsilon_p(1+i), &
                    nb_contam_points(i+1),random_x0=iid_temp)
             END IF
             IF (contam_dist(i+1) .EQ. "NORMAL_") THEN
                iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(i+1))
                iid_temp = iid_temp * SQRT(contam_norm_var(i+1)) + contam_norm_mean(i+1)
                CALL CONTAMINATE(T+DGP_nb_discard,iid(:,i),epsilon_p(1+i), &
                    nb_contam_points(i+1),random_x0=iid_temp)
             END IF
             IF (contam_dist(i+1) .EQ. "STUDENT") THEN
                iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(i+1),contam_df(i+1))
                iid_temp = iid_temp * SQRT(contam_student_var(i+1)) + contam_student_mean(i+1)
                CALL CONTAMINATE(T+DGP_nb_discard,iid(:,i),epsilon_p(1+i), &
                    nb_contam_points(i+1),random_x0=iid_temp)
             END IF
          ELSE !now perform an empty run for comparison with situations of no contamination
             IF (contam_dist(i+1) .EQ. "UNIFORM") THEN
                iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(i+1))
                !iid_temp = iid_temp * (contam_upper_bound(i+1) -  contam_lower_bound(i+1)) + &
                !     contam_lower_bound(i+1)
             END IF
             IF (contam_dist(i+1) .EQ. "NORMAL_") THEN
                iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(i+1))
                !iid_temp = iid_temp * SQRT(contam_norm_var(i+1)) + contam_norm_mean(i+1)
             END IF
             IF (contam_dist(i+1) .EQ. "STUDENT") THEN
                iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(i+1),contam_df(i+1))
                !iid_temp = iid_temp * SQRT(contam_student_var(i+1)) + contam_student_mean(i+1)
             END IF
             CALL CONTAMINATE(T+DGP_nb_discard,iid(:,i),0._realkind, &
                nb_contam_points(i+1),x0=zero_additive_contam_points)
          END IF
       END DO
       !verify if it is necessary to simulate the series again
       is_contaminated = 0
       DO i=1, DGP_nb_iid_series
          IF (contam_model(1+i) .EQ. 1) THEN
             is_contaminated = 1
          END IF
       END DO
       IF ((is_simulated .EQ. 0) .OR. (is_contaminated .EQ. 1)) THEN
          !simulate the contaminated process
          IF (DGP_name .EQ. "ARMA(p,q)_MODEL") THEN
             CALL ARMA_SIMULATION(DGP_nb_par,DGP_par_vec,DGP_order_model_dim(1), DGP_order_model_dim(2), &
                  T, time_series, DGP_nb_discard, DGP_nb_iid_series, iid)
          END IF

          IF (DGP_name .EQ. "STOC_VOLATILITY") THEN
             CALL SVM1_SIMULATION(DGP_nb_par,DGP_par_vec,DGP_order_model_dim(1),DGP_order_model_dim(2), &
                  T,time_series,DGP_nb_discard,DGP_nb_iid_series,iid, simulation_kind, &
                  simulated_not_observable=simulated_not_observable)
          END IF

          IF (DGP_name .EQ. "ARMA_GARCH_MOD_") THEN
             CALL ARMA_GARCH(DGP_nb_par,DGP_par_vec,DGP_order_model_dim(1), DGP_order_model_dim(2), &
                  DGP_order_model_dim(3),DGP_order_model_dim(4),T,time_series,DGP_nb_discard, &
                  DGP_nb_iid_series,iid)
          END IF
          
          IF (DGP_name .EQ. "SWITCH_AR_MODEL") THEN
                 CALL SWITCH_AR(DGP_order_model_dim(1),DGP_order_model_dim(2), DGP_nb_par, DGP_par_vec, &
                 T, time_series, DGP_nb_discard, DGP_nb_iid_series, iid, &
                 states_series=simulated_not_observable)
          END IF
       END IF

       !then for the additive outliers
       IF (contam_model(1) .EQ. 1) THEN
          is_contaminated = 1
          iid(1:DGP_nb_discard,1) = 0._realkind
          iid(DGP_nb_discard+1:T+DGP_nb_discard,1) = time_series(1:T)
          IF (contam_dist(1) .EQ. "_DIRAC_") THEN
             additive_contam_points(:)=contam_x0(1,:)
             CALL CONTAMINATE(T+DGP_nb_discard,iid(:,1),epsilon_p(1),nb_contam_points(1),x0=additive_contam_points)
          END IF
          IF (contam_dist(1) .EQ. "UNIFORM") THEN
             iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(1))
             iid_temp = iid_temp * (contam_upper_bound(1) -  contam_lower_bound(1)) + &
                  contam_lower_bound(1)
             CALL CONTAMINATE(T+DGP_nb_discard,iid(:,1),epsilon_p(1),nb_contam_points(1),random_x0=iid_temp)
          END IF
          IF (contam_dist(1) .EQ. "NORMAL_") THEN
             iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(1))
             iid_temp = iid_temp * SQRT(contam_norm_var(1)) + contam_norm_mean(1)
             CALL CONTAMINATE(T+DGP_nb_discard,iid(:,1),epsilon_p(1),nb_contam_points(1),random_x0=iid_temp)
          END IF
          IF (contam_dist(1) .EQ. "STUDENT") THEN
             iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(1),contam_df(1))
             iid_temp = iid_temp * SQRT(contam_student_var(1)) + contam_student_mean(1)
             CALL CONTAMINATE(T+DGP_nb_discard,iid(:,1),epsilon_p(1),nb_contam_points(1),random_x0=iid_temp)
          END IF
          time_series(1:T) = iid(DGP_nb_discard+1:T+DGP_nb_discard,1)
       ELSE
          IF (contam_dist(1) .EQ. "UNIFORM") THEN
             iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(1))
             !iid_temp = iid_temp * (contam_upper_bound(1) -  contam_lower_bound(1)) + &
             !     contam_lower_bound(1)
          END IF
          IF (contam_dist(1) .EQ. "NORMAL_") THEN
             iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(1))
             !iid_temp = iid_temp * SQRT(contam_norm_var(1)) + contam_norm_mean(1)
          END IF
          IF (contam_dist(1) .EQ. "STUDENT") THEN
             iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(1),contam_df(1))
             !iid_temp = iid_temp * SQRT(contam_student_var(1)) + contam_student_mean(1)
          END IF
          CALL CONTAMINATE(T+DGP_nb_discard,iid(:,1),0._realkind,nb_contam_points(1),x0=zero_additive_contam_points)
       END IF
       !save the contaminated process if necessary, i.e. if nb_rep_monte_carlo=1 and
       !is_contaminated .EQ. 1 and write_output_2 = true
       IF ((nb_rep_monte_carlo .GT. 0) .AND. (is_contaminated .EQ. 1) .AND. write_output_2) THEN 
          OPEN (UNIT=unit_nb,FILE="series_contaminated.dat",FORM="FORMATTED")
          DO i=1, T    
             write (unit=unit_nb,fmt=*) time_series(i) 
          END DO
          CLOSE(UNIT=unit_nb)
          IF (ALLOCATED(simulated_not_observable)) THEN
             !write the inobservable generated series in a file
             OPEN (UNIT=unit_nb,FILE="series_not_observable_contaminated.dat",FORM="FORMATTED")
             DO i=1, T   
                WRITE (unit=unit_nb,fmt=*) simulated_not_observable(i,:) 
             END DO
             CLOSE(UNIT=unit_nb)
          END IF
       END IF
    ELSE !(DGP_pol_dim .LE. 1)
       !This means that the AR_GARCH_ with a SNP_DST with lags has been selected
       !simulate the iid UNIFORM distributed necessary to simulate the SNP_ARCH_LAG model
       DO i=1, DGP_nb_iid_series
          iid_temp = IID_SIM(T+DGP_nb_discard,1,"UNIFORM")
          iid(1:T+DGP_nb_discard,i) = iid_temp(1:T+DGP_nb_discard,1)
       END DO
       !simulate the uncontaminated process
       is_simulated = 1

       CALL AR_GARCH_SNP_LAG(DGP_nb_par, DGP_par_vec, 1+DGP_order_model_dim(1), &
            DGP_order_model_dim(2), DGP_order_model_dim(3), &
            DGP_pol_dim, DGP_pol_deg, DGP_nb_pol_coef, DGP_pol_coef, DGP_epsilon_0, &
            T, time_series, DGP_nb_discard, &
            DGP_nb_iid_series, iid, 0)
       IF ((nb_rep_monte_carlo .GT. 0) .AND. write_output_2) THEN 
          !write the generated series in a file
          OPEN (UNIT=unit_nb,FILE="series_uncontaminated.dat",FORM="FORMATTED") 
          DO i=1, T   
             write (unit=unit_nb,fmt=*) time_series(i) 
          END DO
          CLOSE(UNIT=unit_nb)
       END IF
       !verify if it is necessary to simulate the series again
       is_contaminated = 0
       DO i=1, DGP_nb_iid_series
          IF (contam_model(1+i) .EQ. 1) THEN
             is_contaminated = 1
          END IF
       END DO
       IF ((is_simulated .EQ. 0) .OR. (is_contaminated .EQ. 1)) THEN
          !simulate the contaminated process
          CALL AR_GARCH_SNP_LAG(DGP_nb_par, DGP_par_vec, 1+DGP_order_model_dim(1), &
               DGP_order_model_dim(2), DGP_order_model_dim(3), &
               DGP_pol_dim, DGP_pol_deg, DGP_nb_pol_coef, DGP_pol_coef, DGP_epsilon_0, &
               T, time_series, DGP_nb_discard, &
               DGP_nb_iid_series, iid, 1) !The last 1 is "consider_contamination": 1 -> innov. out.
       END IF  !((is_simulated .EQ. 0) .OR. (is_contaminated .EQ. 1))
       !then for the additive outliers
       IF (contam_model(1) .EQ. 1) THEN
          is_contaminated = 1
          iid(1:DGP_nb_discard,1) = 0._realkind
          iid(DGP_nb_discard+1:T+DGP_nb_discard,1) = time_series(1:T)
          IF (contam_dist(1) .EQ. "_DIRAC_") THEN
             additive_contam_points(:)=contam_x0(1,:)
             CALL CONTAMINATE(T+DGP_nb_discard,iid(:,1),epsilon_p(1),nb_contam_points(1),x0=additive_contam_points)
          END IF
          IF (contam_dist(1) .EQ. "UNIFORM") THEN
             iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(1))
             iid_temp = iid_temp * (contam_upper_bound(1) -  contam_lower_bound(1)) + &
                  contam_lower_bound(1)
             CALL CONTAMINATE(T+DGP_nb_discard,iid(:,1),epsilon_p(1),nb_contam_points(1),random_x0=iid_temp)
          END IF
          IF (contam_dist(1) .EQ. "NORMAL_") THEN
             iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(1))
             iid_temp = iid_temp * SQRT(contam_norm_var(1)) + contam_norm_mean(1)
             CALL CONTAMINATE(T+DGP_nb_discard,iid(:,1),epsilon_p(1),nb_contam_points(1),random_x0=iid_temp)
          END IF
          IF (contam_dist(1) .EQ. "STUDENT") THEN
             iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(1),contam_df(1))
             iid_temp = iid_temp * SQRT(contam_student_var(1)) + contam_student_mean(1)
             CALL CONTAMINATE(T+DGP_nb_discard,iid(:,1),epsilon_p(1),nb_contam_points(1),random_x0=iid_temp)
          END IF
          time_series(1:T) = iid(DGP_nb_discard+1:T+DGP_nb_discard,1)
       ELSE
          iid(1:T+DGP_nb_discard,1) = 0._realkind
          IF (contam_dist(1) .EQ. "UNIFORM") THEN
             iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(1))
             !iid_temp = iid_temp * (contam_upper_bound(1) -  contam_lower_bound(1)) + &
             !     contam_lower_bound(1)
          END IF
          IF (contam_dist(1) .EQ. "NORMAL_") THEN
             iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(1))
             !iid_temp = iid_temp * SQRT(contam_norm_var(i+1)) + contam_norm_mean(i+1)
          END IF
          IF (contam_dist(1) .EQ. "STUDENT") THEN
             iid_temp = IID_SIM(T+DGP_nb_discard,1,contam_dist(1),contam_df(1))
             !iid_temp = iid_temp * SQRT(contam_student_var(1)) + contam_student_mean(1)
          END IF
          CALL CONTAMINATE(T+DGP_nb_discard,iid(:,1),0._realkind,nb_contam_points(1),x0=zero_additive_contam_points)
      END IF
       !save the contaminated process if necessary, i.e. if nb_rep_monte_carlo=1 and
       !is_contaminated .EQ. 1 and write_output_2 = true
       IF ((nb_rep_monte_carlo .GT. 0) .AND. (is_contaminated .EQ. 1) .AND. write_output_2) THEN
          OPEN (UNIT=unit_nb,FILE="series_contaminated.dat",FORM="FORMATTED") 
          DO i=1, T    
             write (unit=unit_nb,fmt=*) time_series(i) 
          END DO
          CLOSE(UNIT=unit_nb)
       END IF
    END IF  !(DGP_pol_dim .LE. 1)

    IF (DGP_name .EQ. "STOC_VOLATILITY") THEN
       DEALLOCATE(simulated_not_observable)
    END IF

    RETURN
  END SUBROUTINE DGP_SIMULATION


  SUBROUTINE STRUCTURAL_MODEL_SIMULATION(n,rho,T,simulated_series,nb_discard,nb_iid_series,random_iid)
    USE S2_COMMON, ONLY: name_SM, SM_dim_parts
    USE MYFUNCTIONS, ONLY: ARMA_SIMULATION, SVM1_SIMULATION,  ARMA_GARCH, SWITCH_AR
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: rho
    INTEGER, INTENT(IN) :: T
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:T) :: simulated_series 
    INTEGER, INTENT(IN) :: nb_discard
    INTEGER, INTENT(IN) :: nb_iid_series
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T+nb_discard,nb_iid_series) :: random_iid

    INTEGER :: nb_mod_par
    
    IF (name_SM .EQ. "ARMA(p,q)_MODEL") THEN
       CALL ARMA_SIMULATION(n,rho,SM_dim_parts(1),SM_dim_parts(2),T,simulated_series, &
            nb_discard,nb_iid_series,random_iid)
    END IF

    IF (name_SM .EQ. "STOC_VOLATILITY") THEN
       CALL SVM1_SIMULATION(n,rho,SM_dim_parts(1),SM_dim_parts(2),T,simulated_series, &
            nb_discard,nb_iid_series,random_iid)
    END IF

    IF (name_SM .EQ. "ARMA_GARCH_MOD_") THEN
       CALL ARMA_GARCH(n,rho, SM_dim_parts(1), SM_dim_parts(2), SM_dim_parts(3), &
            SM_dim_parts(4),T,simulated_series,nb_discard,nb_iid_series,random_iid)
    END IF

    IF (name_SM .EQ. "SWITCH_AR_MODEL") THEN
       nb_mod_par = SM_dim_parts(1)**2-SM_dim_parts(1)+(SM_dim_parts(2)+2)*SM_dim_parts(1)
       CALL SWITCH_AR(SM_dim_parts(1),SM_dim_parts(2),nb_mod_par,rho,T,simulated_series, &
            nb_discard, nb_iid_series, random_iid)
    END IF

    RETURN
  END SUBROUTINE STRUCTURAL_MODEL_SIMULATION
END MODULE SIMULATION




