
MODULE STARTUP
  USE DATATYPES
  IMPLICIT NONE

CONTAINS

  SUBROUTINE MOVE_TO(unit_nb, string)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: unit_nb
    CHARACTER (LEN=*), INTENT(IN) :: string

    CHARACTER (LEN=5) :: text

    DO
       READ (UNIT=unit_nb,FMT="(A5)") text
       IF (text .EQ. string) EXIT 
    END DO

    RETURN
  END SUBROUTINE MOVE_TO


  SUBROUTINE NPSOL_SETUP(index)
    USE S1_COMMON 
    USE S2_COMMON, ONLY: SM_nb_par, SM_nb_free_par, is_fixed_par_SM, is_restricted_par_SM, &
         is_gt_par_SM, is_lt_par_SM, gt_par_SM_value, lt_par_SM_value, nclin2, ncnln2, &
         is_gt_lin_constr_SM, is_lt_lin_constr_SM, A_lin_constr_SM, &
         gt_lin_constr_SM_value, lt_lin_constr_SM_value
    IMPLICIT NONE
    INTERFACE
       SUBROUTINE NPOPTN( string )
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN) :: string
       END SUBROUTINE NPOPTN
    END INTERFACE
    INTEGER, INTENT(IN) :: index

    INTEGER :: i, j, number_parameters
    REAL (KIND=realkind) :: eps
    CHARACTER(LEN=40) :: txt

    eps = 1._realkind

    !Modify the file npsol_option_file.txt and set the number of parameters
    !according to the required step
    IF (Nolist .EQ. 1) THEN
       CALL NPOPTN("Nolist")
    END IF

    write(txt,FMT="(A,F12.9)") "Step limit ", step_limit(index)
    CALL NPOPTN(txt)

    IF (npsol_first_call) THEN
       IF (infinite_bound_size .NE. -1._realkind) THEN
          write(txt,FMT="(A,F12.9)") "Infinite bound size ", infinite_bound_size
          CALL NPOPTN(txt)                                       
       END IF

       IF (feasibility_tolerance .NE. -1._realkind) THEN
          write(txt,FMT="(A,F12.9)") "Feasibility tolerance ", feasibil_toler
          CALL NPOPTN(txt) 
       END IF

       IF (crash_tolerance .NE. 0.01_realkind) THEN
          write(txt,FMT="(A,F12.9)") "Crash Tolerance ", crash_tolerance
          CALL NPOPTN(txt) 
       END IF

       IF (major_iteration_limit .NE. -1) THEN
          write(txt,FMT="(A,I4)") "Major iteration limit ", major_iteration_limit
          CALL NPOPTN(txt)
       END IF

       IF (minor_iteration_limit .NE. -1) THEN
          write(txt,FMT="(A,I4)") "Minor iteration limit ", minor_iteration_limit
          CALL NPOPTN(txt)
       END IF

       IF (Print_File .EQ. 0) THEN
          CALL NPOPTN("Print File = 0")
       ELSE
          IF (Print_File .NE. -1) THEN
             write(txt,FMT="(A,I4)") "Print File ", Print_File
             CALL NPOPTN(txt)
          END IF
       END IF

       IF (major_print_level .NE. 10) THEN
          write(txt,FMT="(A,I4)") "Major print level ", major_print_level
          CALL NPOPTN(txt)
       END IF

       IF (minor_print_level .NE. 0) THEN
          write(txt,FMT="(A,I4)") "Minor print level ", minor_print_level
          CALL NPOPTN(txt)
       END IF
       npsol_first_call = .FALSE.
    END IF

    IF (function_precision(index) .NE. -1._realkind) THEN
       write(txt,FMT="(A,F12.9)") "Function precision ", function_precision(index)
       CALL NPOPTN(txt)
    END IF

    IF (optimality_tolerance(index) .NE. -1._realkind) THEN
       write(txt,FMT="(A,F12.9)") "Optimality tolerance ", optimality_tolerance(index)
       CALL NPOPTN(txt)
    END IF

    IF (line_search_tolerance(index) .NE. 0.9_realkind) THEN
       write(txt,FMT="(A,F12.9)") "Linesearch tolerance ", line_search_tolerance(index)
       CALL NPOPTN(txt)
    END IF

    IF (index .EQ.1) THEN
       IF (derivative_level .EQ. 0) THEN
          CALL NPOPTN("Derivative level 0")
          CALL NPOPTN("Verify No")
       ELSE
          CALL NPOPTN("Derivative level 1")
          IF (verify_gradients .EQ. 1) THEN
             CALL NPOPTN("Verify level 3")
          ELSE
             CALL NPOPTN("Verify No")
          END IF
       END IF
       number_parameters = AM_nb_free_par
       nclin = nclin1
       ncnln = ncnln1
    END IF

    IF (index .EQ. 3) THEN
       CALL NPOPTN("Derivative level 0")
       CALL NPOPTN("Verify No")
       number_parameters = AM_nb_free_par
       nclin = nclin1
       ncnln = ncnln1
    END IF

    IF ((index .EQ. 2) .OR. (index .EQ. 4)) THEN
       CALL NPOPTN("Derivative level 0")
       CALL NPOPTN("Verify No")
       CALL NPOPTN("Difference interval 0.00001")
       number_parameters = SM_nb_free_par
       nclin = nclin2
       ncnln = ncnln2
    END IF

    !need to deallocate first
    DEALLOCATE(bl,bu,A,cJac,istate,clamda,RR,iw,w)

    nctotl = number_parameters + nclin + ncnln
    ALLOCATE(bl(1:nctotl),bu(1:nctotl))
    !Set bounds in bl and bu vectors (see manual of npsol)
    DO i=1,nctotl
       bl(i) = - bigbnd
       bu(i) = bigbnd
    END DO

    ldA = MAX(1,nclin)
    IF (nclin .EQ. 0) THEN
       ALLOCATE(A(1:1,1:1)) 
    ELSE
       ALLOCATE(A(1:ldA,1:number_parameters))
       IF ((index .EQ. 1) .OR. (index .EQ. 3)) THEN
          A = A_lin_constr_AM
       ELSE
          A = A_lin_constr_SM
       END IF
    END IF

    ldJ = MAX(1,ncnln)

    IF (ncnln .EQ. 0) THEN   
       ALLOCATE(cJac(1:ldJ,1:1))
    ELSE
       ALLOCATE(cJac(1:ldJ,1:number_parameters))
    END IF

    ALLOCATE(istate(1:nctotl))
    ldR = number_parameters
    !!!!
    ALLOCATE(clamda(1:nctotl),RR(1:ldR,1:number_parameters+20))

    leniw = 3*number_parameters + nclin + 2*ncnln
    IF ((nclin .EQ. 0) .AND. (ncnln .EQ. 0)) THEN
       lenw = 20*number_parameters
    ELSE
       IF ((ncnln .EQ. 0)) THEN
          lenw = 2*number_parameters**2 + 20*number_parameters + 11*nclin
       ELSE
          lenw = 2*number_parameters**2 + number_parameters*nclin + 2*number_parameters*ncnln + &
               20*number_parameters + 11*nclin + 21*ncnln
       END IF
    END IF
    ALLOCATE(iw(1:leniw),w(1:lenw))

    !Now insert the restriction on the parameters of the process
    IF ((index .EQ. 1) .OR. (index .EQ. 3)) THEN !estimation of the auxiliary model
       j = 0
       !loop over all parameters of the AM
       DO i=1, AM_nb_par
          IF ((.NOT. is_fixed_par_AM(i)) .AND. (.NOT. is_restricted_par_AM(i))) THEN
             j = j + 1
             IF (is_gt_par_AM(i)) THEN
                bl(j) = gt_par_AM_value(i) + restriction_correction_factor
             END IF

             IF (is_lt_par_AM(i)) THEN
                bu(j) = lt_par_AM_value(i) - restriction_correction_factor
             END IF
          END IF
       END DO    !loop over all parameters of the AM
       DO i=1, nclin
          IF (is_gt_lin_constr_AM(i)) THEN
             bl(AM_nb_free_par+i) = gt_lin_constr_AM_value(i) + restriction_correction_factor
          END IF
          IF (is_lt_lin_constr_AM(i)) THEN
             bu(AM_nb_free_par+i) = lt_lin_constr_AM_value(i) - restriction_correction_factor
          END IF
       END DO
    ELSE            !if index .EQ. 2
       j = 0
       !loop over all parameters of the SM
       DO i=1, SM_nb_par
          IF ((.NOT. is_fixed_par_SM(i)) .AND. (.NOT. is_restricted_par_SM(i))) THEN
             j = j + 1
             IF (is_gt_par_SM(i)) THEN
                bl(j) = gt_par_SM_value(i) + restriction_correction_factor
             END IF

             IF (is_lt_par_SM(i)) THEN
                bu(j) = lt_par_SM_value(i) - restriction_correction_factor
             END IF
          END IF
       END DO    !loop over all parameters of the SM
       DO i=1, nclin
          IF (is_gt_lin_constr_SM(i)) THEN
             bl(SM_nb_free_par+i) = gt_lin_constr_SM_value(i) + restriction_correction_factor
          END IF

          IF (is_lt_lin_constr_SM(i)) THEN
             bu(SM_nb_free_par+i) = lt_lin_constr_SM_value(i) - restriction_correction_factor
          END IF
       END DO
    END IF       !if index .EQ. 2

    RETURN
  END SUBROUTINE NPSOL_SETUP

  SUBROUTINE CHECK_SM_PAR_RESTRICTIONS()
    USE S2_COMMON
    IMPLICIT NONE

    INTEGER :: i,j
    REAL (KIND=realkind), DIMENSION(1:SM_nb_par) :: test_real

    DO i=1, SM_nb_par
       IF ((is_fixed_par_SM(i)) .AND. (is_gt_par_SM(i))) THEN
          WRITE(*,*) "Inequality > for parameter", i, " conflicts with strict equality!"
          WRITE(*,*) "Please modify equality/inequality constraints. Program stopped."
          STOP
       END IF
    END DO

    DO i=1, SM_nb_par
       IF ((is_fixed_par_SM(i)) .AND. (is_lt_par_SM(i))) THEN
          WRITE(*,*) "Inequality < for parameter", i, " conflicts with strict equality!"
          WRITE(*,*) "Please modify equality/inequality constraints. Program stopped."
          STOP
       END IF
    END DO

    test_real = lt_par_SM_value - gt_par_SM_value

    DO i=1, SM_nb_par
       IF ((is_gt_par_SM(i)) .AND. (is_lt_par_SM(i))) THEN
          IF (test_real(i) .LE. 0._realkind) THEN
             WRITE(*,*) "The inequality retrictions are wrong. Program stopped."
             STOP
          END IF
       END IF
    END DO
    !count the number of excluded parameters
    j = 0
    DO i=1, SM_nb_par
       IF (is_fixed_par_SM(i)) THEN
          j = j + 1
       END IF
    END DO

    IF (SM_nb_par .EQ. j) THEN
       WRITE(*,*) "All parameters of the SM are excluded. No parameters to estimate. Program stopped"
       STOP
    END IF

    SM_nb_not_fixed_par = SM_nb_par - j
    ALLOCATE(rho_rest_SM_start(1:SM_nb_not_fixed_par))
    ALLOCATE(rho_rest_SM_end(1:SM_nb_not_fixed_par))
    ALLOCATE(rho_rest_SM_old(1:SM_nb_not_fixed_par))
    rho_rest_SM_old = -999999.9_realkind

    RETURN
  END SUBROUTINE CHECK_SM_PAR_RESTRICTIONS


  SUBROUTINE CHECK_AM_PAR_RESTRICTIONS()
    USE S1_COMMON
    IMPLICIT NONE

    INTEGER :: i,j
    REAL (KIND=realkind), DIMENSION(1:AM_nb_par) :: test_real

    DO i=1, AM_nb_par
       IF ((is_restricted_par_AM(i)) .AND. (is_gt_par_AM(i))) THEN
          WRITE(*,*) "AUXILIARY MODEL:"
          WRITE(*,*) "Inequality > for parameter", i, " conflicts with strict equality!"
          WRITE(*,*) "Please modify equality/inequality constraints. Program stopped."
          STOP
       END IF
    END DO

    DO i=1, AM_nb_par
       IF ((is_restricted_par_AM(i)) .AND. (is_lt_par_AM(i))) THEN
          WRITE(*,*) "AUXILIARY MODEL:"
          WRITE(*,*) "Inequality < for parameter", i, " conflicts with strict equality!"
          WRITE(*,*) "Please modify equality/inequality constraints. Program stopped."
          STOP
       END IF
    END DO

    test_real = lt_par_AM_value - gt_par_AM_value

    DO i=1, AM_nb_par
       IF ((is_gt_par_AM(i)) .AND. (is_lt_par_AM(i))) THEN
          IF (test_real(i) .LE. 0._realkind) THEN
             WRITE(*,*) "AUXILIARY MODEL:"
             WRITE(*,*) "The inequality retrictions are wrong. Program stopped."
             STOP
          END IF
       END IF
    END DO
    !count the number of excluded parameters
    j = 0
    DO i=1, AM_nb_par
       IF (is_fixed_par_AM(i)) THEN
          j = j + 1
       END IF
    END DO

    IF (AM_nb_par .EQ. j) THEN
       WRITE(*,*) "All parameters of the AM are excluded. No parameters to estimate. Program stopped"
       STOP
    END IF

    RETURN
  END SUBROUTINE CHECK_AM_PAR_RESTRICTIONS

  SUBROUTINE GRID_POINTS_AM(AM_nb_par,start_vector)
    USE GRID_SEARCH_AM
    IMPLICIT NONE
    INTEGER :: AM_nb_par
    REAL (KIND=realkind), DIMENSION(1:AM_nb_par) :: start_vector
    
    INTEGER :: i,j,counter
    REAL (KIND=realkind) :: interval_lenght

    counter=1
    DO i=1, AM_nb_par
       IF (is_grid_search_AM(i)) THEN
             interval_lenght = (ub_grid_search_AM(i)-lb_grid_search_AM(i))/nb_grid_intervals_AM(i)
             DO j=1,nb_grid_points_AM(i)-1
                grid_points_values_AM(counter) = lb_grid_search_AM(i) + (j-1)*interval_lenght
                counter = counter + 1
             END DO
             grid_points_values_AM(counter) = ub_grid_search_AM(i)
             counter = counter + 1            
       ELSE
             grid_points_values_AM(counter) = start_vector(i)
             counter = counter + 1
       END IF
    END DO

    RETURN
  END SUBROUTINE GRID_POINTS_AM

  SUBROUTINE CHECK_GRID_SEARCH_AM()
    USE S1_COMMON, ONLY: AM_nb_par, val_fix_restr_par_AM
    USE GRID_SEARCH_AM
    IMPLICIT NONE

    INTEGER :: i

    !check the feasibility of the grid search
    do_grid_search_AM = .FALSE.
    DO i=1,AM_nb_par
       IF (is_grid_search_AM(i)) THEN
          do_grid_search_AM = .TRUE.
          EXIT
       END IF
    END DO
    IF (do_grid_search_AM) THEN
       !check the validity of the number of grid intervals
       DO i=1,AM_nb_par
          IF (is_grid_search_AM(i)) THEN
             IF (nb_grid_intervals_AM(i) .LE. 0) THEN
                WRITE(*,*) "The number of grid intervals must be > 0! Program stopped."
                STOP        
             END IF
          END IF
       END DO
       !check the validity of the lower and upper bounds
       DO i=1,AM_nb_par
          IF (is_grid_search_AM(i)) THEN
             IF (lb_grid_search_AM(i) .GT. ub_grid_search_AM(i)) THEN
                WRITE(*,*) "AM grid search: lower bound > upper bound. Program stopped."
                STOP        
             END IF
          END IF
       END DO
       !compute the number of grid points for every coefficients 
       !(fixed and restricted coeff. included!)
       ALLOCATE(nb_grid_points_AM(1:AM_nb_par))
       nb_grid_points_AM = nb_grid_intervals_AM + 1
       !compute the number of grid searches
       nb_grid_searches_AM = 1
       DO i=1,AM_nb_par
          nb_grid_searches_AM = nb_grid_searches_AM*nb_grid_points_AM(i)
       END DO
       !allocate the vector containing the resuts
       ALLOCATE(grid_search_results_AM(1:nb_grid_searches_AM))
       !compute the total number of grid points and the cumulative number of grid points (used
       !to recover the values from the grid_points_values_AM vector)
       ALLOCATE(cum_grid_points_AM(1:AM_nb_par))
       tot_nb_grid_points_AM = 0
       DO i=1,AM_nb_par
          cum_grid_points_AM(i) = tot_nb_grid_points_AM
          tot_nb_grid_points_AM = tot_nb_grid_points_AM+nb_grid_points_AM(i)
       END DO
       !compute the grid points for every coefficient
       ALLOCATE(grid_points_values_AM(1:tot_nb_grid_points_AM))     
       CALL GRID_POINTS_AM(AM_nb_par,val_fix_restr_par_AM)
    END IF !IF (do_grid_search)

    RETURN
  END SUBROUTINE CHECK_GRID_SEARCH_AM


  SUBROUTINE GRID_POINTS_SM(SM_nb_par,start_vector)
    USE GRID_SEARCH_SM
    IMPLICIT NONE
    INTEGER :: SM_nb_par
    REAL (KIND=realkind), DIMENSION(1:SM_nb_par) :: start_vector
    
    INTEGER :: i,j,counter
    REAL (KIND=realkind) :: interval_lenght

    counter=1
    DO i=1, SM_nb_par
       IF (is_grid_search_SM(i)) THEN
             interval_lenght = (ub_grid_search_SM(i)-lb_grid_search_SM(i))/nb_grid_intervals_SM(i)
             DO j=1,nb_grid_points_SM(i)-1
                grid_points_values_SM(counter) = lb_grid_search_SM(i) + (j-1)*interval_lenght
                counter = counter + 1
             END DO
             grid_points_values_SM(counter) = ub_grid_search_SM(i)
             counter = counter + 1            
       ELSE
             grid_points_values_SM(counter) = start_vector(i)
             counter = counter + 1
       END IF
    END DO

    RETURN
  END SUBROUTINE GRID_POINTS_SM

  
  
  SUBROUTINE CHECK_GRID_SEARCH_SM()
    USE S2_COMMON, ONLY: SM_nb_par, val_fix_restr_par_SM
    USE GRID_SEARCH_SM
    IMPLICIT NONE

    INTEGER :: i

    !check the feasibility of the grid search
    do_grid_search_SM = .FALSE.
    DO i=1,SM_nb_par
       IF (is_grid_search_SM(i)) THEN
          do_grid_search_SM = .TRUE.
          EXIT
       END IF
    END DO
    IF (do_grid_search_SM) THEN
       !check the validity of the number of grid intervals
       DO i=1,SM_nb_par
          IF (is_grid_search_SM(i)) THEN
             IF (nb_grid_intervals_SM(i) .LE. 0) THEN
                WRITE(*,*) "The number of grid intervals must be > 0! Program stopped."
                STOP        
             END IF
          END IF
       END DO
       !check the validity of the lower and upper bounds
       DO i=1,SM_nb_par
          IF (is_grid_search_SM(i)) THEN
             IF (lb_grid_search_SM(i) .GT. ub_grid_search_SM(i)) THEN
                WRITE(*,*) "SM grid search: lower bound > upper bound. Program stopped."
                STOP        
             END IF
          END IF
       END DO
       !compute the number of grid points for every coefficients
       !(fixed and restricted coeff. included!)
       ALLOCATE(nb_grid_points_SM(1:SM_nb_par))
       nb_grid_points_SM = nb_grid_intervals_SM + 1
       !compute the number of grid searches
       nb_grid_searches_SM = 1
       DO i=1,SM_nb_par
          nb_grid_searches_SM = nb_grid_searches_SM*nb_grid_points_SM(i)
       END DO
       !allocate the vector containing the resuts
       ALLOCATE(grid_search_results_SM(1:nb_grid_searches_SM))
       !compute the total number of grid points and the cumulative number of grid points (used
       !to recover the values from the grid_points_values_SM vector)
       ALLOCATE(cum_grid_points_SM(1:SM_nb_par))
       tot_nb_grid_points_SM = 0
       DO i=1,SM_nb_par
          cum_grid_points_SM(i) = tot_nb_grid_points_SM
          tot_nb_grid_points_SM = tot_nb_grid_points_SM+nb_grid_points_SM(i)
       END DO
       !compute the grid points for every coefficient
       ALLOCATE(grid_points_values_SM(1:tot_nb_grid_points_SM))
       CALL GRID_POINTS_SM(SM_nb_par,val_fix_restr_par_SM) 
    END IF !IF (do_grid_search)

    RETURN
  END SUBROUTINE CHECK_GRID_SEARCH_SM


  SUBROUTINE INITIATION()
    USE S1_COMMON
    USE S2_COMMON
    USE ROB_COMMON
    USE Z_MATRIX
    USE MONTE_CARLO
    USE DGP_DEFINITION
    USE CONTAMINATION_MODEL
    USE GRID_SEARCH_AM
    USE GRID_SEARCH_SM
    USE MYFUNCTIONS, ONLY: M_MATRIX_CONSTRUCTION
    USE DEVEL
    IMPLICIT NONE

    INTEGER :: i,j
    CHARACTER (LEN=1)  :: text
    CHARACTER (LEN=6)  :: text6
    CHARACTER (LEN=11) :: text11
    CHARACTER (LEN=2)  :: gt
    CHARACTER (LEN=2)  :: lt
    REAL (KIND=realkind) :: eps
    CHARACTER (LEN=space_lenght) :: space

    count = 100

    version = "SNP V3.0_beta"
    !set up the starting time
    CALL date_and_time (date=date_start, time=time_start)
    
    !some other constant    
    normal_c = (2*pi)**(-0.5)
    eps = 1._realkind
    space = " "

    OPEN (UNIT=tmp_unit_nb, FILE="setup.dat",STATUS="OLD",ACTION="READ")

    !#0 SECTION GENERAL
    CALL MOVE_TO(tmp_unit_nb,"0.01)")
    READ (UNIT=tmp_unit_nb,FMT=*) estimation_kind

    CALL MOVE_TO(tmp_unit_nb,"0.02)")
    READ (UNIT=tmp_unit_nb,FMT=*) step

    CALL MOVE_TO(tmp_unit_nb,"0.03)")
    READ (UNIT=tmp_unit_nb,FMT=*) simulate_data

    CALL MOVE_TO(tmp_unit_nb,"0.04)")
    READ (UNIT=tmp_unit_nb,FMT=*) T1

    ! se non simuli
    IF (simulate_data .EQ. 0) THEN
       nb_rep_monte_carlo = 0
       estimate_all_mc_simulations = 1
       nb_selected_simulations = 0
    ELSE

       CALL MOVE_TO(tmp_unit_nb,"0.05)")
       READ (UNIT=tmp_unit_nb,FMT=*) text
       READ (UNIT=tmp_unit_nb,FMT=*) text
       READ (UNIT=tmp_unit_nb,FMT=*) estimate_all_mc_simulations

       IF (estimate_all_mc_simulations .EQ. 1) THEN
          CALL MOVE_TO(tmp_unit_nb,"0.06)")
          READ (UNIT=tmp_unit_nb,FMT="(A20,I7)") text, nb_rep_monte_carlo
          nb_selected_simulations = 0
       ELSE 
          CALL MOVE_TO(tmp_unit_nb,"0.06)")
          READ (UNIT=tmp_unit_nb,FMT="(A20,I7)") text, nb_selected_simulations
          IF (nb_selected_simulations .LT. 1) THEN
              WRITE(*,*) "The number of selected simulations must be at least 1. Program stopped."
              STOP
          ELSE
             ALLOCATE(selected_simulations(1:nb_selected_simulations))
             DO i=1, nb_selected_simulations
               READ (UNIT=tmp_unit_nb, FMT="(A20,I7)") text, selected_simulations(i)
             END DO
             nb_rep_monte_carlo = MAXVAL(selected_simulations)
          END IF
       END IF
       
       IF (nb_rep_monte_carlo .EQ. 0) THEN
          WRITE(*,*) "The number of MC repetitions must be at least 1. Program stopped."
          STOP
       END IF
    END IF
    
    IF (nb_rep_monte_carlo .EQ. 0) THEN
       skip_first_n = 0
    ELSE
       CALL MOVE_TO(tmp_unit_nb,"0.07)")
       READ (UNIT=tmp_unit_nb,FMT=*) skip_first_n
    END IF

    CALL MOVE_TO(tmp_unit_nb,"0.08)")
    READ (UNIT=tmp_unit_nb,FMT=*) write_output_1
    READ (UNIT=tmp_unit_nb,FMT=*) write_output_2
    READ (UNIT=tmp_unit_nb,FMT=*) write_output_3

    CALL MOVE_TO(tmp_unit_nb,"0.09)")
    READ (UNIT=tmp_unit_nb,FMT=*) text
    READ (UNIT=tmp_unit_nb,FMT=*) text
    READ (UNIT=tmp_unit_nb,FMT=*) text
    READ (UNIT=tmp_unit_nb,FMT=*) do_grid_search
    
    CALL MOVE_TO(tmp_unit_nb,"0.10)")
    READ (UNIT=tmp_unit_nb,FMT=*) text
    READ (UNIT=tmp_unit_nb,FMT=*) text
    READ (UNIT=tmp_unit_nb,FMT=*) do_grid_search_first_n_loops


    
    IF (write_output_2) THEN
       OPEN (UNIT=log_unit_nb, FILE="out_log")
       WRITE(UNIT=log_unit_nb,FMT=*) "START READING SECTION #0"
       WRITE(UNIT=log_unit_nb,FMT=*) space, "estimation_kind..........", estimation_kind
       WRITE(UNIT=log_unit_nb,FMT=*) space, "step.....................", step
       WRITE(UNIT=log_unit_nb,FMT=*) space, "simulate_data............", simulate_data
       WRITE(UNIT=log_unit_nb,FMT=*) space, "T1.......................", T1
       WRITE(UNIT=log_unit_nb,FMT=*) space, "estimate_all_mc_simulat..", estimate_all_mc_simulations
       WRITE(UNIT=log_unit_nb,FMT=*) space, "nb_rep_monte_carlo.......", nb_rep_monte_carlo
       WRITE(UNIT=log_unit_nb,FMT=*) space, "nb_selected_simulations..", nb_selected_simulations
       IF (nb_selected_simulations .EQ. 0) THEN
            WRITE(UNIT=log_unit_nb,FMT=*) space, "selected_simulations....."
       ELSE
         DO i=1,nb_selected_simulations
            WRITE(UNIT=log_unit_nb,FMT=*) space, "selected_simulations.....", i, selected_simulations(i)
         END DO
       END IF
       WRITE(UNIT=log_unit_nb,FMT=*) space, "skip_first_n.............", skip_first_n
       WRITE(UNIT=log_unit_nb,FMT=*) space, "write_output_1............", write_output_1
       WRITE(UNIT=log_unit_nb,FMT=*) space, "write_output_2............", write_output_2
       WRITE(UNIT=log_unit_nb,FMT=*) space, "write_output_3............", write_output_3
       WRITE(UNIT=log_unit_nb,FMT=*) space, "do_grid_search...........", do_grid_search
       WRITE(UNIT=log_unit_nb,FMT=*) space, "do_grid_..._first_n_loops", do_grid_search_first_n_loops
       WRITE(UNIT=log_unit_nb,FMT=*)
       WRITE(UNIT=log_unit_nb,FMT=*) "START READING SECTION #1"
    END IF

    !#1 GENERATOR PROCESS OF THE DATA (ONLY IF 0.03 = 1)
    CALL MOVE_TO(tmp_unit_nb,"1.01)")
    READ (UNIT=tmp_unit_nb,FMT="(A15)") DGP_name
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "DGP_NAME................. ", DGP_name
    END IF

    CALL MOVE_TO(tmp_unit_nb,"1.02)")
    READ (UNIT=tmp_unit_nb,FMT=*) DGP_nb_iid_series
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "DGP_nb_iid_series........", DGP_nb_iid_series
    END IF
    !check the consistency between DGP_name and DGP_nb_iid_series
    IF (DGP_name .EQ. "MOVING__AVERAGE") THEN
       IF (DGP_nb_iid_series .NE. 1) THEN
          WRITE(*,*) "Wrong number of iid series! Program stopped."
          STOP
       END IF
    END IF
    IF (DGP_name .EQ. "AR_GARCH_MODEL_") THEN
       IF (DGP_nb_iid_series .NE. 1) THEN
          WRITE(*,*) "Wrong number of iid series! Program stopped."
          STOP
       END IF
    END IF
    IF (DGP_name .EQ. "STOC_VOLATILITY") THEN
       IF (DGP_nb_iid_series .NE. 2) THEN
          WRITE(*,*) "Wrong number of iid series! Program stopped."
          STOP
       END IF
    END IF

    ALLOCATE(DGP_innov_dist(DGP_nb_iid_series))
    ALLOCATE(DGP_deg_freedom(DGP_nb_iid_series))
    CALL MOVE_TO(tmp_unit_nb,"1.03)")
    DO i=1, DGP_nb_iid_series
       READ (UNIT=tmp_unit_nb,FMT="(A7)") DGP_innov_dist(i)
       IF (DGP_innov_dist(i) .EQ. "STUDENT") THEN
          BACKSPACE (UNIT=tmp_unit_nb)
          READ (UNIT=tmp_unit_nb,FMT=*) DGP_innov_dist(i), DGP_deg_freedom(i)
       ELSE
          DGP_deg_freedom(i) = -100
       END IF
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "DGP_innov_dist(",i,")...... ", DGP_innov_dist(i), &
               DGP_deg_freedom(i)
       END IF
    END DO
    !check if at least one innovation distribution is a SNP distribution
    DGP_pol_dim = 0
    j=0
    DO i=1, DGP_nb_iid_series
       IF (DGP_innov_dist(i) .EQ. "SNP_DST") THEN
          j=1
       END IF
    END DO
    IF (j .EQ. 1) THEN 
       CALL MOVE_TO(tmp_unit_nb,"1.04)")
       READ (UNIT=tmp_unit_nb,FMT=*) DGP_pol_dim, text
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "DGP_pol_dim..............", DGP_pol_dim
       END IF
       IF (DGP_pol_dim .LT. 1) THEN
          WRITE(*,*) "The polynomial dimension must be at least 1. Program stopped."
          STOP
       END IF
       ALLOCATE(DGP_pol_deg(1:DGP_pol_dim))
       IF (DGP_pol_dim .EQ. 1) THEN
          READ (UNIT=tmp_unit_nb,FMT=*) DGP_pol_deg(1), text
          DGP_nb_pol_coef = 1 + DGP_pol_deg(1)
          IF (write_output_2) THEN
             WRITE(UNIT=log_unit_nb,FMT=*) space, "DGP_pol_deg(1)...........", DGP_pol_deg(1)
             WRITE(UNIT=log_unit_nb,FMT=*) space, "DGP_nb_pol_coef..........", DGP_nb_pol_coef
          END IF
       END IF
       IF (DGP_pol_dim .GT. 1) THEN
          DO i=1, DGP_pol_dim
             READ (UNIT=tmp_unit_nb,FMT=*) DGP_pol_deg(i), text
             IF (write_output_2) THEN
                WRITE(UNIT=log_unit_nb,FMT=*) space, "DGP_pol_deg(",i,").......", DGP_pol_deg(i)
             END IF
          END DO
          DGP_nb_pol_coef = 1
          DO i=1, DGP_pol_dim
             DGP_nb_pol_coef = DGP_nb_pol_coef*(DGP_pol_deg(i)+1)
          END DO
          IF (write_output_2) THEN
             WRITE(UNIT=log_unit_nb,FMT=*) space, "DGP_nb_pol_coef.......", DGP_nb_pol_coef
          END IF
       END IF

       !read the values of the parameters
       ALLOCATE (DGP_pol_coef(DGP_nb_pol_coef))
       CALL MOVE_TO(tmp_unit_nb,"1.05)")
       READ (UNIT=tmp_unit_nb,FMT=*) DGP_epsilon_0
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "DGP_epsilon_0............", DGP_epsilon_0
       END IF
       DO i=1, DGP_nb_pol_coef
          READ (UNIT=tmp_unit_nb,FMT=*) DGP_pol_coef(i)
          IF (write_output_2) THEN
             WRITE(UNIT=log_unit_nb,FMT=*) space, "DGP_pol_coef(",i,")........", DGP_pol_coef(i)
          END IF
       END DO

       !check that the 1.03 selection table has been respected
       IF (DGP_name .NE. "AR_GARCH_MODEL_") THEN
          IF (DGP_pol_dim .GT. 1) THEN
             WRITE(*,*) "This distribution is not available for the selected model. Program stoped."
             STOP
          END IF
       END IF
    END IF

    CALL MOVE_TO(tmp_unit_nb,"1.06)")
    READ (UNIT=tmp_unit_nb,FMT=*) DGP_model_dim
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "DGP_model_dim............", DGP_model_dim
    END IF
    ALLOCATE(DGP_order_model_dim(1:DGP_model_dim))
    CALL MOVE_TO(tmp_unit_nb,"1.07)")
    DO i=1, DGP_model_dim
       READ (UNIT=tmp_unit_nb,FMT=*) DGP_order_model_dim(i)
    END DO
    IF (write_output_2) THEN
       DO i=1, DGP_model_dim
          WRITE(UNIT=log_unit_nb,FMT=*) space, "DGP_order_model_dim(",i,").", DGP_order_model_dim(i)
       END DO
    END IF

    !compute the total number of parameters in the DGP model
    IF (DGP_name .EQ. "ARMA(p,q)_MODEL") THEN
       DGP_nb_par = 2 + DGP_order_model_dim(1)+DGP_order_model_dim(2)
    END IF
    IF (DGP_name .EQ. "STOC_VOLATILITY") THEN
       DGP_nb_par = 3 + DGP_order_model_dim(1) + DGP_order_model_dim(2)
    END IF
    IF (DGP_name .EQ. "ARMA_GARCH_MOD_") THEN
       DGP_nb_par = 2 + DGP_order_model_dim(1) + DGP_order_model_dim(2) + &
            DGP_order_model_dim(3) + DGP_order_model_dim(4)
    END IF

    IF (DGP_name .EQ. "SWITCH_AR_MODEL") THEN
       DGP_nb_par = DGP_order_model_dim(1)**2-DGP_order_model_dim(1) + &
            (DGP_order_model_dim(2)+2)*DGP_order_model_dim(1)
    END IF

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "DGP_nb_par...............", DGP_nb_par
    END IF

    ALLOCATE(DGP_par_vec(DGP_nb_par))

    CALL MOVE_TO(tmp_unit_nb,"1.08)")
    DO i=1, DGP_nb_par
       READ (UNIT=tmp_unit_nb,FMT=*) DGP_par_vec(i)
    END DO
    IF (write_output_2) THEN
       DO i=1, DGP_nb_par
          WRITE(UNIT=log_unit_nb,FMT=*) space, "DGP_par_vec(",i,").........", DGP_par_vec(i)
       END DO
    END IF

    CALL MOVE_TO(tmp_unit_nb,"1.09)")
    READ (UNIT=tmp_unit_nb,FMT=*) DGP_nb_discard
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "DGP_nb_discard...........", DGP_nb_discard
       WRITE(UNIT=log_unit_nb,FMT=*)
       WRITE(UNIT=log_unit_nb,FMT=*) "START READING SECTION #2"
    END IF


    !#2 CONTAMINATION MODEL OF THE DATA (ONLY IF 0.03 = 1)

    ALLOCATE(contam_model(1+DGP_nb_iid_series))
    ALLOCATE(epsilon_p(1+DGP_nb_iid_series))
    ALLOCATE(contam_dist(1+DGP_nb_iid_series))
    ALLOCATE(contam_lower_bound(1+DGP_nb_iid_series))
    ALLOCATE(contam_upper_bound(1+DGP_nb_iid_series))
    ALLOCATE(contam_norm_mean(1+DGP_nb_iid_series))
    ALLOCATE(contam_norm_var(1+DGP_nb_iid_series))
    ALLOCATE(contam_df(1+DGP_nb_iid_series))
    ALLOCATE(contam_student_mean(1+DGP_nb_iid_series))
    ALLOCATE(contam_student_var(1+DGP_nb_iid_series))

    CALL MOVE_TO(tmp_unit_nb,"2.01)")
    READ(UNIT=tmp_unit_nb,FMT=*) contam_model(1), text
    READ(UNIT=tmp_unit_nb,FMT=*) contam_model(2), text
    DO i=2, DGP_nb_iid_series
       READ(UNIT=tmp_unit_nb,FMT=*) contam_model(1+i)
    END DO
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "contam_model(1)..........", contam_model(1)
       WRITE(UNIT=log_unit_nb,FMT=*) space, "contam_model(2)..........", contam_model(2)
       DO i=2, DGP_nb_iid_series
          WRITE(UNIT=log_unit_nb,FMT=*) space, "contam_model(",1+i,")........", contam_model(1+i)
       END DO
    END IF

    CALL MOVE_TO(tmp_unit_nb,"2.02)")
    READ(UNIT=tmp_unit_nb,FMT=*) epsilon_p(1), text
    READ(UNIT=tmp_unit_nb,FMT=*) epsilon_p(2), text
    DO i=2, DGP_nb_iid_series
       READ(UNIT=tmp_unit_nb,FMT=*) epsilon_p(1+i)
    END DO
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "epsilon_p(1).............", epsilon_p(1)
       WRITE(UNIT=log_unit_nb,FMT=*) space, "epsilon_p(2).............", epsilon_p(2)
       DO i=2, DGP_nb_iid_series
          WRITE(UNIT=log_unit_nb,FMT=*) space, "epsilon_p(",1+i,")...........", epsilon_p(1+i)
       END DO
    END IF

    ALLOCATE(nb_contam_points(DGP_nb_iid_series+1))
    CALL MOVE_TO(tmp_unit_nb,"2.03)")
    DO i=1, DGP_nb_iid_series+1 
       READ(UNIT=tmp_unit_nb,FMT=*) nb_contam_points(i) 
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "nb_contam_points(i)......", nb_contam_points(i)
       END IF
    END DO
    
    i=MAXVAL(nb_contam_points)
    ALLOCATE(contam_x0(1+DGP_nb_iid_series,i))
    CALL MOVE_TO(tmp_unit_nb,"2.04)")
    DO i=1, DGP_nb_iid_series + 1
       READ(UNIT=tmp_unit_nb,FMT="(A7,A2)",ADVANCE="NO") contam_dist(i)
       IF (contam_dist(i) .EQ. "_DIRAC_") THEN
          DO j=1, nb_contam_points(i)-1
            READ(UNIT=tmp_unit_nb,FMT="(A2,F5.1)",ADVANCE="NO") text, contam_x0(i,j)
          END DO
          j=MAX(1,nb_contam_points(i))
          READ(UNIT=tmp_unit_nb,FMT="(A2,F5.1)") text, contam_x0(i,j)
       END IF
       IF (contam_dist(i) .EQ. "UNIFORM") THEN
          READ(UNIT=tmp_unit_nb,FMT=*) contam_dist(i), contam_lower_bound(i), contam_upper_bound(i)
       END IF
       IF (contam_dist(i) .EQ. "NORMAL_") THEN
          READ(UNIT=tmp_unit_nb,FMT=*) contam_dist(i), contam_norm_mean(i), contam_norm_var(i)
       END IF
       IF (contam_dist(i) .EQ. "STUDENT") THEN
          READ(UNIT=tmp_unit_nb,FMT=*) contam_dist(i), contam_df(i), contam_student_mean(i), &
               contam_student_var(i)
       END IF
    END DO

    IF (write_output_2) THEN
       DO i=1, DGP_nb_iid_series+1
          IF (contam_dist(i) .EQ. "_DIRAC_") THEN
             WRITE(UNIT=log_unit_nb,FMT="(A,A,A7,A,F5.1)",ADVANCE="NO") space, &
             " Cont. Dist............... ", contam_dist(i), ", ", contam_x0(i,1)
             DO j=2,nb_contam_points(i)-1
               WRITE(UNIT=log_unit_nb,FMT="(A,F5.1)",ADVANCE="NO") ", ", contam_x0(i,j)
             END DO
             IF (nb_contam_points(i) .GT. 1) THEN
               WRITE(UNIT=log_unit_nb,FMT="(A,F5.1)") ", ", contam_x0(i,nb_contam_points(i))
             END IF
          END IF
          IF (contam_dist(i) .EQ. "UNIFORM") THEN
             WRITE(UNIT=log_unit_nb,FMT=*) space, "Cont. Dist...............", contam_dist(i), &
                  contam_lower_bound(i), contam_upper_bound(i)
          END IF
          IF (contam_dist(i) .EQ. "NORMAL_") THEN
             WRITE(UNIT=log_unit_nb,FMT=*) space, "Cont. Dist...............", contam_dist(i)
             WRITE(UNIT=log_unit_nb,FMT=*) space, "      Mean...............", contam_norm_mean(i)
             WRITE(UNIT=log_unit_nb,FMT=*) space, "      Var................", contam_norm_var(i)
          END IF
          IF (contam_dist(i) .EQ. "STUDENT") THEN
             WRITE(UNIT=log_unit_nb,FMT=*) space, "Cont. Dist...............", contam_dist(i), contam_df(i)
             WRITE(UNIT=log_unit_nb,FMT=*) space, "      Df.................", contam_df(i)
             WRITE(UNIT=log_unit_nb,FMT=*) space, "      Mean...............", contam_student_mean(i)
             WRITE(UNIT=log_unit_nb,FMT=*) space, "      Var................", contam_student_var(i)
          END IF
       END DO
       WRITE(UNIT=log_unit_nb,FMT=*)
       WRITE(UNIT=log_unit_nb,FMT=*) "START READING SECTION #3 (PART I)"
    END IF

    !#3 SECTION FIRST STEP --> ESTIMATION OF THE AUXILIARY MODEL
    !#FIRST STEP IN THE SETUP OF THE AUXILIARY MODEL
    CALL MOVE_TO(tmp_unit_nb,"3.01)")
    READ (UNIT=tmp_unit_nb,FMT=*) run_number
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "run_number...............", run_number
    END IF

    CALL MOVE_TO(tmp_unit_nb,"3.02)")
    READ (UNIT=tmp_unit_nb,FMT=*) nb_lags_mu
    !consider the constant term as part of the cond. mean function
    dim_mu = nb_lags_mu + 1
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "nb_lags_mu...............", nb_lags_mu
       WRITE(UNIT=log_unit_nb,FMT=*) space, "dim_mu (const. included).", dim_mu
    END IF


    CALL MOVE_TO(tmp_unit_nb,"3.03)")
    READ (UNIT=tmp_unit_nb,FMT=*) dim_arch
    READ (UNIT=tmp_unit_nb,FMT=*) dim_garch
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "dim_arch.................", dim_arch
       WRITE(UNIT=log_unit_nb,FMT=*) space, "dim_garch................", dim_garch
    END IF

    CALL MOVE_TO(tmp_unit_nb,"3.04)")
    READ (UNIT=tmp_unit_nb,FMT=*) pol_dim
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "pol_dim..................", pol_dim
    END IF

    IF ((dim_arch .EQ. 0) .AND. (dim_garch .GT. 0)) THEN
       WRITE(*,*) "The arch part must be at least 1 when garch > 0!"
       STOP 
    END IF
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*)
       WRITE(UNIT=log_unit_nb,FMT=*) "START READING SECTION #4"
    END IF


    !#4 SECTION STEP 2 --> ESTIMATION OF THE STRUCTURAL MODEL
    IF ((step .EQ. 2) .OR. (step .EQ. 3)) THEN
       !import the data for the step 2
       CALL MOVE_TO(tmp_unit_nb,"4.01)")
       READ (UNIT=tmp_unit_nb,FMT="(A15)") name_SM
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "name_SM.................."," ",  name_SM
       END IF

       CALL MOVE_TO(tmp_unit_nb,"4.02)")
       READ (UNIT=tmp_unit_nb,FMT=*) SM_model_parts
       ALLOCATE(SM_dim_parts(1:SM_model_parts))
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "SM_model_parts...........", SM_model_parts
       END IF

       CALL MOVE_TO(tmp_unit_nb,"4.03)")
       DO i=1, SM_model_parts
          READ (UNIT=tmp_unit_nb,FMT=*) SM_dim_parts(i)
       END DO
       IF (write_output_2) THEN
          DO i=1, SM_model_parts
             WRITE(UNIT=log_unit_nb,FMT=*) space, "SM_dim_parts(",i,")........", SM_dim_parts(i)
          END DO
       END IF

       !compute the total number of parameters in the structural model
       IF (name_SM .EQ. "ARMA(p,q)_MODEL") THEN
          SM_nb_par = 2 + SM_dim_parts(1)+SM_dim_parts(2)
       END IF
       IF (name_SM .EQ. "STOC_VOLATILITY") THEN
          SM_nb_par = 3 + SM_dim_parts(1) + SM_dim_parts(2)
       END IF
       IF (name_SM .EQ. "ARMA_GARCH_MOD_") THEN
          SM_nb_par = 2 + SM_dim_parts(1) + SM_dim_parts(2) + SM_dim_parts(3) + &
               SM_dim_parts(4)
          IF ((SM_dim_parts(3) .EQ. 0) .AND. (SM_dim_parts(4) .GT. 0)) THEN
             WRITE(*,*) "The arch dimension must be >= 1 when the garch dimension is > 0!"
             STOP
          END IF
       END IF
       IF (name_SM .EQ. "SWITCH_AR_MODEL") THEN
          SM_nb_par = (SM_dim_parts(1)**2-SM_dim_parts(1))+(SM_dim_parts(2)+2)*SM_dim_parts(1)
       END IF


       

       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "SM_nb_par................", SM_nb_par
       END IF

       CALL MOVE_TO(tmp_unit_nb,"4.04)")
       ALLOCATE(rho_SM_start(1:SM_nb_par))
       ALLOCATE(rho_SM_end(1:SM_nb_par))
       ALLOCATE(is_fixed_par_SM(1:SM_nb_par))
       ALLOCATE(is_restricted_par_SM(1:SM_nb_par))
       ALLOCATE(val_fix_restr_par_SM(1:SM_nb_par))
       ALLOCATE(is_gt_par_SM(1:SM_nb_par))
       ALLOCATE(is_lt_par_SM(1:SM_nb_par))
       ALLOCATE(gt_par_SM_value(1:SM_nb_par))
       ALLOCATE(lt_par_SM_value(1:SM_nb_par))

       READ (UNIT=tmp_unit_nb,FMT="(A)") text
       DO i=1, SM_nb_par
          READ (UNIT=tmp_unit_nb,FMT=*) is_fixed_par_SM(i), is_restricted_par_SM(i), &
               rho_SM_start(i), is_gt_par_SM(i), gt, gt_par_SM_value(i), &
               is_lt_par_SM(i), lt, lt_par_SM_value(i)
       END DO

       IF (write_output_2) THEN
100       FORMAT (A,L1,A,L1,A,E11.4,A,L1,A,E11.4,A,L1,A,E11.4)
          WRITE(UNIT=log_unit_nb,FMT=*) space, "exclude restrict     value       LBound           UBound"
          WRITE(UNIT=log_unit_nb,FMT=*) space,"PARAMETERS:"
          DO i=1, SM_nb_par
             WRITE(UNIT=log_unit_nb,FMT=100) "       ", is_fixed_par_SM(i), &
                  "       ", is_restricted_par_SM(i), &
                  "      ", rho_SM_start(i), &
                  "    ", is_gt_par_SM(i),  &
                  " > ", gt_par_SM_value(i), &
                  "  ", is_lt_par_SM(i), &
                  " < ", lt_par_SM_value(i)
          END DO
       END IF

       SM_nb_free_par=0
       DO i=1, SM_nb_par
          IF ((.NOT. is_fixed_par_SM(i)) .AND. (.NOT. is_restricted_par_SM(i))) THEN
             SM_nb_free_par = SM_nb_free_par + 1
          END IF
       END DO

       val_fix_restr_par_SM = rho_SM_start

       CALL CHECK_SM_PAR_RESTRICTIONS()

       !allocate the variable used in the grid_search
       ALLOCATE(is_grid_search_SM(1:SM_nb_par))
       ALLOCATE(lb_grid_search_SM(1:SM_nb_par))
       ALLOCATE(ub_grid_search_SM(1:SM_nb_par))
       ALLOCATE(nb_grid_intervals_SM(1:SM_nb_par))

       CALL MOVE_TO(tmp_unit_nb,"4.05)")
       READ (UNIT=tmp_unit_nb,FMT="(A)") text
       is_grid_search_SM = .FALSE.
       lb_grid_search_SM =  999999.9_realkind
       ub_grid_search_SM = -999999.9_realkind

       DO i=1, SM_nb_par
          IF ((.NOT. is_fixed_par_SM(i)) .AND. (.NOT. is_restricted_par_SM(i))) THEN
             READ (UNIT=tmp_unit_nb,FMT=*) is_grid_search_SM(i)
             IF (is_grid_search_SM(i)) THEN
                BACKSPACE (UNIT=tmp_unit_nb)
                READ (UNIT=tmp_unit_nb,FMT=*) is_grid_search_SM(i), lb_grid_search_SM(i), &
                     ub_grid_search_SM(i), nb_grid_intervals_SM(i)
             ELSE
                lb_grid_search_SM(i) = 0._realkind
                ub_grid_search_SM(i) = 0._realkind
             END IF
          END IF
       END DO

       DO i=1, SM_nb_par
          IF (.NOT. is_grid_search_SM(i)) THEN
             nb_grid_intervals_SM(i) = 0
          END IF
       END DO

       IF (write_output_2) THEN
          DO i=1, SM_nb_par
             WRITE(UNIT=log_unit_nb,FMT=*) space, "grid search_SM(",i,")......", is_grid_search_SM(i), &
                  lb_grid_search_SM(i), ub_grid_search_SM(i), nb_grid_intervals_SM(i)
          END DO
       END IF

       CALL CHECK_GRID_SEARCH_SM()

       CALL MOVE_TO(tmp_unit_nb,"4.06)")
       READ (UNIT=tmp_unit_nb,FMT=*) nclin2
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "nclin2...................", nclin2
       END IF

       IF (nclin2 .GT. 0) THEN
          CALL MOVE_TO(tmp_unit_nb,"4.07)")
          READ (UNIT=tmp_unit_nb,FMT="(A)") text
          ALLOCATE(A_lin_constr_SM(1:nclin2,1:SM_nb_free_par))
          ALLOCATE(is_gt_lin_constr_SM(1:nclin2))
          ALLOCATE(is_lt_lin_constr_SM(1:nclin2))
          ALLOCATE(gt_lin_constr_SM_value(1:nclin2))
          ALLOCATE(lt_lin_constr_SM_value(1:nclin2))
          DO i=1, nclin2
             READ (UNIT=tmp_unit_nb,FMT=*) A_lin_constr_SM(i,:), is_gt_lin_constr_SM(i), &
                  gt_lin_constr_SM_value(i), is_lt_lin_constr_SM(i), lt_lin_constr_SM_value(i)
          END DO
          IF (write_output_2) THEN
             DO i=1,nclin2
                WRITE(UNIT=log_unit_nb,FMT=*) space, "A(",i,",:).etc.............", A_lin_constr_SM(i,:), &
                     is_gt_lin_constr_SM(i), &
                     gt_lin_constr_SM_value(i), is_lt_lin_constr_SM(i), lt_lin_constr_SM_value(i)
             END DO
          END IF
          !verify the values
          DO i=1, nclin2
             IF ((.NOT. is_lt_lin_constr_SM(i)) .AND. (.NOT. is_gt_lin_constr_SM(i))) THEN
                WRITE(*,*) "Problems with the linear constraint", i, "of the SM."
                WRITE(*,*) "At least 1 linear bound must be active. Program stopped."
                STOP
             END IF
             IF (is_lt_lin_constr_SM(i) .AND. is_gt_lin_constr_SM(i)) THEN
                IF (gt_lin_constr_SM_value(i) .GT. lt_lin_constr_SM_value(i)) THEN
                   WRITE(*,*) "Wrong bounds on the linear constraints of the SM. Program stopped."
                   STOP
                END IF
             END IF
          END DO
       END IF

       CALL MOVE_TO(tmp_unit_nb,"4.08)")
       READ (UNIT=tmp_unit_nb,FMT=*) nb_iid_series
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "nb_iid_series............", nb_iid_series
       END IF

       CALL MOVE_TO(tmp_unit_nb,"4.09)")
       READ (UNIT=tmp_unit_nb,FMT=*) T2
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "T2.......................", T2
       END IF

       CALL MOVE_TO(tmp_unit_nb,"4.10)")
       READ (UNIT=tmp_unit_nb,FMT=*) nb_discard_step_2
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "nb_discard_step_2........", nb_discard_step_2
       END IF

       CALL MOVE_TO(tmp_unit_nb,"4.11)")
       READ (UNIT=tmp_unit_nb,FMT=*) max_T2
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "max_T2...................", max_T2
       END IF

       CALL MOVE_TO(tmp_unit_nb,"4.12)")
       READ (UNIT=tmp_unit_nb,FMT=*) text6
       READ (UNIT=tmp_unit_nb,FMT=*) expectation_analysis 
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "expectation_analysis.....", expectation_analysis
       END IF
       IF (expectation_analysis) THEN
          BACKSPACE (UNIT=tmp_unit_nb)
          READ (UNIT=tmp_unit_nb,FMT=*) expectation_analysis, nb_rep_exp_analysis
       ELSE
          nb_rep_exp_analysis = 0
       END IF
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "nb_rep_exp_analysis......", nb_rep_exp_analysis
       END IF

       CALL MOVE_TO(tmp_unit_nb,"4.13)")
       READ (UNIT=tmp_unit_nb,FMT=*) text6
       READ (UNIT=tmp_unit_nb,FMT=*) change_random_iid
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "change_random_iid........", change_random_iid
       END IF

       CALL MOVE_TO(tmp_unit_nb,"4.14)")
       READ (UNIT=tmp_unit_nb,FMT=*) delta_M_rho
       READ (UNIT=tmp_unit_nb,FMT=*) delta_method_M_rho

       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "delta_M_rho..............", delta_M_rho
          WRITE(UNIT=log_unit_nb,FMT=*) space, "delta_method_M_rho.......", delta_method_M_rho         
       END IF

       CALL MOVE_TO(tmp_unit_nb,"4.15)")
       READ (UNIT=tmp_unit_nb,FMT=*) method_derivative_M_rho

       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "method_derivative_M_rho..", method_derivative_M_rho
       END IF

       CALL MOVE_TO(tmp_unit_nb,"4.16)")
       READ (UNIT=tmp_unit_nb,FMT=*) compute_accuracy_M_rho

       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "compute_accuracy_M_rho...", compute_accuracy_M_rho
       END IF

       ncnln2 = 0
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "ncnln2...................", ncnln2
       END IF

    ELSE

       CALL MOVE_TO(tmp_unit_nb,"4.08)")
       READ (UNIT=tmp_unit_nb,FMT=*) nb_iid_series
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "nb_iid_series............", nb_iid_series
       END IF

       CALL MOVE_TO(tmp_unit_nb,"4.09)")
       READ (UNIT=tmp_unit_nb,FMT=*) T2
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "T2.......................", T2
       END IF

       CALL MOVE_TO(tmp_unit_nb,"4.10)")
       READ (UNIT=tmp_unit_nb,FMT=*) nb_discard_step_2
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "nb_discard_step_2........", nb_discard_step_2
       END IF

       CALL MOVE_TO(tmp_unit_nb,"4.11)")
       READ (UNIT=tmp_unit_nb,FMT=*) max_T2
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "max_T2...................", max_T2
       END IF

       expectation_analysis = .FALSE.
       nb_rep_exp_analysis = 0
       change_random_iid = .FALSE.
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "expectation_analysis.....", expectation_analysis
          WRITE(UNIT=log_unit_nb,FMT=*) space, "nb_rep_exp_analysis......", nb_rep_exp_analysis
          WRITE(UNIT=log_unit_nb,FMT=*) space, "change_random_iid........", change_random_iid
       END IF

    END IF

    IF (max_T2 .LT. T2) THEN
       WRITE(*,*) "The value in 4.11 must be greater or equal to that in 4.09"
       STOP
    END IF
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*)
       WRITE(UNIT=log_unit_nb,FMT=*) "START READING SECTION #5"
    END IF


    !#5 SECTION ROBUST ESTIMATION
    IF (estimation_kind .EQ. 1) THEN
       !now the SECTION ROBUST ESTIMATION
       CALL MOVE_TO(tmp_unit_nb,"5.01)")
       READ (UNIT=tmp_unit_nb,FMT=*) bound
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "bound....................", bound
       END IF

       CALL MOVE_TO(tmp_unit_nb,"5.02)")
       READ (UNIT=tmp_unit_nb,FMT=*) weights_metric
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "method for weights.......", weights_metric
       END IF

       CALL MOVE_TO(tmp_unit_nb,"5.03)")
       READ (UNIT=tmp_unit_nb,FMT=*) nb_A_updates_start
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "nb_A_updates_start.......", nb_A_updates_start
       END IF

       READ (UNIT=tmp_unit_nb,FMT=*) nb_A_updates
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "nb_A_updates.............", nb_A_updates
       END IF

       
       CALL MOVE_TO(tmp_unit_nb,"5.04)")
       READ (UNIT=tmp_unit_nb,FMT=*) nb_B_updates
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "nb_B_updates.............", nb_B_updates
       END IF

       CALL MOVE_TO(tmp_unit_nb,"5.05)")
       READ (UNIT=tmp_unit_nb,FMT=*) nb_M_rho_updates
       
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "nb_M_rho_updates.........", nb_M_rho_updates
       END IF


       IF (nb_M_rho_updates .NE. -1) THEN
          CALL MOVE_TO(tmp_unit_nb,"5.06)")
          READ (UNIT=tmp_unit_nb,FMT=*) wait_update_M_rho
       ELSE
          wait_update_M_rho = 1
       END IF
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "wait_update_M_rho........", wait_update_M_rho
       END IF
       
       CALL MOVE_TO(tmp_unit_nb,"5.07)")
       READ (UNIT=tmp_unit_nb,FMT=*) rob_start_value(1)
       READ (UNIT=tmp_unit_nb,FMT=*) rob_start_value(2)
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "rob_start_value(1).......", rob_start_value(1)
          WRITE(UNIT=log_unit_nb,FMT=*) space, "rob_start_value(2).......", rob_start_value(2)
       END IF

       CALL MOVE_TO(tmp_unit_nb,"5.08)")
       READ (UNIT=tmp_unit_nb,FMT=*) rob_iteration_limit
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "rob_iteration_limit......", rob_iteration_limit
       END IF

       CALL MOVE_TO(tmp_unit_nb,"5.09)")
       READ (UNIT=tmp_unit_nb,FMT=*) rob_parameter_tolerance
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "rob_parameter_tolerance..", rob_parameter_tolerance
       END IF

       CALL MOVE_TO(tmp_unit_nb,"5.10)")
       READ (UNIT=tmp_unit_nb,FMT=*) A_updating_kind
       IF (write_output_2) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "A_updating_kind..........", A_updating_kind
       END IF

       IF (step .EQ. 2) THEN
          WRITE(*,*) "Robust estimation of the second step only has no sense!"
          STOP
       END IF
    END IF

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*)
       WRITE(UNIT=log_unit_nb,FMT=*) "START READING SECTION #3 (PART II)"
    END IF


    ALLOCATE (imported_series(1:T1))
    ALLOCATE (modified_series(1:T1))

    !set the maximal length used by vectors and matrices used with modified_series and
    !simulated_series
    T_max = MAX(T1,T2)
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "T_max....................", T_max
    END IF

    !allocate the parameters of the conditional variance used to calculate the likelihood and
    !the gradient (subroutine Z_MATRIX_UPDATE 
    dim_sigma = dim_arch + dim_garch + 1
    sig2_start = MIN(1,dim_mu+dim_arch-dim_garch)
    ALLOCATE(sig2(sig2_start:T_max),inv_sigma(1:T_max))
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "dim_sigma................", dim_sigma
    END IF

    !set the first time period at which it's possible to compute the conditional density
    !of the auxiliary model
    t_start = MAX(dim_mu + dim_arch,pol_dim)
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "t_start..................", t_start
    END IF

    !#SECOND STEP IN THE SETUP OF THE AUXILIARY MODEL
    !set the dimension of the polynomial
    !IF pol_dim = 0 => The process is gaussian.
    IF (pol_dim .LT. 0) THEN
       WRITE(*,*) "The polynomial must have a positive dimension!"
       STOP
    END IF

    REWIND (UNIT=tmp_unit_nb) 
    IF (pol_dim .EQ. 0) THEN
       nb_pol_coef = 0
       !we need at least one degree of z to compute the normal density part
       deg_z_max = 1
       ALLOCATE (U_temp(1:T_max), U_square(1:T_max), S_Z(1:T_max,1:deg_z_max))
       !We need to allocate pol_deg even if we don't use it because of LIKELIHOOD_UPDATE procedure
       ALLOCATE (pol_deg(1:1))
       pol_deg(1) = 0
    ELSE
       IF (pol_dim .EQ. 1) THEN 
          ALLOCATE (pol_deg(1:1))
          !pol_deg(1) = 3 means that the "z" variable has maximal degree equal 3
          !pol_deg(1) must be > 0 otherwise f(epsilon, y_lag)=f(epsilon)*f(y_lag) and
          !the resulting cond. density is normal (the same as the pol_dim = 0)
          CALL MOVE_TO(tmp_unit_nb,"3.05)")
          READ (UNIT=tmp_unit_nb,FMT=*) pol_deg(1)
          IF (write_output_2) THEN
             WRITE(UNIT=log_unit_nb,FMT=*) space, "pol_deg(1)..................", pol_deg(1)
          END IF
          deg_z_max = pol_deg(1)
          IF (pol_deg(1) .LE. 0) THEN
             WRITE(*,*) "In this case you must set pol_dim = 0!"
             STOP
          END IF
          nb_pol_coef = pol_deg(1) + 1
          ALLOCATE (M(1:pol_deg(1)+1,1:pol_deg(1)+1))
          CALL M_MATRIX_CONSTRUCTION(pol_deg(1)+1,M)
          IF (write_output_2) THEN
             WRITE(UNIT=log_unit_nb,FMT=*) space, "M matrix:"
             DO i=1,pol_deg(1)+1
                WRITE(UNIT=log_unit_nb,FMT=*) space, "         ", M(i,:)
             END DO
          END IF
          ALLOCATE (U_temp(1:T_max), U_square(1:T_max), S_Z(1:T_max,1:deg_z_max))
       ELSE
          ALLOCATE (pol_deg(1:pol_dim))
          CALL MOVE_TO(tmp_unit_nb,"3.05)")
          DO i=1,pol_dim
             READ (UNIT=tmp_unit_nb,FMT=*) pol_deg(i)
          END DO

          deg_z_max = pol_deg(1)

          IF (pol_deg(1) .LE. 0) THEN
             WRITE(*,*) "In this case you must set pol_dim = 0!"
             STOP
          END IF

          !Check if the degrees are correct
          DO i=2, pol_dim
             IF (pol_deg(i) .LT. 0) THEN
                WRITE(*,*) "Error by definion of polynomial's degrees"
                STOP
             END IF
          END DO
          !Compute the total number of polinomial coefficients
          nb_pol_coef = 1
          DO i=1, pol_dim
             nb_pol_coef = nb_pol_coef*(pol_deg(i)+1)
          END DO

          ALLOCATE (M(1:pol_deg(1)+1,1:pol_deg(1)+1))
          CALL M_MATRIX_CONSTRUCTION(pol_deg(1)+1,M)
          IF (write_output_2) THEN
             WRITE(UNIT=log_unit_nb,FMT=*) space, "M matrix:"
             DO i=1,pol_deg(1)+1
                WRITE(UNIT=log_unit_nb,FMT=*) space, "         ", M(i,:)
             END DO
          END IF
          ALLOCATE (U_temp(1:T_max), U_square(1:T_max), S_Z(1:T_max,1:deg_z_max))

          deg_lag_max = MAXVAL(pol_deg(2:pol_dim))
          IF (deg_lag_max .EQ. 0) THEN
             WRITE(*,*) "The polynomial is bad defined. In this case set pol_dim = 1!"
             STOP
          END IF
          dim_rows_X = 0
          DO i=2, pol_dim
             dim_rows_X = dim_rows_X+pol_deg(i)
          END DO

          !In the X matrix the constant row corresponding to the potence 0 is omitted.
          ALLOCATE (X(1:T_max,1:dim_rows_X))

          !Initialisation of the nb_rep vector
          ALLOCATE (nb_rep(1:pol_dim)) 
          nb_rep = 1
          DO i=1, pol_dim-1
             nb_rep(pol_dim-i) = nb_rep(pol_dim-i+1)*(pol_deg(pol_dim-i+1)+1)
          END DO

          DO i=2,pol_dim
             IF (pol_deg(i) .EQ. 0) THEN
                nb_rep(i) = 0
             END IF
          END DO

          !Initialisation of the nb_jumps vector
          ALLOCATE (nb_jumps(1:pol_dim))
          nb_jumps = 1
          DO i=2, pol_dim
             nb_jumps(i) = nb_jumps(i-1)*(pol_deg(i-1)+1)
          END DO
          DO i=2,pol_dim
             IF (pol_deg(i) .EQ. 0) THEN
                nb_jumps(i) = 0
             END IF
          END DO

          !Initialisation of the high_jumps vector
          ALLOCATE (high_jumps(1:pol_dim))
          DO i=2, pol_dim
             high_jumps(i) = nb_rep(i)*(pol_deg(i)+1)
          END DO
       END IF
    END IF

    AM_nb_par = dim_mu + dim_sigma + nb_pol_coef

    ALLOCATE(beta_AM_start(1:AM_nb_par),beta_old(1:AM_nb_par),beta_AM_end(1:AM_nb_par))
    ALLOCATE(is_fixed_par_AM(1:AM_nb_par))
    ALLOCATE(is_restricted_par_AM(1:AM_nb_par))
    ALLOCATE(is_gt_par_AM(1:AM_nb_par))
    ALLOCATE(is_lt_par_AM(1:AM_nb_par))
    ALLOCATE(val_fix_restr_par_AM(1:AM_nb_par))
    ALLOCATE(gt_par_AM_value(1:AM_nb_par))
    ALLOCATE(lt_par_AM_value(1:AM_nb_par))

    !Now import the Hessian and score covariance matrix computation methods
    CALL MOVE_TO(tmp_unit_nb,"3.06)")
    READ (UNIT=tmp_unit_nb,FMT=*) sel_score_cov_mat_s1
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "sel_score_cov_mat_s1.....", sel_score_cov_mat_s1
    END IF

    CALL MOVE_TO(tmp_unit_nb,"3.07)")
    READ (UNIT=tmp_unit_nb,FMT=*) lags_score_cov_mat_s1
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "lags_score_cov_mat_s1....", lags_score_cov_mat_s1
    END IF

    CALL MOVE_TO(tmp_unit_nb,"3.08)")
    READ (UNIT=tmp_unit_nb,FMT=*) sel_cov_estimates_s1
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "sel_cov_estimates_s1.....", sel_cov_estimates_s1
    END IF

    CALL MOVE_TO(tmp_unit_nb,"3.09)")
    READ (UNIT=tmp_unit_nb,FMT=*) sel_hessian_mat_s1
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "sel_hessian_mat_s1.......", sel_hessian_mat_s1
    END IF

    CALL MOVE_TO(tmp_unit_nb,"3.10)")
    READ (UNIT=tmp_unit_nb,FMT=*)  delta_M_beta
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "delta_M_beta.............", delta_M_beta
    END IF

    !Now import the starting values of the conditional mean, conditional variance and
    !if dim_pol>0 the values of epsilon_0 and the polynomial coefficients.
    CALL MOVE_TO(tmp_unit_nb,"3.11)")
    READ (UNIT=tmp_unit_nb,FMT="(A)") text
    DO i=1, dim_mu
       READ (UNIT=tmp_unit_nb,FMT=*) is_fixed_par_AM(i), is_restricted_par_AM(i), &
            beta_AM_start(i), is_gt_par_AM(i), gt, gt_par_AM_value(i), &
            is_lt_par_AM(i), lt, lt_par_AM_value(i)
    END DO

    CALL MOVE_TO(tmp_unit_nb,"3.12)")
    READ (UNIT=tmp_unit_nb,FMT="(A)") text
    DO i=1, dim_sigma
       READ (UNIT=tmp_unit_nb,FMT=*) is_fixed_par_AM(dim_mu+i), is_restricted_par_AM(dim_mu+i), &
            beta_AM_start(dim_mu+i), is_gt_par_AM(dim_mu+i), gt, gt_par_AM_value(dim_mu+i), &
            is_lt_par_AM(dim_mu+i), lt, lt_par_AM_value(dim_mu+i)
    END DO

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "exclude restrict     value       LBound           UBound"
       WRITE(UNIT=log_unit_nb,FMT=*) space,"PARAMETERS:"
       DO i=1, dim_mu+dim_sigma
          WRITE(UNIT=log_unit_nb,FMT=100) "       ", is_fixed_par_AM(i), &
               "       ", is_restricted_par_AM(i), &
               "      ", beta_AM_start(i), &
               "    ", is_gt_par_AM(i),  &
               " > ", gt_par_AM_value(i), &
               "  ", is_lt_par_AM(i), &
               " < ", lt_par_AM_value(i)
       END DO
    END IF

    IF (nb_pol_coef .GT. 0) THEN
       CALL MOVE_TO(tmp_unit_nb,"3.13)")
       READ (UNIT=tmp_unit_nb,FMT="(A)") text
       DO i=1, nb_pol_coef
          READ (UNIT=tmp_unit_nb,FMT=*) is_fixed_par_AM(dim_mu+dim_sigma+i), &
               is_restricted_par_AM(dim_mu+dim_sigma+i), beta_AM_start(dim_mu+dim_sigma+i), &
               is_gt_par_AM(dim_mu+dim_sigma+i), gt, gt_par_AM_value(dim_mu+dim_sigma+i), &
               is_lt_par_AM(dim_mu+dim_sigma+i), lt, lt_par_AM_value(dim_mu+dim_sigma+i)
       END DO
       IF (write_output_2) THEN
          j = dim_mu+dim_sigma
          DO i=1, nb_pol_coef
             WRITE(UNIT=log_unit_nb,FMT=100) "       ", is_fixed_par_AM(j+i), &
                  "       ", is_restricted_par_AM(j+i), &
                  "      ", beta_AM_start(j+i), &
                  "    ", is_gt_par_AM(j+i),  &
                  " > ", gt_par_AM_value(j+i), &
                  "  ", is_lt_par_AM(j+i), &
                  " < ", lt_par_AM_value(j+i)
          END DO
       END IF
    END IF

    val_fix_restr_par_AM = beta_AM_start

    CALL CHECK_AM_PAR_RESTRICTIONS()

    !set the total number of parameter to be estimated and how many of them are in
    !the mean / variance / polynomial function
    j=0
    DO i=1, dim_mu
       IF (is_fixed_par_AM(i)) THEN 
          j=j+1
       END IF
    END DO
    nb_not_fixed_mu = dim_mu-j

    j=0
    DO i=dim_mu+1, dim_mu+dim_sigma
       IF (is_fixed_par_AM(i)) THEN
          j=j+1
       END IF
    END DO
    nb_not_fixed_sigma = dim_sigma-j

    j=0
    DO i=dim_mu+dim_sigma+1, AM_nb_par
       IF (is_fixed_par_AM(i)) THEN
          j=j+1
       END IF
    END DO
    nb_not_fixed_pol = AM_nb_par-dim_mu-dim_sigma-j 

    AM_nb_not_fixed_par = nb_not_fixed_mu + nb_not_fixed_sigma + nb_not_fixed_pol 

    DO i=1, AM_nb_par
       IF (is_fixed_par_AM(i) .AND. is_restricted_par_AM(i)) THEN
          WRITE(*,*) "Problem in the AM specification: parameter nb.", i,"."
          WRITE(*,*) "The parameter can not be excluded and restricted at the same time."
          STOP
       END IF
    END DO

    IF (AM_nb_not_fixed_par .EQ. 0) THEN
       WRITE(*,*) "AUXILIARY MODEL:"
       WRITE(*,*) "All parameters are excluded! Stop."
       STOP
    END IF

    i=MAX(1,nb_not_fixed_mu)
    ALLOCATE(j1(1:i))
    i=MAX(1,nb_not_fixed_sigma)
    ALLOCATE(j2(1:i))
    i=MAX(1,nb_not_fixed_pol )
    ALLOCATE(j3(1:i))

    j1=0
    j2=0
    j3=0

    j=0
    DO i=1, dim_mu
       IF (.NOT. is_fixed_par_AM(i)) THEN 
          j=j+1
          j1(j)=i
       END IF
    END DO

    j=0
    DO i=dim_mu+1, dim_mu+dim_sigma
       IF (.NOT. is_fixed_par_AM(i)) THEN
          j=j+1
          j2(j)=i-dim_mu
       END IF
    END DO

    j=0
    DO i=dim_mu+dim_sigma+1, AM_nb_par
       IF (.NOT. is_fixed_par_AM(i)) THEN
          j=j+1
          j3(j)=i-dim_mu-dim_sigma
       END IF
    END DO

    AM_nb_free_par=0
    DO i=1, AM_nb_par
       IF ((.NOT. is_fixed_par_AM(i)) .AND. (.NOT. is_restricted_par_AM(i))) THEN
          AM_nb_free_par = AM_nb_free_par + 1
       END IF
    END DO

    IF (AM_nb_free_par .EQ. 0) THEN
       WRITE(*,*) "The number of free parameters in the AM is 0!"
       STOP
    END IF

    !allocate the variable used in the grid_search
    ALLOCATE(is_grid_search_AM(1:AM_nb_par))
    ALLOCATE(lb_grid_search_AM(1:AM_nb_par))
    ALLOCATE(ub_grid_search_AM(1:AM_nb_par))
    ALLOCATE(nb_grid_intervals_AM(1:AM_nb_par))

    CALL MOVE_TO(tmp_unit_nb,"3.14)")
    READ (UNIT=tmp_unit_nb,FMT="(A)") text
    is_grid_search_AM = .FALSE.
    lb_grid_search_AM =  999999.9_realkind
    ub_grid_search_AM = -999999.9_realkind

    !only for the AM_nb_free_par is checked if a grid search is necessary or not
    DO i=1, AM_nb_par
       IF ((.NOT. is_fixed_par_AM(i)) .AND. (.NOT. is_restricted_par_AM(i))) THEN
          READ (UNIT=tmp_unit_nb,FMT=*) is_grid_search_AM(i)
          IF (is_grid_search_AM(i)) THEN
             BACKSPACE (UNIT=tmp_unit_nb)
             READ (UNIT=tmp_unit_nb,FMT=*) is_grid_search_AM(i), lb_grid_search_AM(i), &
                  ub_grid_search_AM(i), nb_grid_intervals_AM(i)
          ELSE
             lb_grid_search_AM(i) = 0._realkind
             ub_grid_search_AM(i) = 0._realkind
          END IF
       END IF
    END DO

    DO i=1, AM_nb_par
       IF (.NOT. is_grid_search_AM(i)) THEN
          nb_grid_intervals_AM(i) = 0
       END IF
    END DO

    IF (write_output_2) THEN
       DO i=1, AM_nb_par
          WRITE(UNIT=log_unit_nb,FMT=*) space, "grid search_AM(",i,")......", is_grid_search_AM(i), &
               lb_grid_search_AM(i), ub_grid_search_AM(i), nb_grid_intervals_AM(i)
       END DO
    END IF

    CALL CHECK_GRID_SEARCH_AM()

    CALL MOVE_TO(tmp_unit_nb,"3.15)")
    READ (UNIT=tmp_unit_nb,FMT=*) epsilon_0
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "epsilon_0................", epsilon_0
    END IF

    CALL MOVE_TO(tmp_unit_nb,"3.16)")
    READ (UNIT=tmp_unit_nb,FMT=*) nclin1
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "nclin1...................", nclin1
    END IF

    IF (nclin1 .GT. 0) THEN
       CALL MOVE_TO(tmp_unit_nb,"3.17)")
       READ (UNIT=tmp_unit_nb,FMT="(A)") text
       ALLOCATE(A_lin_constr_AM(1:nclin1,1:AM_nb_free_par))
       ALLOCATE(is_gt_lin_constr_AM(1:nclin1))
       ALLOCATE(is_lt_lin_constr_AM(1:nclin1))
       ALLOCATE(gt_lin_constr_AM_value(1:nclin1))
       ALLOCATE(lt_lin_constr_AM_value(1:nclin1))
       DO i=1, nclin1
          READ (UNIT=tmp_unit_nb,FMT=*) A_lin_constr_AM(i,:), is_gt_lin_constr_AM(i), &
               gt_lin_constr_AM_value(i), is_lt_lin_constr_AM(i), lt_lin_constr_AM_value(i)
       END DO
       IF (write_output_2) THEN
          DO i=1,nclin1
             WRITE(UNIT=log_unit_nb,FMT=*) space, "A(",i,",:).etc.............", A_lin_constr_AM(i,:), &
                  is_gt_lin_constr_AM(i), &
                  gt_lin_constr_AM_value(i), is_lt_lin_constr_AM(i), lt_lin_constr_AM_value(i)
          END DO
       END IF
       !verify the values
       DO i=1, nclin1
          IF ((.NOT. is_lt_lin_constr_AM(i)) .AND. (.NOT. is_gt_lin_constr_AM(i))) THEN
             WRITE(*,*) "Problems with the linear constraint", i, "of the AM."
             WRITE(*,*) "At least 1 linear bound must be active. Program stopped."
             STOP
          END IF
          IF (is_lt_lin_constr_AM(i) .AND. is_gt_lin_constr_AM(i)) THEN
             IF (gt_lin_constr_AM_value(i) .GT. lt_lin_constr_AM_value(i)) THEN
                WRITE(*,*) "Wrong bounds on the linear constraints of the AM. Program stopped."
                STOP
             END IF
          END IF
       END DO
    END IF

    CALL MOVE_TO(tmp_unit_nb,"3.18)")
    READ (UNIT=tmp_unit_nb,FMT=*) ncnln1
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "ncnln1...................", ncnln1
    END IF

    !Now read the integer used to determine if you want center and scale the lags.
    CALL MOVE_TO(tmp_unit_nb,"3.19)")
    READ (UNIT=tmp_unit_nb,FMT=*) center
    READ (UNIT=tmp_unit_nb,FMT=*) scale
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "center...................", center
       WRITE(UNIT=log_unit_nb,FMT=*) space, "scale....................", scale
    END IF

    !Now select if you want to compute the exact gradient (1) or not (0)
    CALL MOVE_TO(tmp_unit_nb,"3.20)")
    READ (UNIT=tmp_unit_nb,FMT=*) derivative_level
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "derivative_level.........", derivative_level
    END IF
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*)
       WRITE(UNIT=log_unit_nb,FMT=*) "START READING SECTION #6"
    END IF


    !#6 SECTION NPSOL SETUP AND PRINT OPTIONS
    npsol_first_call = .TRUE.

    CALL MOVE_TO(tmp_unit_nb,"6.01)")
    READ (UNIT=tmp_unit_nb,FMT=*) infinite_bound_size
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "infinite_bound_size......", infinite_bound_size
    END IF

    CALL MOVE_TO(tmp_unit_nb,"6.02)")
    READ (UNIT=tmp_unit_nb,FMT=*) feasibility_tolerance
    IF (feasibility_tolerance .EQ. -1._realkind) THEN
       feasibil_toler = sqrt(EPSILON(eps))
    ELSE
       feasibil_toler = feasibility_tolerance
    END IF

    !chech if the selected feasibility_toler value is compatible with the
    !restriction_correction_factor
    IF (feasibil_toler .GT. restriction_correction_factor) THEN
       WRITE(*,*) "The feasibility tolerance parameter in 6.03 must be 0 < feas.tol. < 1.E-07 !"
       STOP
    END IF
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "feasibility_tolerance....", feasibility_tolerance
       IF (feasibility_tolerance .EQ. -1._realkind) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "        default value....", feasibil_toler
       END IF
    END IF

    CALL MOVE_TO(tmp_unit_nb,"6.03)")
    READ (UNIT=tmp_unit_nb,FMT="(A)") text
    READ (UNIT=tmp_unit_nb,FMT=*) text11, function_precision(1)
    READ (UNIT=tmp_unit_nb,FMT=*) text11, function_precision(2)
    READ (UNIT=tmp_unit_nb,FMT=*) text11, function_precision(3)
    READ (UNIT=tmp_unit_nb,FMT=*) text11, function_precision(4)

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "function_precision.....s1", function_precision(1)
       WRITE(UNIT=log_unit_nb,FMT=*) space, "function_precision.....s2", function_precision(2)
       WRITE(UNIT=log_unit_nb,FMT=*) space, "function_precision s1 rob", function_precision(3)
       WRITE(UNIT=log_unit_nb,FMT=*) space, "function_precision s2 rob", function_precision(4)
       IF (function_precision(1) .EQ. -1._realkind) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "         default value s1", (EPSILON(eps))**0.9
       END IF
       IF (function_precision(2) .EQ. -1._realkind) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "         default value s2", (EPSILON(eps))**0.9
       END IF
       IF (function_precision(3) .EQ. -1._realkind) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "     default value s1 rob", (EPSILON(eps))**0.9
       END IF
       IF (function_precision(4) .EQ. -1._realkind) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "     default value s2 rob", (EPSILON(eps))**0.9
       END IF
    END IF

    CALL MOVE_TO(tmp_unit_nb,"6.04)")
    READ (UNIT=tmp_unit_nb,FMT="(A)") text
    READ (UNIT=tmp_unit_nb,FMT=*) text11, optimality_tolerance(1)
    READ (UNIT=tmp_unit_nb,FMT=*) text11, optimality_tolerance(2)
    READ (UNIT=tmp_unit_nb,FMT=*) text11, optimality_tolerance(3)
    READ (UNIT=tmp_unit_nb,FMT=*) text11, optimality_tolerance(4)

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "optimality_tolerance   s1", optimality_tolerance(1)
       WRITE(UNIT=log_unit_nb,FMT=*) space, "optimality_tolerance   s2", optimality_tolerance(2)
       WRITE(UNIT=log_unit_nb,FMT=*) space, "optimal._tolerance s1 rob", optimality_tolerance(3)
       WRITE(UNIT=log_unit_nb,FMT=*) space, "optimal._tolerance s2 rob", optimality_tolerance(4)

       IF (optimality_tolerance(1) .EQ. -1._realkind) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "         default value s1", (EPSILON(eps))**0.72
       END IF
       IF (optimality_tolerance(2) .EQ. -1._realkind) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "         default value s2", (EPSILON(eps))**0.72
       END IF
       IF (optimality_tolerance(3) .EQ. -1._realkind) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "     default value s1 rob", (EPSILON(eps))**0.72
       END IF
       IF (optimality_tolerance(4) .EQ. -1._realkind) THEN
          WRITE(UNIT=log_unit_nb,FMT=*) space, "     default value s2 rob", (EPSILON(eps))**0.72
       END IF
    END IF

    CALL MOVE_TO(tmp_unit_nb,"6.05)")
    READ (UNIT=tmp_unit_nb,FMT="(A)") text
    READ (UNIT=tmp_unit_nb,FMT=*) text11, line_search_tolerance(1)
    READ (UNIT=tmp_unit_nb,FMT=*) text11, line_search_tolerance(2)
    READ (UNIT=tmp_unit_nb,FMT=*) text11, line_search_tolerance(3)
    READ (UNIT=tmp_unit_nb,FMT=*) text11, line_search_tolerance(4)
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "line_search_tolerance  s1", line_search_tolerance(1)
       WRITE(UNIT=log_unit_nb,FMT=*) space, "line_search_tolerance  s2", line_search_tolerance(2)
       WRITE(UNIT=log_unit_nb,FMT=*) space, "line_s._tolerance  s1 rob", line_search_tolerance(3)
       WRITE(UNIT=log_unit_nb,FMT=*) space, "line_s._tolerance  s2 rob", line_search_tolerance(4)
    END IF
    
    CALL MOVE_TO(tmp_unit_nb,"6.06)")
    READ (UNIT=tmp_unit_nb,FMT="(A)") text
    READ (UNIT=tmp_unit_nb,FMT=*) text11, step_limit(1)
    READ (UNIT=tmp_unit_nb,FMT=*) text11, step_limit(2)
    READ (UNIT=tmp_unit_nb,FMT=*) text11, step_limit(3)
    READ (UNIT=tmp_unit_nb,FMT=*) text11, step_limit(4)
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "step_limit for step 1....", step_limit(1)
       WRITE(UNIT=log_unit_nb,FMT=*) space, "step_limit for step 3....", step_limit(2)
       WRITE(UNIT=log_unit_nb,FMT=*) space, "step_limit for step 1 rob", step_limit(3)
       WRITE(UNIT=log_unit_nb,FMT=*) space, "step_limit for step 2 rob", step_limit(4)
    END IF

    CALL MOVE_TO(tmp_unit_nb,"6.07)")
    READ (UNIT=tmp_unit_nb,FMT=*) crash_tolerance
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "crash_tolerance..........", crash_tolerance
    END IF

    CALL MOVE_TO(tmp_unit_nb,"6.08)")
    READ (UNIT=tmp_unit_nb,FMT=*) major_iteration_limit
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "major_iteration_limit....", major_iteration_limit
    END IF

    CALL MOVE_TO(tmp_unit_nb,"6.09)")
    READ (UNIT=tmp_unit_nb,FMT=*) minor_iteration_limit
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "minor_iteration_limit....", minor_iteration_limit
    END IF

    CALL MOVE_TO(tmp_unit_nb,"6.10)")
    READ (UNIT=tmp_unit_nb,FMT="(A)") text
    READ (UNIT=tmp_unit_nb,FMT="(A)") text
    READ (UNIT=tmp_unit_nb,FMT=*) Nolist
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "Nolist...................", Nolist
    END IF

    CALL MOVE_TO(tmp_unit_nb,"6.11)")
    READ (UNIT=tmp_unit_nb,FMT="(A)") text
    READ (UNIT=tmp_unit_nb,FMT=*) Print_File
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "Print_File...............", Print_File
    END IF

    CALL MOVE_TO(tmp_unit_nb,"6.12)")
    READ (UNIT=tmp_unit_nb,FMT="(A)") text
    READ (UNIT=tmp_unit_nb,FMT="(A)") text
    READ (UNIT=tmp_unit_nb,FMT="(A)") text
    READ (UNIT=tmp_unit_nb,FMT=*) major_print_level
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "major_print_level........", major_print_level
    END IF

    CALL MOVE_TO(tmp_unit_nb,"6.13)")
    READ (UNIT=tmp_unit_nb,FMT="(A)") text
    READ (UNIT=tmp_unit_nb,FMT="(A)") text
    READ (UNIT=tmp_unit_nb,FMT="(A)") text
    READ (UNIT=tmp_unit_nb,FMT=*) minor_print_level
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "minor_print_level........", minor_print_level
    END IF

    CALL MOVE_TO(tmp_unit_nb,"6.14)")
    READ (UNIT=tmp_unit_nb,FMT=*) verify_gradients
    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*) space, "verify_gradients.........", verify_gradients
    END IF

    IF (write_output_2) THEN
       WRITE(UNIT=log_unit_nb,FMT=*)
    END IF
    CLOSE (UNIT=tmp_unit_nb)

    !allocate the vectors containing the minimum of the robust steps
    
    i=MAX(1,nb_rep_monte_carlo)
    ALLOCATE(minimum_step_1_vector(i))
    ALLOCATE(minimum_step_2_vector(i))
    ALLOCATE(minimum_step_1_rob_vector(i))
    ALLOCATE(minimum_step_2_rob_vector(i))
    
    !allocate the vectors necessary to store the Monte Carlo results
    IF ((nb_rep_monte_carlo .GT. 0) .AND. (step .EQ. 2)) THEN
       WRITE(*,*) "Is not possible to do monte carlo only with step2"
       STOP
    ELSE
       IF (nb_rep_monte_carlo .GT. 0) THEN
          ALLOCATE(mc_results_step_1(AM_nb_par,nb_rep_monte_carlo))
          ALLOCATE(npsol_inform_step_1(1:nb_rep_monte_carlo))
          ALLOCATE(npsol_iter_step_1(1:nb_rep_monte_carlo))
          mc_results_step_1 = -999999.9_realkind
          npsol_inform_step_1 = 9
          npsol_iter_step_1 = 0
          IF (step .EQ. 3) THEN
             ALLOCATE(mc_results_step_2(SM_nb_par,nb_rep_monte_carlo,1:(1+nb_rep_exp_analysis)))
             ALLOCATE(hansen_test_results(1:nb_rep_monte_carlo,1:(1+nb_rep_exp_analysis)))
             ALLOCATE(npsol_inform_step_2(1:nb_rep_monte_carlo,1:(1+nb_rep_exp_analysis)))
             ALLOCATE(npsol_iter_step_2(1:nb_rep_monte_carlo,1:(1+nb_rep_exp_analysis)))
             mc_results_step_2 = -999999.9_realkind
             npsol_inform_step_2 = 9
             npsol_iter_step_2 = 0
          END IF
          IF (estimation_kind .EQ. 1) THEN
             ALLOCATE(hansen_test_results_rob(1:nb_rep_monte_carlo,1:(1+nb_rep_exp_analysis)))
             ALLOCATE(mc_results_step_1_rob(AM_nb_par,nb_rep_monte_carlo))
             ALLOCATE(npsol_inform_step_1_rob(1:nb_rep_monte_carlo))
             ALLOCATE(npsol_iter_step_1_rob(1:nb_rep_monte_carlo))
             ALLOCATE(inform_rob(1:nb_rep_monte_carlo))
             ALLOCATE(iter_rob(1:nb_rep_monte_carlo))
             mc_results_step_1_rob = -999999.9_realkind
             npsol_inform_step_1_rob = 9
             npsol_iter_step_1_rob = 0
             IF (step .EQ. 3) THEN
                ALLOCATE(mc_results_step_2_rob(SM_nb_par,nb_rep_monte_carlo))
                ALLOCATE(npsol_inform_step_2_rob(1:nb_rep_monte_carlo))
                ALLOCATE(npsol_iter_step_2_rob(1:nb_rep_monte_carlo))
                mc_results_step_2_rob = -999999.9_realkind
                npsol_inform_step_2_rob = 9
                npsol_iter_step_2_rob = 0
             END IF
          END IF
       END IF
    END IF

    !dimension variables used to compute the likelihood of the auxiliary model
    IF (pol_dim .GT. 1) THEN
       ALLOCATE(X_product(1:T_max,1:nb_pol_coef))
       ALLOCATE(b(1:T_max,1:pol_deg(1)+1)) 
    END IF

    !now dimension the matrices and vectors used to compute and store the gradient
    !of the auxiliary model
    ALLOCATE(grad_beta_U2(1:T_max,1:dim_mu))
    ALLOCATE(grad_beta_sig2(1:T_max,1:dim_mu))
    ALLOCATE(grad_sigma_sig2(1:T_max,1:dim_sigma))
    ALLOCATE(grad_const(1:T_max))
    ALLOCATE(grad_l_beta(1:dim_mu))
    ALLOCATE(grad_l_sigma(1:dim_sigma))
    ALLOCATE(sum_grad_beta_U2(1:dim_mu))

    IF ((pol_dim .GT.0)) THEN
       ALLOCATE(numerator(1:T_max))
       ALLOCATE(denumerator(1:T_max))
       ALLOCATE(pol(1:T_max))
       ALLOCATE(grad_beta_Z(1:T_max,1:dim_mu))
       ALLOCATE(grad_beta_pol(1:T_max,1:dim_mu))
       ALLOCATE(grad_sigma_pol(1:T_max,1:dim_sigma))
       ALLOCATE(grad_theta_pol(1:T_max,1:nb_pol_coef))
       ALLOCATE(fixed_grad_P(1:T_max))
       ALLOCATE(grad_l_theta(1:nb_pol_coef))
       ALLOCATE(grad_theta_denumerator(1:nb_pol_coef))        
       ALLOCATE(grad_theta_qMq(1:T_max,1:nb_pol_coef))
    END IF

    ALLOCATE(grad_l(1:T_max,1:AM_nb_par))
    i = T_max-t_start+1
    ALLOCATE(grad_i(1:i,1:AM_nb_not_fixed_par))
    ALLOCATE(score_cov(1:AM_nb_not_fixed_par,1:AM_nb_not_fixed_par))
    ALLOCATE(inv_score_cov(1:AM_nb_not_fixed_par,1:AM_nb_not_fixed_par))
    IF ((step .NE.1) .AND. (expectation_analysis)) THEN
       ALLOCATE(inv_score_cov_mc(1:AM_nb_not_fixed_par,1:AM_nb_not_fixed_par,1:nb_rep_monte_carlo))
    END IF
    ALLOCATE(M_beta(1:AM_nb_not_fixed_par,1:AM_nb_free_par))
    ALLOCATE(cov_estimates_s1(1:AM_nb_free_par,1:AM_nb_free_par))
    ALLOCATE(random_iid_step_2(1:(T2+nb_discard_step_2),nb_iid_series))
    ALLOCATE(simulated_series_step_2(1:T2))
    ALLOCATE(random_iid_rob_update(1:(T2+nb_discard_step_2),nb_iid_series))

    !allocate the variable for the robust estimation
    N1 = T1 - t_start + 1
    N2 = T2 - t_start + 1
    IF (estimation_kind .EQ. 1) THEN
       IF (weights_metric .EQ. 0) THEN
          ALLOCATE(rob_A(1:AM_nb_not_fixed_par,1:AM_nb_not_fixed_par))
       ELSE
          ALLOCATE(rob_A_tilde(1:SM_nb_free_par,1:SM_nb_free_par))
          ALLOCATE(M_rho(1:AM_nb_not_fixed_par,1:SM_nb_free_par))
          ALLOCATE(B_0(1:AM_nb_not_fixed_par,1:AM_nb_not_fixed_par))
          ALLOCATE(S(1:AM_nb_not_fixed_par,1:AM_nb_not_fixed_par))
          ALLOCATE(psi_N1_tilde(1:N1,1:SM_nb_free_par))
          ALLOCATE(psi_N2_tilde(1:N2,1:SM_nb_free_par))
       END IF
       ALLOCATE(rob_weights_N1(1:N1,1))
       ALLOCATE(rob_weights_N2(1:N2,1))
       ALLOCATE(psi_N1(1:N1,1:AM_nb_not_fixed_par))
       ALLOCATE(psi_N2(1:N2,1:AM_nb_not_fixed_par))
       ALLOCATE(score_N1(1:N1,1:AM_nb_not_fixed_par))
       ALLOCATE(score_N2(1:N2,1:AM_nb_not_fixed_par))
    END IF

    ALLOCATE(bl(1),bu(1),A(1,1),cJac(1,1),istate(1),clamda(1),RR(1,1),iw(1),w(1))
    bigbnd_default = (10._realkind)**20
    IF (infinite_bound_size .EQ. -1._realkind) THEN
       bigbnd = bigbnd_default
    ELSE
       bigbnd = infinite_bound_size
    END IF
    !the random seed
    CALL random_seed(size=size_seed)
    ALLOCATE (seed(size_seed))
    CALL RANDOM_SEED(get=seed)

    RETURN
  END SUBROUTINE INITIATION


  SUBROUTINE IMPORT_SERIES(unit_nb)
    USE S1_COMMON, ONLY: T1, imported_series, empirical_is_mean, empirical_is_variance, simulate_data
    USE DGP_DEFINITION
    USE CONTAMINATION_MODEL
    USE SIMULATION, ONLY: DGP_SIMULATION
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: unit_nb

    INTEGER :: i

    IF (simulate_data .EQ. 1) THEN
       CALL DGP_SIMULATION(T1,imported_series,unit_nb)      
    ELSE
       OPEN (UNIT=unit_nb, FILE="data.txt",ACTION="READ")
       !OPEN (UNIT=unit_nb, FILE="data.txt",STATUS="OLD",ACTION="READ")
       DO i=1, T1
          READ (UNIT=unit_nb,FMT=*) imported_series(i)
       END DO
       CLOSE (UNIT=unit_nb)
    END IF

    empirical_is_mean = 0.0_realkind
    DO i=1, T1
       empirical_is_mean = empirical_is_mean + imported_series(i)
    END DO
    empirical_is_mean = empirical_is_mean / T1

    !Computation of the empirical variance of the process
    empirical_is_variance = 0.0_realkind
    DO i=1, T1
       empirical_is_variance = empirical_is_variance + imported_series(i)**2
    END DO
    empirical_is_variance = empirical_is_variance/T1
    empirical_is_variance = empirical_is_variance - empirical_is_mean**2

    RETURN
  END SUBROUTINE IMPORT_SERIES


  SUBROUTINE MODIFY_SERIES()
    USE S1_COMMON, ONLY: T1, imported_series, modified_series, empirical_is_mean, &
         empirical_is_variance, empirical_ms_mean, empirical_ms_variance, &
         center, scale 
    IMPLICIT NONE

    INTEGER :: i
    REAL (KIND=realkind) :: inv_stdv

    empirical_ms_mean = empirical_is_mean
    empirical_ms_variance = empirical_is_variance
    modified_series = imported_series

    IF (center .EQ. 1) THEN 
       DO i=1, T1
          modified_series(i) = modified_series(i) - empirical_ms_mean    
       END DO
       empirical_ms_mean = 0._realkind
    END IF

    IF (scale .EQ. 1) THEN   
       inv_stdv = 1._realkind/SQRT(empirical_ms_variance)
       DO i=1, T1
          modified_series(i) = inv_stdv*modified_series(i)
       END DO
       empirical_ms_variance = 1._realkind
    END IF

    RETURN
  END SUBROUTINE MODIFY_SERIES

END MODULE STARTUP



MODULE LIKELIHOOD_UTILITIES
  USE DATATYPES
  IMPLICIT NONE

CONTAINS

  SUBROUTINE X_MATRIX_CONSTRUCTION(T, working_series)
    USE S1_COMMON, ONLY:  pol_deg, pol_dim, t_start, nb_pol_coef 
    USE Z_MATRIX, ONLY: X, X_product, nb_jumps, nb_rep, high_jumps
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: T
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T) :: working_series

    INTEGER :: i, j, k, s, count, time, nb_variable, deg_variable

    k = pol_dim - 1 !Dimension of the X vector
    !Computation of the powers of x at every date t starting from k+1
    count = 1
    DO j = 1, k
       IF (pol_deg(j+1) .GT. 0) THEN
          DO time=k+1,T
             X(time,count) = working_series(time-j)
          END DO
          count = count + 1
          DO i=2, pol_deg(j+1)
             DO time=k+1,T
                X(time,count) = X(time,count-1)*working_series(time-j)
             END DO
             count = count + 1
          END DO
       END IF
    END DO


    !computation of the product of the powers of all lagged X_t according to pol
    X_product = 1._realkind    
    DO time=t_start, T
       i = 1
       DO nb_variable=2, pol_dim
          DO deg_variable=1, pol_deg(nb_variable)
             DO j=1, nb_jumps(nb_variable)
                DO s=1, nb_rep(nb_variable)
                   X_product(time,nb_rep(nb_variable)*deg_variable + (j-1)*high_jumps(nb_variable) + s) = &
                        X_product(time,nb_rep(nb_variable)*deg_variable + (j-1)*high_jumps(nb_variable) + s) * &
                        X(time,i)                        
                END DO
             END DO
             i = i+1
          END DO
       END DO
    END DO
    !end computation of the product of the powers of all lagged X_t accordind to pol

    RETURN
  END SUBROUTINE X_MATRIX_CONSTRUCTION


  SUBROUTINE ARCH(n,parameters,sigma_2,U_square,T,dim_mu,dim_arch)
    USE S1_COMMON, ONLY: T_max
    USE Z_MATRIX,  ONLY: sig2_start 
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, T, dim_mu, dim_arch
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: parameters
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(sig2_start:T_max) :: sigma_2
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T_max) :: U_square

    INTEGER :: i,j

    sigma_2(1:T) = parameters(1)
    DO j=1, dim_arch
       DO i=dim_mu+dim_arch, T
          sigma_2(i) = sigma_2(i) + parameters(j+1)*U_square(i-j) 
       END DO
    END DO

    RETURN
  END SUBROUTINE ARCH


  SUBROUTINE GARCH(dim_sigma,parameters,sigma_2,U_square,T,dim_mu,dim_arch,dim_garch)
    USE S1_COMMON, ONLY: T_max
    USE Z_MATRIX,  ONLY: sig2_start, U_var_uncond
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: dim_sigma,T,dim_mu,dim_arch,dim_garch  
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:dim_sigma) :: parameters
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(sig2_start:T_max) :: sigma_2
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T_max) :: U_square

    INTEGER :: i,j

    DO i=dim_mu+dim_arch-dim_garch,dim_mu+dim_arch-1
       sigma_2(i) = U_var_uncond 
    END DO

    DO i=dim_mu+dim_arch, T
       sigma_2(i) = parameters(1)
    END DO

    DO j=2, dim_arch+1
       DO i=dim_mu+dim_arch, T
          sigma_2(i) = sigma_2(i) + parameters(j)*U_square(i-j+1) 
       END DO
    END DO

    DO i=dim_mu+dim_arch, T
       DO j=1, dim_garch
          sigma_2(i) = sigma_2(i) + parameters(1+dim_arch+j)*sigma_2(i-j)
       END DO
    END DO

    RETURN
  END SUBROUTINE GARCH


  SUBROUTINE LIKELIHOOD_UPDATE(beta_restricted,T,working_series,which_series,mode)
    USE S1_COMMON, ONLY: dim_mu, dim_sigma, dim_arch, dim_garch, AM_nb_par, &
         AM_nb_free_par, t_start, pol_dim, pol_deg, nb_pol_coef, &
         deg_z_max, epsilon_0, M, derivative_level, write_output_3, dev_unit_nb
    USE Z_MATRIX
    USE UTILITIES, ONLY: COMPLETE_BETA_RESTRICTED
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_free_par) :: beta_restricted
    INTEGER, INTENT(IN) :: T
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T) :: working_series
    INTEGER, INTENT(IN) :: which_series
    INTEGER, INTENT(IN OUT) :: mode

    LOGICAL :: modify, modify_U_temp, modify_sig2
    INTEGER :: time, i, j, k, s, store_mode
    REAL (KIND=realkind), DIMENSION(1:AM_nb_par) :: beta
    REAL (KIND=realkind), DIMENSION(t_start:T) :: sum_rec

    CALL COMPLETE_BETA_RESTRICTED(beta_restricted,beta)

    !IF ((write_output_3) .AND. (which_series .EQ. 1)) THEN
    !   WRITE(UNIT=dev_unit_nb,FMT=*) "chiamata a LIKELIHOOD_UPDATE"
    !END IF

    modify = .FALSE.
    modify_U_temp = .FALSE.
    modify_sig2 = .FALSE.

    !store the initial value of mode if the subroutine is used with which_series 2 or 3
    IF (which_series .NE. 1) THEN
       store_mode = mode
       mode = 1
    END IF

    IF (which_series .NE. series_old) THEN 
       modify = .TRUE.
       modify_U_temp = .TRUE.
       modify_sig2 = .TRUE.
       series_old = which_series
       beta_old(:) = beta(:)  

       IF (which_series .EQ. 1) THEN
          is_new1 = 0
       END IF
       IF (which_series .EQ. 2) THEN
          is_new2 = 0
       END IF
       IF (which_series .EQ. 3) THEN
          is_new3 = 0
       END IF
    ELSE
       IF (which_series .EQ. 1) THEN
          IF (is_new1 .EQ. 1) THEN
             modify = .TRUE.
             modify_U_temp = .TRUE.
             modify_sig2 = .TRUE.
             is_new1 = 0
             beta_old(:) = beta(:)
          END IF
       END IF
       IF (which_series .EQ. 2) THEN
          IF (is_new2 .EQ. 1) THEN
             modify = .TRUE.
             modify_U_temp = .TRUE.
             modify_sig2 = .TRUE.
             is_new2 = 0
             beta_old(:) = beta(:)
          END IF
       END IF
       IF (which_series .EQ. 3) THEN
          IF (is_new3 .EQ. 1) THEN
             modify = .TRUE.
             modify_U_temp = .TRUE.
             modify_sig2 = .TRUE.
             is_new3 = 0
             beta_old(:) = beta(:)
          END IF
       END IF
    END IF

    !if some modification is needed at this point, this means that we have to recalculate X
    IF (modify .AND. (pol_dim .GT. 1)) THEN
       CALL X_MATRIX_CONSTRUCTION(T, working_series)
    END IF

    IF (.NOT. modify) THEN
       DO i=1, dim_mu
          IF (beta(i) .NE. beta_old(i)) THEN
             modify = .TRUE.
             modify_U_temp = .TRUE.
             modify_sig2 = .TRUE.
             beta_old(:) = beta(:)
             EXIT
          END IF
       END DO
    END IF

    IF (modify_U_temp) THEN
       !Modifica il vettore U_temp
       DO j=1, T
          U_temp(j) = working_series(j) - beta(1) 
       END DO
       DO j=2, dim_mu 
          DO k=dim_mu, T
             U_temp(k) = U_temp(k) - beta(j)*working_series(k-j+1)
          END DO
       END DO

       !now compute the square of U_temp e U_var_uncond da usare nel garch
       !quale valore di partenza
       U_var_uncond = 0._realkind
       DO k=dim_mu, T
          U_square(k) = U_temp(k)**2
          U_var_uncond = U_var_uncond + U_square(k)
       END DO
       U_var_uncond = U_var_uncond/(T-dim_mu+1)
    END IF

    IF (.NOT. modify_U_temp) THEN
       DO i=1, dim_sigma
          IF (beta(dim_mu+i) .NE. beta_old(dim_mu+i)) THEN
             modify = .TRUE.
             modify_sig2 = .TRUE.
             beta_old(dim_mu+1:dim_mu+dim_sigma) = beta(dim_mu+1:dim_mu+dim_sigma)
             EXIT
          END IF
       END DO
    END IF

    IF ((modify_U_temp) .OR. (modify_sig2)) THEN
       IF (dim_sigma .EQ. 1) THEN
          DO j=1, T
             sig2(j) = beta(dim_mu+1)
             inv_sigma(j) = 1._realkind/SQRT(sig2(j))
          END DO
       ELSE
          IF (dim_garch .EQ. 0) THEN
             CALL ARCH(dim_sigma,beta(dim_mu+1:dim_mu+dim_sigma),sig2,U_square,T,dim_mu,dim_arch)
          ELSE
             CALL GARCH(dim_sigma,beta(dim_mu+1:dim_mu+dim_sigma),sig2,U_square,T,dim_mu,dim_arch,dim_garch)
          END IF
          DO j=dim_mu+dim_arch, T
             inv_sigma(j) = 1._realkind/SQRT(sig2(j))
          END DO
       END IF

       !remember that in dim_mu is included the constant term of the condit. mean
       DO i=dim_mu+dim_arch, T
          S_Z(i,1) = U_temp(i)*inv_sigma(i)
       END DO
       IF (deg_z_max .GT. 1) THEN
          !remember that in dim_mu is included the constant term of the condit. mean 
          DO i=2, deg_z_max
             DO j=dim_mu+dim_arch, T         
                S_Z(j,i)=S_Z(j,i-1)*S_Z(j,1)
             END DO
          END DO
       END IF
    END IF

    IF (pol_dim .EQ. 1) THEN
       !computation of the constant denumerator (because pol_dim .EQ. 1)
       constant_denumerator = epsilon_0
       DO j=1, pol_deg(1)+1
          constant_denumerator = constant_denumerator + M(j,j)*beta(dim_mu+dim_sigma+j)**2
          DO i=j+1, pol_deg(1)+1
             constant_denumerator = constant_denumerator + 2*M(i,j)*beta(dim_mu+dim_sigma+j) & 
                  *beta(dim_mu+dim_sigma+i)
          END DO
       END DO
       !end computation of the constant denumerator (because pol_dim .EQ. 1)

       !computation of the polynomial pol and of pol^2 + epsilon_0
       numerator(t_start:T) = beta(dim_mu+dim_sigma+1)
       DO j=1, deg_z_max
          DO time=t_start,T
             numerator(time) = numerator(time) + beta(dim_mu+dim_sigma+j+1)*S_Z(time,j)
          END DO
       END DO
       DO time=t_start,T
          pol(time)=numerator(time)
          numerator(time) = numerator(time)**2 + epsilon_0
       END DO
       !end computation of the polynomial pol and of pol^2 + epsilon_0
    END IF !(pol_dim .EQ. 1)

    IF (pol_dim .GT. 1) THEN
       !compute the necessary variables to compute the log-likelihood
       !Computation of the vector b used to perform the quadratic form b'*M*b
       b = 0._realkind
       DO i=1, pol_deg(1)+1
          DO j=(i-1)*nb_rep(1)+1, i*nb_rep(1)
             DO time=t_start, T
                b(time,i) = b(time,i) + beta(dim_mu+dim_sigma+j)*X_product(time,j)
             END DO
          END DO
       END DO
       !end computation of the vector b used to perform the quadratic form b'*M*b

       !computation of the numerator of the conditional density function
       DO time=t_start,T
          numerator(time) = b(time,1)
       END DO

       DO i=1, pol_deg(1)
          DO time=t_start,T
             numerator(time) = numerator(time) + b(time,i+1)*S_Z(time,i)
          END DO
       END DO
       DO time=t_start,T
          pol(time) = numerator(time)
          numerator(time) = numerator(time)**2 + epsilon_0
       END DO
       !end computation of the numerator of the conditional density function

       !Computation of the denumerator of the conditional density function      
       denumerator = 0._realkind
       DO j=1, pol_deg(1)+1
          DO time=t_start,T
             denumerator(time) = denumerator(time) + M(j,j)*b(time,j)**2
          END DO
          DO i=j+1, pol_deg(1)+1
             DO time=t_start,T
                denumerator(time) = denumerator(time) + 2*M(i,j)*b(time,i)*b(time,j)
             END DO
          END DO
       END DO
       DO time=t_start, T
          denumerator(time) = denumerator(time) + epsilon_0
       END DO
       !end computation of the denumerator of the conditional density function

       !end compute the necessary variables to compute the log-likelihood
    END IF   !end IF (pol_dim .GT. 1)

    IF ((mode .GT. 0) .AND. (derivative_level .EQ. 1)) THEN
       !computation of the gradient of the normal part!
       !computation of the gradient of U2 with respect to beta
       !       IF (modify_U_temp) THEN
       grad_beta_U2 = 0._realkind
       sum_grad_beta_U2 = 0._realkind
       DO time=dim_mu,T
          grad_beta_U2(time,1)=-2._realkind*U_temp(time)
       END DO
       DO j=2,dim_mu
          DO time=dim_mu,T       
             grad_beta_U2(time,j)=grad_beta_U2(time,1)*working_series(time-j+1)
          END DO
       END DO

       DO j=1, dim_mu
          DO time=dim_mu,T
             sum_grad_beta_U2(j) = sum_grad_beta_U2(j) + &
                  grad_beta_U2(time,j)
          END DO
       END DO
       sum_grad_beta_U2(1:dim_mu) = sum_grad_beta_U2(1:dim_mu)/(T-dim_mu+1)
       !END IF
       !end computation of the gradient of U2 with respect to beta  

       ! IF ((modify_U_temp) .OR. (modify_sig2)) THEN
       !computation of the gradient of sig2 with respect to beta
       grad_beta_sig2 = 0._realkind
       DO j=1, dim_arch
          DO s=1, dim_mu
             DO time=dim_mu+dim_arch,dim_mu+dim_arch+dim_garch
                grad_beta_sig2(time,s)=grad_beta_sig2(time,s) + &
                     beta(dim_mu+1+j)*grad_beta_U2(time-j,s)
             END DO
          END DO
       END DO
       DO time=dim_mu+dim_arch,dim_mu+dim_arch+dim_garch
          DO j=1, MIN(dim_garch,time-dim_mu-dim_arch)
             grad_beta_sig2(time,1:dim_mu)=grad_beta_sig2(time,1:dim_mu) + &
                  beta(dim_mu+1+dim_arch+j)*grad_beta_sig2(time-j,1:dim_mu)
          END DO

          DO j=MIN(dim_garch,time-dim_mu-dim_arch)+1,dim_garch
             grad_beta_sig2(time,1:dim_mu)=grad_beta_sig2(time,1:dim_mu) + &
                  beta(dim_mu+1+dim_arch+j)*sum_grad_beta_U2(1:dim_mu)
          END DO
       END DO

       DO j=1, dim_arch
          DO s=1, dim_mu
             DO time=dim_mu+dim_arch+dim_garch+1,T
                grad_beta_sig2(time,s)=grad_beta_sig2(time,s) + &
                     beta(dim_mu+1+j)*grad_beta_U2(time-j,s)
             END DO
          END DO
       END DO

       DO time=dim_mu+dim_arch+dim_garch+1,T
          DO j=1, dim_garch
             grad_beta_sig2(time,1:dim_mu)=grad_beta_sig2(time,1:dim_mu) + &
                  beta(dim_mu+1+dim_arch+j)*grad_beta_sig2(time-j,1:dim_mu)
          END DO
       END DO
       !end of computation of the gradient of sig2 with respect to beta

       !computation of the gradient of sig2 with respect to the param. of the ARCH/GARCH
       grad_sigma_sig2(:,1) = 1._realkind
       DO j=1, dim_arch
          DO time=dim_mu+dim_arch,dim_mu+dim_arch+dim_garch
             grad_sigma_sig2(time,1+j) = U_square(time-j) 
          END DO
       END DO

       DO time=dim_mu+dim_arch,dim_mu+dim_arch+dim_garch
          DO j=1, MIN(dim_garch,time-dim_mu-dim_arch)
             grad_sigma_sig2(time,1+dim_arch+j) = sig2(time-j)
          END DO

          DO j=MIN(dim_garch,time-dim_mu-dim_arch)+1,dim_garch
             grad_sigma_sig2(time,1+dim_arch+j) = U_var_uncond 
          END DO
       END DO

       DO time=dim_mu+dim_arch,dim_mu+dim_arch+dim_garch
          DO j=1, MIN(dim_garch,time-dim_mu-dim_arch)
             grad_sigma_sig2(time,1:dim_sigma) = grad_sigma_sig2(time,1:dim_sigma)+ & 
                  beta(dim_mu+1+dim_arch+j)*grad_sigma_sig2(time-j,1:dim_sigma)
          END DO

          !DO j=MIN(dim_garch,i-dim_mu-dim_arch)+1,dim_garch
          !   grad_sigma_sig2(i,1:dim_sigma) = grad_sigma_sig2(i,1:dim_sigma) + & 
          !        0._realkind
          !END DO
       END DO

       DO j=1, dim_arch
          DO time=dim_mu+dim_arch+dim_garch+1,T 
             grad_sigma_sig2(time,1+j) = U_square(time-j)
          END DO
       END DO

       DO j=1, dim_garch
          DO time=dim_mu+dim_arch+dim_garch+1,T
             grad_sigma_sig2(time,1+dim_arch+j) = sig2(time-j)
          END DO
       END DO

       DO time=dim_mu+dim_arch+dim_garch+1,T
          DO j=1, dim_garch
             grad_sigma_sig2(time,1:dim_sigma) = grad_sigma_sig2(time,1:dim_sigma) + &
                  beta(dim_mu+1+dim_arch+j)*grad_sigma_sig2(time-j,1:dim_sigma)
          END DO
       END DO
       !end computation of the gradient of sig2 with respect to the param. of the ARCH/GARCH


       !computation of the costant common to all components of the second summand of the gradient
       DO time=t_start,T
          grad_const(time) = 1._realkind-U_square(time)/sig2(time)
       END DO
       !end computation of the costant common to all components of the second summand of the gradient
       !END IF

       grad_l = 0._realkind
       DO j=1, dim_mu
          DO time=t_start,T
             grad_l(time,j)=(-0.5_realkind)/sig2(time)*(grad_beta_U2(time,j) + &
                  grad_const(time)*grad_beta_sig2(time,j))
          END DO
       END DO
       DO j=1, dim_sigma
          DO time=t_start,T
             grad_l(time,dim_mu+j) =(-0.5_realkind)/sig2(time)*grad_const(time) * &
                  grad_sigma_sig2(time,j)
          END DO
       END DO

       IF (pol_dim .EQ. 1) THEN
          !compute the gradient of the fix denumerator:
          !-log(theta'*M*theta+epsilon_0)^(-1)
          grad_theta_denumerator(1:nb_pol_coef) = 0._realkind
          DO i=1, nb_pol_coef
             DO j=1, nb_pol_coef
                grad_theta_denumerator(j) = grad_theta_denumerator(j) + & 
                     M(j,i)*beta(dim_mu+dim_sigma+i)
             END DO
          END DO
          grad_theta_denumerator(1:nb_pol_coef) = -2._realkind / constant_denumerator * &
               grad_theta_denumerator(1:nb_pol_coef)
          !end compute the gradient of the fixed denumerator

          !compute the constant common to all components of the gradient of
          !log(pol^2+epsilon_0)
          DO time=t_start, T
             fixed_grad_P(time) = 2._realkind*pol(time)/numerator(time)
          END DO
          !end compute the constant common to all components of ...


          !computation of the gradient of Z with respect to beta
          DO i=dim_mu+dim_arch,T
             grad_beta_Z(i,1)=-1._realkind
          END DO
          DO j=2,dim_mu
             DO i=dim_mu+dim_arch,T
                grad_beta_Z(i,j)=-working_series(i-j+1)
             END DO
          END DO
          DO j=1, dim_mu
             DO i=dim_mu+dim_arch,T
                grad_beta_Z(i,j) = grad_beta_Z(i,j)*inv_sigma(i) + &
                     (-0.5_realkind)*S_Z(i,1)/sig2(i)*grad_beta_sig2(i,j)
             END DO
          END DO
          !end computation of the gradient of Z with respect to beta

          !compute the derivative of pol with respect to beta and the param. of sigma
          sum_rec = beta(dim_mu+dim_sigma+2)
          DO j=2, pol_deg(1)
             DO time=t_start, T
                sum_rec(time) = sum_rec(time) + beta(dim_mu+dim_sigma+1+j)*j*S_Z(time,j-1) 
             END DO
          END DO

          DO j=1,dim_mu
             DO time=t_start, T
                grad_beta_pol(time,j)= sum_rec(time) * fixed_grad_P(time)* &
                     grad_beta_Z(time,j)
             END DO
          END DO
          DO j=1,dim_sigma
             DO time=t_start, T
                grad_sigma_pol(time,j)= sum_rec(time)*fixed_grad_P(time)*(-0.5_realkind)* &
                     S_Z(time,1)/sig2(time)*grad_sigma_sig2(time,j)
             END DO
          END DO
          !end compute the derivative of pol with respect to beta and the param. of sigma

          !compute the derivative with respect to the polynomial param.
          DO time=t_start, T
             grad_theta_pol(time,1)=1._realkind
          END DO

          DO j=1, deg_Z_max
             DO time=t_start, T
                grad_theta_pol(time,j+1)=S_Z(time,j)
             END DO
          END DO

          DO j=1, nb_pol_coef
             DO time=t_start, T
                grad_theta_pol(time,j)=grad_theta_pol(time,j) * &
                     fixed_grad_P(time)
             END DO
          END DO
          !end compute the derivative with respect to the polynomial param.

          !add the results to the gradient of the likelihood
          DO j=1, dim_mu
             DO time=t_start,T
                grad_l(time,j) = grad_l(time,j) + grad_beta_pol(time,j)
             END DO
          END DO

          DO j=1, dim_sigma
             DO time=t_start,T
                grad_l(time,dim_mu+j) = grad_l(time,dim_mu+j) + &
                     grad_sigma_pol(time,j)
             END DO
          END DO
          DO j=1, nb_pol_coef
             DO time=t_start,T
                grad_l(time,dim_mu+dim_sigma+j) = grad_theta_pol(time,j) + &
                     grad_theta_denumerator(j)
             END DO
          END DO
          !end add the results to the gradient of the likelihood
       END IF    !with respect to IF (pol_dim .EQ. 1)

       IF (pol_dim .GT. 1) THEN
          !compute the exact gradient of the likelihood if necessary
          !compute the constant common to all components of the gradient of
          !log(pol^2+epsilon_0)
          DO time=t_start, T
             fixed_grad_P(time) = 2._realkind*pol(time)/numerator(time)
          END DO
          !end compute the constant common to all components of ...

          !computation of the gradient of Z with respect to beta
          DO time=dim_mu+dim_arch,T
             grad_beta_Z(time,1)=-1._realkind 
          END DO

          DO j=2,dim_mu
             DO time=dim_mu+dim_arch,T       
                grad_beta_Z(time,j)=-working_series(time-j+1)
             END DO
          END DO
          DO j=1, dim_mu
             DO time=dim_mu+dim_arch,T  
                grad_beta_Z(time,j) = grad_beta_Z(time,j)*inv_sigma(time) + &
                     (-0.5_realkind)*S_Z(time,1)/sig2(time)*grad_beta_sig2(time,j)
             END DO
          END DO
          !end computation of the gradient of Z with respect to beta

          !compute the derivative of pol with respect to beta and the param. of sigma
          DO time=t_start, T
             sum_rec(time) = b(time,2)          
          END DO
          DO j=2, pol_deg(1)
             DO time=t_start, T
                sum_rec(time) = sum_rec(time) + b(time,1+j)*j*S_Z(time,j-1) 
             END DO
          END DO

          DO j=1, dim_mu 
             DO time=t_start, T
                grad_beta_pol(time,j)= sum_rec(time) * fixed_grad_P(time)* &
                     grad_beta_Z(time,j)
             END DO
          END DO

          DO j=1, dim_sigma 
             DO time=t_start, T
                grad_sigma_pol(time,j)= sum_rec(time)*fixed_grad_P(time)*(-0.5_realkind)* &
                     S_Z(time,1)/sig2(time)*grad_sigma_sig2(time,j)
             END DO
          END DO
          !end compute the derivative of pol with respect to beta and the param. of sigma

          !compute the derivative with respect to the polynomial param.
          DO i=1, nb_rep(1)
             DO time=t_start, T
                grad_theta_pol(time,i)= X_product(time,i)
             END DO
          END DO

          DO i=2, pol_deg(1)+1
             DO j=(i-1)*nb_rep(1)+1, i*nb_rep(1)
                DO time=t_start, T
                   grad_theta_pol(time,j)= S_Z(time,i-1)*X_product(time,j)
                END DO
             END DO
          END DO

          DO j=1,nb_pol_coef
             DO time=t_start, T
                grad_theta_pol(time,j)=grad_theta_pol(time,j) * &
                     fixed_grad_P(time)
             END DO
          END DO
          !end compute the derivative with respect to the polynomial param.

          !compute the gradient of the denumerator:
          !-log(b'*M*b+epsilon_0)
          grad_theta_denumerator = 0._realkind
          grad_theta_qMq = 0._realkind
          DO time=t_start,T
             DO i=1, pol_deg(1)+1
                sum_rec(time)=0._realkind
                DO j=1, pol_deg(1)+1
                   sum_rec(time) = sum_rec(time) + M(j,i)*b(time,j)
                END DO

                DO j=(i-1)*nb_rep(1)+1, i*nb_rep(1)
                   grad_theta_qMq(time,j) = grad_theta_qMq(time,j) + &
                        2._realkind*sum_rec(time)*X_product(time,j)
                END DO
             END DO
          END DO
          !end compute the gradient of the denumerator


          !add the results to the gradient of the likelihood
          DO j=1, dim_mu
             DO time=t_start,T
                grad_l(time,j) = grad_l(time,j) + grad_beta_pol(time,j)
             END DO
          END DO

          DO j=1, dim_sigma
             DO time=t_start,T
                grad_l(time,dim_mu+j) = grad_l(time,dim_mu+j) + &
                     grad_sigma_pol(time,j)
             END DO
          END DO
          DO j=1, nb_pol_coef
             DO time=t_start,T
                grad_l(time,dim_mu+dim_sigma+j) = grad_theta_pol(time,j) - &
                     grad_theta_qMq(time,j)/denumerator(time)
             END DO
          END DO
          !end add the results to the gradient of the likelihood
       END IF       !end IF (pol_dim .GT. 1)
    END IF          !end IF ((mode .GT. 0) .AND. (derivative_level .EQ. 1))

    !restore the initial value of mode
    IF (which_series .NE. 1) THEN
       mode = store_mode
    END IF

    RETURN
  END SUBROUTINE LIKELIHOOD_UPDATE


  SUBROUTINE SNP_ORT_FUNC(T)
    USE S1_COMMON, ONLY: grad_i, t_start, dim_mu, nb_not_fixed_mu, dim_sigma, &
         nb_not_fixed_sigma, nb_not_fixed_pol , j1, j2, j3, &
         dev_unit_nb
      !rimuovi
    USE Z_MATRIX,  ONLY: grad_l
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: T

    INTEGER :: time, j

    !This subroutine copies the gradient elements of the vector
    !grad_l into the vector grad_i. T represents the last element in the grad_l vector to be
    !copied.

    DO j=1, nb_not_fixed_mu 
       DO time=t_start,T
          grad_i(time-t_start+1,j) = grad_l(time,j1(j))
       END DO
    END DO
    DO j=1, nb_not_fixed_sigma
       DO time=t_start,T 
          grad_i(time-t_start+1,nb_not_fixed_mu+j) = grad_l(time,dim_mu+j2(j))
       END DO
    END DO

    DO j=1, nb_not_fixed_pol 
       DO time=t_start,T 
          grad_i(time-t_start+1,nb_not_fixed_mu+nb_not_fixed_sigma+j) = &
               grad_l(time,dim_mu+dim_sigma+j3(j))
       END DO
    END DO

    RETURN
  END SUBROUTINE SNP_ORT_FUNC


  SUBROUTINE GRADIENT_STEP_1(mode,n,g,T)
    USE S1_COMMON, ONLY: t_start, dim_mu, nb_not_fixed_mu, dim_sigma, &
         nb_not_fixed_sigma, nb_not_fixed_pol , j1, j2, j3, derivative_level, &
         is_restricted_par_AM
    USE Z_MATRIX, ONLY: grad_l
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: mode, n, T
    REAL (KIND=realkind), INTENT(OUT), DIMENSION(1:n) :: g

    INTEGER :: i,k,s,v,time

    IF ((mode .GT. 0) .AND. (derivative_level .EQ. 1)) THEN
       g = 0._realkind
       DO i=1, nb_not_fixed_mu
          k=j1(i)
          IF (.NOT. is_restricted_par_AM(k)) THEN
             DO time=t_start,T
                g(i) = g(i) - grad_l(time,k)
             END DO
          END IF
       END DO
       DO i=1, nb_not_fixed_sigma
          k=j2(i)
          IF (.NOT. is_restricted_par_AM(dim_mu+k)) THEN
             s = nb_not_fixed_mu+i
             v = dim_mu+k
             DO time=t_start,T
                g(s) = g(s) - grad_l(time,v)
             END DO
          END IF
       END DO
       DO i=1, nb_not_fixed_pol 
          k=j3(i)
          IF (.NOT. is_restricted_par_AM(dim_mu+dim_sigma+k)) THEN
             s = nb_not_fixed_mu+nb_not_fixed_sigma+i
             v = dim_mu+dim_sigma+k
             DO time=t_start,T
                g(s) = g(s) - grad_l(time,v)
             END DO
          END IF
       END DO

       g = g/REAL(T-t_start+1)

    END IF
    RETURN
  END SUBROUTINE GRADIENT_STEP_1

  SUBROUTINE CDF_CONSTRUCTION(T,likelihood)
    USE S1_COMMON, ONLY: t_start, normal_c, pol_dim
    USE Z_MATRIX,  ONLY: inv_sigma, numerator, constant_denumerator, &
         denumerator, S_Z 
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: T
    REAL (KIND=realkind), INTENT(OUT) :: likelihood

    INTEGER :: time
    REAL (KIND=realkind), DIMENSION(1:T) :: log_cond_density
    !write(*,*) "inizio subroutine CDF_CONSTRUCTION"

    IF (pol_dim .EQ. 0) THEN
       DO time=t_start, T
          log_cond_density(time) = LOG(normal_c*inv_sigma(time)) - 0.5_realkind*S_Z(time,1)**2
       END DO
    ELSE
       IF (pol_dim .EQ. 1) THEN
          DO time=t_start,T
             log_cond_density(time) = LOG(numerator(time)*normal_c*inv_sigma(time)/constant_denumerator) &
                  - 0.5_realkind*S_Z(time,1)**2
          END DO
       ELSE
          DO time=t_start,T
             log_cond_density(time) = LOG(numerator(time)*normal_c*inv_sigma(time)/denumerator(time)) &
                  - 0.5_realkind*S_Z(time,1)**2
          END DO
       END IF
    END IF

    likelihood = 0.0_realkind

    DO time=t_start, T 
       likelihood = likelihood - log_cond_density(time)
    END DO
    likelihood = likelihood/REAL(T-t_start+1)

    RETURN
  END SUBROUTINE CDF_CONSTRUCTION
END MODULE LIKELIHOOD_UTILITIES


MODULE COVARIANCE_ESTIMATES
  USE DATATYPES
  IMPLICIT NONE

CONTAINS
  SUBROUTINE EST_COV_SCORE_AM_STEP1(beta_restricted,T,working_series,which_series)
    USE S1_COMMON,   ONLY: AM_nb_not_fixed_par, score_cov, sel_score_cov_mat_s1, lags_score_cov_mat_s1,&
         t_start, grad_i, derivative_level, write_output_3, dev_unit_nb
    USE Z_MATRIX
    USE MYFUNCTIONS, ONLY:  MAT_CROSS_LAG, EYE, INV, MM
    USE LIKELIHOOD_UTILITIES, ONLY:  SNP_ORT_FUNC, LIKELIHOOD_UPDATE
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_not_fixed_par) :: beta_restricted
    INTEGER, INTENT(IN) :: T
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T) :: working_series
    INTEGER, INTENT(IN) :: which_series

    !this subroutine computes the covariance matrix estimate of the score of the
    !auxiliary model. As input we need the vector of free model's parameters.
    INTEGER :: j, T_tmp

    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par,1:AM_nb_not_fixed_par) :: gamma
    REAL (KIND=realkind), DIMENSION(1:T-t_start+1,1:AM_nb_not_fixed_par) :: H


    T_tmp = T-t_start+1
    IF (derivative_level .NE. 1) THEN
       WRITE(*,*) "Impossible to compute the weighting matrix for the quadratic form h'Sh: exact score needed!"
       WRITE(*,*) "The weighting matrix is set equal to the identity matrix!"
       score_cov=EYE(AM_nb_not_fixed_par)
    ELSE
       j = 1
       CALL LIKELIHOOD_UPDATE(beta_restricted,T,working_series,which_series,j)
       !compute the matrix grad_i for the orthogonality cond.
       CALL SNP_ORT_FUNC(T)

       H(1:T_tmp,:) = grad_i(1:T_tmp,:)

       !estimate the covariance matrix of the score
       score_cov=MM(TRANSPOSE(H),H,AM_nb_not_fixed_par,T_tmp,AM_nb_not_fixed_par)

       IF (sel_score_cov_mat_s1 .EQ. 2) THEN
          !use the Newey-West method
          j=1
          DO
             IF (j .GT. lags_score_cov_mat_s1) EXIT
             gamma = MAT_CROSS_LAG(H,T_tmp,AM_nb_not_fixed_par,j)
             score_cov = score_cov+(1.0_realkind-REAL(j)/REAL(lags_score_cov_mat_s1+1)) * &
                  (gamma+TRANSPOSE(gamma))
             j=j+1
          END DO
          !end use the Newey-West method
       END IF
       IF (sel_score_cov_mat_s1 .EQ. 3) THEN
          !use the Gallant method
       END IF
       score_cov = score_cov/REAL(T_tmp)
    END IF !(derivative_level)

    RETURN
  END SUBROUTINE EST_COV_SCORE_AM_STEP1

  SUBROUTINE SUM_AND_UPDATE_PSI_STEP_1(beta_restricted, T, working_series, which_series, T_tmp, H)
    USE S1_COMMON,   ONLY: AM_nb_free_par, AM_nb_not_fixed_par, grad_i
    USE LIKELIHOOD_UTILITIES, ONLY: SNP_ORT_FUNC, LIKELIHOOD_UPDATE

    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_free_par) :: beta_restricted
    INTEGER, INTENT(IN) :: T
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T) :: working_series
    INTEGER, INTENT(IN) :: which_series
    INTEGER, INTENT(IN) :: T_tmp
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:AM_nb_not_fixed_par) :: H

    INTEGER :: i, j

    !compute the ortogonality function of the modified_series at the beta
    !and save the grad_i in psi_matrix_tmp for later use

    i = 1
    CALL LIKELIHOOD_UPDATE(beta_restricted,T,working_series,which_series,i)
    !compute the matrix grad_i for the orthogonality cond.
    CALL SNP_ORT_FUNC(T)

    H = 0._realkind
    DO j=1, AM_nb_not_fixed_par
       DO i=1, T_tmp
          H(j) = H(j) + grad_i(i,j)
       END DO
    END DO

    H = H/REAL(T_tmp)
    
    RETURN
  END SUBROUTINE SUM_AND_UPDATE_PSI_STEP_1


  SUBROUTINE APPROX_M_BETA(beta_restricted, M_beta, T, working_series, which_series, T_tmp)
    USE DATATYPES
    USE S1_COMMON,   ONLY: AM_nb_par, AM_nb_not_fixed_par, AM_nb_free_par, grad_i, &
         write_output_3, dev_unit_nb, delta_M_beta
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_free_par) :: beta_restricted
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:AM_nb_not_fixed_par,1:AM_nb_free_par) :: M_beta
    INTEGER, INTENT(IN) :: T
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T) :: working_series
    INTEGER, INTENT(IN) :: which_series
    INTEGER, INTENT(IN) :: T_tmp

    REAL (KIND=realkind), DIMENSION(1:AM_nb_free_par)      :: beta_restricted_1
    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par) :: H_1
    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par) :: H_2

    INTEGER :: i_AM_nb_free_par

100 FORMAT(25(ES10.3,TR1))

    IF (write_output_3) THEN
       WRITE(UNIT=dev_unit_nb,FMT=*) "calcolo di M_beta"
       WRITE(UNIT=dev_unit_nb,FMT=*) " delta:", delta_M_beta      
    END IF

    !loop over all AM_nb_free_par (compute the directional derivatives)
    DO i_AM_nb_free_par=1, AM_nb_free_par

       IF (write_output_3) THEN
          WRITE(UNIT=dev_unit_nb,FMT=*) "chiamata a APPROX_M_BETA, coefficiente nr. ", i_AM_nb_free_par
          WRITE(UNIT=dev_unit_nb,FMT=*) "delta +"
       END IF

       beta_restricted_1 = beta_restricted
       beta_restricted_1(i_AM_nb_free_par) = beta_restricted_1(i_AM_nb_free_par) + delta_M_beta

       CALL SUM_AND_UPDATE_PSI_STEP_1(beta_restricted_1, T, working_series, which_series, T_tmp, H_1)

       !simulate the series
       IF (write_output_3) THEN
          WRITE(UNIT=dev_unit_nb,FMT=*) "delta -"       
       END IF

       beta_restricted_1 = beta_restricted
       beta_restricted_1(i_AM_nb_free_par) = beta_restricted_1(i_AM_nb_free_par) - delta_M_beta

       CALL SUM_AND_UPDATE_PSI_STEP_1(beta_restricted_1, T, working_series, which_series, T_tmp, H_2)

       M_beta(:,i_AM_nb_free_par) = (H_1(:) - H_2(:) ) / (2._realkind * delta_M_beta)

       IF (write_output_3) THEN
          WRITE(UNIT=dev_unit_nb,FMT=*) "derivative:"
          WRITE(UNIT=dev_unit_nb,FMT=100) M_beta(:,i_AM_nb_free_par)
          WRITE(UNIT=dev_unit_nb,FMT=*) 
       END IF

    END DO  !i_AM_nb_free_par

       IF (AM_nb_free_par .EQ. AM_nb_not_fixed_par) THEN
            M_beta = (M_beta + TRANSPOSE(M_beta))*0.5_realkind
       END IF
            RETURN
  END SUBROUTINE APPROX_M_BETA


  SUBROUTINE EST_COV_AM_FREE_PAR_STEP1(beta_restricted,T,working_series,which_series)
    USE DATATYPES
    USE S1_COMMON,   ONLY: AM_nb_par, AM_nb_not_fixed_par, AM_nb_free_par, score_cov, sel_hessian_mat_s1, &
         t_start, grad_i, sel_cov_estimates_s1, cov_estimates_s1, derivative_level, M_beta, &
         is_restricted_par_AM, is_fixed_par_AM
    USE Z_MATRIX
    USE MYFUNCTIONS, ONLY: INV, MM, EYE
    USE LIKELIHOOD_UTILITIES, ONLY:  SNP_ORT_FUNC, LIKELIHOOD_UPDATE
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_not_fixed_par) :: beta_restricted
    INTEGER, INTENT(IN) :: T
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T) :: working_series
    INTEGER, INTENT(IN) :: which_series

    INTEGER :: i, j, T_tmp
    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par,1:AM_nb_free_par) :: gamma
    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par,1:AM_nb_not_fixed_par) :: V
    REAL (KIND=realkind), DIMENSION(1,1:AM_nb_not_fixed_par) :: gam_temp

    IF (derivative_level .NE. 1) THEN
       WRITE(*,*) "Impossible to compute the covariance of the AM free parameters."
       WRITE(*,*) "The gradient is not available. Cov = I"
       cov_estimates_s1 = EYE(AM_nb_not_fixed_par)
    END IF !(derivative_level)

    !estimate the covariance matrix of the free parameters in the auxiliary model

    !first we need to estimate the inverse of the Hessian matrix: select the method
    IF (sel_cov_estimates_s1 .EQ. 1) THEN

    CALL EST_COV_SCORE_AM_STEP1(beta_restricted,T,working_series,which_series)

       IF (AM_nb_not_fixed_par .EQ. AM_nb_free_par) THEN

          cov_estimates_s1 = INV(score_cov,AM_nb_free_par)

       ELSE

          V = INV(score_cov,AM_nb_not_fixed_par)
          i=1
          DO j=1, AM_nb_par
             IF ((.NOT. is_fixed_par_AM(j)) .AND. (.NOT. is_restricted_par_AM(j))) THEN
                gamma(:,i) = V(:,j)
                i=i+1
             END IF
          END DO

          cov_estimates_s1 = MM(score_cov,gamma,AM_nb_not_fixed_par,AM_nb_not_fixed_par,AM_nb_free_par)
          cov_estimates_s1 = MM(TRANSPOSE(gamma),cov_estimates_s1,AM_nb_free_par,AM_nb_not_fixed_par,AM_nb_free_par)

       END IF
 
    END IF

    IF (sel_cov_estimates_s1 .EQ. 2) THEN
       T_tmp = T-t_start+1
       IF (sel_hessian_mat_s1 .EQ. 1) THEN

          j = 1
          CALL LIKELIHOOD_UPDATE(beta_restricted,T,working_series,which_series,j)
          !compute the matrix grad_i for the orthogonality cond.
          CALL SNP_ORT_FUNC(T)
          !guarda che qui puoi scrivere piu' velocemente in quanto gamma e' uguale a X'X!! 
          !estimates the inverse of the Hessian using the outproduct of the gradient
          V = 0._realkind
          DO j=1, T_tmp
             gam_temp(1,:) = grad_i(j,:)
             V =  V + MM(TRANSPOSE(gam_temp),gam_temp,AM_nb_not_fixed_par,1,AM_nb_not_fixed_par)
          END DO
          V = V / REAL(T_tmp)

          V = INV(V,AM_nb_not_fixed_par)
          !end stimates the inverse of the Hessian using the outproduct of the gradient

          i=1
          DO j=1, AM_nb_par
             IF ((.NOT. is_fixed_par_AM(j)) .AND. (.NOT. is_restricted_par_AM(j))) THEN
                gamma(:,i) = V(:,j)
                i=i+1
             END IF
          END DO

          !end first we need to estimate the inverse of the Hessian matrix: select the method
          cov_estimates_s1 = MM(score_cov,gamma,AM_nb_not_fixed_par,AM_nb_not_fixed_par,AM_nb_free_par)
          cov_estimates_s1 = MM(gamma,cov_estimates_s1,AM_nb_free_par,AM_nb_not_fixed_par,AM_nb_free_par)

       ELSE

          CALL APPROX_M_BETA(beta_restricted, M_beta, T, working_series, which_series, T_tmp)
          cov_estimates_s1 = MM(INV(score_cov,AM_nb_not_fixed_par),M_beta,AM_nb_not_fixed_par,AM_nb_not_fixed_par,AM_nb_free_par)
          cov_estimates_s1 = MM(TRANSPOSE(M_beta),cov_estimates_s1,AM_nb_free_par,AM_nb_not_fixed_par,AM_nb_free_par)
          cov_estimates_s1 = INV(cov_estimates_s1,AM_nb_free_par)

       END IF

    END IF

    RETURN
  END SUBROUTINE EST_COV_AM_FREE_PAR_STEP1
END MODULE COVARIANCE_ESTIMATES



MODULE OUTPUT_UTILITIES
  USE DATATYPES
  IMPLICIT NONE

CONTAINS

  SUBROUTINE WRITE_EMPTY(unit_nb,nb_lines)
    IMPLICIT NONE
    !this subroutine writes nb_lines empty lines in unit unit_nb
    INTEGER, INTENT(IN) :: unit_nb
    INTEGER, INTENT(IN) :: nb_lines

    INTEGER :: i

    DO i=1,nb_lines
       WRITE(UNIT=unit_nb,FMT=*)
    END DO
    RETURN
  END SUBROUTINE WRITE_EMPTY

  SUBROUTINE NEW_RUN(i)
    USE S1_COMMON
    USE S2_COMMON
    USE ROB_COMMON
    USE DGP_DEFINITION
    USE CONTAMINATION_MODEL
    USE MONTE_CARLO, ONLY: nb_rep_monte_carlo
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: i

    INTEGER :: j
    WRITE(UNIT=i,FMT="(A)") "******************"
    WRITE(UNIT=i,FMT="(A)") "#0 SECTION GENERAL"
    WRITE(UNIT=i,FMT="(A)") "******************"
    WRITE(UNIT=i,FMT=*)

    WRITE(UNIT=i,FMT="(A)") "0.01) Select the kind of estimation: 0=classic, 1=robust:"
    WRITE(UNIT=i,FMT=*) estimation_kind
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "0.02) Select the kind of job: 1=only first step, 2=only second step, 3=both steps:"
    WRITE(UNIT=i,FMT=*) step
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "0.03) Simulate the dataset or use real data: 0=real data, 1=simulate"
    WRITE(UNIT=i,FMT=*) simulate_data
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "0.04) Number of observations:"
    WRITE(UNIT=i,FMT=*) T1
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "0.05) Number of monte carlo repetitions:"
    WRITE(UNIT=i,FMT=*) nb_rep_monte_carlo
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "0.06) Select the output you want, 0=no output or 1=output:"
    WRITE(UNIT=i,FMT=*) write_output_1
    WRITE(UNIT=i,FMT=*) write_output_2
    WRITE(UNIT=i,FMT=*) write_output_3
    CALL WRITE_EMPTY(i,5)





    WRITE(UNIT=i,FMT="(A)") "***************************************************"
    WRITE(UNIT=i,FMT="(A)") "#1 GENERATOR PROCESS OF THE DATA (ONLY IF 0.03 = 1)"
    WRITE(UNIT=i,FMT="(A)") "***************************************************"
    WRITE(UNIT=i,FMT=*)

    WRITE(UNIT=i,FMT="(A)") "1.01) Name of the DGP (A15)"
    WRITE(UNIT=i,FMT="(A)") DGP_name
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "Select between: MOVING__AVERAGE - STOC_VOLATILITY - AR_GARCH_MODEL_"
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "1.02) Number of iid series necessary for the simulation:"
    WRITE(UNIT=i,FMT=*) DGP_nb_iid_series
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "Selection table:"
    WRITE(UNIT=i,FMT="(A)") "MOVING__AVERAGE     1"
    WRITE(UNIT=i,FMT="(A)") "STOC_VOLATILITY     2"
    WRITE(UNIT=i,FMT="(A)") "AR_GARCH_MODEL_     1"
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "1.03) Distribution of the innovations of every iid series (A7 or A7 I for STUDENT):"
    DO j=1, DGP_nb_iid_series
       IF (DGP_innov_dist(j) .EQ. "STUDENT") THEN
          WRITE(UNIT=i,FMT="(A)") DGP_innov_dist(j), DGP_deg_freedom(j)
       ELSE
          WRITE(UNIT=i,FMT="(A)") DGP_innov_dist(j)
       END IF
    END DO
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "Selection table: MOVING__AVERAGE   STOC_VOLATILITY   AR_GARCH_MODEL_"
    WRITE(UNIT=i,FMT="(A)") "UNIFORM ......         YES               YES              YES"
    WRITE(UNIT=i,FMT="(A)") "NORMAL_ ......         YES               YES              YES"
    WRITE(UNIT=i,FMT="(A)") "STUDENT df....         YES               YES              YES"
    WRITE(UNIT=i,FMT="(A)") "SNP_DST ......     WITHOUT LAGS      WITHOUT LAGS         YES"
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "1.04) SNP_DST dimensions (only if the innov. distribution = SNP_DST)"
    WRITE(UNIT=i,FMT=*) DGP_pol_dim, &
         "             -> dimension of the polynomial, i.e. 1 if no lags in the polynomial part (I)"
    IF (DGP_pol_dim .LT.1) THEN  
       WRITE(UNIT=i,FMT=*) &
            "0             -> order of dimension 1 (I)"
       WRITE(UNIT=i,FMT=*) &
            "0             -> order of dimension 2 (I), and so on for every dimension ..."
    ELSE
       IF(DGP_pol_dim .EQ. 1) THEN
          WRITE(UNIT=i,FMT=*) DGP_pol_deg(1), &
               " -> dimension of the polynomial, i.e. 1 if no lags in the polynomial part (I)"
          WRITE(UNIT=i,FMT=*) &
               "0             -> order of dimension 2 (I), and so on for every dimension ..."
       ELSE
          WRITE(UNIT=i,FMT=*) DGP_pol_deg(1), &
               " -> dimension of the polynomial, i.e. 1 if no lags in the polynomial part (I)"
          WRITE(UNIT=i,FMT=*) DGP_pol_deg(2), &
               " -> order of dimension 2 (I), and so on for every dimension ..."
          DO j=3, DGP_pol_dim
             WRITE(UNIT=i,FMT=*) DGP_pol_deg(j)
          END DO
       END IF
    END IF
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "1.05) SNP_DST parameter's values (only if the innov. distribution = SNP_DST)"
    IF (DGP_pol_dim .EQ. 0) THEN
       WRITE(UNIT=i,FMT=*) "0             -> epsilon_0 of the polynomial (R)"
       WRITE(UNIT=i,FMT=*) "0             -> first polynomial coefficient"
       WRITE(UNIT=i,FMT=*) "0             -> second polynomial coefficient, and so on ..."
    ELSE
       WRITE(UNIT=i,FMT=*) DGP_epsilon_0,    "-> epsilon_0 of the polynomial (R)"
       WRITE(UNIT=i,FMT=*) DGP_pol_coef(1),  "-> first polynomial coefficient"
       WRITE(UNIT=i,FMT=*) DGP_pol_coef(2),  "-> second polynomial coefficient, and so on ..."
       DO j=3,DGP_nb_pol_coef 
          WRITE(UNIT=i,FMT=*) DGP_pol_coef(j)
       END DO
    END IF
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "1.06) Number of parts of the DGP model (I):"
    WRITE(UNIT=i,FMT=*) DGP_model_dim
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "Selection table:"
    WRITE(UNIT=i,FMT="(A)") "MOVING__AVERAGE     1   (only MA-part)"
    WRITE(UNIT=i,FMT="(A)") "STOC_VOLATILITY     2   (AR-part for the mean, AR-part for the volatility)"
    WRITE(UNIT=i,FMT="(A)") "AR_GARCH_MODEL_     3   (AR-part for the mean, ARCH-part and GARCH-part)"
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "1.07) Order of every part (I):"
    DO j=1, DGP_model_dim
       WRITE(UNIT=i,FMT=*) DGP_order_model_dim(j)
    END DO
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "Selection table:"
    WRITE(UNIT=i,FMT="(A)") "MOVING__AVERAGE     MA-order"
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "STOC_VOLATILITY     AR-order for the mean function"
    WRITE(UNIT=i,FMT="(A)") "                    AR-order for the var. function"
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "AR_GARCH_MODEL_     AR-order for the mean function"
    WRITE(UNIT=i,FMT="(A)") "                    ARCH-order"
    WRITE(UNIT=i,FMT="(A)") "                    GARCH-order"
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "1.08) Parameters values (R):"
    DO j=1, DGP_nb_par
       WRITE(UNIT=i,FMT=*) DGP_par_vec(j)
    END DO
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") &
         "The parameters must be written in vertical order. If you don't want the constant, then set its value = 0."
    WRITE(UNIT=i,FMT="(A)") "MOVING__AVERAGE     constant, lags(1) ...  lags(n), variance"
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "STOC_VOLATILITY     MEAN: constant, lags(1) ... lags(n1)"
    WRITE(UNIT=i,FMT="(A)") "                    VAR.: constant, lags(1) ... lags(n2), sigma"
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "AR_GARCH_MODEL_     MEAN: constant, lags(1) ... lags(n1)"
    WRITE(UNIT=i,FMT="(A)") "                    ARCH: constant, lags(1) ... lags(n2)"
    WRITE(UNIT=i,FMT="(A)") "                    GARCH: lags(1) ... lags(n3)"
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "1.09) The number of observations to discard in the simulation:"
    WRITE(UNIT=i,FMT=*) DGP_nb_discard
    CALL WRITE_EMPTY(i,5)





    WRITE(UNIT=i,FMT="(A)") "*****************************************************"
    WRITE(UNIT=i,FMT="(A)") "#2 CONTAMINATION MODEL OF THE DATA (ONLY IF 0.03 = 1)"
    WRITE(UNIT=i,FMT="(A)") "*****************************************************"
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "2.01) Contamination model: 0=no contamination, 1= contamination"
    WRITE(UNIT=i,FMT=*) contam_model(1), "         -> Additive outliers"
    WRITE(UNIT=i,FMT=*) contam_model(2), "         -> Innovation (as many as number of IID series!)"
    DO j=3, DGP_nb_iid_series
       WRITE(UNIT=i,FMT=*) contam_model(j)
    END DO
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "2.02) Epsilon contamination:"
    WRITE(UNIT=i,FMT=*) epsilon_p(1), "         -> Additive outliers"
    WRITE(UNIT=i,FMT=*) epsilon_p(2), "         -> Innovation (as many as number of IID series!)"
    DO j=3, DGP_nb_iid_series
       WRITE(UNIT=i,FMT=*) epsilon_p(j)
    END DO
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "2.03) Single point for the contamination:"
!    WRITE(UNIT=i,FMT=*) x0(1), "          -> Additive outliers"
!    WRITE(UNIT=i,FMT=*) x0(2), "          -> Innovation (as many as number of IID series!)"
    DO j=3, DGP_nb_iid_series
 !      WRITE(UNIT=i,FMT=*) x0(j)
    END DO
    CALL WRITE_EMPTY(i,5)





    WRITE(UNIT=i,FMT="(A)") "***********************************************************"
    WRITE(UNIT=i,FMT="(A)") "#3 SECTION FIRST STEP --> ESTIMATION OF THE AUXILIARY MODEL"
    WRITE(UNIT=i,FMT="(A)") "***********************************************************"
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "3.01) Run number (I):"
    WRITE(UNIT=i,FMT=*) run_number
    CALL WRITE_EMPTY(i,3)    

    WRITE(UNIT=i,FMT="(A)") "3.02) Number of lags in the conditional mean, >=0:"
    WRITE(UNIT=i,FMT=*) dim_mu-1
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "3.03) Volatility estimation: first arch then garch orders:"
    WRITE(UNIT=i,FMT=*) dim_arch
    WRITE(UNIT=i,FMT=*) dim_garch
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "3.04) Dimension of the polynomial, pol_dim >=0:"
    WRITE(UNIT=i,FMT=*) pol_dim
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "3.05) Maximal degree of every polynomial, pol_deg(i) >=0:"
    IF (pol_dim .GT. 0) THEN
       DO j=1, pol_dim
          WRITE(UNIT=i,FMT=*) pol_deg(j)
       END DO
    END IF
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "3.06) Score covariance matrix: 1=Outproduct, 2=Newey-West, 3=Gallant"
    WRITE(UNIT=i,FMT=*) sel_score_cov_mat_s1
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "3.07) Number of lags for the Newey-West method, >=0:"
    WRITE(UNIT=i,FMT=*) lags_score_cov_mat_s1
    CALL WRITE_EMPTY(i,3) 

    WRITE(UNIT=i,FMT="(A)") "3.08) Estimates covariance matrix: 1=Identity Hessian, 2=Estimate the Hessian"
    WRITE(UNIT=i,FMT=*) sel_cov_estimates_s1
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "3.09) Hessian matrix estimation: 1=Outproduct, 2=second derivatives...2not implemented!"
    WRITE(UNIT=i,FMT=*) sel_hessian_mat_s1
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "3.10) Starting values for the cond. mean, first the constant then AR-parameters (I R):"
    DO j=1, dim_mu
       WRITE(UNIT=i,FMT=*) is_fixed_par_AM(j), beta_AM_end(j)
    END DO
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") &
         "3.11) Starting values for the cond. variance. First the constant, arch and garch param. (I R):"
    DO j=1, dim_sigma
       WRITE(UNIT=i,FMT=*) is_fixed_par_AM(dim_mu+j), beta_AM_end(dim_mu+j)
    END DO
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "3.12) Lower bound parameter in the SNP density, epsilon_0:"
    WRITE(UNIT=i,FMT=*) epsilon_0
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "3.13) Starting values for the polynomial coefficient (I R):"
    IF (nb_pol_coef .GT. 0) THEN
       DO j=1, nb_pol_coef
          WRITE(UNIT=i,FMT=*) is_fixed_par_AM(dim_mu+dim_sigma+j), beta_AM_end(dim_mu+dim_sigma+j)    
       END DO
    END IF
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "3.14) The number of linear constraints on the auxiliary parameters:" 
    WRITE(UNIT=i,FMT=*) nclin
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "3.15) The number of non linear constraints on the auxiliary parameters:"
    WRITE(UNIT=i,FMT=*) ncnln
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "3.16) Center and scale: 1 (0) means (not) to center the series. The same for scale."
    WRITE(UNIT=i,FMT=*) center
    WRITE(UNIT=i,FMT=*) scale
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "3.17) Compute the gradient: 0=numerical approximation, 1=exact computation"
    WRITE(UNIT=i,FMT=*) derivative_level
    CALL WRITE_EMPTY(i,5)





    WRITE(UNIT=i,FMT="(A)") "********************************************************"
    WRITE(UNIT=i,FMT="(A)") "#4 SECTION STEP 2 --> ESTIMATION OF THE STRUCTURAL MODEL"
    WRITE(UNIT=i,FMT="(A)") "********************************************************"
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "4.01) Name of the structural model (CHARACTER 15):"
    WRITE(UNIT=i,FMT="(A)") name_SM
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "Select between: MOVING__AVERAGE - STOC_VOLATILITY - AUTOREG_PROCESS"
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "4.02) Number of parts of the structural model (I):"
    WRITE(UNIT=i,FMT=*) SM_model_parts
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "Selection table:"
    WRITE(UNIT=i,FMT="(A)") "MOVING__AVERAGE     1   (only MA-part)"
    WRITE(UNIT=i,FMT="(A)") "STOC_VOLATILITY     2   (AR-part for the mean, AR-part for the volatility)"
    WRITE(UNIT=i,FMT="(A)") "AUTOREG_PROCESS     1   (AR-part for the mean)"
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "4.03) Order of every part (I):"
    DO j=1, SM_model_parts
       WRITE(UNIT=i,FMT=*) SM_dim_parts(j)
    END DO
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "Selection table:"
    WRITE(UNIT=i,FMT="(A)") "MOVING__AVERAGE     MA-order"
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "STOC_VOLATILITY     AR-order for the mean function"
    WRITE(UNIT=i,FMT="(A)") "                    AR-order for the var. function"
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "AUTOREG_PROCESS     AR-part for the mean"

    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "4.04) Starting values for the structural parameters (I R):"
    DO j=1, SM_nb_par
       WRITE(UNIT=i,FMT=*) is_fixed_par_SM(j), rho_SM_start(j)
    END DO
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") &
         "The parameters must be written in vertical order. If you don't want the constant, then restrict its value = 0."
    WRITE(UNIT=i,FMT="(A)") "MOVING__AVERAGE     constant, lags(1) ...  lags(n), variance"
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "STOC_VOLATILITY     MEAN: constant, lags(1) ... lags(n1)"
    WRITE(UNIT=i,FMT="(A)") "                    VAR.: constant, lags(1) ... lags(n2), sigma"
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "AUTOREG_PROCESS     MEAN: constant, lags(1) ... lags(n), variance"
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "4.05) parameter > inequality conditions (I R):"
    DO j=1, SM_nb_par
       WRITE(UNIT=i,FMT=*) is_gt_par_SM(j), gt_par_SM_value(j)
    END DO
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "4.06) parameter < inequality conditions (I R):"
    DO j=1, SM_nb_par
       WRITE(UNIT=i,FMT=*) is_lt_par_SM(j), lt_par_SM_value(j)
    END DO
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "4.07) Number of iid series necessary for the simulation of the SM:" 
    WRITE(UNIT=i,FMT=*) nb_iid_series
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "4.08) The length of the simulation used to approximate the expectation:"
    WRITE(UNIT=i,FMT=*) T2
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "4.09) The number of observations to discard in the simulation:"
    WRITE(UNIT=i,FMT=*) nb_discard_step_2
    CALL WRITE_EMPTY(i,5)





    WRITE(UNIT=i,FMT="(A)") "**********************************************************"
    WRITE(UNIT=i,FMT="(A)") "#5 SECTION ROBUST ESTIMATION --> ONLY ACTIVE IF (0.01) = 1"
    WRITE(UNIT=i,FMT="(A)") "**********************************************************"
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "5.01) Set the bound c of the robustness degree:"
    WRITE(UNIT=i,FMT=*) bound
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") "5.02) Set the metric used to calculate the weights: 0=simple, 1=efficient"
    WRITE(UNIT=i,FMT=*) weights_metric
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "5.03) Classic estimates as starting values: t=yes, f=no"
    WRITE(UNIT=i,FMT=*) rob_start_value(1)
    WRITE(UNIT=i,FMT=*) rob_start_value(2)
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "5.04) Select the maximal nb. of loops in the robust estimation:"
    WRITE(UNIT=i,FMT=*) rob_iteration_limit
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "5.05) Select the epsilon convergence criteria to stop the loops:"
    WRITE(UNIT=i,FMT=*) rob_parameter_tolerance
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "5.06) Select how the A matrix must be computed: 1=use observed data, 2=by simulation"
    WRITE(UNIT=i,FMT=*) A_updating_kind
    CALL WRITE_EMPTY(i,5)

    WRITE(UNIT=i,FMT="(A)") "5.07) Number of updates in the calculation of A and alpha"
    WRITE(UNIT=i,FMT=*) nb_A_updates
    CALL WRITE_EMPTY(i,5)



    WRITE(UNIT=i,FMT="(A)") "***********************************************************************"
    WRITE(UNIT=i,FMT="(A)") "#6 SECTION NPSOL SETUP AND PRINT OPTIONS: USE -1 TO ACCEPT THE DEFAULTS"
    WRITE(UNIT=i,FMT="(A)") "***********************************************************************"
    WRITE(UNIT=i,FMT=*)
    WRITE(UNIT=i,FMT="(A)") &
         "6.01) Set the bound used to define the infinity when bounding the parameters (R): default 10^20."
    WRITE(UNIT=i,FMT=*) infinite_bound_size
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") &
         "6.02) Set the max. acceptable absolute violation in the constraints (R): default e^0.5."
    WRITE(UNIT=i,FMT=*) feasibility_tolerance
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "6.03) Set the function precision, denoted eR, i.e. the measure of the accuracy"
    WRITE(UNIT=i,FMT="(A)") "      in the objective function computation (R): default e^0.9."
    WRITE(UNIT=i,FMT=*) function_precision
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "6.04) Set the optimality tolerance accuracy in the approximation of the final"
    WRITE(UNIT=i,FMT="(A)") "      solution of the problem (R), default eR^0.8."
    WRITE(UNIT=i,FMT=*) optimality_tolerance
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "6.05) Set the maximal number of major iterations permitted before termination:"
    WRITE(UNIT=i,FMT=*) major_iteration_limit
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") &
         "6.06) Set the maximal number of iterations for the optimization phase of each QP subproblem:"
    WRITE(UNIT=i,FMT=*) minor_iteration_limit
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "6.07) Set the Nolist option: 0=print the option listing, 1=suppress it."
    WRITE(UNIT=i,FMT=*) Nolist
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") &
         "6.08) Print File option:  0=suppresses any printing not associated with the optional parameters,"
    WRITE(UNIT=i,FMT="(A)") "                         -1=default, i=unit number to which the output must be sent."
    WRITE(UNIT=i,FMT=*) Print_File
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "6.09) Set the major print level:  0=is the default, i.e. no output exept error messages"
    WRITE(UNIT=i,FMT="(A)") "                                  1=the final solution"
    WRITE(UNIT=i,FMT="(A)") &
         "                                  5=one line of output for each iteration (but no final solution)"
    WRITE(UNIT=i,FMT="(A)") &
         "                                 10=the final solution and one line of output for each iteration"
    WRITE(UNIT=i,FMT=*) major_print_level
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "6.10) Set the minor print level:  0=is the default, i.e. no output exept error messages"
    WRITE(UNIT=i,FMT="(A)") "                                  1=The final QP solution"
    WRITE(UNIT=i,FMT="(A)") &
         "                                  5=One line of output for each minor iter. (but no final QP solution)"
    WRITE(UNIT=i,FMT="(A)") &
         "                                 10=The final QP solution and one brief line of output for each min. iter."
    WRITE(UNIT=i,FMT=*) minor_print_level
    CALL WRITE_EMPTY(i,3)

    WRITE(UNIT=i,FMT="(A)") "6.11) Verify gradients options: 0=don't verify, 1=verify."
    WRITE(UNIT=i,FMT=*) verify_gradients

    RETURN
  END SUBROUTINE NEW_RUN


  SUBROUTINE WRITE_RESULTS_SNP_STEP_1(f,start_log_lik)
    USE S1_COMMON, ONLY: T1, t_start, dim_mu, dim_sigma, dim_arch, dim_garch, pol_dim, &
         pol_deg, epsilon_0, center, scale, empirical_is_mean, is_fixed_par_AM, is_restricted_par_AM, &
         empirical_is_variance, inform, iter, cov_estimates_s1, AM_nb_not_fixed_par, &
         beta_AM_end, nb_not_fixed_mu, nb_not_fixed_sigma, nb_pol_coef, version, tmp_unit_nb, &
         AM_nb_par, AM_nb_free_par
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN) :: f
    REAL (KIND=realkind), INTENT(IN) :: start_log_lik

    INTEGER :: i,j,s
    REAL (KIND=realkind) :: stdv, t_value
    CHARACTER (LEN=21) :: lab
    CHARACTER (LEN=21) :: lab_fixed
    CHARACTER (LEN=21) :: lab_restr


100 FORMAT (A,I12)
101 FORMAT (A,E12.6)
102 FORMAT (A,I2,A,E12.6,A,E12.6,A,E12.6)
104 FORMAT (A,A12)
105 FORMAT (A,E22.10)

    OPEN (UNIT=tmp_unit_nb,FILE="results_step_1.txt",FORM="FORMATTED")

    WRITE(UNIT=tmp_unit_nb,FMT="(A)") version
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Total number of obs. T1  ", T1
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Start at period........  ", t_start
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Nb. initial lags.......  ", t_start-1
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "T1- initial lags.......  ", T1-t_start+1
    WRITE(UNIT=tmp_unit_nb,FMT=*)

    IF (.NOT. is_fixed_par_AM(1)) THEN
       WRITE(UNIT=tmp_unit_nb,FMT=104)    "Est. with free constant  ", "         yes"
    ELSE
       WRITE(UNIT=tmp_unit_nb,FMT=104)    "Est. with free constant  ", "          no"
       WRITE(UNIT=tmp_unit_nb,FMT=101)    "Value of the constant..  ", beta_AM_end(1)
    END IF
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "AR.....................  ", dim_mu-1
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "ARCH...................  ", dim_arch
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "GARCH..................  ", dim_garch
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Polynomial Dimension...  ", pol_dim
    IF (pol_dim .GT. 0) THEN
       WRITE(UNIT=tmp_unit_nb,FMT="(A)")  "           Degrees.....  "
       DO j=1, pol_dim
          WRITE(UNIT=tmp_unit_nb,FMT=100) &
               "Degree.................  ", pol_deg(j)
       END DO
    END IF
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Nb. total parameters...  ", AM_nb_par
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Nb. fixed parameters...  ", AM_nb_par - AM_nb_not_fixed_par
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Nb. estim parameters...  ", AM_nb_free_par
    WRITE(UNIT=tmp_unit_nb,FMT=*)

    WRITE(UNIT=tmp_unit_nb,FMT=101)       "epsilon_0..............  ", epsilon_0
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Center option..........  ", center
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Scale option...........  ", scale
    WRITE(UNIT=tmp_unit_nb,FMT=*)

    WRITE(UNIT=tmp_unit_nb,FMT=101)       "Mean of the series.....  ", empirical_is_mean
    WRITE(UNIT=tmp_unit_nb,FMT=101)       "Stdv of the series.....  ", SQRT(empirical_is_variance)

    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Convergence code.......  ", inform
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Numb. of maj. iter.....  ", iter
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT=*)

    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)")     "         RESULTS"
    WRITE(UNIT=tmp_unit_nb,FMT=*)   
    WRITE(UNIT=tmp_unit_nb,FMT="(A,A,A,A)") "AR                      ", "   ESTIMATE  ", "          STDV  ", &
         "          t-TEST "

    lab       = ")..................  "
    lab_fixed = ") fixed at.........  "
    lab_restr = ") restricted at....  "

    i=1
    s=0
    DO j=1, dim_mu
       IF ((.NOT. is_fixed_par_AM(j))) THEN
          IF ((.NOT. is_restricted_par_AM(j))) THEN
             stdv = SQRT(cov_estimates_s1(i,i)/REAL(T1-t_start+1))
             t_value = beta_AM_end(j)/stdv
          END IF
          i=i+1
       ELSE
          stdv = 0._realkind
          t_value = 0._realkind
       END IF

       IF (is_fixed_par_AM(j)) THEN
          WRITE(UNIT=tmp_unit_nb,FMT=102)    "C(", s, lab_fixed, beta_AM_end(j), "      ", stdv, &
               "      ", t_value
       ELSE
          IF (is_restricted_par_AM(j)) THEN
             WRITE(UNIT=tmp_unit_nb,FMT=102) "C(", s, lab_restr, beta_AM_end(j), "      ", stdv, &
                  "      ", t_value
          ELSE
             WRITE(UNIT=tmp_unit_nb,FMT=102) "C(", s, lab      , beta_AM_end(j), "      ", stdv, &
                  "      ", t_value
          END IF
       END IF
       s=s+1
    END DO

    WRITE(UNIT=tmp_unit_nb,FMT=*) 
    WRITE(UNIT=tmp_unit_nb,FMT="(A)")     "ARCH"
    s=0
    DO j=dim_mu+1, dim_mu+dim_arch+1
       IF ((.NOT. is_fixed_par_AM(j))) THEN
          IF ((.NOT. is_restricted_par_AM(j))) THEN
             stdv = SQRT(cov_estimates_s1(i,i)/REAL(T1-t_start+1))
             t_value = beta_AM_end(j)/stdv
          END IF
          i=i+1
       ELSE
          stdv = 0._realkind
          t_value = 0._realkind
       END IF

       IF (is_fixed_par_AM(j)) THEN
          WRITE(UNIT=tmp_unit_nb,FMT=102)    "A(", s, lab_fixed, beta_AM_end(j), "      ", stdv, &
               "      ", t_value
       ELSE
          IF (is_restricted_par_AM(j)) THEN
             WRITE(UNIT=tmp_unit_nb,FMT=102) "A(", s, lab_restr, beta_AM_end(j), "      ", stdv, &
                  "      ", t_value
          ELSE
             WRITE(UNIT=tmp_unit_nb,FMT=102) "A(", s, lab      , beta_AM_end(j), "      ", stdv, &
                  "      ", t_value
          END IF
       END IF
       s=s+1
    END DO
    s=1
    IF (dim_garch .GT. 0) THEN
       WRITE(UNIT=tmp_unit_nb,FMT=*) 
       WRITE(UNIT=tmp_unit_nb,FMT="(A)")  "GARCH"
       DO j=dim_mu+dim_arch+2, dim_mu+dim_arch+1+dim_garch
          IF ((.NOT. is_fixed_par_AM(j))) THEN
             IF ((.NOT. is_restricted_par_AM(j))) THEN
                stdv = SQRT(cov_estimates_s1(i,i)/REAL(T1-t_start+1))
                t_value = beta_AM_end(j)/stdv
             END IF
             i=i+1
          ELSE
             stdv = 0._realkind
             t_value = 0._realkind
          END IF

          IF (is_fixed_par_AM(j)) THEN
             WRITE(UNIT=tmp_unit_nb,FMT=102)    "G(", s, lab_fixed, beta_AM_end(j), "      ", stdv, &
                  "      ", t_value
          ELSE
             IF (is_restricted_par_AM(j)) THEN
                WRITE(UNIT=tmp_unit_nb,FMT=102) "G(", s, lab_restr, beta_AM_end(j), "      ", stdv, &
                     "      ", t_value
             ELSE
                WRITE(UNIT=tmp_unit_nb,FMT=102) "G(", s, lab      , beta_AM_end(j), "      ", stdv, &
                     "      ", t_value
             END IF
          END IF
          s=s+1
       END DO
    END IF
    s=0
    IF (pol_dim .GT. 0) THEN
       WRITE(UNIT=tmp_unit_nb,FMT=*)
       WRITE(UNIT=tmp_unit_nb,FMT="(A)")  "Pol. coeff."
       DO j=dim_mu+dim_sigma+1,dim_mu+dim_sigma+nb_pol_coef
          IF ((.NOT. is_fixed_par_AM(j))) THEN
             IF ((.NOT. is_restricted_par_AM(j))) THEN
                stdv = SQRT(cov_estimates_s1(i,i)/REAL(T1-t_start+1))
                t_value = beta_AM_end(j)/stdv
             END IF
             i=i+1
          ELSE
             stdv = 0._realkind
             t_value = 0._realkind
          END IF

          IF (is_fixed_par_AM(j)) THEN
             WRITE(UNIT=tmp_unit_nb,FMT=102)    "P(", s, lab_fixed, beta_AM_end(j), "      ", stdv, &
                  "      ", t_value
          ELSE
             IF (is_restricted_par_AM(j)) THEN
                WRITE(UNIT=tmp_unit_nb,FMT=102) "P(", s, lab_restr, beta_AM_end(j), "      ", stdv, &
                     "      ", t_value
             ELSE
                WRITE(UNIT=tmp_unit_nb,FMT=102) "P(", s, lab      , beta_AM_end(j), "      ", stdv, &
                     "      ", t_value
             END IF
          END IF
          s=s+1
       END DO
    END IF

    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT=105)       "-2 ln likelihood.......  ", 2*T1*f
    WRITE(UNIT=tmp_unit_nb,FMT=105)       "-2 ln likel. at start..  ", 2*T1*start_log_lik
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT=*)


    WRITE(UNIT=tmp_unit_nb,FMT="(A)")     "                 DESCRIPTION OF THE CONVERGENCE CODE"
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "0   The iterates have converged to a point x that satisfies the optimality con-"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    ditions to the accuracy requested by the Linear feasibility tolerance,"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    the Nonlinear feasibility tolerance, and the Optimality tolerance."
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    That is, the active constraint residuals and the reduced gradient are negligible"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    at x."
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "1   The final iterate x satisfies the optimality conditions to the accuracy re"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    quested, but the sequence of iterates has not yet converged. NPSOL was"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    terminated because no further improvement could be made in the merit func-"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    tion."
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "2   The linear constraints and bounds could not be satisfied. The problem has"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    no feasible solution."
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "3   The nonlinear constraints could not be satisfied. The problem may have no"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    feasible solution."
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "4   The Major iteration limit was reached."
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "6   x does not satisfy the first-order optimality conditions to the required accu-"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    racy, and no improved point for the merit function could be found during the"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    final linesearch."
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "7   The function derivatives returned by funcon or funobj appear to be incor-"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    rect."
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "9   An input parameter was invalid."

    CLOSE (UNIT=tmp_unit_nb)

    RETURN
  END SUBROUTINE WRITE_RESULTS_SNP_STEP_1

  SUBROUTINE WRITE_RESULTS_SNP_STEP_2(f)
    USE S1_COMMON, ONLY: T1, version, inform, iter, AM_nb_par, AM_nb_not_fixed_par, tmp_unit_nb
    USE S2_COMMON, ONLY: T2, nb_discard_step_2, rho_SM_end, SM_nb_par, SM_nb_not_fixed_par, &
         is_fixed_par_SM
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN) :: f

    INTEGER :: j

    OPEN (UNIT=tmp_unit_nb,FILE="results_step_2.txt",FORM="FORMATTED")
100 FORMAT (A,I12)
101 FORMAT (A,E12.6)
102 FORMAT (A,I1,A,E12.6)
103 FORMAT (A,I2,A,E12.6)

    WRITE(UNIT=tmp_unit_nb,FMT="(A)") version
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Nb. structural parameters  ", SM_nb_par
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Nb. auxiliary  parameters  ", AM_nb_par
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Nb. free struct. paramet.  ", SM_nb_not_fixed_par
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Nb. free auxil.  paramet.  ", AM_nb_not_fixed_par
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Total number of obs. T1..  ", T1
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Simulation length    T2..  ", T2
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Nb. discarded sim. obs...  ", nb_discard_step_2
    WRITE(UNIT=tmp_unit_nb,FMT=*)

    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Kind of convergence......  ", inform
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Numb. of maj. iter.......  ", iter
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT=*)



    WRITE(UNIT=tmp_unit_nb,FMT="(A)")     "         RESULTS"
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A,A)")   "RHO                        ", "  ESTIMATE  "
    WRITE(UNIT=tmp_unit_nb,FMT=*)

    DO j=1, SM_nb_par
       IF (.NOT. is_fixed_par_SM(j)) THEN
          IF (j .LE. 9) THEN
             WRITE(UNIT=tmp_unit_nb,FMT=102) "C(", j, ").....................  ", rho_SM_end(j)
          ELSE
             WRITE(UNIT=tmp_unit_nb,FMT=103)  "C(", j, ").................  ", rho_SM_end(j)
          END IF
       ELSE
          IF (j .LE. 9) THEN
             WRITE(UNIT=tmp_unit_nb,FMT=102) "C(", j, ") fixed at............  ", rho_SM_end(j)
          ELSE
             WRITE(UNIT=tmp_unit_nb,FMT=103)  "C(", j, ") fixed at........  ", rho_SM_end(j)
          END IF
       END IF
    END DO

    WRITE(UNIT=tmp_unit_nb,FMT=*)

    WRITE(UNIT=tmp_unit_nb,FMT=101)       "objective function.......  ", f
    WRITE(UNIT=tmp_unit_nb,FMT=101)       "T1*(objective function)..  ", T1*f
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT=*)


    WRITE(UNIT=tmp_unit_nb,FMT="(A)")     "           DESCRIPTION OF THE CONVERGENCE CODE"
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "0   The iterates have converged to a point x that satisfies the optimality con-"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    ditions to the accuracy requested by the Linear feasibility tolerance,"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    the Nonlinear feasibility tolerance, and the Optimality tolerance."
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    That is, the active constraint residuals and the reduced gradient are negligible"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    at x."
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "1   The final iterate x satisfies the optimality conditions to the accuracy re"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    quested, but the sequence of iterates has not yet converged. NPSOL was"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    terminated because no further improvement could be made in the merit func-"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    tion."
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "2   The linear constraints and bounds could not be satisfied. The problem has"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    no feasible solution."
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "3   The nonlinear constraints could not be satisfied. The problem may have no"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    feasible solution."
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "4   The Major iteration limit was reached."
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "6   x does not satisfy the first-order optimality conditions to the required accu-"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    racy, and no improved point for the merit function conld be found during the"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    final linesearch."
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "7   The function derivatives returned by funcon or funobj appear to be incor-"
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "    rect."
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A)") "9   An input parameter was invalid."

    CLOSE (UNIT=tmp_unit_nb)

    RETURN
  END SUBROUTINE WRITE_RESULTS_SNP_STEP_2


  SUBROUTINE WRITE_RESULTS_SNP_ROB_STEP_1()
    USE S1_COMMON, ONLY: dim_mu, dim_sigma, dim_arch, dim_garch, pol_dim, &
         beta_AM_end, j1, j2, j3, nb_pol_coef, tmp_unit_nb
    IMPLICIT NONE

    INTEGER :: j

102 FORMAT (A,I1,A,E12.6)
103 FORMAT (A,I2,A,E12.6)


    !write the estimated parameters of the auxiliary model
    OPEN (UNIT=tmp_unit_nb,FILE="results_step_1_rob.txt",FORM="FORMATTED")
    WRITE(UNIT=tmp_unit_nb,FMT=*) "ESTIMATED ROBUST PARAMETERS: AUXILIARY MODEL"
    WRITE(UNIT=tmp_unit_nb,FMT=*)   
    WRITE(UNIT=tmp_unit_nb,FMT="(A,A)") "AR                      ", "   ESTIMATE  "

    DO j=1, dim_mu
       IF (j .LT. 11) THEN
          WRITE(UNIT=tmp_unit_nb,FMT=102)    "C(", j-1, ")...................  ", beta_AM_end(j)
       ELSE
          WRITE(UNIT=tmp_unit_nb,FMT=103)    "C(", j-1, ")..................  ", beta_AM_end(j)
       END IF
    END DO

    WRITE(UNIT=tmp_unit_nb,FMT=*) 
    WRITE(UNIT=tmp_unit_nb,FMT="(A)")     "ARCH"

    DO j=1, dim_arch+1
       IF (j .LT. 11) THEN
          WRITE(UNIT=tmp_unit_nb,FMT=102)    "A(", j-1, ")...................  ", &
               beta_AM_end(dim_mu+j)
       ELSE
          WRITE(UNIT=tmp_unit_nb,FMT=103)    "A(", j-1,  ")..................  ", &
               beta_AM_end(dim_mu+j)
       END IF
    END DO

    IF (dim_garch .GT. 0) THEN
       WRITE(UNIT=tmp_unit_nb,FMT=*) 
       WRITE(UNIT=tmp_unit_nb,FMT="(A)")  "GARCH"
       DO j=1, dim_garch
          IF (j .LE. 9) THEN
             WRITE(UNIT=tmp_unit_nb,FMT=102) "G(", j, ")..................   ", &
                  beta_AM_end(dim_mu+dim_arch+j+1)

          ELSE
             WRITE(UNIT=tmp_unit_nb,FMT=103) "G(", j, ").................   ", &
                  beta_AM_end(dim_mu+dim_arch+j+1)
          END IF
       END DO
    END IF

    IF (pol_dim .GT. 0) THEN
       WRITE(UNIT=tmp_unit_nb,FMT=*)
       WRITE(UNIT=tmp_unit_nb,FMT="(A)")  "Pol. coeff."
       DO j=1,nb_pol_coef
          IF (j .LT. 11) THEN
             WRITE(UNIT=tmp_unit_nb,FMT=102) "P(", j-1,")...................  ", &
                  beta_AM_end(dim_mu+dim_sigma+j)

          ELSE
             WRITE(UNIT=tmp_unit_nb,FMT=103) "P(", j-1, ")..................  ", &
                  beta_AM_end(dim_mu+dim_sigma+j)
          END IF
       END DO
    END IF
    CLOSE(UNIT=tmp_unit_nb)
    RETURN
  END SUBROUTINE WRITE_RESULTS_SNP_ROB_STEP_1


  SUBROUTINE WRITE_RESULTS_SNP_ROB_STEP_2()
    USE S1_COMMON,  ONLY: T1, N1, t_start, version, AM_nb_par, AM_nb_not_fixed_par, tmp_unit_nb
    USE S2_COMMON,  ONLY: T2, nb_discard_step_2, rho_SM_end, SM_nb_par, SM_nb_not_fixed_par, &
         is_fixed_par_SM
    USE ROB_COMMON, ONLY: rob_weights_N1
    IMPLICIT NONE

    INTEGER :: i,j

    OPEN (UNIT=tmp_unit_nb,FILE="results_step_2_rob.txt",FORM="FORMATTED")
100 FORMAT (A,I12)
102 FORMAT (A,I1,A,E12.6)
103 FORMAT (A,I2,A,E12.6)

    WRITE(UNIT=tmp_unit_nb,FMT="(A)") version
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Nb. structural parameters  ", SM_nb_par
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Nb. auxiliary  parameters  ", AM_nb_par
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Nb. free struct. paramet.  ", SM_nb_not_fixed_par
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Nb. free auxil.  paramet.  ", AM_nb_not_fixed_par
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Total number of obs. T1..  ", T1
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Simulation length    T2..  ", T2
    WRITE(UNIT=tmp_unit_nb,FMT=100)       "Nb. discarded sim. obs...  ", nb_discard_step_2
    WRITE(UNIT=tmp_unit_nb,FMT=*)

    !WRITE(UNIT=tmp_unit_nb,FMT=100)       "Kind of convergence......  ", inform
    !WRITE(UNIT=tmp_unit_nb,FMT=100)       "Numb. of maj. iter.......  ", iter
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT=*)



    WRITE(UNIT=tmp_unit_nb,FMT="(A)")     "         RESULTS"
    WRITE(UNIT=tmp_unit_nb,FMT=*)
    WRITE(UNIT=tmp_unit_nb,FMT="(A,A)")   "RHO                        ", "  ESTIMATE  "
    WRITE(UNIT=tmp_unit_nb,FMT=*)

    DO j=1, SM_nb_par
       IF (.NOT. is_fixed_par_SM(j)) THEN
          IF (j .LE. 9) THEN
             WRITE(UNIT=tmp_unit_nb,FMT=102) "C(", j, ").....................  ", rho_SM_end(j)
          ELSE
             WRITE(UNIT=tmp_unit_nb,FMT=103)  "C(", j, ").................  ", rho_SM_end(j)
          END IF
       ELSE
          IF (j .LE. 9) THEN
             WRITE(UNIT=tmp_unit_nb,FMT=102) "C(", j, ") fixed at............  ", rho_SM_end(j)
          ELSE
             WRITE(UNIT=tmp_unit_nb,FMT=103)  "C(", j, ") fixed at........  ", rho_SM_end(j)
          END IF
       END IF
    END DO

    !write the weights of the robust estimation in a file
    OPEN (UNIT=tmp_unit_nb,FILE="series_weights.dat",FORM="FORMATTED")
    DO i=1, t_start-1
       WRITE(UNIT=tmp_unit_nb,FMT=*) 1._realkind
    END DO
    DO i=1,N1
       WRITE(UNIT=tmp_unit_nb,FMT=*) rob_weights_N1(i,1)
    END DO
    CLOSE(UNIT=tmp_unit_nb)

    RETURN
  END SUBROUTINE WRITE_RESULTS_SNP_ROB_STEP_2


  SUBROUTINE OUTPUT_SNP_STEP_1(f, start_log_lik)
    USE S1_COMMON,   ONLY: T1, AM_nb_par, AM_nb_not_fixed_par, beta_AM_end, &
         run_number, pol_dim, dim_mu, dim_sigma, modified_series, nb_pol_coef, &
         epsilon_0, score_cov,cov_estimates_s1, normal_c, tmp_unit_nb
    USE UTILITIES,   ONLY: RESTRICT_BETA
    USE MONTE_CARLO, ONLY: nb_rep_monte_carlo
    USE MYFUNCTIONS, ONLY: SNP_DENSITY_FUNC
    USE COVARIANCE_ESTIMATES, ONLY: EST_COV_SCORE_AM_STEP1, EST_COV_AM_FREE_PAR_STEP1
    IMPLICIT NONE
    REAL (KIND=realkind), INTENT(IN) :: f, start_log_lik

    INTEGER :: i,j
    REAL (KIND=realkind) :: ir
    REAL (KIND=realkind), DIMENSION(1:AM_nb_not_fixed_par) :: beta_restricted

    CALL RESTRICT_BETA(beta_AM_end,beta_restricted)

    !compute the covariance of the score
    CALL EST_COV_SCORE_AM_STEP1(beta_restricted,T1,modified_series,1)

    !write the estimated covariance matrix of the score
    OPEN (UNIT=tmp_unit_nb,FILE="AM_score_cov.txt",FORM="FORMATTED")
    DO i=1, AM_nb_not_fixed_par
       DO  j=1, AM_nb_not_fixed_par
          WRITE(UNIT=tmp_unit_nb,FMT=*) i , j, score_cov(i,j)
       END DO
    END DO
    CLOSE (UNIT=tmp_unit_nb)

    !compute the covariance of the estimates
    CALL EST_COV_AM_FREE_PAR_STEP1(beta_restricted,T1,modified_series,1)

    !write the estimated covariance matrix of the estimates
    OPEN (UNIT=tmp_unit_nb,FILE="AM_param_cov.txt",FORM="FORMATTED")
    DO i=1, AM_nb_not_fixed_par
       DO  j=1, AM_nb_not_fixed_par
          WRITE(UNIT=tmp_unit_nb,FMT=*) i,j,cov_estimates_s1(i,j)
       END DO
    END DO
    CLOSE (UNIT=tmp_unit_nb)

    CALL WRITE_RESULTS_SNP_STEP_1(f,start_log_lik)

    !write the value of the log-likelihood to the likelihood.txt file
    OPEN (UNIT=tmp_unit_nb,FILE="likelihood.txt",STATUS="UNKNOWN",POSITION="APPEND",FORM="FORMATTED")
    WRITE(UNIT=tmp_unit_nb,FMT="(I4,A,E12.6)") run_number, " Log-likelihood.........", f 
    CLOSE (UNIT=tmp_unit_nb)

    !Write the new setup.new file
    run_number = run_number + 1 

    OPEN (UNIT=tmp_unit_nb,FILE="setup.new",FORM="FORMATTED")
    !    CALL NEW_RUN(tmp_unit_nb)
    CLOSE (UNIT=tmp_unit_nb)

    !if the polynomial dimension is 1 write the values of the estimated density between [-4.5,4.5]
    IF ((pol_dim .EQ. 1) .AND. (nb_rep_monte_carlo .EQ. 0)) THEN  
       !write the estimated snp density function
       OPEN (UNIT=tmp_unit_nb,FILE="density_snp.dat",FORM="FORMATTED")
       ir = 0.0_realkind
       DO
          IF (ir .GT. 9._realkind) EXIT
          WRITE(UNIT=tmp_unit_nb,FMT=*) -4.5_realkind+ir, SNP_DENSITY_FUNC(-4.5_realkind+ir, &
               beta_AM_end(dim_mu+dim_sigma+1:AM_nb_par), nb_pol_coef, epsilon_0)
          ir = ir + .01_realkind
       END DO
       CLOSE (UNIT=tmp_unit_nb)
       !write the standard normal density
       OPEN (UNIT=tmp_unit_nb,FILE="density_standard_normal.dat",FORM="FORMATTED")
       ir = 0.0_realkind
       DO
          IF (ir .GT. 9._realkind) EXIT
          WRITE(UNIT=tmp_unit_nb,FMT=*) -4.5_realkind+ir, normal_c*exp(-0.5_realkind*(-4.5_realkind+ir)**2)
          ir = ir + .01_realkind
       END DO
       CLOSE (UNIT=tmp_unit_nb)
    END IF

    RETURN
  END SUBROUTINE OUTPUT_SNP_STEP_1

  SUBROUTINE OUTPUT_SNP_STEP_2(f)
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN) :: f

    CALL WRITE_RESULTS_SNP_STEP_2(f)  

    RETURN
  END SUBROUTINE OUTPUT_SNP_STEP_2


  SUBROUTINE OUTPUT_SNP_MC()
    USE S1_COMMON,   ONLY: beta_AM_end, epsilon_0, pol_dim, dim_mu, dim_sigma, step, &
         AM_nb_par, nb_pol_coef, normal_c, is_fixed_par_AM, val_fix_restr_par_AM,&
         tmp_unit_nb, npsol_inform_step_1, npsol_iter_step_1
    USE S2_COMMON,   ONLY: npsol_inform_step_2, npsol_iter_step_2, hansen_test_results
    USE MONTE_CARLO, ONLY: mc_results_step_1, mc_results_step_2, nb_rep_monte_carlo, &
         skip_first_n, estimate_all_mc_simulations, selected_simulations
    USE MYFUNCTIONS, ONLY: SNP_DENSITY_FUNC           
    IMPLICIT NONE

    INTEGER :: i,j,sub_i
    LOGICAL :: ok
    REAL (KIND=realkind) :: ir

100 FORMAT (I5,TR1,I1,TR1,I3,TR1,25(ES12.4,TR1))

    OPEN (UNIT=tmp_unit_nb,FILE="mc_results_step_1.dat",FORM="FORMATTED")
    sub_i=1
    DO i=1,nb_rep_monte_carlo
       IF (estimate_all_mc_simulations .EQ. 1) THEN
          IF (skip_first_n .LT. i) THEN
             ok=.true.
          ELSE
             ok=.false.
          END IF
       ELSE
          IF (i .EQ. selected_simulations(sub_i)) THEN
             sub_i = sub_i + 1
             ok=.true.
          ELSE
             ok=.false.
          END IF
       END IF

       IF (ok) THEN
         WRITE(UNIT=tmp_unit_nb,FMT=100) i, npsol_inform_step_1(i), npsol_iter_step_1(i), mc_results_step_1(:,i)
      END IF
    END DO
    CLOSE(UNIT=tmp_unit_nb)

    IF (pol_dim .EQ. 1) THEN
       OPEN (UNIT=tmp_unit_nb,FILE="density_snp_based_on_average_mc.dat",FORM="FORMATTED")
       DO i=1, AM_nb_par
          IF (.NOT. is_fixed_par_AM(i)) THEN
             ir = 0._realkind
             DO j=1, nb_rep_monte_carlo
                ir = ir + mc_results_step_1(i,j)
             END DO
             ir = ir / nb_rep_monte_carlo
             beta_AM_end(i) = ir
          ELSE
             beta_AM_end(i)=val_fix_restr_par_AM(i)
          END IF
       END DO

       ir = 0.0_realkind
       DO
          IF (ir .GT. 14._realkind) EXIT
          WRITE(UNIT=tmp_unit_nb,FMT=*) -7._realkind+ir, &
               SNP_DENSITY_FUNC(-7._realkind+ir, beta_AM_end(dim_mu+dim_sigma+1:AM_nb_par), &
               nb_pol_coef, epsilon_0)
          ir = ir + .01_realkind
       END DO
       CLOSE (UNIT=tmp_unit_nb)
       !write the standard normal density
       OPEN (UNIT=tmp_unit_nb,FILE="density_standard_normal.dat",FORM="FORMATTED")
       ir = 0.0_realkind
       DO
          IF (ir .GT. 14._realkind) EXIT
          WRITE(UNIT=tmp_unit_nb,FMT=*) -7._realkind+ir, normal_c*exp(-0.5_realkind*(-4.5_realkind+ir)**2)
          ir = ir + .01_realkind
       END DO
       CLOSE (UNIT=tmp_unit_nb)
       !write the values of the average theta vector 
       OPEN (UNIT=tmp_unit_nb,FILE="average_snp_theta_vector.dat",FORM="FORMATTED")
       WRITE(UNIT=tmp_unit_nb,FMT=*)  beta_AM_end(dim_mu+dim_sigma+1:AM_nb_par)
       CLOSE (UNIT=tmp_unit_nb)

    END IF

    IF (step .EQ. 3) THEN
       OPEN (UNIT=tmp_unit_nb,FILE="mc_results_step_2.dat",FORM="FORMATTED")
       sub_i=1
       DO i=1,nb_rep_monte_carlo
          IF (estimate_all_mc_simulations .EQ. 1) THEN
             IF (skip_first_n .LT. i) THEN
                ok=.true.
             ELSE
                ok=.false.
             END IF
          ELSE
             IF (i .EQ. selected_simulations(sub_i)) THEN
                sub_i = sub_i + 1
                ok=.true.
             ELSE
                ok=.false.
             END IF
          END IF

          IF (ok) THEN
             WRITE(UNIT=tmp_unit_nb,FMT=100) i, npsol_inform_step_2(i,1), npsol_iter_step_2(i,1), mc_results_step_2(:,i,1)
         END IF
       END DO
       CLOSE(UNIT=tmp_unit_nb)

       OPEN (UNIT=tmp_unit_nb,FILE="hansen_test_step_2.dat",FORM="FORMATTED")
       DO i=1,nb_rep_monte_carlo
          WRITE(UNIT=tmp_unit_nb,FMT=*) i, hansen_test_results(i,1)
       END DO
       CLOSE(UNIT=tmp_unit_nb)
    END IF

    RETURN
  END SUBROUTINE OUTPUT_SNP_MC


  SUBROUTINE OUTPUT_SNP_ROB_STEP_1()
    USE S1_COMMON,   ONLY: AM_nb_par, beta_AM_end, pol_dim, dim_mu, dim_sigma, nb_pol_coef, &
         epsilon_0, normal_c, tmp_unit_nb
    USE MYFUNCTIONS, ONLY: SNP_DENSITY_FUNC
    USE MONTE_CARLO, ONLY: nb_rep_monte_carlo
    IMPLICIT NONE

    REAL (KIND=realkind) :: ir

    CALL WRITE_RESULTS_SNP_ROB_STEP_1()

    !if the polynomial dimension is 1 write the values of the estimated density between [-4.5,4.5]
    IF ((pol_dim .EQ. 1) .AND. (nb_rep_monte_carlo .EQ. 0)) THEN  
       !write the estimated snp density function
       OPEN (UNIT=tmp_unit_nb,FILE="density_snp_rob.dat",FORM="FORMATTED")
       ir = 0.0_realkind
       DO
          IF (ir .GT. 9._realkind) EXIT
          WRITE(UNIT=tmp_unit_nb,FMT=*) -4.5_realkind+ir, SNP_DENSITY_FUNC(-4.5_realkind+ir, &
               beta_AM_end(dim_mu+dim_sigma+1:AM_nb_par), nb_pol_coef, epsilon_0)
          ir = ir + .01_realkind
       END DO
       CLOSE (UNIT=tmp_unit_nb)
       !write the standard normal density
       OPEN (UNIT=tmp_unit_nb,FILE="density_standard_normal.dat",FORM="FORMATTED")
       ir = 0.0_realkind
       DO
          IF (ir .GT. 9._realkind) EXIT
          WRITE(UNIT=tmp_unit_nb,FMT=*) -4.5_realkind+ir, normal_c*exp(-0.5_realkind*(-4.5_realkind+ir)**2)
          ir = ir + .01_realkind
       END DO
       CLOSE (UNIT=tmp_unit_nb)
    END IF

    RETURN
  END SUBROUTINE OUTPUT_SNP_ROB_STEP_1



  SUBROUTINE OUTPUT_SNP_ROB_STEP_2()
    IMPLICIT NONE

    CALL WRITE_RESULTS_SNP_ROB_STEP_2()

    RETURN
  END SUBROUTINE OUTPUT_SNP_ROB_STEP_2



  SUBROUTINE OUTPUT_SNP_ROB_MC()
    USE S1_COMMON,   ONLY: N1, step, AM_nb_par, beta_AM_end, pol_dim, is_fixed_par_AM, &
         dim_mu, dim_sigma, val_fix_restr_par_AM, nb_pol_coef, epsilon_0, tmp_unit_nb, &
         npsol_inform_step_1, npsol_iter_step_1, t_start, minimum_step_1_vector
    USE S2_COMMON,   ONLY: npsol_inform_step_2, npsol_iter_step_2, &
         minimum_step_2_vector
    USE ROB_COMMON
    USE MONTE_CARLO
    USE MYFUNCTIONS, ONLY: SNP_DENSITY_FUNC
    IMPLICIT NONE

    INTEGER :: i, j, sub_i
    REAL (KIND=realkind) :: ir
    LOGICAL :: ok

100 FORMAT (I5,TR1,I1,TR1,I3,TR1,25(ES12.4,TR1))
101 FORMAT (I5,TR1,I1,TR1,I3,TR1,I1,TR1,I3,TR1,25(ES12.4,TR1))

    IF (rob_start_value(1)) THEN
       !first write the non robust estimates of the auxiliary parameters
       OPEN (UNIT=tmp_unit_nb,FILE="mc_results_step_1.dat",FORM="FORMATTED")
       sub_i=1
       DO i=1,nb_rep_monte_carlo
          IF (estimate_all_mc_simulations .EQ. 1) THEN
             IF (skip_first_n .LT. i) THEN
                ok=.true.
             ELSE
                ok=.false.
             END IF
          ELSE
             IF (i .EQ. selected_simulations(sub_i)) THEN
                sub_i = sub_i + 1
                ok=.true.
             ELSE
                ok=.false.
             END IF
          END IF

          IF (ok) THEN
             WRITE(UNIT=tmp_unit_nb,FMT=100) i, npsol_inform_step_1(i), npsol_iter_step_1(i), &
                  minimum_step_1_vector(i), mc_results_step_1(:,i)
          END IF
       END DO
       CLOSE(UNIT=tmp_unit_nb)

       !write density calculated at the average value of the density parameters (non robust estimation)
       IF (pol_dim .EQ. 1) THEN
          OPEN (UNIT=tmp_unit_nb,FILE="density_snp_based_on_average_mc_step_1.dat",FORM="FORMATTED")
          DO i=1, AM_nb_par
             IF (.NOT. is_fixed_par_AM(i)) THEN
                ir = 0._realkind
                DO j=1, nb_rep_monte_carlo
                   ir = ir + mc_results_step_1(i,j)
                END DO
                ir = ir / nb_rep_monte_carlo
                beta_AM_end(i) = ir
             ELSE
                beta_AM_end(i)=val_fix_restr_par_AM(i)
             END IF
          END DO

          ir = 0.0_realkind
          DO
             IF (ir .GT. 9._realkind) EXIT
             WRITE(UNIT=tmp_unit_nb,FMT=*) -4.5_realkind+ir, &
                  SNP_DENSITY_FUNC(-4.5_realkind+ir, beta_AM_end(dim_mu+dim_sigma+1:AM_nb_par), &
                  nb_pol_coef, epsilon_0)
             ir = ir + .01_realkind
          END DO
          CLOSE (UNIT=tmp_unit_nb)
       END IF
    END IF

    !write the robust estimates of the auxiliary parameters
    OPEN (UNIT=tmp_unit_nb,FILE="mc_results_step_1_rob.dat",FORM="FORMATTED")
    sub_i=1
    DO i=1,nb_rep_monte_carlo
       IF (estimate_all_mc_simulations .EQ. 1) THEN
          IF (skip_first_n .LT. i) THEN
             ok=.true.
          ELSE
             ok=.false.
          END IF
       ELSE
          IF (i .EQ. selected_simulations(sub_i)) THEN
             sub_i = sub_i + 1
             ok=.true.
          ELSE
             ok=.false.
          END IF
       END IF

       IF (ok) THEN
          WRITE(UNIT=tmp_unit_nb,FMT=101) i, inform_rob(i), iter_rob(i), npsol_inform_step_1_rob(i), &
               npsol_iter_step_1_rob(i), minimum_step_1_rob_vector(i), mc_results_step_1_rob(:,i)
       END IF
    END DO
    CLOSE(UNIT=tmp_unit_nb)

    !write density calculated at the average value of the density parameters (robust estimation)
    IF (pol_dim .EQ. 1) THEN
       OPEN (UNIT=tmp_unit_nb,FILE="density_snp_based_on_average_mc_rob_step_1.dat",FORM="FORMATTED")
       DO i=1, AM_nb_par
          IF (.NOT. is_fixed_par_AM(i)) THEN
             ir = 0._realkind
             DO j=1, nb_rep_monte_carlo
                ir = ir + mc_results_step_1_rob(i,j)
             END DO
             ir = ir / nb_rep_monte_carlo
             beta_AM_end(i) = ir
          ELSE
             beta_AM_end(i) = val_fix_restr_par_AM(i)
          END IF
       END DO

       ir = 0.0_realkind
       DO
          IF (ir .GT. 9._realkind) EXIT
          WRITE(UNIT=tmp_unit_nb,FMT=*) -4.5_realkind+ir, &
               SNP_DENSITY_FUNC(-4.5_realkind+ir, beta_AM_end(dim_mu+dim_sigma+1:AM_nb_par), &
               nb_pol_coef, epsilon_0)
          ir = ir + .01_realkind
       END DO
       CLOSE (UNIT=tmp_unit_nb)
    END IF
    !write the non robust estimates of the structural parameters
    IF ((step .EQ. 3) .AND. (rob_start_value(2))) THEN
       OPEN (UNIT=tmp_unit_nb,FILE="mc_results_step_2.dat",FORM="FORMATTED")
       sub_i=1
       DO i=1,nb_rep_monte_carlo
          IF (estimate_all_mc_simulations .EQ. 1) THEN
             IF (skip_first_n .LT. i) THEN
                ok=.true.
             ELSE
                ok=.false.
             END IF
          ELSE
             IF (i .EQ. selected_simulations(sub_i)) THEN
                sub_i = sub_i + 1
                ok=.true.
             ELSE
                ok=.false.
             END IF
          END IF

          IF (ok) THEN
             WRITE(UNIT=tmp_unit_nb,FMT=100)  i, npsol_inform_step_2(i,1), npsol_iter_step_2(i,1), &
                  minimum_step_2_vector(i), mc_results_step_2(:,i,1)
          END IF
       END DO
       CLOSE(UNIT=tmp_unit_nb)
    END IF

    !write the robust estimates of the structural parameters
    IF (step .EQ. 3) THEN
       OPEN (UNIT=tmp_unit_nb,FILE="mc_results_step_2_rob.dat",FORM="FORMATTED")
       sub_i=1
       DO i=1,nb_rep_monte_carlo
          IF (estimate_all_mc_simulations .EQ. 1) THEN
             IF (skip_first_n .LT. i) THEN
                ok=.true.
             ELSE
                ok=.false.
             END IF
          ELSE
             IF (i .EQ. selected_simulations(sub_i)) THEN
                sub_i = sub_i + 1
                ok=.true.
             ELSE
                ok=.false.
             END IF
          END IF

          IF (ok) THEN
             WRITE(UNIT=tmp_unit_nb,FMT=101) i, inform_rob(i), iter_rob(i), &
                  npsol_inform_step_2_rob(i), npsol_iter_step_2_rob(i), &
                  minimum_step_2_rob_vector(i), mc_results_step_2_rob(:,i)
          END IF
       END DO
       CLOSE(UNIT=tmp_unit_nb)

       !write the weights of the robust estimation in a file
       IF (nb_rep_monte_carlo .GE. 1) THEN
          OPEN (UNIT=tmp_unit_nb,FILE="series_weights.dat",FORM="FORMATTED")
          DO i=1, t_start-1
             WRITE(UNIT=tmp_unit_nb,FMT=*) 1._realkind
          END DO
          DO i=1,N1
             WRITE(UNIT=tmp_unit_nb,FMT=*) rob_weights_N1(i,1)
          END DO
          CLOSE(UNIT=tmp_unit_nb)
       END IF

       IF (nb_rep_monte_carlo .GT. 0) THEN
          OPEN (UNIT=tmp_unit_nb,FILE="hansen_test_step_2_rob.dat",FORM="FORMATTED")
          sub_i=1
          DO i=1,nb_rep_monte_carlo
             IF (estimate_all_mc_simulations .EQ. 1) THEN
                IF (skip_first_n .LT. i) THEN
                   ok=.true.
                ELSE
                   ok=.false.
                END IF
             ELSE
                IF (i .EQ. selected_simulations(sub_i)) THEN
                   sub_i = sub_i + 1
                   ok=.true.
                ELSE
                   ok=.false.
                END IF
             END IF

             IF (ok) THEN
                WRITE(UNIT=tmp_unit_nb,FMT=*) i, hansen_test_results_rob(i,1)
             END IF
          END DO
          CLOSE(UNIT=tmp_unit_nb)
       END IF
    END IF
    RETURN
  END SUBROUTINE OUTPUT_SNP_ROB_MC



  SUBROUTINE OUTPUT_EXPECTATION_MC_STEP_2()
    USE S1_COMMON,   ONLY: tmp_unit_nb
    USE S2_COMMON,   ONLY: npsol_inform_step_2, npsol_iter_step_2, hansen_test_results
    USE MONTE_CARLO, ONLY: mc_results_step_1, mc_results_step_2, nb_rep_monte_carlo, &
         nb_rep_exp_analysis
    IMPLICIT NONE

    INTEGER :: mc_counter, exp_counter

    !write the results of the hansen test
    OPEN (UNIT=tmp_unit_nb,FILE="mc_hansen_test_exp.dat",FORM="FORMATTED")
    DO mc_counter=1, nb_rep_monte_carlo
       DO exp_counter=1, nb_rep_exp_analysis+1
          WRITE(UNIT=tmp_unit_nb,FMT=*) mc_counter, exp_counter, hansen_test_results(mc_counter,exp_counter)
       END DO
    END DO
    !write the results of the structural parameters
    OPEN (UNIT=tmp_unit_nb,FILE="mc_results_step_2_exp.dat",FORM="FORMATTED")
    DO mc_counter=1, nb_rep_monte_carlo
       DO exp_counter=1, nb_rep_exp_analysis+1
          WRITE(UNIT=tmp_unit_nb,FMT=*) mc_counter, exp_counter-1, npsol_inform_step_2(mc_counter, exp_counter), &
               npsol_iter_step_2(mc_counter,exp_counter), mc_results_step_2(:,mc_counter,exp_counter)
       END DO
    END DO

    RETURN
  END SUBROUTINE OUTPUT_EXPECTATION_MC_STEP_2
END MODULE OUTPUT_UTILITIES














