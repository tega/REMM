MODULE DATATYPES
  IMPLICIT  NONE

  INTRINSIC SELECTED_INT_KIND, SELECTED_REAL_KIND
  !Per Manno utilizzare questi parametri.
  !INTEGER , PARAMETER :: intrkind = SELECTED_INT_KIND (R = 15)
  !INTEGER , PARAMETER :: realkind = SELECTED_REAL_KIND(P = 15)

  !Per casa utilizzare questi parametri
  INTEGER , PARAMETER :: intrkind = SELECTED_INT_KIND (R = 9)
  INTEGER , PARAMETER :: realkind = SELECTED_REAL_KIND(P = 9)  

  PRIVATE
  PUBLIC :: intrkind, realkind
END MODULE DATATYPES

MODULE DEVEL
   USE DATATYPES
   IMPLICIT NONE
   INTEGER :: count !count viene inizializzato in setup
    
END MODULE DEVEL

MODULE S1_COMMON
  USE DATATYPES
  IMPLICIT NONE

  !common variables defining the nb of observations and the time series containing them
  INTEGER :: T1
  INTEGER :: N1      !the nb. of observations in step1 without the starting values (lags)
  INTEGER :: T_max
  REAL (KIND=realkind), DIMENSION(:), ALLOCATABLE :: imported_series
  REAL (KIND=realkind), DIMENSION(:), ALLOCATABLE :: modified_series

  !common variables defining the conditional mean mu of the markov process
  INTEGER :: nb_lags_mu   !number of lags in the conditional mean
  INTEGER :: dim_mu       !number of lags + 1 (for the constant)
  INTEGER :: do_grid_search !variable used to indicate for which method the grid
                            !search is done
  INTEGER :: do_grid_search_first_n_loops
  
  !common variables defining the joint density function
  INTEGER :: pol_dim !Dimension of the polynomial
  INTEGER :: nb_pol_coef
  INTEGER :: dim_sigma
  INTEGER :: dim_arch
  INTEGER :: dim_garch
  INTEGER, ALLOCATABLE, DIMENSION(:) :: pol_deg !=(Maximal exponent in the i-th dimension)
  REAL (KIND=realkind) :: epsilon_0             !the value of the polynomial epsilon
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: M

  !common variables used to compute the conditional density function
  INTEGER :: deg_lag_max !Maximal exponent between all exponents from 
  !the 2-nd to the pol_dim-th dimension
  INTEGER :: deg_z_max   !Maximal exponent between all exponents from 
  !the 1-st pol dimension
  INTEGER :: dim_rows_X  !Number of rows of the X matrix containing the 
  !trasformed "lag y" and its powers
  INTEGER :: t_start     !First time at which the conditional density is computed

  !common variables containing all the estimation parameters
  INTEGER :: AM_nb_par
  INTEGER :: AM_nb_not_fixed_par
  INTEGER :: AM_nb_free_par
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: beta_AM_start
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: beta_AM_end

  !common variables used to compute the matrix of variance and covariance of the gradient,
  !the pseudo_hessian matrix
  INTEGER :: sel_score_cov_mat_s1  !=2 is Newey-West: =3 is Gallant 
  INTEGER :: sel_cov_estimates_s1
  INTEGER :: sel_hessian_mat_s1
  INTEGER :: lags_score_cov_mat_s1
  REAL (KIND=realkind) :: delta_M_beta
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: M_beta
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: grad_i
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: score_cov
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: inv_score_cov
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: cov_estimates_s1

  !common variables used to compute the number of free parameters and their location
  !among the mu, sigma2 and theta vectors which contains the parameters of the mean,
  !variance and polynomial functions
  INTEGER :: nb_not_fixed_mu       !number of free parameters in the conditional mean part
  INTEGER :: nb_not_fixed_sigma    !number of free parameters in the conditional variance part
  INTEGER :: nb_not_fixed_pol      !number of free parameters in the polynomial part
  INTEGER, ALLOCATABLE, DIMENSION(:) :: j1 !vector used to store the indices of the not excluded parameters
  !                                        !in the conditional mean part
  INTEGER, ALLOCATABLE, DIMENSION(:) :: j2 !idem for the condit. variance
  INTEGER, ALLOCATABLE, DIMENSION(:) :: j3 !idem for the polynomial part
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_fixed_par_AM 
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_restricted_par_AM
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_gt_par_AM
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_lt_par_AM
  REAL (KIND=realkind), PARAMETER :: restriction_correction_factor = 0.0000001_realkind
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: val_fix_restr_par_AM
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: gt_par_AM_value
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: lt_par_AM_value

  !variables used for the linear constraints
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_gt_lin_constr_AM
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_lt_lin_constr_AM
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: A_lin_constr_AM
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: gt_lin_constr_AM_value
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: lt_lin_constr_AM_value

  !common variables used by the NPSOL routine
  LOGICAL, DIMENSION(4) :: use_npsol
  INTEGER :: nclin1
  INTEGER :: ncnln1
  INTEGER :: nclin
  INTEGER :: ncnln
  INTEGER :: ldA
  INTEGER :: ldJ
  INTEGER :: ldR
  INTEGER :: nctotl
  INTEGER :: leniw
  INTEGER :: lenw
  INTEGER, ALLOCATABLE, DIMENSION(:) :: istate
  INTEGER, ALLOCATABLE, DIMENSION(:) :: iw
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: bl 
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: bu
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: clamda
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: w 
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: A
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: cJac
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: RR
  INTEGER :: inform   ! not necessary to be initialized in the setup
  INTEGER :: iter     ! not necessary to be initialized in the setup
  INTEGER, ALLOCATABLE, DIMENSION(:) :: npsol_inform_step_1
  INTEGER, ALLOCATABLE, DIMENSION(:) :: npsol_iter_step_1
  LOGICAL :: npsol_first_call

  !optional input parameters in the NPSOL routine: e indicates the machine epsilon!
  REAL (KIND=realkind) :: infinite_bound_size     !bound used to bound the parameters. 10^20
  REAL (KIND=realkind) :: bigbnd_default
  REAL (KIND=realkind) :: bigbnd
  REAL (KIND=realkind) :: feasibility_tolerance   !max. acceptable absolute violation in lin.
  !                                                and nonlin. constraints. e^0.5
  REAL (KIND=realkind) :: crash_tolerance
  REAL (KIND=realkind), DIMENSION(1:4) :: &
                          step_limit              !the max. change in variables at the first step of line search
  REAL (KIND=realkind), DIMENSION(1:4) :: &
                          function_precision      !=:eR --> measure of the accuracy in the
  !                                                objective function computation. e^0.9
  REAL (KIND=realkind), DIMENSION(1:4) :: &
                          optimality_tolerance    !accuracy in the approximation of the final
  !                                                solution of the problem. eR^0.8
  REAL (KIND=realkind), DIMENSION(1:4) :: &
                          line_search_tolerance
  INTEGER :: minor_iteration_limit                !max. number of iterations for the optim.
  !                                                phase of each QP subproblem
  INTEGER :: major_iteration_limit                !max. number of maj. iterations permitted
  !                                                before termination
  INTEGER :: major_print_level
  INTEGER :: minor_print_level
  INTEGER :: verify_gradients
  INTEGER :: derivative_level                     !which derivatives are provided by the user
  INTEGER :: Nolist                               !used to suppress the option listing
  INTEGER :: Print_File                           !used to suppress the option listing

  !other common variables and parameters
  CHARACTER (LEN=13) :: version
  INTEGER :: center
  INTEGER :: scale
  INTEGER :: run_number
  INTEGER :: step     
  INTEGER :: estimation_kind
  INTEGER :: simulate_data
  INTEGER :: space_lenght
  INTEGER :: size_seed
  INTEGER, ALLOCATABLE, DIMENSION(:) :: seed
  LOGICAL :: expectation_analysis
  LOGICAL :: write_output_1   !often associated with tmp_unit_nb, i.e. the results of a monte carlo simulation
  LOGICAL :: write_output_2   !often associated with log_unit_nb, i.e. setup.dat file
  LOGICAL :: write_output_3   !often associated with dev_unit_nb, generally used in the developpement phase 
  INTEGER, PARAMETER :: tmp_unit_nb = 14    !unit number used for temporary input output operations
  INTEGER, PARAMETER :: log_unit_nb = 15    !unit number used for the log file
  INTEGER, PARAMETER :: dev_unit_nb = 17    !unit number used for the development output file
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: minimum_step_1_vector
  CHARACTER (LEN=10) :: time_start, time_end
  CHARACTER (LEN=8)  :: date_start, date_end

  REAL (KIND=realkind) :: normal_c                   
  REAL (KIND=realkind) :: empirical_is_mean
  REAL (KIND=realkind) :: empirical_is_variance
  REAL (KIND=realkind) :: empirical_ms_mean
  REAL (KIND=realkind) :: empirical_ms_variance
  REAL (KIND=realkind) :: feasibil_toler
  REAL (KIND=realkind), PARAMETER :: pi = 3.14159265358979323846_realkind
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:,:) :: exp_analysis_results

  SAVE  
END MODULE S1_COMMON

MODULE Z_MATRIX
  USE DATATYPES
  IMPLICIT NONE
  !in this module are saved the last values of the arguments used in the
  !Z_MATRIX_CONSTRUCTION subroutine

  INTEGER :: series_old  !1=first  step (used time_series -> T1) 
  !                       2=second step (used simulated_series -> T2)
  !                       3=robust step (used simulated_series -> T3)
  INTEGER :: is_new1     !=1 if the time_series has been modified, else 0
  INTEGER :: is_new2     !=1 if the simulated_series has been modified, else 0
  INTEGER :: is_new3     !=1 if the simulated_series for rob. has been modified, else 0
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: beta_old

  !common variables used to compute the arch/garch
  REAL (KIND=realkind) :: U_var_uncond 
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: U_square

  !some variables used to compute the density
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nb_rep
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nb_jumps
  INTEGER, ALLOCATABLE, DIMENSION(:) :: high_jumps
  REAL (KIND=realkind) :: constant_denumerator
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:)   :: numerator
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:)   :: denumerator
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:)   :: U_temp
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:)   :: sig2
  INTEGER :: sig2_start !the first component of the dimension of sig2 
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:)   :: inv_sigma
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: S_Z
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: X
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: X_product
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: b

  !common variables used to compute and store the gradient
  !fist without polynomial part
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: grad_l
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:)   :: grad_l_beta
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:)   :: grad_l_sigma
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:)   :: grad_const
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:)   :: sum_grad_beta_U2
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: grad_beta_U2 
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: grad_beta_sig2
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: grad_sigma_sig2
  !then with polynomial part
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:)   :: grad_l_theta
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:)   :: pol
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:)   :: grad_theta_denumerator
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:)   :: fixed_grad_P
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: grad_beta_pol
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: grad_sigma_pol  
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: grad_theta_pol
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: grad_beta_Z
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: grad_theta_qMq

  SAVE
END MODULE Z_MATRIX

MODULE S2_COMMON
  USE DATATYPES
  IMPLICIT NONE

  !The common variables of the structural model or for short SM
  INTEGER :: T2
  INTEGER :: N2 !the number of observation in the simulation step without the starting values (lags)
  INTEGER :: max_T2
  INTEGER :: SM_model_parts
  INTEGER :: nb_iid_series
  INTEGER :: SM_nb_par
  INTEGER :: SM_nb_not_fixed_par
  INTEGER :: SM_nb_free_par
  INTEGER, ALLOCATABLE, DIMENSION(:) :: SM_dim_parts
  INTEGER :: nb_discard_step_2
  LOGICAL :: change_random_iid
  CHARACTER (LEN=15) :: name_SM
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: rho_SM_start
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: rho_SM_end
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: rho_rest_SM_start
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: rho_rest_SM_end
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: rho_rest_SM_old
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: minimum_step_2_vector

  LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_fixed_par_SM
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_restricted_par_SM
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_gt_par_SM
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_lt_par_SM
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: val_fix_restr_par_SM
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: gt_par_SM_value
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: lt_par_SM_value

  !variables used for the linear constraints
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_gt_lin_constr_SM
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_lt_lin_constr_SM
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: A_lin_constr_SM
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: gt_lin_constr_SM_value
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: lt_lin_constr_SM_value

  !the variables for the NPSOL
  INTEGER :: nclin2
  INTEGER :: ncnln2

  !the variables used for the monte carlo integral approximation 
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: random_iid_step_2
  REAL (KIND=realkind), DIMENSION(:), ALLOCATABLE :: simulated_series_step_2

  !the vector containing the hansen specification test results
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: hansen_test_results
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: npsol_inform_step_2
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: npsol_iter_step_2

  !the matrix M_rho and the epsilon used to approximate it (numerical derivative)
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: M_rho !the M_rho matrix used to define the metric for the weights
  REAL (KIND=realkind) :: delta_M_rho
  INTEGER :: delta_method_M_rho
  INTEGER :: method_derivative_M_rho
  LOGICAL :: compute_accuracy_M_rho

  SAVE
END MODULE S2_COMMON

MODULE ROB_COMMON
  USE DATATYPES
  IMPLICIT NONE

  INTEGER :: weights_metric
  LOGICAL, DIMENSION(1:2) :: rob_start_value
  INTEGER :: rob_iteration_limit  !the maximal number of iterations before the robust algorithm
  INTEGER :: A_updating_kind      !it determines if the A_matrix must be estimated by simulation (=2)
                                  !or with the oberved data (=1)
  INTEGER :: nb_A_updates
  INTEGER :: nb_A_updates_start
  INTEGER :: nb_B_updates
  INTEGER :: nb_M_rho_updates
  INTEGER :: wait_update_M_rho
  INTEGER :: count_update_M_rho   !it counts the number of times the procedure updating M_rho
                                  !is called
  INTEGER :: M_rho_delta_selection

  
  REAL (KIND=realkind) :: rob_parameter_tolerance !accuracy in the approximation of the paramet. estimates.
  REAL (KIND=realkind) :: bound
  REAL (KIND=realkind) :: minimum_step_1_rob, minimum_step_2_rob
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: minimum_step_1_rob_vector
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: minimum_step_2_rob_vector
 
  !the scaling matrix of the psi function
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: rob_A       !for the method with orthogonal projection
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: rob_A_tilde !for the correct metric
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: B_0 !the asymptotic variance_covariance matrix of the mean of the truncated scores 
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: S !the weighting matrix S used in the second step of the EMM 
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: rob_weights_N1 !the weights for the real data
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: rob_weights_N2 !the weights of the simulated data
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: psi_N1 !matrix used in computations
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: psi_N2 !matrix used in computations
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: psi_N1_tilde !matrix used in computations
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: psi_N2_tilde !matrix used in computations
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: score_N1 !matrix used in computations
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: score_N2 !matrix used in computations
  !the variable used for the update of alpha
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: random_iid_rob_update

  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: hansen_test_results_rob
  INTEGER, ALLOCATABLE, DIMENSION(:) :: npsol_inform_step_1_rob
  INTEGER, ALLOCATABLE, DIMENSION(:) :: npsol_iter_step_1_rob
  INTEGER, ALLOCATABLE, DIMENSION(:) :: npsol_inform_step_2_rob
  INTEGER, ALLOCATABLE, DIMENSION(:) :: npsol_iter_step_2_rob
  INTEGER, ALLOCATABLE, DIMENSION(:) :: inform_rob
  INTEGER, ALLOCATABLE, DIMENSION(:) :: iter_rob

  SAVE
END MODULE ROB_COMMON

MODULE DGP_DEFINITION
  USE DATATYPES
  IMPLICIT NONE

  CHARACTER (LEN=15) :: DGP_name              !name of the selected model
  INTEGER            :: DGP_nb_par            !total nb of parameters of the model
  INTEGER            :: DGP_model_dim         !dimensions of the model, i.e. AR-MA is 2
  INTEGER            :: DGP_nb_iid_series     !nb of iid series used to simulate the model
  INTEGER            :: DGP_nb_discard        !nb of initial observations discarded

  CHARACTER (LEN=7), ALLOCATABLE, DIMENSION(:)    :: DGP_innov_dist
  INTEGER, ALLOCATABLE, DIMENSION(:)              :: DGP_deg_freedom
  INTEGER, ALLOCATABLE, DIMENSION(:)              :: DGP_order_model_dim
  INTEGER, ALLOCATABLE, DIMENSION(:)              :: DGP_nb_par_vec
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: DGP_par_vec

  !variables used for the SNP density
  INTEGER :: DGP_pol_dim
  INTEGER, ALLOCATABLE, DIMENSION(:) :: DGP_pol_deg
  INTEGER :: DGP_nb_pol_coef
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: DGP_pol_coef
  REAL (KIND=realkind) :: DGP_epsilon_0

  SAVE
END MODULE DGP_DEFINITION

MODULE CONTAMINATION_MODEL
  USE DATATYPES
  IMPLICIT NONE

  INTEGER, ALLOCATABLE, DIMENSION(:) :: contam_model
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: epsilon_p
  CHARACTER (LEN=7), ALLOCATABLE, DIMENSION(:)    :: contam_dist
  INTEGER, ALLOCATABLE, DIMENSION(:)              :: nb_contam_points  
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: contam_x0
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: contam_lower_bound
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: contam_upper_bound
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: contam_norm_mean
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: contam_norm_var
  INTEGER, ALLOCATABLE, DIMENSION(:) :: contam_df
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: contam_student_mean
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: contam_student_var

  SAVE
END MODULE CONTAMINATION_MODEL

MODULE MONTE_CARLO
  USE DATATYPES
  IMPLICIT NONE

  !common variables used in monte carlo
  INTEGER :: mc_counter
  INTEGER :: estimate_all_mc_simulations
  INTEGER :: nb_rep_monte_carlo
  INTEGER :: nb_selected_simulations
  INTEGER :: nb_rep_exp_analysis
  INTEGER :: skip_first_n
  INTEGER, ALLOCATABLE, DIMENSION(:) :: selected_simulations

  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: mc_results_step_1
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: mc_results_step_1_rob
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: mc_results_step_2_rob
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:,:) :: mc_results_step_2
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:,:) :: inv_score_cov_mc

  SAVE
END MODULE MONTE_CARLO


MODULE GRID_SEARCH_AM
   USE DATATYPES
   IMPLICIT NONE

   INTEGER :: nb_grid_searches_AM
   INTEGER :: tot_nb_grid_points_AM
   INTEGER, ALLOCATABLE, DIMENSION(:) :: nb_grid_points_AM
   INTEGER, ALLOCATABLE, DIMENSION(:) :: nb_grid_intervals_AM
   INTEGER, ALLOCATABLE, DIMENSION(:) :: cum_grid_points_AM
   REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: lb_grid_search_AM
   REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: ub_grid_search_AM
   REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: grid_points_values_AM
   REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: grid_search_results_AM
   LOGICAL :: do_grid_search_AM
   LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_grid_search_AM
   LOGICAL :: grid_search_AM_successfull

   SAVE
END MODULE GRID_SEARCH_AM


MODULE GRID_SEARCH_SM
   USE DATATYPES
   IMPLICIT NONE

   INTEGER :: nb_grid_searches_SM
   INTEGER :: tot_nb_grid_points_SM
   INTEGER, ALLOCATABLE, DIMENSION(:) :: nb_grid_points_SM
   INTEGER, ALLOCATABLE, DIMENSION(:) :: nb_grid_intervals_SM
   INTEGER, ALLOCATABLE, DIMENSION(:) :: cum_grid_points_SM
   REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: lb_grid_search_SM
   REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: ub_grid_search_SM
   REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: grid_points_values_SM
   REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: grid_search_results_SM
   LOGICAL :: do_grid_search_SM
   LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_grid_search_SM
   LOGICAL :: grid_search_SM_successfull

   SAVE
END MODULE GRID_SEARCH_SM


MODULE MYFUNCTIONS
  USE DATATYPES
  IMPLICIT NONE

CONTAINS

  FUNCTION CHOL(A,n1) RESULT (L)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n1
    REAL (KIND=realkind), INTENT(IN), DIMENSION(:,:) :: A
    REAL (KIND=realkind), DIMENSION(1:n1,1:n1) :: L

    REAL (KIND=realkind) :: temp
    INTEGER :: i,k,s


    !This procedure returns a lower diagonal matrix L such that
    !LL'= A
    L = 0._realkind
    L(1,1) = SQRT(A(1,1))
    DO i=2,n1
       L(i,1) = A(i,1)/L(1,1)
    END DO

    DO k=2,n1
       temp=0._realkind
       DO s=1,k-1
          temp = temp + L(k,s)*L(k,s)
       END DO
       L(k,k)=SQRT(A(k,k)-temp)
       IF (k .EQ. n1) THEN
          EXIT
       END IF

       DO i=k+1,n1
          temp = 0.0_realkind
          DO s=1, k-1
             temp = temp + L(i,s)*L(k,s)
          END DO
          L(i,k)=(A(i,k)-temp)/L(k,k)
       END DO
    END DO

    RETURN
  END FUNCTION CHOL


  FUNCTION EYE(k) RESULT (I)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: k

    REAL (KIND=realkind), DIMENSION(1:k,1:k) :: I
    INTEGER :: m, n

    DO m=1,k
       DO n=1,k
          I(n,m)=0._realkind
       END DO
       I(m,m)=1._realkind
    END DO

    RETURN
  END FUNCTION EYE


  FUNCTION MM(m1,m2,n1,n2,n3) RESULT (matrix)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n1,n2,n3
    REAL (KIND=realkind), INTENT(IN), DIMENSION(:,:) :: m1, m2
    REAL (KIND=realkind), DIMENSION(1:n1,1:n3) :: matrix

    INTEGER :: i, j, l

    matrix = 0._realkind

    if (n1 >= n2) then
       DO l=1,n2
          DO j=1,n3
             DO i=1,n1
                matrix(i,j) = matrix(i,j)+m1(i,l)*m2(l,j)
             END DO
          END DO
       END DO
    ELSE
       !CDIR NOLOOPCHG
       DO i=1,n1
          DO j=1,n3
             DO l=1,n2
                matrix(i,j) = matrix(i,j)+m1(i,l)*m2(l,j)
             END DO
          END DO
       END DO
    END IF

    RETURN
  END FUNCTION MM


  FUNCTION MMM(m1,m2,m3,n1,n2,n3,n4) RESULT (matrix)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    REAL (KIND=realkind), INTENT(IN), DIMENSION(:,:) :: m1, m2, m3
    REAL (KIND=realkind), DIMENSION(1:n1,1:n4) :: matrix

    REAL (KIND=realkind), DIMENSION(1:n1,1:n3) :: matrix_temp

    matrix_temp = MM(m1,m2,n1,n2,n3)
    matrix = MM(matrix_temp,m3,n1,n3,n4)

    RETURN
  END FUNCTION MMM


  FUNCTION PROD_C(array_1,array_2,n1,n2) RESULT (array_res)
    !Calcola il prodotto componente per componente fra due matrici della medesima dimensione.
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n1,n2
    REAL (KIND=realkind), INTENT(IN), DIMENSION(:,:) :: array_1, array_2
    REAL (KIND=realkind), DIMENSION(1:n1,1:n2) :: array_res

    INTEGER :: i,j

    !Inserire un controllo della dimensione!
    DO j=1,n2
       DO i=1,n1
          array_res(i,j) = array_1(i,j)*array_2(i,j)
       END DO
    END DO

    RETURN
  END FUNCTION PROD_C


  !MB\
  SUBROUTINE preset_rnd_numbers(t, v1_v, v2_v, r_v)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: t
    REAL (KIND=realkind), DIMENSION(1:t), INTENT(OUT) :: v1_v, v2_v, r_v

    REAL (KIND=realkind), DIMENSION(1:t) :: v1_tmp, v2_tmp
    INTEGER, DIMENSION(1:t) :: done
    INTEGER :: cnt, i

    DONE = 0
    cnt = 0

    DO WHILE (cnt < t)
       CALL RANDOM_NUMBER(v1_tmp(1:t))
       CALL RANDOM_NUMBER(v2_tmp(1:t))
       DO i = 1, t
          IF(done(i) == 0) THEN
             v1_v(i) = 2.0_realkind*v1_tmp(i)-1.0_realkind
             v2_v(i) = 2.0_realkind*v2_tmp(i)-1.0_realkind
             r_v(i) = v1_v(i)**2 + v2_v(i)**2
             IF ((r_v(i) .GT. 0.0_realkind) .AND. (r_v(i) .LT. 1.0_realkind)) THEN
                done(i) = 1
             END IF
          END IF
       END DO
       cnt = SUM(done)
    END DO

    RETURN
  END SUBROUTINE preset_rnd_numbers
  !MB/


  FUNCTION RND_N(t,n) RESULT (z)
    ! Normal (Gaussian) Deviates
    ! p.202 Numerical Recipes
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: t,n

    INTEGER ::  i,j,q,s,t_is_odd
    REAL (KIND=realkind), DIMENSION(1:t,1:n) ::  z
    REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: z_tmp,v1_v, v2_v, r_v, fac

    IF (t/2.0_realkind-t/2 .EQ. 0.0_realkind) THEN
       t_is_odd = 0
    ELSE
       t_is_odd = 1
    END IF

    s=t+t_is_odd
    q=t-s/2    
    ALLOCATE(v1_v(1:s),v2_v(1:s),r_v(1:s),z_tmp(1:q),fac(1:q))

    DO j = 1, n
       CALL preset_rnd_numbers(s, v1_v, v2_v, r_v)

       IF (t_is_odd .EQ. 0) THEN
          DO i = 1, q
             fac(i) = SQRT(-2.0_realkind*LOG(r_v(i))/r_v(i))
             z(i,j)   = v1_v(i)*fac(i)
             z_tmp(i) = v2_v(i)*fac(i)
          END DO
          z((1+q):t,j)=z_tmp(1:q)
       ELSE
          DO i = 1, q
             fac(i) = SQRT(-2.0_realkind*LOG(r_v(i))/r_v(i))
             z(i,j)   = v1_v(i)*fac(i)
             z_tmp(i) = v2_v(i)*fac(i)
          END DO
          fac(i) = SQRT(-2.0_realkind*LOG(r_v(i))/r_v(i))
          z(i,j)   = v1_v(i)*fac(i)
          z((2+q):t,j)=z_tmp(1:q)
       END IF
    END DO

    DEALLOCATE(v1_v,v2_v,r_v,z_tmp,fac)
    RETURN
  END FUNCTION RND_N


  FUNCTION RND_CHI2(m,n,gl) RESULT (chi2)
    ! Ritorna una matrice di dimensione (m,n) di variabili aleatorie
    ! indipendenti distribuite secondo una chi2 a gl gradi di liberta'.
    ! m = number of rows, n = number of columns
    ! gl = degree of freedom
    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: m,n,gl
    REAL (KIND=realkind), DIMENSION(1:m,1:n) :: chi2

    REAL (KIND=realkind), DIMENSION(1:m,1:n) :: y
    INTEGER :: i

    chi2=0._realkind
    DO i=1,gl
       y = RND_N(m,n)
       y = PROD_C(y,y,m,n)
       chi2=chi2+y
    END DO

    RETURN
  END FUNCTION RND_CHI2


  FUNCTION RND_STD(m,n,gl) RESULT (std)
    ! Ritorna una matrice di dimensione (m,n) di variabili aleatorie
    ! indipendenti distribuite secondo una student a gl gradi di liberta'.
    ! m = number of rows, n = number of columns
    ! gl = degree of freedom
    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: m,n,gl
    REAL (KIND=realkind), DIMENSION(1:m,1:n) :: std

    REAL (KIND=realkind), DIMENSION(1:m,1:n) :: chi2
    REAL (KIND=realkind), DIMENSION(1:m,1:n) :: normal
    INTEGER :: i,j

    chi2 = RND_CHI2(m,n,gl)
    normal = RND_N(m,n)
    DO j=1,n
       DO i=1,m
          std(i,j)=normal(i,j)/sqrt(chi2(i,j)/gl)
       END DO
    END DO

    RETURN
  END FUNCTION RND_STD

  FUNCTION FACTORIAL(n) RESULT(fact_result)
    IMPLICIT NONE

    REAL (kind=realkind) ::  fact_result, real_i
    INTEGER, INTENT(IN)  ::  n

    INTEGER :: i

    fact_result = 1._realkind
    real_i = 2._realkind  

    DO i=2,n
       fact_result = fact_result * real_i 
       real_i = real_i + 1._realkind
    END DO
    RETURN
  END FUNCTION FACTORIAL


  SUBROUTINE M_MATRIX_CONSTRUCTION(n,M)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL (kind=realkind),  INTENT(IN OUT), DIMENSION(1:n,1:n) :: M

    INTEGER :: i, j, is_zero 

    M=0._realkind

    DO j=1, n
       M(j,j) = FACTORIAL(j+j-2)/(2._realkind**(j-1)*FACTORIAL(j-1))
       is_zero = 0
       DO i=j+1, n
          IF (is_zero .EQ. 1) THEN
             M(i,j) = FACTORIAL(i+j-2)/(2._realkind**((i+j)/2-1)*FACTORIAL((i+j)/2-1))
             M(j,i) = FACTORIAL(i+j-2)/(2._realkind**((i+j)/2-1)*FACTORIAL((i+j)/2-1))
             is_zero = 0
          ELSE
             M(i,j) = 0._realkind
             M(j,i) = 0._realkind 
             is_zero = 1
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE M_MATRIX_CONSTRUCTION


  FUNCTION INV_NORMAL_DIST_FUNC_APPROX_1(uniform, mu, sigma) RESULT (F_inv)
    !This function computes the inverse of the Normal distribution function
    !with a an error < .45E-03.See Fishman, page 151.
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN) :: uniform, mu, sigma
    REAL (KIND=realkind) :: F_inv

    REAL (KIND=realkind), PARAMETER :: c0 = 2.515517_realkind
    REAL (KIND=realkind), PARAMETER :: c1 = 0.802853_realkind
    REAL (KIND=realkind), PARAMETER :: c2 = 0.010328_realkind
    REAL (KIND=realkind), PARAMETER :: d1 = 1.432788_realkind
    REAL (KIND=realkind), PARAMETER :: d2 = 0.189269_realkind
    REAL (KIND=realkind), PARAMETER :: d3 = 0.001308_realkind

    REAL (KIND=realkind) :: t

    t = (-LOG((1.0_realkind-ABS(1.0_realkind-2.0*REAL(uniform))/2.0_realkind)**2))**0.5
    F_inv = mu + (uniform-0.5_realkind) / ABS(uniform-0.5_realkind) * sigma * &
         (t - (c0 + c1*t + c2*t**2)/(1 + d1*t + d2*t**2 + d3*t**3))
    RETURN
  END FUNCTION INV_NORMAL_DIST_FUNC_APPROX_1


  FUNCTION NORMAL_DIST_FUNC_APPROX_1(x, mu, sigma) RESULT (F)
    !This function computes the distribution function of the normal distribution
    !with a an error < 1.E-07.See Lamberton-Lapeyre, page 168.
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN) :: x, mu, sigma
    REAL (KIND=realkind) :: F

    REAL (KIND=realkind), PARAMETER :: p = 0.2316419_realkind
    REAL (KIND=realkind), PARAMETER :: b1 = 0.31938153_realkind
    REAL (KIND=realkind), PARAMETER :: b2 = -0.356563782_realkind
    REAL (KIND=realkind), PARAMETER :: b3 = 1.781477937_realkind
    REAL (KIND=realkind), PARAMETER :: b4 = -1.821255978_realkind
    REAL (KIND=realkind), PARAMETER :: b5 = 1.330274429_realkind
    REAL (KIND=realkind), PARAMETER :: pi = 3.14159265358979323846_realkind
    REAL (KIND=realkind) :: t, y

    y = (x-mu)/sigma
    t = 1._realkind/(1._realkind + p*y)
    IF (y .GE. 0) THEN
       t = 1._realkind/(1._realkind + p*y)
       F = 1.0_realkind - 1.0_realkind/SQRT(2.0_realkind*pi)*EXP(-0.5*y**2)* &
            (b1*t + b2*t**2 + b3*t**3 +b4*t**4 + b5*t**5)
    ELSE
       t = 1._realkind/(1._realkind - p*y)
       F = 1.0_realkind/SQRT(2.0_realkind*pi)* &
            EXP(-0.5*y**2)*(b1*t + b2*t**2 + b3*t**3 +b4*t**4 + b5*t**5)
    END IF
    RETURN
  END FUNCTION NORMAL_DIST_FUNC_APPROX_1


  FUNCTION SNP_DIST_FUNC_APPROX(x, theta, n, epsilon_0) RESULT (F)
    !This function computes the distribution function of the univariate snp 
    !distribution as given by Gallant & Nichka. theta is the vector of coefficient
    !given in P(x|theta)**2 + espilon_0.
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN) :: x, epsilon_0
    INTEGER, INTENT(IN) :: n !n-1 is the degree of the polynomial so that n is the nb
    !of polynomial coefficients.
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: theta
    REAL (KIND=realkind) :: F

    INTEGER :: i,j
    REAL (KIND=realkind) :: temp, constant
    REAL (KIND=realkind), DIMENSION(1:2*n-1) :: a, b
    REAL (KIND=realkind), PARAMETER :: pi = 3.14159265358979323846_realkind
    REAL (KIND=realkind), DIMENSION(1:n,1:n) :: M

    CALL M_MATRIX_CONSTRUCTION(n,M)
    constant = epsilon_0
    DO j=1, n
       constant = constant + M(j,j)*theta(j)**2
       DO i=j+1, n
          constant = constant + 2*M(i,j)*theta(j)*theta(i)
       END DO
    END DO

    !The vector a contains the coefficients of the polinomial P(x|theta)**2    
    a = 0.0_realkind
    DO i=1, n
       DO j=1, n
          a(i+j-1) = a(i+j-1) + theta(i)*theta(j)
       END DO
    END DO

    !a(1) contains the constant of the polynomial. We have to add epsilon_0.
    a(1) = a(1) + epsilon_0

    b = a
    temp = 0.0_realkind
    IF (n .GT. 1) THEN
       DO i=2*(n-1), 2, -1
          temp = temp - b(i+1)*x**(i-1)
          b(i-1) = b(i-1) + b(i+1)*(i-1)
       END DO
       temp = temp-b(2)
       temp = temp*(2*pi)**(-0.5)*EXP(-0.5*x**2)
    END IF
    temp = temp + b(1)*NORMAL_DIST_FUNC_APPROX_1(x, 0.0_realkind, 1.0_realkind)

    F = temp/constant
    RETURN
  END FUNCTION SNP_DIST_FUNC_APPROX

  FUNCTION SNP_DENSITY_FUNC(x, theta, n, epsilon_0) RESULT (f)
    !This function computes the density function of the univariate snp 
    !distribution as given by Gallant & Nichka. theta is the vector of coefficient
    !given in P(x|theta)**2 + espilon_0.
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN) :: x, epsilon_0
    INTEGER, INTENT(IN) :: n !n-1 is the degree of the polynomial so that n is the nb
    !of polynomial coefficients.
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: theta
    REAL (KIND=realkind) :: f

    INTEGER :: i,j
    REAL (KIND=realkind) :: constant
    REAL (KIND=realkind), PARAMETER :: pi = 3.14159265358979323846_realkind
    REAL (KIND=realkind), DIMENSION(1:n,1:n) :: M

    CALL M_MATRIX_CONSTRUCTION(n,M)
    constant = epsilon_0
    DO j=1, n
       constant = constant + M(j,j)*theta(j)**2
       DO i=j+1, n
          constant = constant + 2*M(i,j)*theta(j)*theta(i)
       END DO
    END DO

    !computation of P(x|theta)    
    f = 0.0_realkind
    DO i=1, n
       f = f + theta(i)*x**(i-1)
    END DO
    f = f**2 + epsilon_0
    f = f*(2*pi)**(-0.5)*EXP(-0.5*x**2)   
    f = f/constant
    RETURN
  END FUNCTION SNP_DENSITY_FUNC


  FUNCTION INV_NORMAL_DIST_FUNC_APPROX_2(uniform, mu, sigma) RESULT (F_inv)
    !This function computes the inverse of the Normal distribution function
    !by approximation using NORMAL_DIST_FUNC_APPROX_1.
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN) :: uniform, mu, sigma
    REAL (KIND=realkind) :: F_inv

    REAL (KIND=realkind) :: epsilon, stepp, y, temp
    INTEGER :: index
    epsilon = 0.00000001_realkind
    stepp = 0.1_realkind   

    y = mu
    IF ((NORMAL_DIST_FUNC_APPROX_1(y,mu,sigma)-uniform) .GT. 0.0_realkind) THEN
       index = 1
       y = y - sigma*stepp 
    ELSE
       index = 0
       y = y + sigma*stepp
    END IF

    DO
       temp = NORMAL_DIST_FUNC_APPROX_1(y,mu,sigma)
       IF (ABS(temp-uniform) .LT. epsilon)  EXIT
       IF ((temp-uniform) .GT. 0.0_realkind) THEN
          IF (index .EQ. 1) THEN
             y = y - sigma*stepp
          ELSE
             stepp = stepp/10.0_realkind
             y = y - sigma*stepp
             index = 1
          END IF
       ELSE
          IF (index .EQ. 0) THEN
             y = y + sigma*stepp
          ELSE
             stepp = stepp/10.0_realkind
             y = y + sigma*stepp
             index = 0
          END IF
       END IF
    END DO
    F_inv = y
    RETURN
  END FUNCTION INV_NORMAL_DIST_FUNC_APPROX_2

  FUNCTION INV_SNP_DIST_FUNC(uniform, theta, n, epsilon_0) RESULT (F_inv)
    !This function computes the inverse of the SNP distribution function
    !by approximation using SNP_DIST_FUNC_APPROX.
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL (KIND=realkind), INTENT(IN) :: uniform, epsilon_0
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: theta
    REAL (KIND=realkind) :: F_inv


    REAL (KIND=realkind) :: epsilon, stepp, y, temp
    INTEGER :: index
    epsilon = 1.0E-9_realkind
    stepp = 1.0_realkind   

    y = 0.0_realkind
    IF ((SNP_DIST_FUNC_APPROX(y, theta, n, epsilon_0)-uniform) .GT. 0.0_realkind) THEN
       index = 1
       y = y - stepp 
    ELSE
       index = 0
       y = y + stepp
    END IF

    DO
       temp = SNP_DIST_FUNC_APPROX(y, theta, n, epsilon_0)
       IF (ABS(temp-uniform) .LT. epsilon)  EXIT
       IF ((temp-uniform) .GT. 0.0_realkind) THEN
          IF (index .EQ. 1) THEN
             y = y - stepp
          ELSE
             stepp = stepp/10.0_realkind
             y = y - stepp
             index = 1
          END IF
       ELSE
          IF (index .EQ. 0) THEN
             y = y + stepp
          ELSE
             stepp = stepp/10.0_realkind
             y = y + stepp
             index = 0
          END IF
       END IF
    END DO
    F_inv = y
    RETURN
  END FUNCTION INV_SNP_DIST_FUNC



  FUNCTION IID_SIM(nb_row,nb_col,distribution_kind,df,theta, &
       nb_pol_coef,epsilon) RESULT (simulated_iid)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nb_row
    INTEGER, INTENT(IN) :: nb_col
    CHARACTER(LEN=*), INTENT(IN) :: distribution_kind
    INTEGER, INTENT(IN), OPTIONAL :: df
    REAL (KIND=realkind), INTENT(IN), DIMENSION(:), OPTIONAL :: theta
    INTEGER, INTENT(IN), OPTIONAL :: nb_pol_coef
    REAL (KIND=realkind), INTENT(IN), OPTIONAL :: epsilon
    REAL (KIND=realkind), DIMENSION(1:nb_row,1:nb_col) :: simulated_iid

    INTEGER :: i, j

    IF (distribution_kind .EQ. "NORMAL_") THEN
       simulated_iid = RND_N(nb_row,nb_col)
    END IF

    IF (distribution_kind .EQ. "STUDENT") THEN
       IF (df .GT. 2) THEN 
          simulated_iid = RND_STD(nb_row,nb_col,df)
          simulated_iid = simulated_iid/SQRT((1._realkind*df)/(df-2))
       END IF
       IF (df .LE. 0) THEN 
          WRITE(*,*) "Error in IID_SIM: distribution_kind is STUDENT but the deg. of"
          WRITE(*,*) "freedom are less or equal 0! Program stopped."
          STOP          
       END IF
       IF (df .LE. 2) THEN
          simulated_iid = RND_STD(nb_row,nb_col,df)
       END IF
    END IF

    IF (distribution_kind .EQ. "UNIFORM") THEN
       CALL RANDOM_NUMBER(simulated_iid(1:nb_row,1:nb_col))
    END IF

    IF (distribution_kind .EQ. "SNP_DST") THEN
       CALL RANDOM_NUMBER(simulated_iid(1:nb_row,1:nb_col))
       DO j=1, nb_col
          DO i=1, nb_row
             simulated_iid(i,j) = INV_SNP_DIST_FUNC(simulated_iid(i,j), & 
                  theta,nb_pol_coef,epsilon) 
          END DO
       END DO
    END IF

    RETURN
  END FUNCTION IID_SIM


  FUNCTION MAT_CROSS_LAG(A,n1,n2,tau) RESULT(C)
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(:,:) :: A
    INTEGER, INTENT(IN) :: n1,n2,tau
    REAL (KIND=realkind), DIMENSION(1:n2,1:n2) :: C
    !This procedure computes the sum of the cross product between the columns of A
    !and those of A with lag tau, i.e.
    !C=(A(t))'*A(t-tau). If the dimension of A is Txr, C will be a rxr matrix.

    INTEGER :: i,j,k

    C=0._realkind
    DO j=1,n2   
       DO k=1,n2
          DO i=1,n1-tau 
             C(j,k)=C(j,k)+A(i+tau,j)*A(i,k)
          END DO
       END DO
    END DO

    RETURN
  END FUNCTION MAT_CROSS_LAG


  FUNCTION INV(A,n) RESULT (B)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL (KIND=realkind), INTENT(IN), DIMENSION(:,:) :: A

    REAL (KIND=realkind), DIMENSION(1:n,1:n) :: B

    INTEGER :: i,j,m
    REAL (KIND=realkind), DIMENSION(1:n,1:2*n) :: AB
    REAL (KIND=realkind), DIMENSION(1:2*n) :: temp

    AB(1:n,1:n)=A
    AB(1:n,n+1:2*n)=EYE(n)
    DO j=1,n !Loop over first n columns of matrix
       !Find row index m of element with largest magnitude.
       m=j         !Start on diagonal
       DO i=j+1,n  !Loop over rows below diagonal
          IF (ABS(AB(i,j))>ABS(AB(m,j))) m=i
       END DO
       !Exchange row m with row j
       temp=AB(j,:)
       AB(j,:)=AB(m,:)
       AB(m,:)=temp
       AB(j,:)=AB(j,:)/AB(j,j) !Divide row j by a(j,j)
       DO i=1,n !Loop over each row of matrix.
          !Subtract a(i,j) times row j from row i.
          IF (i/=j) AB(i,:)=AB(i,:)-AB(i,j)*AB(j,:)
       END DO
    END DO
    B=AB(1:n,n+1:2*n)

    RETURN
  END FUNCTION INV


  SUBROUTINE ZERO_ONE(n,pb,x,m)
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: n 
    REAL (KIND=realkind), INTENT(IN) :: pb 
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:n) :: x 
    INTEGER, INTENT(OUT), OPTIONAL :: m

    INTEGER :: i,s 
    REAL (KIND=realkind), DIMENSION(1:n) :: z 

    CALL RANDOM_NUMBER(z(1:n))
    s = 0
    DO i=1,n 
       IF (z(i) .LE. pb) THEN 
          x(i)=1._realkind
          s = s + 1  
       ELSE
          x(i)=0._realkind 
       END IF
    END DO

    IF (PRESENT(m)) THEN 
       m=s
    END IF

    RETURN 
  END SUBROUTINE ZERO_ONE


  SUBROUTINE CONTAMINATE(T,time_series,epsilon,nb_points,x0,random_x0)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: T
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:T) :: time_series
    REAL (KIND=realkind), INTENT(IN) :: epsilon
    INTEGER, INTENT(IN)  :: nb_points
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:nb_points), OPTIONAL :: x0
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T,1), OPTIONAL :: random_x0

    INTEGER :: i,j
    REAL (KIND=realkind), DIMENSION(1:T) :: ts_temp_1
    REAL (KIND=realkind), DIMENSION(1:T) :: ts_temp_2
    REAL (KIND=realkind), DIMENSION(1:nb_points+1) :: interval

    CALL ZERO_ONE(T,epsilon,ts_temp_1)

    DO i=1, nb_points+1
      interval(i) = REAL(i-1)/REAL(nb_points)
    END DO
    
    CALL RANDOM_NUMBER(ts_temp_2(1:T))

    IF (epsilon .EQ. 0._realkind) THEN

    ELSE
       IF (PRESENT(x0)) THEN
          DO i=1,T
            IF (ts_temp_1(i) .EQ. 1) THEN
               DO j=1,nb_points
                  IF ((ts_temp_2(i) .GE. interval(j)) .AND. &
                     (ts_temp_2(i) .LT. interval(j+1))) THEN
                     time_series(i) = x0(j)
                     EXIT
                  END IF
               END DO
            END IF
          END DO
       ELSE
          DO i=1,T
            IF (ts_temp_1(i) .EQ. 1) THEN
               time_series(i) = random_x0(i,1)
            END IF
          END DO
       END IF
    END IF

    RETURN 
  END SUBROUTINE CONTAMINATE


  SUBROUTINE ARMA_SIMULATION(nb_par,theta,order_ar,order_ma,T,simulated_series, &
       nb_discard, nb_iid_series, iid)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nb_par
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:nb_par) :: theta
    INTEGER, INTENT(IN) :: order_ar,order_ma,T
    REAL (KIND=realkind), INTENT(OUT), DIMENSION(1:T) :: simulated_series
    INTEGER, INTENT(IN) :: nb_discard, nb_iid_series
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T+nb_discard,nb_iid_series) :: iid

    !theta(1) contains the constant term of the ar part
    !theta(2) ... theta(1+order_ar) the AR coeff.,
    !theta(2+order_ar) the variance of the innovation,
    !theta(3+order_ar) ... theta(2+order_ar+order_ma) the MA coeff.

    REAL (KIND=realkind) :: c, stdv, y_start
    REAL (KIND=realkind), DIMENSION((1-order_ma):T+nb_discard) :: epsilon
    REAL (KIND=realkind), DIMENSION(1:T+nb_discard) :: y_sim_temp
    INTEGER :: start, time, i, j

    !compute the mean of the process
    c = 0._realkind
    DO i=2,order_ar+1 
       c = c + theta(i)
    END DO

    y_start = theta(1)/(1._realkind-c)

    !simulate the MA part
    epsilon((1-order_ma):0) = 0._realkind
    stdv = SQRT(theta(2+order_ar))

    DO time=1, T+nb_discard
       epsilon(time) = stdv*iid(time,1)
    END DO

    DO time=1, T+nb_discard
       y_sim_temp(time) = epsilon(time)
    END DO

    DO j=1, order_ma
       DO time=1, T+nb_discard
          y_sim_temp(time) = y_sim_temp(time) + theta(2+order_ar+j)*epsilon(time-j)
       END DO
    END DO

    !now the AR part
    y_sim_temp(1:order_ar) = y_start + y_sim_temp(1:order_ar)
    start = order_ar+1

    IF (start .GT. nb_discard + 1) THEN
       WRITE(*,*) "Fatal error in ARMA_SIMULATION."
       WRITE(*,*) "The starting point of the simulation must be at least nb_discard+1."
       STOP
    END IF
    DO i=start, T+nb_discard
       y_sim_temp(i) = theta(1) + y_sim_temp(i)
       DO j=1, order_ar
          y_sim_temp(i) = y_sim_temp(i) + theta(1+j)*y_sim_temp(i-j)
       END DO
    END DO

    DO time=1,T
       simulated_series(time) = y_sim_temp(nb_discard+time)
    END DO

    RETURN
  END SUBROUTINE ARMA_SIMULATION

  SUBROUTINE AR_AS_MA(n, beta, T, ts_tmp, nb_discard, iid)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, T, nb_discard
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: beta
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:T+nb_discard) :: ts_tmp
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T+nb_discard,1) :: iid
    !beta must contains:
    !beta(1) = the constant of the ar process
    !beta(2) ... beta(n-1) = the lag coefficients
    !beta(n) = the stdv of the innovation errors

    INTEGER :: i, j, p     !p is the order of the ar-process
    REAL (KIND=realkind) :: c, mean
    REAL (KIND=realkind), DIMENSION(1:n-2) :: last_p_coef
    REAL (KIND=realkind), DIMENSION(1:T+nb_discard) :: sigma_vector

    mean = 0._realkind

    DO i=2, n-1
       mean = mean + beta(i)
    END DO

    mean = beta(1)/(1._realkind - mean)

    sigma_vector(1:T+nb_discard) = beta(n)*iid(1:T+nb_discard,1)    

    ts_tmp(1:T+nb_discard) = sigma_vector(1:T+nb_discard)

    p = n-2

    last_p_coef = 0
    last_p_coef(1) = 1._realkind

    DO i=2, T+nb_discard-4
       c=0._realkind

       DO j=1, p
          c = c + beta(1+j)*last_p_coef(j)
       END DO

       DO j=i, T+nb_discard
          ts_tmp(j) = ts_tmp(j) + c*sigma_vector(j-i+1)
       END DO

       IF (c .LT. 1.E-8_realkind) THEN
          EXIT
       ELSE
          IF (p .GE. 2) THEN
             last_p_coef(2:p) = last_p_coef(1:p-1)
          END IF
          last_p_coef(1) = c
       END IF

    END DO

    IF (c .GE. 1.E-8_realkind) THEN

       DO i=T+nb_discard-3, T+nb_discard
          c=0._realkind
          DO j=1, p
             c = c + beta(1+j)*last_p_coef(j)
          END DO

          DO j=i, T+nb_discard
             ts_tmp(j) = ts_tmp(j) + c*sigma_vector(j-i+1)
          END DO

          IF (p .GE. 2) THEN
             last_p_coef(2:p) = last_p_coef(1:p-1)
          END IF
          last_p_coef(1) = c
       END DO

    END IF

    ts_tmp = ts_tmp + mean    

    RETURN
  END SUBROUTINE AR_AS_MA


  SUBROUTINE SWITCH_AR(n,m,nb_par,beta,T,simulated_series,nb_discard, nb_iid_series, iid,states_series)
    IMPLICIT NONE
    !n: the number of states of the Markov chain
    !m: the order of the autoregressiv model
    !nb_par: the total number of parameters in the model, i.e. n*n-n+(m+2)*n
    !beta: the vector of parameters beta=(transition probab.|autoregressive
    !      parameters state 1|autoregressive parameters state 2 ...)
    !T: the desired lenght of the time series
    !simulated_series: the vector to be returned with the simulated series
    !nb_discard: the number of observation to be discarded
    !nb_iid_series: the number of iid series used to simulated the process (2)
    !iid: the vector containing the simulated iid observations. The first
    !column contains the innovation for the construction of the y_t process
    !while the second column contains the uniform random variables used to
    !construct the sequence of regimes, i.e. unobserved states variables.
    
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: m
    INTEGER, INTENT(IN) :: nb_par
    INTEGER, INTENT(IN) :: T
    INTEGER, INTENT(IN) :: nb_discard
    INTEGER, INTENT(IN) :: nb_iid_series
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:nb_par) :: beta
    REAL (KIND=realkind), INTENT(OUT), DIMENSION(1:T) :: simulated_series
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T+nb_discard,1:nb_iid_series) :: iid
    REAL (KIND=realkind), INTENT(OUT), DIMENSION(1:T,1), OPTIONAL :: states_series
    
    INTEGER :: i,j,s
    INTEGER :: nb_cond_par
    INTEGER :: nb_tran_par
    INTEGER, DIMENSION(T+nb_discard) :: states
    REAL (KIND=realkind) :: q
    REAL (KIND=realkind), DIMENSION(1:n) :: ones
    REAL (KIND=realkind), DIMENSION(1:(n+1),1) :: base
    REAL (KIND=realkind), DIMENSION(1:n,1) :: pi_0 
    REAL (KIND=realkind), DIMENSION(1:n,1:n) :: C
    REAL (KIND=realkind), DIMENSION(1:n,1:n) :: Id
    REAL (KIND=realkind), DIMENSION(1:n,1:n) :: P
    REAL (KIND=realkind), DIMENSION(1:(n+1),1:n) :: A
    REAL (KIND=realkind), DIMENSION(1:n,1) :: cum_pi_0
    REAL (KIND=realkind), DIMENSION(1:n,1:n) :: cum_P
    REAL (KIND=realkind), DIMENSION(1:(m+2),n) :: beta_ar
    REAL (KIND=realkind), DIMENSION(1:(n+1)) :: mean
    REAL (KIND=realkind), DIMENSION(1:(T+nb_discard)) :: y_sim
    
    ones = 1._realkind
    base = 0._realkind
    base(n+1,1) = 1._realkind
    nb_tran_par=n*n-n
    DO j=1, n
      q=0._realkind
      DO i=1, n-1
         q=q+beta(i+(j-1)*(n-1))
         P(i,j)=beta(i+(j-1)*(n-1))
      END DO
      P(n,j) = 1._realkind-q
    END DO

    nb_cond_par=(m+2)*n
    DO i=1,n
       beta_ar(1:(m+2),i)=beta((nb_tran_par+(m+2)*(i-1)+1):(nb_tran_par+(m+2)*i))
    END DO
    beta_ar(m+2,:) = sqrt(beta_ar(m+2,:))
    
    Id=EYE(n)
    A(1:n,1:n) = Id(1:n,1:n) - P(1:n,1:n)
    A(n+1,1:n) = ones(1:n)
    !compute pi_0=(A'A)^(-1)A'e
    C = INV(MM(transpose(A),A,n,n+1,n),n)
    pi_0 = MM(MM(C,transpose(A),n,n,n+1),base,n,n+1,1)
    !compute the cumulated probabilities used for the simulation
    CALL CUMULATE_PROBABILITIES(n,1,pi_0,cum_pi_0)
    CALL CUMULATE_PROBABILITIES(n,n,P,cum_P)
    
    !simulate the states of the markov process 
    !first from the ergodic probability pi_0
    DO j=1,n
        IF (iid(1,2) .LE. cum_pi_0(j,1)) THEN
           states(1) = j
           EXIT
        END IF
    END DO
    
    !simulate the states
    DO s=2,T+nb_discard
        i=states(s-1)
        DO j=1,n
           IF (IID(s,2) .LE. cum_P(j,i)) THEN
                 states(s) = j
              EXIT
           END IF
        END DO
    END DO
    IF (PRESENT(states_series)) THEN
       DO i=1,T
          states_series(i,1) = states(i+nb_discard)
       END DO
    END IF
    
    !computation of the conditional mean of every state
    mean=0._realkind
    DO j=1,n
       DO i=1,m
          mean(j) = mean(j) + beta_ar(1+i,j)
       END DO
       mean(j) = beta_ar(1,j)/(1._realkind-mean(j))
       mean(n+1) = mean(n+1) + mean(j)*pi_0(j,1)
    END DO

    y_sim = 0
    y_sim(1:m) = mean(n+1)
    s=m+2
    DO i=m+1,T+nb_discard
       y_sim(i) = beta_ar(1,states(i))
       DO j=1,m
          y_sim(i) = y_sim(i) + beta_ar(1+j,states(i))*y_sim(i-j)
       END DO
       y_sim(i) = y_sim(i) + beta_ar(s,states(i))*iid(i,1)
    END DO
    DO i=1, T
       simulated_series(i) = y_sim(i+nb_discard)
    END DO
    
    RETURN

  CONTAINS

     SUBROUTINE CUMULATE_PROBABILITIES(n,m,P,cum_P)
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: n  !the number of rows (states) in the transition matrix P
     INTEGER, INTENT(IN) :: m !the number of columns (states) in the transition matrix P
     REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n,1:m) :: P
     REAL (KIND=realkind), INTENT(OUT), DIMENSION(1:n,1:m) :: cum_P

     INTEGER :: i, j

     DO j=1,m
        cum_P(1,j) = P(1,j)
        DO i=2,n-1
           cum_P(i,j) = cum_P(i-1,j) + P(i,j)
        END DO
        cum_P(n,j) = 1._realkind
     END DO
     
     RETURN
     END SUBROUTINE CUMULATE_PROBABILITIES     
  END SUBROUTINE SWITCH_AR




  SUBROUTINE SVM1_SIMULATION(n, beta, n1, n2, T, simulated_series, nb_discard, nb_iid_series, iid, simulation_kind, &
       simulated_not_observable)
    IMPLICIT NONE

    !This subroutine simulates an univariate lognormal stochastic autoregressive
    !volatility model (see e.g. T. G. Anderson et al., Journal of Econometrics 91, 
    !(1999) (p. 61-87)
    !MEAN: constant, lags(1) ... lags(n1)
    !VAR.: constant, lags(1) ... lags(n2), sigma   <sigma is the stdv of the innovation error>
    !n = |1 + n1| + |1 + n2 + 1| = 3 + n1 + n2

    INTEGER, INTENT(IN) :: n
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: beta
    INTEGER, INTENT(IN) :: n1
    INTEGER, INTENT(IN) :: n2
    INTEGER, INTENT(IN) :: T
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:T) :: simulated_series
    INTEGER, INTENT(IN) :: nb_discard
    INTEGER, INTENT(IN) :: nb_iid_series
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T+nb_discard,nb_iid_series) :: iid
    CHARACTER (LEN=2), INTENT(IN), OPTIONAL :: simulation_kind
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:T,1), OPTIONAL :: simulated_not_observable

    REAL (KIND=realkind), DIMENSION(1:T+nb_discard) :: y_sim, sigma_sim
    REAL (KIND=realkind), DIMENSION(1:n) :: beta_modif
    REAL (KIND=realkind) :: y_start, sigma_start
    INTEGER :: i, j, start
    LOGICAL :: is_MA

    IF (nb_iid_series .LT. 2) THEN
       WRITE(*,*) "The number of iid series must be at least 2. Fatal error in SVM1_SIMULATION."
       STOP
    END IF
    IF (n2 .LE. 0) THEN
       WRITE(*,*) "This is a degenerate Stoc. vol. model. Plese select another model. Program stopped"
       STOP
    END IF
    IF (beta(n) .LE. 0) THEN
       WRITE(*,*) "The variance of the innovation of the volatility equation must be .GT. 0!"
       WRITE(*,*) "Fatal error in SVM1_SIMULATION."
       STOP
    END IF

    !because we simulate the stdv and not the variance we have to divite beta(n)/2 and
    !beta(2+n1)/2
    beta_modif = beta
    beta_modif(2+n1) = beta(2+n1)/2._realkind
    beta_modif(n) = beta(n)/2._realkind

    !first simulate the volatility model
    IF (PRESENT(simulation_kind)) THEN
       IF (simulation_kind .EQ. "MA") THEN
          is_MA = .TRUE.
       ELSE
          is_MA = .FALSE.
       END IF
    ELSE
       is_MA = .FALSE.
    END IF


    IF (is_MA) THEN

       CALL AR_AS_MA(2+n2, beta_modif(2+n1:n), T, sigma_sim, nb_discard, iid(:,2))

    ELSE !simulate using the autoregressiv scheme
       !compute the mean of the volatility of the process
       sigma_start = 0._realkind
       DO i=1, n2
          sigma_start = sigma_start + beta_modif(2+n1+i)
       END DO

       sigma_start = beta_modif(2+n1)/(1._realkind - sigma_start)

       !compute the starting period for the simulation of the volatility
       start = n2+1
       IF (start .GE. nb_discard + 1) THEN
          WRITE(*,*) "Fatal error in SVM1_SIMULATION."
          WRITE(*,*) "The starting point of the simulation must be at least nb_discard+1."
          STOP
       END IF

       sigma_sim(1:(start-1)) = sigma_start

       sigma_sim(start:nb_discard+T) = beta_modif(2+n1)

       DO i=start,nb_discard+T
          DO j=1, n2
             sigma_sim(i) = sigma_sim(i) + beta_modif(2+n1+j)*sigma_sim(i-j)
          END DO
          sigma_sim(i) = sigma_sim(i) + beta_modif(n)*iid(i,2)
       END DO


    END IF

    sigma_sim = exp(sigma_sim)

    !compute the mean of the process
    y_start = 0._realkind
    DO i=1, n1
       y_start = y_start + beta_modif(1+i)
    END DO
    y_start = beta_modif(1)/(1._realkind - y_start)


    !compute the starting period for the simulation
    start = n1+1

    IF (start .GT. nb_discard + 1) THEN
       WRITE(*,*) "Fatal error in SVM1_SIMULATION."
       WRITE(*,*) "The starting point of the simulation must be at least nb_discard+1."
       STOP
    END IF

    y_sim(1:(start-1)) = y_start
    y_sim(start:nb_discard+T)=beta_modif(1)

    DO i=start,nb_discard+T
       DO j=1, n1
          y_sim(i) = y_sim(i) + beta_modif(1+j)*y_sim(i-j)
       END DO
       y_sim(i) = y_sim(i) + sigma_sim(i)*iid(i,1)
    END DO

    DO i=1,T
       simulated_series(i) = y_sim(nb_discard+i)
    END DO

    IF (PRESENT(simulated_not_observable)) THEN
       DO i=1,T
          simulated_not_observable(i,1) = sigma_sim(nb_discard+i)
       END DO
    END IF

    RETURN
  END SUBROUTINE SVM1_SIMULATION


  SUBROUTINE AR_SIMULATION(n,beta,T,simulated_series,nb_discard,nb_iid_series, iid)
    USE DATATYPES 
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: beta
    INTEGER, INTENT(IN) :: T
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:T) :: simulated_series
    INTEGER, INTENT(IN) :: nb_discard
    INTEGER, INTENT(IN) :: nb_iid_series
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T+nb_discard,nb_iid_series) :: iid

    !beta(1) contains the constant of the process,
    !beta(2) ... beta(n-1) the AR coeff.,
    !beta(n) the sigma2 coeff.

    REAL (KIND=realkind) :: stdv, y_start, c
    REAL (KIND=realkind), DIMENSION(1:T+nb_discard) :: y_sim_temp
    INTEGER :: i, j, start

    !simulate an AR(2)
    y_sim_temp = 0._realkind
    c = 0._realkind
    DO i=2, n-1
       c = c + beta(i)
    END DO

    y_start = beta(1)/(1._realkind-c)

    stdv = SQRT(beta(n))

    y_sim_temp(1:n-2) = y_start
    start = n-1
    IF (start .GT. nb_discard + 1) THEN
       WRITE(*,*) "Fatal error in AR_SIMULATION."
       WRITE(*,*) "The starting point of the simulation must be at least nb_discard+1."
       STOP
    END IF
    DO i=start, T+nb_discard
       y_sim_temp(i) = beta(1) + stdv*iid(i,1)
       DO j=1, n-2
          y_sim_temp(i) = y_sim_temp(i) + beta(1+j)*y_sim_temp(i-j)
       END DO
    END DO
    DO i=1,T
       simulated_series(i) = y_sim_temp(nb_discard+i)
    END DO

    RETURN
  END SUBROUTINE AR_SIMULATION


  SUBROUTINE AR_GARCH(n,beta_true,dim_mu,dim_arch,dim_garch,T,simulated_series, &
       nb_discard,nb_iid_series,iid)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: beta_true
    INTEGER, INTENT(IN) :: dim_mu
    INTEGER, INTENT(IN) :: dim_arch
    INTEGER, INTENT(IN) :: dim_garch
    INTEGER, INTENT(IN) :: T
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:T) :: simulated_series
    INTEGER, INTENT(IN) :: nb_discard
    INTEGER, INTENT(IN) :: nb_iid_series
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T+nb_discard,nb_iid_series) :: iid

    INTEGER :: dim_sigma, i, j, start
    REAL (KIND=realkind), DIMENSION(1:T+nb_discard) :: ts_temp
    REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: sigma2
    REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: U2
    REAL (KIND=realkind) :: mean, variance, root_sigma2  

    dim_sigma = 1 + dim_arch + dim_garch

    !now start simulate the process
    !Computation of the mean of the process 
    mean = 1._realkind 
    DO i=2, dim_mu 
       mean = mean - beta_true(i)  
    END DO
    mean = beta_true(1)/mean 
    root_sigma2 = SQRT(beta_true(dim_mu+1)) 

    IF (dim_sigma .EQ. 1) THEN 
       !First the case when there is no ARCH/GARCH effect 
       ts_temp(1:dim_mu) = mean 
       DO i=dim_mu, T+nb_discard 
          ts_temp(i) = beta_true(1) 
          DO j=2, dim_mu 
             ts_temp(i) = ts_temp(i) + beta_true(j)*ts_temp(i-j+1) 
          END DO
          ts_temp(i) = ts_temp(i) + root_sigma2*iid(i,1)            
       END DO
       simulated_series(1:T) = ts_temp(nb_discard+1:T+nb_discard) 
    ELSE 
       !Then the case of an ARCH process
       ALLOCATE (sigma2(1:T+nb_discard))
       ALLOCATE (U2(1:T+nb_discard))  
       IF (dim_garch .EQ. 0) THEN  
          !Computation of the variance of the process under the requirement 
          !that the cond. exp. of the innovation is equal 0 
          variance = 1._realkind 
          DO i=2, dim_sigma 
             variance = variance - beta_true(dim_mu+i)  
          END DO
          variance = beta_true(dim_mu+1)/variance 
          !Now start the construction of the time series. Replace the missing lags  
          !with the unconditional moments 
          start = MAX(dim_mu,dim_sigma) 
          ts_temp(1:start) = mean 
          U2(1:start) = variance 
          DO i=start, T+nb_discard 
             ts_temp(i) = beta_true(1) 
             DO j=2, dim_mu 
                ts_temp(i) = ts_temp(i) + beta_true(j)*ts_temp(i-j+1) 
             END DO
             sigma2(i) = beta_true(dim_mu+1) 
             DO j=2, dim_sigma 
                sigma2(i) = sigma2(i) + beta_true(dim_mu+j)*U2(i-j+1) 
             END DO
             U2(i) = sigma2(i)*iid(i,1)**2 
             ts_temp(i) = ts_temp(i) + SQRT(sigma2(i))*iid(i,1) 
          END DO
          simulated_series(1:T) = ts_temp(nb_discard+1:T+nb_discard) 
       ELSE 
          !The last case is GARCH process 
          !Computation of the variance of the process under the requirement 
          !that the cond. exp. of the innovation is equal 0 
          variance = 1._realkind 
          DO i=2, dim_sigma 
             variance = variance - beta_true(dim_mu+i)  
          END DO
          variance = beta_true(dim_mu+1)/variance 
          start = MAX(dim_mu,dim_arch+1) 
          start = MAX(start,dim_garch+1) 
          ts_temp(1:start) = mean 
          U2(1:start) = variance 
          sigma2(1:start) = variance 
          DO i=start, T+nb_discard 
             ts_temp(i) = beta_true(1) 
             DO j=2, dim_mu 
                ts_temp(i) = ts_temp(i) + beta_true(j)*ts_temp(i-j+1) 
             END DO
             sigma2(i) = beta_true(dim_mu+1) 
             DO j=1, dim_arch 
                sigma2(i) = sigma2(i) + beta_true(dim_mu+1+j)*U2(i-j) 
             END DO
             DO j=1, dim_garch 
                sigma2(i) = sigma2(i) + beta_true(dim_mu+dim_arch+1+j)*sigma2(i-j) 
             END DO
             U2(i) = sigma2(i)*iid(i,1)**2 
             ts_temp(i) = ts_temp(i) + SQRT(sigma2(i))*iid(i,1) 
          END DO
          simulated_series(1:T) = ts_temp(nb_discard+1:T+nb_discard) 
       END IF

       IF (ALLOCATED(U2)) DEALLOCATE(U2)
       IF (ALLOCATED(sigma2)) DEALLOCATE(sigma2)

    END IF


    RETURN
  END SUBROUTINE AR_GARCH


  SUBROUTINE ARMA_GARCH(n,beta_true,order_ar,order_ma,dim_arch,dim_garch,T,simulated_series, &
       nb_discard,nb_iid_series,iid)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: beta_true
    INTEGER, INTENT(IN) :: order_ar
    INTEGER, INTENT(IN) :: order_ma
    INTEGER, INTENT(IN) :: dim_arch
    INTEGER, INTENT(IN) :: dim_garch
    INTEGER, INTENT(IN) :: T
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:T) :: simulated_series
    INTEGER, INTENT(IN) :: nb_discard
    INTEGER, INTENT(IN) :: nb_iid_series
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T+nb_discard,nb_iid_series) :: iid

    INTEGER :: dim_mu, dim_sigma, i, j, start
    REAL (KIND=realkind), DIMENSION(1:T+nb_discard) :: ts_temp
    REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: sigma2
    REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: U2
    REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: eps
    REAL (KIND=realkind) :: mean, variance, root_sigma2  

    dim_mu    = 1 + order_ar + order_ma
    dim_sigma = 1 + dim_arch + dim_garch

    !now start simulate the process
    !Computation of the mean of the process

    mean = 1._realkind 
    DO i=2, order_ar+1
       mean = mean - beta_true(i)  
    END DO
    mean = beta_true(1)/mean 
    ALLOCATE (eps(1:T+nb_discard))

    IF (dim_sigma .EQ. 1) THEN 
       !First the case when there is no ARCH/GARCH effect 
       root_sigma2 = SQRT(beta_true(dim_mu+1))
       eps(1:T+nb_discard) = root_sigma2*iid(1:T+nb_discard,1)  
       start = MAX(order_ar+1,order_ma+1)
       ts_temp(1:start) = mean
       DO i=start, T+nb_discard
          ts_temp(i) = beta_true(1) 
          DO j=2, order_ar+1  
             ts_temp(i) = ts_temp(i) + beta_true(j)*ts_temp(i-j+1) 
          END DO
          ts_temp(i) = ts_temp(i) + eps(i)
          DO j=1, order_ma
             ts_temp(i) = ts_temp(i) + beta_true(1+order_ar+j)*eps(i-j)
          END DO
       END DO
       simulated_series(1:T) = ts_temp(nb_discard+1:T+nb_discard)
    ELSE 
       !Then the case of an ARCH process
       ALLOCATE (sigma2(1:T+nb_discard))
       ALLOCATE (U2(1:T+nb_discard))  
       IF (dim_garch .EQ. 0) THEN  
          !Computation of the variance of the process under the requirement 
          !that the cond. exp. of the innovation is equal 0 
          variance = 1._realkind 
          DO i=2, dim_sigma 
             variance = variance - beta_true(dim_mu+i)  
          END DO
          variance = beta_true(dim_mu+1)/variance 
          !Now start the construction of the time series. Replace the missing lags  
          !with the unconditional moments
          start = MAX(order_ma+1,dim_arch+1)
          start = MAX(order_ar+1,start) 
          ts_temp(1:start) = mean
          U2(1:start) = variance
          eps(1:start) = SQRT(variance)
          !first simulate the innovations 
          DO i=start, T+nb_discard 
             sigma2(i) = beta_true(dim_mu+1) 
             DO j=2, dim_sigma 
                sigma2(i) = sigma2(i) + beta_true(dim_mu+j)*U2(i-j+1) 
             END DO
             eps(i) = SQRT(sigma2(i))*iid(i,1)
             U2(i) = eps(i)**2 
          END DO
       ELSE 
          !The last case is GARCH process 
          !Computation of the variance of the process under the requirement 
          !that the cond. exp. of the innovation is equal 0 
          variance = 1._realkind 
          DO i=2, dim_sigma 
             variance = variance - beta_true(dim_mu+i)  
          END DO
          variance = beta_true(dim_mu+1)/variance 
          start = MAX(order_ma+1,dim_arch+1)
          start = MAX(start,dim_garch+1)
          start = MAX(order_ar+1,start)
          ts_temp(1:start) = mean 
          U2(1:start) = variance 
          sigma2(1:start) = variance 
          eps(1:start) = SQRT(variance)
          !first simulate the variance
          DO i=start, T+nb_discard 
             sigma2(i) = beta_true(dim_mu+1) 
             DO j=1, dim_arch 
                sigma2(i) = sigma2(i) + beta_true(dim_mu+1+j)*U2(i-j) 
             END DO
             DO j=1, dim_garch 
                sigma2(i) = sigma2(i) + beta_true(dim_mu+dim_arch+1+j)*sigma2(i-j) 
             END DO
             eps(i) = SQRT(sigma2(i))*iid(i,1)
             U2(i) = eps(i)**2
          END DO
       END IF

       !then simulate the observed process
       DO i=start, T+nb_discard 
          ts_temp(i) = beta_true(1) 
          DO j=2, order_ar+1
             ts_temp(i) = ts_temp(i) + beta_true(j)*ts_temp(i-j+1) 
          END DO
          ts_temp(i) = ts_temp(i) + eps(i)
          DO j=1, order_ma
             ts_temp(i) = ts_temp(i) + beta_true(1+order_ar+j)*eps(i-j)
          END DO
       END DO
       simulated_series(1:T) = ts_temp(nb_discard+1:T+nb_discard)

       IF (ALLOCATED(U2)) DEALLOCATE(U2)
       IF (ALLOCATED(sigma2)) DEALLOCATE(sigma2)

    END IF
    DEALLOCATE (eps)

    RETURN
  END SUBROUTINE ARMA_GARCH


  SUBROUTINE AR_GARCH_SNP_LAG(n, beta_true, dim_mu, dim_arch, dim_garch, &
       pol_dim, pol_deg, nb_pol_coef, pol_coef, epsilon_0, &
       T, simulated_series, nb_discard, &
       nb_iid_series, random_unif, consider_contamination)
    USE CONTAMINATION_MODEL
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: beta_true
    !common variables defining the conditional mean mu of the markov process 
    INTEGER, INTENT(IN) :: dim_mu                 !number of lags + 1 (for the constant) 
    INTEGER, INTENT(IN) :: dim_arch
    INTEGER, INTENT(IN) :: dim_garch
    !common variables defining the joint density function
    INTEGER, INTENT(IN) :: pol_dim                !Dimension of the polynomial
    INTEGER, INTENT(IN), DIMENSION(1:pol_dim) :: pol_deg !=(Maximal exponent in the i-th dimension) 
    INTEGER, INTENT(IN) :: nb_pol_coef
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:nb_pol_coef) :: pol_coef
    REAL (KIND=realkind), INTENT(IN) :: epsilon_0
    INTEGER, INTENT(IN) :: T
    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:T) :: simulated_series
    INTEGER, INTENT(IN) :: nb_discard
    INTEGER, INTENT(IN) :: nb_iid_series
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:T+nb_discard,nb_iid_series) :: random_unif
    INTEGER, INTENT(IN) :: consider_contamination    !0=no consider, 1=consider


    !common variables used to compute the conditional density function 
    INTEGER :: deg_lag_max !Maximal exponent between all exponents from the 2-nd to the pol_dim-th dimension 
    REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: X 
    INTEGER, ALLOCATABLE, DIMENSION(:) :: nb_rep, nb_jumps, high_jumps 
    INTEGER :: dim_sigma, dim_rows_X

    !other variables used
    REAL (KIND=realkind), DIMENSION(1:T+nb_discard)    :: ts_temp
    REAL (KIND=realkind), DIMENSION(1:nb_pol_coef)     :: theta
    REAL (KIND=realkind), DIMENSION(1:pol_deg(1)+1)    :: theta1
    REAL (KIND=realkind), DIMENSION(1:T+nb_discard,1)  :: iid
    REAL (KIND=realkind), DIMENSION(1,1)               :: iid_tmp_bis
    REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:)    :: U2
    REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:)    :: sigma2

    REAL (KIND=realkind) :: mean, variance, root_sigma2  
    REAL (KIND=realkind), DIMENSION(1:nb_contam_points(1)) :: additive_contam_points
    REAL (KIND=realkind), DIMENSION(1:nb_contam_points(1)) :: zero_additive_contam_points

    INTEGER :: i, j, r, s, start, count, nb_variable, deg_variable

    !here starts the procedure
    zero_additive_contam_points=0._realkind
    dim_sigma = 1 + dim_arch + dim_garch

    !now start simulate the process
    !Computation of the mean of the process 
    mean = 1._realkind 
    DO i=2, dim_mu 
       mean = mean - beta_true(i)  
    END DO
    mean = beta_true(1)/mean 
    root_sigma2 = SQRT(beta_true(dim_mu+1)) 

    !Check the dimension of the polynomial 
    !pol_dim = 0 means the polynomial doesn't exist, "it's equal 0" 
    IF (pol_dim .LE. 1) THEN 
       WRITE(*,*) "Error in subroutine AR_GARCH_SNP_LAG. Program stopped."
       WRITE(*,*) "The polynomial must have a dimension > 1" 
       STOP 
    END IF

    IF (pol_deg(1) .LE. 0) THEN 
       WRITE(*,*) "Error in subroutine AR_GARCH_SNP_LAG. Program stopped."
       WRITE(*,*) "The degree of the first dimension must be > 0." 
       STOP 
    END IF

    !Check if the degrees of the lags are correct 
    DO i=2, pol_dim 
       IF (pol_deg(i) .LT. 0) THEN 
          WRITE(*,*) "Error in subroutine AR_GARCH_SNP_LAG. Program stopped."
          WRITE(*,*) "Error by definion of polynomial's degrees: negative value found." 
          STOP 
       END IF
    END DO

    deg_lag_max = MAXVAL(pol_deg(2:pol_dim)) 
    IF (deg_lag_max .EQ. 0) THEN 
       WRITE(*,*) "Error in subroutine AR_GARCH_SNP_LAG. Program stopped."
       WRITE(*,*) "At least one lag dimension must have a degree > 0." 
       STOP 
    END IF
    dim_rows_X = 0 
    DO i=2, pol_dim 
       dim_rows_X = dim_rows_X+pol_deg(i) 
    END DO

    ALLOCATE (X(1:dim_rows_X,1)) !Nella matrice X la riga di costanti
    !corrispondente alla potenza 0 e' omesso.
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

    IF (dim_sigma .EQ. 1) THEN
       !First the case when there is no ARCH/GARCH effect
       start = MAX(dim_mu,pol_dim)
       IF (start .GT. nb_discard + 1) THEN
          WRITE(*,*) "Fatal error in AR_GARCH_SNP_LAG."
          WRITE(*,*) "The starting point of the simulation must be at least nb_discard+1."
          STOP
       END IF
       ts_temp(1:start) = mean
       DO s=start, T+nb_discard
          !Now construct the X vector 
          count = 1 
          DO j = 1, pol_dim - 1 
             IF (pol_deg(j+1) .GT. 0) THEN 
                X(count,1) = ts_temp(s-j) 
                count = count + 1 
                DO i=2, pol_deg(j+1) 
                   X(count,1) = X(count-1,1)*ts_temp(s-j) 
                   count = count + 1 
                END DO
             END IF
          END DO
          !Now compute the theta1 vector to be used by INV_SNP_DIST_FUNC 
          theta = pol_coef 
          theta1 = 0.0_realkind 
          i = 1 
          DO nb_variable=2, pol_dim 
             DO deg_variable=1, pol_deg(nb_variable) 
                DO j=1, nb_jumps(nb_variable) 
                   DO r=1, nb_rep(nb_variable) 
                      theta(nb_rep(nb_variable)*deg_variable + (j-1)*high_jumps(nb_variable) + r) = & 
                           theta(nb_rep(nb_variable)*deg_variable + (j-1)*high_jumps(nb_variable) + r) * & 
                           X(i,1)                         
                   END DO
                END DO
                i = i+1 
             END DO
          END DO
          DO i=1, pol_deg(1)+1 
             DO j=(i-1)*nb_rep(1)+1, i*nb_rep(1) 
                theta1(i) = theta1(i) + theta(j)     
             END DO
          END DO

          iid(s,1) = INV_SNP_DIST_FUNC(random_unif(s,1),theta1,nb_pol_coef,epsilon_0) 
          IF (consider_contamination .EQ. 1) THEN
             IF (contam_model(2) .EQ. 1) THEN
                IF (contam_dist(2) .EQ. "_DIRAC_") THEN
                   additive_contam_points(:)=contam_x0(2,:)
                   CALL CONTAMINATE(1,iid(s,1),epsilon_p(2), nb_contam_points(2),x0=additive_contam_points)
                END IF
                IF (contam_dist(2) .EQ. "UNIFORM") THEN
                   iid_tmp_bis = IID_SIM(1,1,contam_dist(2))
                   iid_tmp_bis = iid_tmp_bis + (contam_upper_bound(2) -  contam_lower_bound(2)) + &
                        contam_lower_bound(2)
                   CALL CONTAMINATE(1,iid(s,1),epsilon_p(2),nb_contam_points(2),random_x0=iid_tmp_bis)
                END IF
                IF (contam_dist(2) .EQ. "NORMAL_") THEN
                   iid_tmp_bis = IID_SIM(1,1,contam_dist(2))
                   iid_tmp_bis = iid_tmp_bis * contam_norm_var(2) + contam_norm_mean(2)
                   CALL CONTAMINATE(1,iid(s,1),epsilon_p(2),nb_contam_points(2),random_x0=iid_tmp_bis)
                END IF
                IF (contam_dist(2) .EQ. "STUDENT") THEN
                   iid_tmp_bis = IID_SIM(1,1,contam_dist(2),contam_df(2))
                   iid_tmp_bis = iid_tmp_bis * contam_student_var(2) + contam_student_mean(2)
                   CALL CONTAMINATE(1,iid(s,1),epsilon_p(2),nb_contam_points(2),random_x0=iid_tmp_bis)
                END IF
             ELSE ! perform an empty simulation for comparison with the no contamination
                IF (contam_dist(2) .EQ. "UNIFORM") THEN
                   iid_tmp_bis = IID_SIM(1,1,contam_dist(1))
                   !iid_tmp_bis = iid_tmp_bis + (contam_upper_bound(1) -  contam_lower_bound(1)) + &
                   !     contam_lower_bound(1)
                END IF
                IF (contam_dist(2) .EQ. "NORMAL_") THEN
                   iid_tmp_bis = IID_SIM(1,1,contam_dist(1))
                   !iid_tmp_bis = iid_tmp_bis * contam_norm_var(1) + contam_norm_mean(1)
                END IF
                IF (contam_dist(2) .EQ. "STUDENT") THEN
                   iid_tmp_bis = IID_SIM(1,1,contam_dist(1),contam_df(1))
                   !iid_tmp_bis = iid_tmp_bis * contam_student_var(1) + contam_student_mean(1)
                END IF
                CALL CONTAMINATE(1,iid(s,1),0._realkind,nb_contam_points(2),x0=zero_additive_contam_points)
             END IF
          END IF
          ts_temp(s) = beta_true(1) 
          DO j=2, dim_mu 
             ts_temp(s) = ts_temp(s) + beta_true(j)*ts_temp(s-j+1) 
          END DO
          ts_temp(s) = ts_temp(s) + root_sigma2*iid(s,1)            
       END DO
       simulated_series(1:T) = ts_temp(nb_discard+1:T+nb_discard)
    ELSE 
       !Then the case of an ARCH process 
       IF (dim_garch .EQ. 0) THEN
          ALLOCATE (sigma2(1:T+nb_discard))
          ALLOCATE (U2(1:T+nb_discard)) 
          !Computation of the variance of the process under the requirement 
          !that the cond. exp. of the innovation is equal 0 
          variance = 1._realkind 
          DO i=2, dim_sigma 
             variance = variance - beta_true(dim_mu+i)  
          END DO
          !The so computed variance is not correct! 
          variance = beta_true(dim_mu+1)/variance 
          !Now start the construction of the time series. Replace the missing lags  
          !with the unconditional moments 
          start = MAX(dim_mu,dim_sigma) 
          start = MAX(start,pol_dim)
          IF (start .GT. nb_discard + 1) THEN
             WRITE(*,*) "Fatal error in AR_GARCH_SNP_LAG."
             WRITE(*,*) "The starting point of the simulation must be at least nb_discard+1."
             STOP
          END IF
          ts_temp(1:start) = mean 
          U2(1:start) = variance 
          DO s=start, T+nb_discard 
             !Now construct the X vector 
             count = 1 
             DO j = 1, pol_dim - 1 
                IF (pol_deg(j+1) .GT. 0) THEN 
                   X(count,1) = ts_temp(s-j) 
                   count = count + 1 
                   DO i=2, pol_deg(j+1) 
                      X(count,1) = X(count-1,1)*ts_temp(s-j) 
                      count = count + 1 
                   END DO
                END IF
             END DO
             !Now compute the theta1 vector to be used by INV_SNP_DIST_FUNC 
             theta = pol_coef 
             theta1 = 0.0_realkind 
             i = 1 
             DO nb_variable=2, pol_dim 
                DO deg_variable=1, pol_deg(nb_variable) 
                   DO j=1, nb_jumps(nb_variable) 
                      DO r=1, nb_rep(nb_variable) 
                         theta(nb_rep(nb_variable)*deg_variable + (j-1)*high_jumps(nb_variable) + r) = & 
                              theta(nb_rep(nb_variable)*deg_variable + (j-1)*high_jumps(nb_variable) + r) * & 
                              X(i,1)                         
                      END DO
                   END DO
                   i = i+1 
                END DO
             END DO
             DO i=1, pol_deg(1)+1 
                DO j=(i-1)*nb_rep(1)+1, i*nb_rep(1) 
                   theta1(i) = theta1(i) + theta(j)     
                END DO
             END DO
             iid(s,1) = INV_SNP_DIST_FUNC(random_unif(s,1),theta1,nb_pol_coef,epsilon_0)
             IF (consider_contamination .EQ. 1) THEN
                IF (contam_model(2) .EQ. 1) THEN
                   IF (contam_dist(2) .EQ. "_DIRAC_") THEN
                      additive_contam_points(:)=contam_x0(2,:)
                      CALL CONTAMINATE(1,iid(s,1),epsilon_p(2),nb_contam_points(2), x0=additive_contam_points)
                   END IF
                   IF (contam_dist(2) .EQ. "UNIFORM") THEN
                      iid_tmp_bis = IID_SIM(1,1,contam_dist(2))
                      iid_tmp_bis = iid_tmp_bis + (contam_upper_bound(2) -  contam_lower_bound(2)) + &
                           contam_lower_bound(2)
                      CALL CONTAMINATE(1,iid(s,1),epsilon_p(2),nb_contam_points(2),random_x0=iid_tmp_bis)
                   END IF
                   IF (contam_dist(2) .EQ. "NORMAL_") THEN
                      iid_tmp_bis = IID_SIM(1,1,contam_dist(2))
                      iid_tmp_bis = iid_tmp_bis * contam_norm_var(2) + contam_norm_mean(2)
                      CALL CONTAMINATE(1,iid(s,1),epsilon_p(2),nb_contam_points(2),random_x0=iid_tmp_bis)
                   END IF
                   IF (contam_dist(2) .EQ. "STUDENT") THEN
                      iid_tmp_bis = IID_SIM(1,1,contam_dist(2),contam_df(2))
                      iid_tmp_bis = iid_tmp_bis * contam_student_var(2) + contam_student_mean(2)
                      CALL CONTAMINATE(1,iid(s,1),epsilon_p(2),nb_contam_points(2),random_x0=iid_tmp_bis)
                   END IF
                ELSE
                   IF (contam_dist(2) .EQ. "UNIFORM") THEN
                      iid_tmp_bis = IID_SIM(1,1,contam_dist(1))
                      !iid_tmp_bis = iid_tmp_bis + (contam_upper_bound(1) -  contam_lower_bound(1)) + &
                      !     contam_lower_bound(1)
                   END IF
                   IF (contam_dist(2) .EQ. "NORMAL_") THEN
                      iid_tmp_bis = IID_SIM(1,1,contam_dist(1))
                      !iid_tmp_bis = iid_tmp_bis * contam_norm_var(1) + contam_norm_mean(1)
                   END IF
                   IF (contam_dist(2) .EQ. "STUDENT") THEN
                      iid_tmp_bis = IID_SIM(1,1,contam_dist(1),contam_df(1))
                      !iid_tmp_bis = iid_tmp_bis * contam_student_var(1) + contam_student_mean(1)
                   END IF
                   CALL CONTAMINATE(1,iid(s,1),0._realkind,nb_contam_points(2),x0=zero_additive_contam_points)
                END IF
             END IF
             ts_temp(s) = beta_true(1) 
             DO j=2, dim_mu 
                ts_temp(s) = ts_temp(s) + beta_true(j)*ts_temp(s-j+1) 
             END DO
             sigma2(s) = beta_true(dim_mu+1) 
             DO j=2, dim_sigma 
                sigma2(s) = sigma2(s) + beta_true(dim_mu+j)*U2(s-j+1) 
             END DO
             U2(s) = sigma2(s)*iid(s,1)**2  
             ts_temp(s) = ts_temp(s) + sqrt(sigma2(s))*iid(s,1) 
          END DO
          simulated_series(1:T) = ts_temp(nb_discard+1:T+nb_discard)
       ELSE 
          !The last case is GARCH process 
          ALLOCATE (sigma2(1:T+nb_discard))
          ALLOCATE (U2(1:T+nb_discard)) 
          !Computation of the variance of the process under the requirement 
          !that the cond. exp. of the innovation is equal 0 
          variance = 1._realkind 
          DO i=2, dim_sigma 
             variance = variance - beta_true(dim_mu+i)  
          END DO
          !The so computed variance is not correct!               
          variance = beta_true(dim_mu+1)/variance 
          start = MAX(dim_mu,dim_arch+1) 
          start = MAX(start,dim_garch+1) 
          start = MAX(start,pol_dim)
          IF (start .GT. nb_discard + 1) THEN
             WRITE(*,*) "Fatal error in AR_GARCH_SNP_LAG."
             WRITE(*,*) "The starting point of the simulation must be at least nb_discard+1."
             STOP
          END IF
          ts_temp(1:start) = mean 
          U2(1:start) = variance 
          sigma2(1:start) = SQRT(variance) 
          DO s=start, T+nb_discard 
             !Now construct the X vector 
             count = 1 
             DO j = 1, pol_dim - 1 
                IF (pol_deg(j+1) .GT. 0) THEN 
                   X(count,1) = ts_temp(s-j)  
                   count = count + 1 
                   DO i=2, pol_deg(j+1)
                      X(count,1) = X(count-1,1)*ts_temp(s-j)  
                      count = count + 1 
                   END DO
                END IF
             END DO
             !Now compute the theta1 vector to be used by INV_SNP_DIST_FUNC 
             theta = pol_coef
             theta1 = 0.0_realkind 
             i = 1 
             DO nb_variable=2, pol_dim 
                DO deg_variable=1, pol_deg(nb_variable) 
                   DO j=1, nb_jumps(nb_variable) 
                      DO r=1, nb_rep(nb_variable) 
                         theta(nb_rep(nb_variable)*deg_variable + (j-1)*high_jumps(nb_variable) + r) = & 
                              theta(nb_rep(nb_variable)*deg_variable + (j-1)*high_jumps(nb_variable) + r) * & 
                              X(i,1)                         
                      END DO
                   END DO
                   i = i+1 
                END DO
             END DO
             DO i=1, pol_deg(1)+1 
                DO j=(i-1)*nb_rep(1)+1, i*nb_rep(1) 
                   theta1(i) = theta1(i) + theta(j)     
                END DO
             END DO
             iid(s,1) = INV_SNP_DIST_FUNC(random_unif(s,1),theta1,nb_pol_coef,epsilon_0) 
             IF (consider_contamination .EQ. 1) THEN
                IF (contam_model(2) .EQ. 1) THEN
                   IF (contam_dist(2) .EQ. "_DIRAC_") THEN
                      additive_contam_points(:)=contam_x0(2,:)
                      CALL CONTAMINATE(1,iid(s,1),epsilon_p(2),nb_contam_points(2),x0=additive_contam_points)
                   END IF
                   IF (contam_dist(2) .EQ. "UNIFORM") THEN
                      iid_tmp_bis = IID_SIM(1,1,contam_dist(2))
                      iid_tmp_bis = iid_tmp_bis + (contam_upper_bound(2) -  contam_lower_bound(2)) + &
                           contam_lower_bound(2)
                      CALL CONTAMINATE(1,iid(s,1),epsilon_p(2),nb_contam_points(2),random_x0=iid_tmp_bis)
                   END IF
                   IF (contam_dist(2) .EQ. "NORMAL_") THEN
                      iid_tmp_bis = IID_SIM(1,1,contam_dist(2))
                      iid_tmp_bis = iid_tmp_bis * contam_norm_var(2) + contam_norm_mean(2)
                      CALL CONTAMINATE(1,iid(s,1),epsilon_p(2),nb_contam_points(2),random_x0=iid_tmp_bis)
                   END IF
                   IF (contam_dist(2) .EQ. "STUDENT") THEN
                      iid_tmp_bis = IID_SIM(1,1,contam_dist(2),contam_df(2))
                      iid_tmp_bis = iid_tmp_bis * contam_student_var(2) + contam_student_mean(2)
                      CALL CONTAMINATE(1,iid(s,1),epsilon_p(2),nb_contam_points(2),random_x0=iid_tmp_bis)
                   END IF
                ELSE
                   IF (contam_dist(2) .EQ. "UNIFORM") THEN
                      iid_tmp_bis = IID_SIM(1,1,contam_dist(1))
                      !iid_tmp_bis = iid_tmp_bis + (contam_upper_bound(1) -  contam_lower_bound(1)) + &
                   !     contam_lower_bound(1)
                   END IF
                   IF (contam_dist(2) .EQ. "NORMAL_") THEN
                      iid_tmp_bis = IID_SIM(1,1,contam_dist(1))
                      !iid_tmp_bis = iid_tmp_bis * contam_norm_var(1) + contam_norm_mean(1)
                   END IF
                   IF (contam_dist(2) .EQ. "STUDENT") THEN
                      iid_tmp_bis = IID_SIM(1,1,contam_dist(1),contam_df(1))
                   !iid_tmp_bis = iid_tmp_bis * contam_student_var(1) + contam_student_mean(1)
                   END IF
                   CALL CONTAMINATE(1,iid(s,1),0._realkind,nb_contam_points(2),x0=zero_additive_contam_points)
                END IF
             END IF
             ts_temp(s) = beta_true(1) 
             DO j=2, dim_mu 
                ts_temp(s) = ts_temp(s) + beta_true(j)*ts_temp(s-j+1) 
             END DO
             sigma2(s) = beta_true(dim_mu+1) 
             DO j=1, dim_arch 
                sigma2(s) = sigma2(s) + beta_true(dim_mu+1+j)*U2(s-j) 
             END DO
             DO j=1, dim_garch 
                sigma2(s) = sigma2(s) + beta_true(dim_mu+dim_arch+1+j)*sigma2(s-j)
             END DO
             U2(s) = sigma2(s)*iid(s,1)**2  
             ts_temp(s) = ts_temp(s) + sqrt(sigma2(s))*iid(s,1) 
          END DO
          simulated_series(1:T) = ts_temp(nb_discard+1:T+nb_discard)
       END IF
    END IF

    IF (ALLOCATED(U2)) DEALLOCATE(U2)
    IF (ALLOCATED(sigma2)) DEALLOCATE(sigma2)
    DEALLOCATE(X)
    DEALLOCATE(nb_rep)
    DEALLOCATE(nb_jumps)
    DEALLOCATE(high_jumps)

    RETURN
  END SUBROUTINE AR_GARCH_SNP_LAG

END MODULE MYFUNCTIONS


MODULE MYINTERFACE
  INTERFACE
     SUBROUTINE MODIFY_SERIES()
       USE DATATYPES
       USE S1_COMMON, ONLY: T1, imported_series, modified_series, empirical_is_mean, &
            empirical_is_variance, empirical_ms_mean, empirical_ms_variance, &
            center, scale 
       IMPLICIT NONE
     END SUBROUTINE MODIFY_SERIES
  END INTERFACE

  INTERFACE
     SUBROUTINE npsol(n, nclin, ncnln, ldA, ldJ, ldR, &
          A, bl, bu, &
          funcon, funobj, &
          inform, iter, istate, &
          c, cJac, clamda, f, g, RR, x, &
          iw, leniw, w, lenw)
       USE datatypes
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: n, nclin, ncnln, ldA, ldJ, ldR, &
            leniw, lenw
       INTEGER, INTENT(OUT) :: inform, iter
       INTEGER, INTENT(IN), DIMENSION(1:leniw) :: iw
       INTEGER, INTENT(IN OUT), DIMENSION(1:n+nclin+ncnln) :: istate
       REAL (KIND=realkind), INTENT(OUT) :: f
       REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n+nclin+ncnln) :: bl, bu
       REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:n+nclin+ncnln) :: clamda
       REAL (KIND=realkind), INTENT(IN), DIMENSION(1:lenw) :: w
       REAL (KIND=realkind), INTENT(IN), DIMENSION(1:ldA,1:n) :: A
       REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:n) :: x
       REAL (KIND=realkind), INTENT(OUT), DIMENSION(1:n) :: g
       REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:ncnln) :: c
       REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:ldR,1:n) :: RR
       REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:ldJ,1:n) :: cJac
       INTERFACE
          SUBROUTINE funobj(mode, n, x, f, g, nstate)
            USE DATATYPES
            IMPLICIT NONE

            INTEGER, INTENT(IN OUT) :: mode
            INTEGER, INTENT(IN) :: n, nstate
            REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: x
            REAL (KIND=realkind), INTENT(OUT) :: f
            REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:n) :: g

          END SUBROUTINE funobj

          SUBROUTINE funcon(mode, ncnln, n, ldJ, needc, x, c, cJac, nstate)
            USE DATATYPES
            IMPLICIT NONE

            INTEGER, INTENT(IN OUT) :: mode
            INTEGER, INTENT(IN) :: ncnln, n, ldJ, nstate
            INTEGER, INTENT(IN), DIMENSION(1:n) :: needc
            REAL (KIND=realkind), INTENT(IN), DIMENSION(1:n) :: x
            REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:ncnln) :: c
            REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(1:ldJ,1:n) :: cJac

          END SUBROUTINE funcon
       END INTERFACE
     END SUBROUTINE npsol
  END INTERFACE
END MODULE MYINTERFACE

MODULE OPTIMUM
  USE DATATYPES

  IMPLICIT NONE

  INTEGER :: NCOM
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:) :: PCOM,XICOM

CONTAINS

  FUNCTION F1DIM(X,MY_FUNC)
    USE DATATYPES

    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN) :: X
    REAL (KIND=realkind) :: F1DIM

    INTEGER :: J
    REAL (KIND=realkind), DIMENSION(NCOM) :: XT

    INTERFACE
       FUNCTION MY_FUNC(beta) RESULT(function_value)
         USE DATATYPES

         IMPLICIT NONE

         REAL (KIND=realkind), INTENT(IN), DIMENSION(:) :: beta
         REAL (KIND=realkind) :: function_value
       END FUNCTION MY_FUNC
    END INTERFACE

    DO J=1,NCOM
       XT(J)=PCOM(J)+X*XICOM(J)
    END DO
    F1DIM=MY_FUNC(XT)

    RETURN
  END FUNCTION F1DIM
END MODULE OPTIMUM


MODULE MY_NPSOL
  USE DATATYPES
  USE S1_COMMON, ONLY: major_iteration_limit
  IMPLICIT NONE

CONTAINS
  SUBROUTINE MINIMUM(FUNC,beta_start,FRET,ITER)
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(:) :: beta_start
    REAL (KIND=realkind), INTENT( OUT) :: FRET
    INTEGER, INTENT(OUT) :: ITER

    INTERFACE

       FUNCTION FUNC(beta) RESULT(function_value)
         USE DATATYPES
         IMPLICIT NONE

         REAL (KIND=realkind), INTENT(IN), DIMENSION(:) :: beta
         REAL (KIND=realkind) :: function_value
       END FUNCTION FUNC

    END INTERFACE

    REAL (KIND=realkind) :: FTOL=0.000001_realkind

    CALL DFPMIN(beta_start,FTOL,ITER,FRET)

    RETURN

  CONTAINS

    SUBROUTINE DFPMIN(P,FTOL,ITER,FRET)
      IMPLICIT NONE

      REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(:) :: P
      REAL (KIND=realkind), INTENT(IN) :: FTOL
      INTEGER, INTENT(OUT) :: ITER
      REAL (KIND=realkind), INTENT(OUT) :: FRET
      !P is the starting point (usually a vector)
      !FTOL is the convergence requirement on the function value (the "epsilon")
      !ITER is the number of iterations that were performed
      !FRET is the minimum value of the function

      INTEGER :: N,I,J,ITS
      INTEGER :: ITMAX   !200
      REAL, PARAMETER :: EPS=0.000001_realkind
      REAL (KIND=realkind), DIMENSION(SIZE(P),SIZE(P)) :: HESSIN
      REAL (KIND=realkind), DIMENSION(SIZE(P)) :: XI, G, DG, HDG
      REAL (KIND=realkind) :: FAC,FAD,FAE,FP

      ITMAX = major_iteration_limit
      N=SIZE(P)
      FP=FUNC(P)
      CALL DFUNC(P,G)
      DO I=1,N
         DO J=1,N
            HESSIN(I,J)=0.0_realkind
         END DO
         HESSIN(I,I)=1.0_realkind
         XI(I)=-G(I)
      END DO
      DO ITS=1,ITMAX
         ITER=ITS
         CALL LINMIN(P,XI,FRET)
         IF(2.0_realkind*ABS(FRET-FP).LE.FTOL*(ABS(FRET)+ABS(FP)+EPS)) RETURN
         FP=FRET
         DO I=1,N
            DG(I)=G(I)
         END DO
         FRET=FUNC(P)
         CALL DFUNC(P,G)
         DO I=1,N
            DG(I)=G(I)-DG(I)
         END DO
         DO I=1,N
            HDG(I)=0.0_realkind
            DO J=1,N
               HDG(I)=HDG(I)+HESSIN(I,J)*DG(J)
            END DO
         END DO
         FAC=0.0_realkind
         FAE=0.0_realkind
         DO I=1,N
            FAC=FAC+DG(I)*XI(I)
            FAE=FAE+DG(I)*HDG(I)
         END DO
         FAC=1./FAC
         FAD=1./FAE
         DO I=1,N
            DG(I)=FAC*XI(I)-FAD*HDG(I)
         END DO
         DO I=1,N
            DO J=1,N
               HESSIN(I,J)=HESSIN(I,J)+FAC*XI(I)*XI(J)-FAD*HDG(I)*HDG(J)+FAE*DG(I)*DG(J)
            END DO
         END DO
         DO I=1,N
            XI(I)=0.0_realkind
            DO J=1,N
               XI(I)=XI(I)-HESSIN(I,J)*G(J)
            END DO
         END DO
      END DO
      WRITE(*,*) "too many iterations in DFPMIN"

      RETURN
    END SUBROUTINE DFPMIN


    SUBROUTINE DFUNC(P,G)
      IMPLICIT NONE

      REAL (KIND=realkind), INTENT(IN), DIMENSION(:) :: P
      REAL (KIND=realkind), INTENT(OUT), DIMENSION(SIZE(P)) :: G

      INTEGER :: J,N
      REAL (KIND=realkind), PARAMETER :: eps=6.0E-6_realkind
      REAL (KIND=realkind), DIMENSION(SIZE(P)) :: ax0,dax0,dh,xdh
      REAL (KIND=realkind), DIMENSION(SIZE(P),2) :: max_matrix
      REAL (KIND=realkind), DIMENSION(SIZE(P),SIZE(P)) :: argplus,argminus

      N=SIZE(P)
      ax0=ABS(P)
      WHERE (ax0 /= 0.0_realkind)
         dax0=P/ax0
      ELSEWHERE
         dax0=1.0_realkind
      END WHERE
      DO J=1,N
         max_matrix(J,2)=1.0E-2_realkind
         max_matrix(J,1)=ax0(J)
      END DO
      dh=eps*MAXVAL(max_matrix,2)*dax0

      xdh=P+dh
      dh=xdh-P       !This increases precision slightly
      DO J=1,N
         argplus(:,J)=P(:)
         argplus(J,J)=xdh(J)
         argminus(:,J)=P(:)
         argminus(J,J)=P(J)-dh(J)
      END DO
      DO J=1,N
         G(J)=(FUNC(argplus(:,J))-FUNC(argminus(:,J)))/(2.0_realkind*dh(J))
      END DO

      RETURN
    END SUBROUTINE DFUNC


    SUBROUTINE LINMIN(P,XI,FRET)
      USE DATATYPES
      USE OPTIMUM

      IMPLICIT NONE

      REAL (KIND=realkind), INTENT(IN OUT), DIMENSION(:) :: P,XI
      REAL (KIND=realkind), INTENT(OUT) :: FRET
      !Given an N dim. point P and an N dim. direction XI, moves and resets
      !P to where the function FUNC(P) takes on a minimum along the direction XI
      !from P, and replaces XI by the actual vector displacement that P was moved.
      !Alse returns as FRET the value of FUNC at the returned location P. This is
      !actually all accomplished by calling the routines MNBRAK and BRENT.

      INTEGER :: N,J
      REAL (KIND=realkind), PARAMETER :: TOL=0.0001_realkind
      REAL (KIND=realkind) :: AX,XX,BX,FA,FB,FX,XMIN

      ALLOCATE (PCOM(SIZE(P)),XICOM(SIZE(P)))
      N=SIZE(P)
      NCOM=N
      DO J=1,N
         PCOM(J)=P(J)
         XICOM(J)=XI(J)
      END DO
      AX=0.0_realkind
      XX=1.0_realkind
      BX=2.0_realkind
      CALL MNBRAK(AX,XX,BX,FA,FX,FB)
      CALL BRENT(FRET,AX,XX,BX,TOL,XMIN)
      DO J=1,N
         XI(J)=XMIN*XI(J)
         P(J)=P(J)+XI(J)
      END DO
      DEALLOCATE(PCOM,XICOM)

      RETURN
    END SUBROUTINE LINMIN


    SUBROUTINE MNBRAK(AX,BX,CX,FA,FB,FC)
      USE OPTIMUM

      IMPLICIT NONE

      REAL (KIND=realkind), INTENT(IN OUT) :: AX,BX
      REAL (KIND=realkind), INTENT(OUT) :: CX,FA,FB,FC

      REAL, PARAMETER :: GOLD=1.618034_realkind,GLIMIT=100.0_realkind,TINY=1.E-20_realkind
      REAL (KIND=realkind) :: FU,Q,R,U,ULIM,DUM

      FA=F1DIM(AX,FUNC)
      FB=F1DIM(BX,FUNC)
      IF(FB.GT.FA)THEN
         DUM=AX
         AX=BX
         BX=DUM
         DUM=FB
         FB=FA
         FA=DUM
      END IF
      CX=BX+GOLD*(BX-AX)
      FC=F1DIM(CX,FUNC)
      DO
         IF (FB .LT. FC) RETURN
         R=(BX-AX)*(FB-FC)
         Q=(BX-CX)*(FB-FA)
         U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.0_realkind*SIGN(MAX(ABS(Q-R),TINY),Q-R))
         ULIM=BX+GLIMIT*(CX-BX)
         IF((BX-U)*(U-CX).GT.0.0_realkind)THEN
            FU=F1DIM(U,FUNC)
            IF(FU.LT.FC)THEN
               AX=BX
               FA=FB
               BX=U
               FB=FU
               RETURN
            ELSE IF(FU.GT.FB)THEN
               CX=U
               FC=FU
               RETURN
            END IF
            U=CX+GOLD*(CX-BX)
            FU=F1DIM(U,FUNC)
         ELSE IF((CX-U)*(U-ULIM).GT.0.0_realkind)THEN
            FU=F1DIM(U,FUNC)
            IF(FU.LT.FC)THEN
               BX=CX
               CX=U
               U=CX+GOLD*(CX-BX)
               FB=FC
               FC=FU
               FU=F1DIM(U,FUNC)
            END IF
         ELSE IF((U-ULIM)*(ULIM-CX).GE.0.0_realkind)THEN
            U=ULIM
            FU=F1DIM(U,FUNC)
         ELSE
            U=CX+GOLD*(CX-BX)
            FU=F1DIM(U,FUNC)
         END IF
         AX=BX
         BX=CX
         CX=U
         FA=FB
         FB=FC
         FC=FU
      END DO

      RETURN
    END SUBROUTINE MNBRAK


    SUBROUTINE BRENT(FRET,AX,BX,CX,TOL,XMIN)
      USE DATATYPES
      USE OPTIMUM

      IMPLICIT NONE

      REAL (KIND=realkind), INTENT(OUT) :: FRET
      REAL (KIND=realkind), INTENT(IN) :: AX,BX,CX,TOL
      REAL (KIND=realkind), INTENT(OUT) :: XMIN

      INTEGER, PARAMETER :: ITMAX=100
      REAL (KIND=realkind), PARAMETER :: CGOLD=.3819660_realkind,ZEPS=1.0E-10_realkind
      INTEGER :: ITER
      REAL (KIND=realkind) :: A,B,D,E,ETEMP,FU,FV,FW,FX,P,Q,R,TOL1,TOL2,U,V,W,X,XM

      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.0_realkind
      FX=F1DIM(X,FUNC)
      FV=FX
      FW=FX
      DO ITER=1,ITMAX
         XM=0.5_realkind*(A+B)
         TOL1=TOL*ABS(X)+ZEPS
         TOL2=2.0_realkind*TOL1
         IF(ABS(X-XM).LE.(TOL2-0.5_realkind*(B-A))) THEN
            XMIN=X
            FRET=FX
            RETURN
         END IF
         IF(ABS(E).GT.TOL1) THEN
            R=(X-W)*(FX-FV)
            Q=(X-V)*(FX-FW)
            P=(X-V)*Q-(X-W)*R
            Q=2.0_realkind*(Q-R)
            IF(Q .GT. 0.0_realkind) P=-P
            Q=ABS(Q)
            ETEMP=E
            E=D
            IF(ABS(P).GE.ABS(0.5_realkind*Q*ETEMP).OR.P.LE.Q*(A-X) .OR. P .GE. Q*(B-X)) THEN
               E=MERGE(A-X,B-X,X >= XM)
               D=CGOLD*E
            ELSE
               D=P/Q
               U=X+D
               IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
            END IF
         ELSE
            E=MERGE(A-X,B-X,X >=XM)
            D=CGOLD*E
         END IF
         U=MERGE(X+D,X+SIGN(TOL1,D), ABS(D) >= TOL1)
         FU=F1DIM(U,FUNC)
         IF(FU.LE.FX) THEN
            IF(U.GE.X) THEN
               A=X
            ELSE
               B=X
            END IF
            V=W
            FV=FW
            W=X
            FW=FX
            X=U
            FX=FU
         ELSE
            IF(U.LT.X) THEN
               A=U
            ELSE
               B=U
            END IF
            IF((FU .LE. FW) .OR. (W .EQ. X)) THEN
               V=W
               FV=FW
               W=U
               FW=FU
            ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
               V=U
               FV=FU
            END IF
         END IF
      END DO
      WRITE(*,*) "Brent exceed maximum iterations."

      RETURN
    END SUBROUTINE BRENT
  END SUBROUTINE MINIMUM
END MODULE MY_NPSOL


MODULE UTILITIES
  USE DATATYPES
  IMPLICIT NONE

CONTAINS

  SUBROUTINE RESTRICT_BETA(beta, beta_restricted)
    USE S1_COMMON, ONLY: AM_nb_par, AM_nb_free_par, is_fixed_par_AM, &
         is_restricted_par_AM
    IMPLICIT NONE
    !this predure removes the restricted components of the AM_parameters from the
    !parameter vector beta
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_par) :: beta
    REAL (KIND=realkind), INTENT(OUT), DIMENSION(1:AM_nb_free_par) :: beta_restricted

    INTEGER :: i,j

    j = 1
    DO i=1, AM_nb_par
       IF ((.NOT. is_fixed_par_AM(i)) .AND. (.NOT. is_restricted_par_AM(i))) THEN
          beta_restricted(j) = beta(i)
          j = j + 1
       END IF
    END DO

    RETURN
  END SUBROUTINE RESTRICT_BETA

  SUBROUTINE COMPLETE_BETA_RESTRICTED(beta_restricted, beta)
    USE S1_COMMON, ONLY: AM_nb_par, AM_nb_free_par, is_fixed_par_AM, &
         is_restricted_par_AM, val_fix_restr_par_AM
    IMPLICIT NONE
    !this predure adds the restricted components of the AM_parameters to the complete
    !vector beta
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_free_par) :: beta_restricted
    REAL (KIND=realkind), INTENT(OUT), DIMENSION(1:AM_nb_par) :: beta


    INTEGER :: i,j

    j = 1
    DO i=1, AM_nb_par
       IF ((.NOT. is_fixed_par_AM(i)) .AND. (.NOT. is_restricted_par_AM(i))) THEN
          beta(i) = beta_restricted(j)
          j = j + 1
       ELSE
          beta(i) = val_fix_restr_par_AM(i)
       END IF
    END DO

    RETURN
  END SUBROUTINE COMPLETE_BETA_RESTRICTED

  SUBROUTINE RESTRICT_RHO(rho, rho_restricted)
    USE S2_COMMON, ONLY: SM_nb_par, SM_nb_free_par, is_fixed_par_SM, &
         is_restricted_par_SM
    IMPLICIT NONE
    !this predure removes the restricted components of the SM_parameters from the
    !parameter vector rho
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:SM_nb_par) :: rho
    REAL (KIND=realkind), INTENT(OUT), DIMENSION(1:SM_nb_free_par) :: rho_restricted

    INTEGER :: i,j

    j = 1
    DO i=1, SM_nb_par
       IF ((.NOT. is_fixed_par_SM(i)) .AND. (.NOT. is_restricted_par_SM(i))) THEN
          rho_restricted(j) = rho(i)
          j = j + 1
       END IF
    END DO

    RETURN
  END SUBROUTINE RESTRICT_RHO

  SUBROUTINE COMPLETE_RHO_RESTRICTED(rho_restricted, rho)
    USE S2_COMMON, ONLY: SM_nb_par, SM_nb_free_par, is_fixed_par_SM, is_restricted_par_SM, &
         val_fix_restr_par_SM
    IMPLICIT NONE
    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:SM_nb_free_par) :: rho_restricted
    REAL (KIND=realkind), INTENT(OUT), DIMENSION(1:SM_nb_par) :: rho
    !this predure adds the restricted components of the SM_parameters to the complete
    !vector rho

    INTEGER :: i,j

    j = 1
    DO i=1, SM_nb_par
       IF ((.NOT. is_fixed_par_SM(i)) .AND. (.NOT. is_restricted_par_SM(i))) THEN
          rho(i) = rho_restricted(j)
          j = j + 1
       ELSE
          rho(i) = val_fix_restr_par_SM(i)
       END IF
    END DO

    RETURN
  END SUBROUTINE COMPLETE_RHO_RESTRICTED


  SUBROUTINE CHECK_CONSTRAINTS_STEP_1(beta,fullfilled)
    USE S1_COMMON, ONLY: bl, bu, nclin,  AM_nb_par, AM_nb_free_par, &
         A_lin_constr_AM, is_lt_lin_constr_AM, is_gt_lin_constr_AM, is_fixed_par_AM, &
         is_restricted_par_AM, is_gt_par_AM, is_lt_par_AM
    IMPLICIT NONE

    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:AM_nb_par) :: beta
    LOGICAL, INTENT(OUT) :: fullfilled

    INTEGER :: i,j,s
    REAL (KIND=realkind) :: a

    fullfilled = .TRUE.
    !check that bl < beta
    j=1
    DO i=1, AM_nb_par
       IF ((.NOT. is_fixed_par_AM(i)) .AND. (.NOT. is_restricted_par_AM(i))) THEN
          IF (is_gt_par_AM(i)) THEN
             IF (bl(j) .GT. beta(i)) THEN
                fullfilled = .FALSE.
                RETURN
             END IF
          END IF
          j=j+1
       END IF
    END DO
    !check that beta < bu
    j=1
    DO i=1, AM_nb_par
       IF ((.NOT. is_fixed_par_AM(i)) .AND. (.NOT. is_restricted_par_AM(i)))  THEN
          IF (is_lt_par_AM(i)) THEN
             IF (bu(j) .LT. beta(i)) THEN
                fullfilled = .FALSE.
                RETURN
             END IF
          END IF
          j=j+1
       END IF
    END DO
    !check that bl < A*beta < bu
    DO i=1, nclin
       a = 0._realkind
       j = 1
       DO s=1,AM_nb_par
          IF ((.NOT. is_fixed_par_AM(s)) .AND. (.NOT. is_restricted_par_AM(s)))  THEN
             a = a + A_lin_constr_AM(i,j)*beta(s)
             j=j+1
          END IF
       END DO
       !check that bl < A*beta
       IF (is_gt_lin_constr_AM(i)) THEN
          IF (bl(AM_nb_free_par+i) .GT. a) THEN
             fullfilled = .FALSE.
             RETURN
          END IF
       END IF
       !check that A*beta < bu
       IF (is_lt_lin_constr_AM(i)) THEN
          IF (bu(AM_nb_free_par+i) .LT. a) THEN
             fullfilled = .FALSE.
             RETURN
          END IF
       END IF
    END DO

    RETURN
  END SUBROUTINE CHECK_CONSTRAINTS_STEP_1


  SUBROUTINE CHECK_CONSTRAINTS_STEP_2(rho,fullfilled)
    USE S1_COMMON, ONLY: bl, bu, nclin
    USE S2_COMMON, ONLY: SM_nb_par, SM_nb_free_par, A_lin_constr_SM, &
         is_lt_lin_constr_SM, is_gt_lin_constr_SM, is_fixed_par_SM, &
         is_gt_par_SM, is_restricted_par_SM, is_lt_par_SM
    IMPLICIT NONE


    REAL (KIND=realkind), INTENT(IN), DIMENSION(1:SM_nb_par) :: rho
    LOGICAL, INTENT(OUT) :: fullfilled

    INTEGER :: i,j,s
    REAL (KIND=realkind) :: a

    fullfilled = .TRUE.
    !check that bl < beta
    j=1
    DO i=1, SM_nb_par
       IF ((.NOT. is_fixed_par_SM(i)) .AND. (.NOT. is_restricted_par_SM(i))) THEN
          IF (is_gt_par_SM(i)) THEN
             IF (bl(j) .GT. rho(i)) THEN
                fullfilled = .FALSE.
                RETURN
             END IF
          END IF
          j=j+1
       END IF
    END DO
    !check that beta < bu
    j=1
    DO i=1, SM_nb_par
       IF ((.NOT. is_fixed_par_SM(i)) .AND. (.NOT. is_restricted_par_SM(i))) THEN
          IF (is_lt_par_SM(i)) THEN
             IF (bu(j) .LT. rho(i)) THEN
                fullfilled = .FALSE.
                RETURN
             END IF
          END IF
          j=j+1
       END IF
    END DO
    !check that bl < A*beta < bu
    DO i=1, nclin
       a = 0._realkind
       j = 1
       DO s=1,SM_nb_par
          IF ((.NOT. is_fixed_par_SM(s)) .AND. (.NOT. is_restricted_par_SM(s))) THEN
             a = a + A_lin_constr_SM(i,j)*rho(s)
             j=j+1
          END IF
       END DO
       !check that bl < A*beta
       IF (is_gt_lin_constr_SM(i)) THEN
          IF (bl(SM_nb_free_par+i) .GT. a) THEN
             fullfilled = .FALSE.
             RETURN
          END IF
       END IF
       !check that A*beta < bu
       IF (is_lt_lin_constr_SM(i)) THEN
          IF (bu(SM_nb_free_par+i) .LT. a) THEN
             fullfilled = .FALSE.
             RETURN
          END IF
       END IF
    END DO

    RETURN
  END SUBROUTINE CHECK_CONSTRAINTS_STEP_2

  SUBROUTINE UPDATE_RANDOM_IID()
    USE S2_COMMON,   ONLY: T2, max_T2, nb_discard_step_2, nb_iid_series, random_iid_step_2
    USE MYFUNCTIONS, ONLY: IID_SIM
    USE ROB_COMMON,  ONLY: random_iid_rob_update
    IMPLICIT NONE

    INTEGER :: i, j
    REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: iid_tmp

    ALLOCATE(iid_tmp(1:(max_T2+nb_discard_step_2),1:nb_iid_series))
    iid_tmp = IID_SIM(max_T2+nb_discard_step_2,nb_iid_series,"NORMAL_")
    DO j=1, nb_iid_series
       DO i=1, T2+nb_discard_step_2
          random_iid_step_2(i,j) = iid_tmp(i,j)
       END DO
    END DO
    iid_tmp = IID_SIM(max_T2+nb_discard_step_2,nb_iid_series,"NORMAL_")
    DO j=1, nb_iid_series
       DO i=1, T2+nb_discard_step_2
          random_iid_rob_update(i,j) = iid_tmp(i,j)
       END DO
    END DO
    DEALLOCATE(iid_tmp)


    RETURN
  END SUBROUTINE UPDATE_RANDOM_IID

END MODULE UTILITIES

