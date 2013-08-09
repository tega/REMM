PROGRAM SNP_PROGRAM
  USE DATATYPES
  USE STEPS,       ONLY: EXPECTATION_MC
  USE STARTUP,     ONLY: INITIATION
  USE S1_COMMON,   ONLY: estimation_kind, write_output_2, write_output_3, tmp_unit_nb, log_unit_nb, dev_unit_nb, &
       space_lenght, expectation_analysis, step, date_start, date_end, time_start, time_end
  USE S2_COMMON,   ONLY: T2, max_T2, nb_iid_series, nb_discard_step_2, random_iid_step_2
  USE ROB_COMMON,  ONLY: random_iid_rob_update
  USE MYFUNCTIONS, ONLY: IID_SIM
  USE DGP_DEFINITION, ONLY: DGP_name

  IMPLICIT NONE

  INTERFACE
     SUBROUTINE SNP()
       IMPLICIT NONE
     END SUBROUTINE SNP

     SUBROUTINE SNP_ROBUST()
       IMPLICIT NONE
     END SUBROUTINE SNP_ROBUST
  END INTERFACE

  INTEGER :: i,j
  REAL (KIND=realkind), ALLOCATABLE, DIMENSION(:,:) :: iid_tmp

  !read the input file "input.dat"
  space_lenght = 3
  CALL INITIATION()

  !open the file which contains all the development output
  IF (write_output_3) THEN
     OPEN (UNIT=dev_unit_nb,FILE="out_dev",FORM="FORMATTED")
  END IF

  !Generate the iid series used for the second step and the robust estimation (even if you
  !don't need it!). max_T2 is the maximal lenght of the simulation (not neccessary = T2)
  IF (write_output_2) THEN
     WRITE(UNIT=log_unit_nb,FMT=*) "START (dd:mm:yy - hh:mm:ss)  ", date_start(7:8), ":", &
          date_start(5:6), ":", date_start(3:4), " - ", &
          time_start(1:2), ":", time_start(3:4), ":", time_start(5:6)
     WRITE(UNIT=log_unit_nb,FMT=*) "SIMULATE IID FOR STEP2 AND ROBUST ESTIMATION"
  END IF

  IF (DGP_name .NE. "SWITCH_AR_MODEL") THEN
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
  ELSE
     ALLOCATE(iid_tmp(1:(max_T2+nb_discard_step_2),1))
     iid_tmp = IID_SIM(max_T2+nb_discard_step_2,1,"NORMAL_")
     DO i=1, T2+nb_discard_step_2
        random_iid_step_2(i,1) = iid_tmp(i,1)
     END DO
     iid_tmp = IID_SIM(max_T2+nb_discard_step_2,1,"UNIFORM")
     DO i=1, T2+nb_discard_step_2
        random_iid_step_2(i,2) = iid_tmp(i,1)
     END DO
     iid_tmp = IID_SIM(max_T2+nb_discard_step_2,1,"NORMAL_")
     DO i=1, T2+nb_discard_step_2
        random_iid_rob_update(i,1) = iid_tmp(i,1)
     END DO
     iid_tmp = IID_SIM(max_T2+nb_discard_step_2,1,"UNIFORM")
     DO i=1, T2+nb_discard_step_2
        random_iid_rob_update(i,2) = iid_tmp(i,1)
     END DO
  END IF
  DEALLOCATE(iid_tmp)

  IF (write_output_2) THEN
     OPEN (UNIT=tmp_unit_nb,FILE="series_iid_for_step_2.dat",FORM="FORMATTED")
     do i=1,T2+nb_discard_step_2
        write(unit=tmp_unit_nb,fmt=*) random_iid_step_2(i,:)
     end do
     close (unit=tmp_unit_nb)
  END IF

  IF (estimation_kind .EQ. 0) THEN
     IF (write_output_2) THEN
        WRITE(UNIT=log_unit_nb,FMT=*) "START SNP()"
     END IF
     CALL SNP()

     IF (write_output_2) THEN
        WRITE(UNIT=log_unit_nb,FMT=*) "EXIT SNP()"
        WRITE(UNIT=log_unit_nb,FMT=*)
     END IF
  ELSE
     IF (write_output_2) THEN
        WRITE(UNIT=log_unit_nb,FMT=*) "START SNP_ROBUST()"
     END IF
     CALL SNP_ROBUST()
     IF (write_output_2) THEN
        WRITE(UNIT=log_unit_nb,FMT=*) "EXIT SNP_ROBUST()"
        WRITE(UNIT=log_unit_nb,FMT=*)
     END IF
  END IF

  IF ((expectation_analysis) .AND. (step .EQ. 3)) THEN
     IF (write_output_2) THEN
        WRITE(UNIT=log_unit_nb,FMT=*) "START THE EXPECTATION ANALYSIS OF STEP_2"
     END IF
     CALL EXPECTATION_MC()
     IF (write_output_2) THEN
        WRITE(UNIT=log_unit_nb,FMT=*) "EXIT THE EXPECTATION ANALYSIS OF STEP_2"
        WRITE(UNIT=log_unit_nb,FMT=*)
     END IF
     WRITE(UNIT=log_unit_nb,FMT=*)
  END IF

  IF (write_output_2) THEN
     CALL date_and_time (date=date_end, time=time_end)
     WRITE(UNIT=log_unit_nb,FMT=*) "STOP (dd:mm:yy - hh:mm:ss)  ", date_end(7:8), ":", &
          date_end(5:6), ":", date_end(3:4), " - ", &
          time_end(1:2), ":", time_end(3:4), ":", time_end(5:6)
     WRITE(UNIT=log_unit_nb,FMT=*) "STOP  SNP_PROGRAM"
  END IF

  !if write_output_2 close the unit "log_unit_nb" which has been opened in INITIATION()
  IF (write_output_2) THEN
     CLOSE(UNIT=log_unit_nb)
  END IF


  IF (write_output_3) THEN
     CLOSE(UNIT=dev_unit_nb)
  END IF

  STOP
END PROGRAM SNP_PROGRAM
