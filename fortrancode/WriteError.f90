      SUBROUTINE WriteError(errtype,problem,problemtype,solution,extraflag,extra1int,extra1ch)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!      This subroutine writes error messages to standard output and
!      to the error file
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!     VARIABLES
!      
!     errtype - 0 = warning message, 1 = fatal error (terminate in WriteError), 
!               2 = fatal error with extra info to provide (return to calling subroutine to terminate)   
!     extra1int - extra integer value to print out
!     extraflag - flag to indicate use of extra integer or character value
!     extra1ch - extra character value to print out 
!     problem - source of error message
!     problemtype - general source of error
!     solution - suggestion for eliminating error
!     suname - name of subroutine where error was generated

      IMPLICIT NONE

      INTEGER, INTENT(in) :: errtype,extra1int
      CHARACTER*3, INTENT(in) :: extraflag
      CHARACTER*220, INTENT(in) :: extra1ch 
      CHARACTER*70, INTENT(in) :: problemtype
      CHARACTER*120, INTENT(in) :: problem,solution

      SELECT CASE (errtype)
       CASE (0)     !  warning message
           WRITE (*,105)
           WRITE (39,105)           
           IF (extraflag(1:2).eq.'no') THEN
             WRITE (*,110) problem
             WRITE (39,111) problem
           ELSE
             IF (extraflag.eq.'int') THEN
               WRITE (*,113) problem,extra1int
               WRITE (39,113) problem,extra1int
             ELSE
               IF (extraflag(1:2).eq.'ch') THEN
                 WRITE (*,110) problem
                 WRITE (*,115) extra1ch
                 WRITE (39,111) problem               
                 WRITE (39,115) extra1ch   
               END IF   
             END IF
           END IF        
           WRITE (20,905) problemtype
           WRITE (20,910)
       CASE (1,2) ! fatal error - terminate execution of SCOOPS3D unless there is extra info to provide
           WRITE (*,100)
           WRITE (39,100)
           IF (extraflag(1:2).eq.'no') THEN
             WRITE (*,110) problem
             WRITE (39,111) problem
           ELSE
             IF (extraflag.eq.'int') THEN
               WRITE (*,113) problem,extra1int
               WRITE (39,113) problem,extra1int
             ELSE
               IF (extraflag(1:2).eq.'ch') THEN
                 WRITE (*,110) problem
                 WRITE (*,115) extra1ch
                 WRITE (39,111) problem               
                 WRITE (39,115) extra1ch   
               END IF   
             END IF
           END IF
           IF (solution(1:2).ne.'no') THEN
             WRITE (*,130) solution 
             WRITE (39,130) solution 
           END IF  
           WRITE (20,900) problemtype
           WRITE (20,910)           
           IF (errtype.eq.1) THEN
             CLOSE (20)
             CLOSE (39)
             STOP
           END IF   
       END SELECT
 
100    FORMAT ('*** ERROR - Scoops3D  execution terminated')
105    FORMAT ('*** WARNING')
110    FORMAT ('problem: ',A)
111    FORMAT ('problem: ',A)
113    FORMAT ('problem: ',A,I6)
115    FORMAT ('   ',A)
130    FORMAT ('possible solution: ',A)
900    FORMAT ('*** ERROR ',A)
905    FORMAT ('*** WARNING ',A)
910    FORMAT ('    Details in "filename_errors_out.txt" ')      

       END SUBROUTINE WriteError
