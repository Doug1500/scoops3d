       LOGICAL FUNCTION isupported(heading)
        
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!        This subroutine checks to see if a feature is  available in the current version of SCOOPS
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!

        CHARACTER*2, INTENT(in) :: heading
        
        CHARACTER*25 :: subname 
        CHARACTER*70 :: problemtype
        CHARACTER*100 :: errmessage,solution
      
        errmessage = ' '
        solution = ' '
        subname = 'isupported.f90'
        problemtype = 'reading input file'  

        SELECT CASE (heading)
                       
            CASE ('ro','RO','fl','FL','ii','II','fa','FA','in','IN')
              isupported = .FALSE.
 
            CASE ('if','IF','nr','NR','su','SU')             
              isupported = .FALSE.
    
            CASE DEFAULT
              isupported = .TRUE.
              
        END SELECT                                                      

        END FUNCTION isupported
