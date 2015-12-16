        SUBROUTINE checksets(i,j,i1,j1,nset,ierr) 
        
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!       This subroutine determines set membership of columns
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!        Called by Subsets
!
!        VARIABLES
!
!        i,j -- DEM column number
!        ii,jj -- nodes of DEM column
!        ierr -- error flag, 1 indicates more subsets found than nnset 
!        in(nx,ny) -- number of nodes of each column bounded by sphere
!        insphere(nx,ny) -- indicates whether DEM column node is bounded by sphere
!        m -- set number
!        mini(nnset),maxi(nnset) -- minimum and maximum i bounds of each subset
!        minj(nnset),maxj(nnset) -- minimum and maximum j bounds of each subset
!        nnset -- maximum number of subsets
!        nset -- number of subsets found
!        nsetmax -- flag for exceeding number of subsets allowed nnset
!        set(nnset) -- number of columns in a contiguous set
!        setmethod -- flag for location of contiguous column in set
!        subset(nx,ny) -- indicates set membership of DEM columns
!        tempm -- variable used in determining inclusion in set
!        xcen,ycen,zcen -- center of search sphere relative to DEM origin
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        
        USE CommonData, ONLY: bfile,filin
        USE GridData, ONLY: xcen,ycen,zcen,nbdy
        USE SetData, ONLY: insphere,mini,minj,maxi,maxj,in,set,subset,&
                           nnset,nsetmax
        
        IMPLICIT NONE

        INTEGER, INTENT(in) :: i,j,i1,j1
        INTEGER, INTENT(out) :: ierr
        INTEGER, INTENT(inout) :: nset
        INTEGER :: ii,jj,m,mm,setmethod,tempm,numset,node,noset

        CHARACTER*70 :: problemtype
        CHARACTER*120 :: errmessage,solution
      
        errmessage = ' '
        solution = ' '
        problemtype = 'determining subsets'                  
        
!     Check how many nodes of column are intersected by trial sphere.
         node = 0 
         noset = 0
         ierr = 0 
         DO jj = j,j+1
           DO ii = i,i+1
             IF (insphere(ii,jj).ne.0) THEN
               node = node + 1
!     If truncated surface node, set flag.
               IF (insphere(ii,jj).eq.-1.and.noset.ne.2) noset = 1
!     If boundary node, set flag.
               IF (nbdy(ii,jj).eq.1) noset = 2
             END IF
           END DO
         END DO
 
!    Column must have two or more adjacent nodes contained to be counted.
         IF (node.gt.1) THEN
!    If one of each set of opposite nodes is intersected by surface,
!    that means two adjacent nodes are intersected, as required for column
!    to be contained in surface.
           IF (((insphere(i,j).ne.0).or.(insphere(i+1,j+1).ne.0)).and.&
              ((insphere(i+1,j).ne.0).or.(insphere(i,j+1).ne.0))) THEN
!    Set array in(i,j) to represent number of nodes in slip surface
!    for future calculation of partial columns.
             in(i,j) = node 
!    If column contains truncating surface or boundary node, flag column. 
!    Keep column in for now until all columns are sorted into sets so
!    the whole set can be flagged.
             IF (noset.eq.1) in(i,j) = -1 
             IF (noset.eq.2) in(i,j) = -2
!    Loop through subsets to find any contiguous existing set with
!    lowest m. Only need to check adjacent columns to the left and below
!    current column because of search order.
             setmethod = 0
             tempm = 0
             mm = 0
             IF (nset.gt.0) THEN ! at least one column already placed in a set. 
               numset = nset          
               DO m = 1,numset               
!    If (i-1,j) or (i,j-1) adjacent column belongs to set m, then current 
!    column may also be a part of set m if columns have a shared node.
                 IF (j.gt.j1) THEN 
                   IF ((subset(i,j-1).eq.m).and. & 
!    Check if shared border nodes are in slip surface.
                      (insphere(i,j).ne.0.or.insphere(i+1,j).ne.0)) THEN
!    Adjacent column below belongs to set m, so this column will also.
                     tempm = m
                     setmethod = 2
                     IF (i.gt.i1) THEN
                       IF ((subset(i-1,j).gt.m).and. & 
                          (insphere(i,j).ne.0.or.insphere(i,j+1).ne.0)) THEN
!    Column to left is in higher set and needs to be changed to set m.
                         mm = subset(i-1,j)
!    If set number is current high, reduce nset.
                         IF (mm.eq.nset) nset = nset-1
                         setmethod = 1
                       END IF
                     END IF
                   END IF
                 END IF
                
                 IF (i.gt.i1) THEN
                   IF ((subset(i-1,j).eq.m).and. &
                     (insphere(i,j).ne.0.or.insphere(i,j+1).ne.0)) THEN
!    Adjacent column to left belongs to set m, so this column will also.
                     tempm = m
                     setmethod = 2
                     IF (j.gt.j1) THEN
                       IF ((subset(i,j-1).gt.m).and. &
                          (insphere(i,j).ne.0.or.insphere(i+1,j).ne.0)) THEN
!    Adjacent column below is in higher set and needs to be changed to set m.
                         mm = subset(i,j-1)
!    If set number is current high, reduce nset.
                         IF (mm.eq.nset) nset = nset-1
                         setmethod = 1
                       END IF
                     END IF 
                   END IF
                 END IF
               
                IF (setmethod.ne.0) EXIT ! If any contiguous columns found.
                                 
               END DO  ! loop on m            
             END IF  ! (nset.ne.0...)
                   
                    
!     If first column found, or if no contiguous columns have been found.
             IF (setmethod.eq.0) THEN
!     If first set found or need to define new set, increment nset
               IF (nset.eq.0) THEN
                 nset = nset + 1
               ELSE
                 IF (set(nset).gt.0) nset = nset + 1
               END IF
               IF (nset.gt.nnset) THEN
                 errmessage = 'too many subsets'
                 solution = 'increase nnset in modules.f90'
                 Call WriteError(0,errmessage,problemtype,solution,'no ',0,' ')
                 PRINT *,'** greater than ',nnset,' trial surfaces at search node, i = ',i,'j = ',j
                 PRINT *,'Scoop3D will continue to next search node'
                 WRITE (39,*) '** greater than ',nnset,' trial surfaces at search node, i = ',i,'j = ',j         
                 WRITE (39,*) 'Scoop3D will continue to next search node'
                 nsetmax = nnset
                 ierr = 1
                 RETURN
               END IF    
!    Assign column set membership.
               subset(i,j) = nset
!    Start column count for set membership.
               set(nset) = 1
!    Keep track of min and max array bounds for columns in set.
               mini(nset) = i
               maxi(nset) = i
               minj(nset) = j
               maxj(nset) = j
             ELSE
                                
!     If found contiguous columns in two different subsets, need to
!     change all the columns and nodes in columns in the higher  
!     numbered subset (mm) to new value of m.
               m = tempm
               IF (setmethod.eq.1) THEN               
                 DO jj = minj(mm),maxj(mm)
                   DO ii = mini(mm),maxi(mm)
                     IF (set(mm).eq.0) EXIT !  No more columns in set.
                     IF (subset(ii,jj).eq.mm) THEN
                       subset(ii,jj) = m
!    Increment column count for set m and decrease count for set mm.
                       set(m) = set(m) + 1
                       set(mm) = set(mm) - 1                          
                     END IF         
                   END DO
                 END DO
                 IF (set(mm).ne.0) THEN
                    errmessage = 'subset not emptied '
                    Call WriteError(1,errmessage,problemtype,'no','no ',0,' ')  
                 END IF

!    Find new set array boundaries from merged sets.
                 mini(m) = MIN(mini(m),mini(mm))
                 maxi(m) = MAX(maxi(m),maxi(mm))
                 minj(m) = MIN(minj(m),minj(mm))
                 maxj(m) = MAX(maxj(m),maxj(mm))   
                 mini(mm) = 0
                 minj(mm) = 0
                 maxi(mm) = 0
                 maxj(mm) = 0       
               END IF 

!    Assign set membership to current column and nodes and update
!    array boundaries.                                                        
               subset(i,j) = m
               set(m) = set(m) + 1
               mini(m) = MIN(mini(m),i)
               maxi(m) = MAX(maxi(m),i)
               minj(m) = MIN(minj(m),j)
               maxj(m) = MAX(maxj(m),j)  
             END IF  ! (setmethod.eq.0)
           END IF  ! (insphere...)
         END IF  ! (node.gt.1)
         
         END SUBROUTINE checksets
