        Subroutine Criteria(nset,newrad,zavg,crit1,crit2,c1min,c1max,c2min,c2max,&
                            va1max,nstate)
        
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!       This subroutine determines if intersected subsets meet volume and area 
!       constraints specified by user.  If constraints are met, factors of safety 
!       will be calculated next. Otherwise, the radius will be adjusted in a way
!       that depends on the current circumstances.  
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!       Called by Scoops3D
!
!       VARIABLES
!
!        c1min -- min. area or volume for search, primary criterion
!        c1max -- max. area or volume for search, primary criterion
!        c2min -- min. area or volume for search, not primary criterion
!        c2max -- max. area or volume for search, not primary criterion
!        crit1(nset) -- variable for volume or area, whichever is 
!           primary criterion
!        crit2(nset) -- variable for volume or area, whichever is 
!           not primary criterion
!        m -- set number
!        mmax -- set number of subset with maximum primary criterion
!        newrad -- indicates whether valid initial range subset has been found.
!           (1=no, 0=yes, 2=no valid sets after 10 radius adjustments.)
!        nin -- number of subsets that fall within criteria range
!        nrange(nnset) -- array of flags for whether set fits volume or area criteria 
!           range
!        nnset -- maximum number of subsets
!        nset -- number of subsets found
!        nstate -- flag used to indicate whether valid sets
!           were found at a particular search node.
!           100 = recalculate radius to find first good set at a search node.
!           300 = increase radius by dr and continue. This flag occurs if good
!                 sets were found, or only sets not grown out of range are too
!                 small to be in range.
!           500 = go to next search grid point because no good sets possible.
!        setflag(nnset) -- indicates valid subset (0), or invalid by containing
!           truncated node (1), or adjacent to DEM boundary (2), or empty (-1)
!        tol -- tolerance on primary control criterion for calculating initial 
!           radius.  (i.e. initial volume must fall between volume and volume+tol)
!        va1max -- maximum volume or area of subsets, whichever is 
!           primary criterion
!        va2max -- max area or volume of subsets, whichever is 
!           not primary criterion
!        zavg(nnset) -- distance from center of subset to sphere center
!        zavgmax -- distance from center of subset to sphere center
!           of subset with maximum primary criterion
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        USE CommonData
        USE GridData, ONLY: xcen,ycen,zcen,zdem,delxy,nx,ny
        USE SetData, ONLY: nrange,nnset,setflag,mini,minj,maxi,maxj,insphere
        USE SearchData, ONLY: zavgmax,tol,iincflag,allin,minclude
        
        IMPLICIT NONE
        
        INTEGER, INTENT(in) :: nset
        INTEGER, INTENT(out) :: nstate
        INTEGER, INTENT(inout) :: newrad
        INTEGER :: mmax,m,nin
        
        REAL(pr), INTENT(in) :: zavg(nnset),crit1(nnset),crit2(nnset)
        REAL(pr), INTENT(in) :: c1min,c1max,c2min,c2max
        REAL(pr), INTENT(out) :: va1max
        REAL(pr) :: va2max,xrel,yrel,zrel
        
        nstate = 0

!     If at a new search point, or no valid initial range radius
!     yet found, find subset with maximum primary criterion. This set
!     will be used to control the initial radius value for this
!     search node (unless another set outgrows this one during the 
!     initial radius adjustment cycle.)
        IF (newrad.eq.1) THEN
          va1max = 0.0_pr
          DO m = 1,nset
            IF (va1max.le.crit1(m)) THEN        
              mmax = m
              va1max = crit1(m)
            END IF           
          END DO
          va2max = crit2(mmax)

!     If criteria 1 is too big and criteria 2 not too small,
!     recalculate smaller radius.
          IF ((va1max.ge.c1min+tol).and.(va2max.ge.c2min)) nstate = 100

!     If criteria 1 is too small and criteria 2 not too big,
!     and set does not contain DEM boundary node,
!     recalculate larger radius.
          IF ((va1max.lt.c1min).and.(setflag(mmax).ne.2))  THEN
            IF ((va2max.le.c2max).or.(c2max.eq.0.0_pr)) nstate = 100
          END IF
              
          IF (nstate.eq.100) THEN
           xrel = (0.5_pr*(mini(mmax)+maxi(mmax))-0.5_pr)*delxy - xcen
           yrel = (0.5_pr*(minj(mmax)+maxj(mmax))-0.5_pr)*delxy - ycen
           zrel = zavg(mmax) - zcen
           zavgmax = SQRT(xrel*xrel + yrel*yrel + zrel*zrel)
         END IF
                 
        END IF
        
!     If initial criteria ranges are just right at this radius, or
!     can't be adjusted to make them right, or
!     prior good radius has been incremented by delr, or
!     radius has been adjusted > 10 times and still can't hit
!     initial ranges, check for any and all valid sets between min                        
!     and max ranges. Set nrange=1 for all valid sets.
        IF (nstate.eq.0) THEN
            nin = 0
            DO m = 1,nset
              nrange(m) = 0
              IF ((crit1(m).ge.c1min).and.(crit1(m).le.c1max).and.&
                  (setflag(m).eq.0)) THEN
                IF ((crit2(m).ge.c2min).and.((crit2(m).le.c2max).or.(c2max.eq.0.0_pr))) THEN
!    Radius will now only change by increments of delr until out of range
!    or criteria cannot be met.
                  newrad = 0  
! check to make sure all cells in include grid are contained in failure surface
                  IF (iincflag.eq.1.and.allin.eq.1.and.m.eq.minclude) nrange(m) = 1
                  IF (iincflag.eq.0) nrange(m) = 1
                  IF (nrange(m).eq.1)  nin = nin + 1                
                END IF
              END IF
            END DO
 
          
!     If no sets fit range, check to see if there are any sets
!     that are too small for crit1 or crit2 and don't contain
!     boundary nodes. If so, will continue to increase 
!     radius by delr until no more possible good sets (nstate=300). 
!     If not, will move to next search node (nstate=500).
            IF (nin.eq.0) THEN
              nstate = 500
              DO m = 1,nset
                IF (iincflag.eq.0) THEN
                  IF ((crit1(m).le.c1max).and.((crit2(m).le.c2max).or.(c2max.eq.0.0_pr)).and.&
                      (setflag(m).ne.2)) nstate = 300
                ELSE
                  IF ((crit1(m).le.c1max).and.((crit2(m).le.c2max).or.(c2max.eq.0.0_pr)).and.&
                      (allin.ne.-2)) nstate = 300
                END IF
              END DO
            END IF           

        END IF
                        
        END SUBROUTINE criteria
            
