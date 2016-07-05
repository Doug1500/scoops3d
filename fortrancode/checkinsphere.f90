        SUBROUTINE checkinsphere(ii,jj,nodecount) 
        
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!       This subroutine determines whether node i,j is contained in
!       current slip surface.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!        Called by Subsets
!
!        VARIABLES
!
!        colcount -- number of columns intersected by sphere
!        delxy -- DEM grid resolution (delta x, delta y)
!        ierr -- error flag, 1 indicates more subsets found than nnset 
!        insphere(nx,ny) -- indicates whether DEM column node is bounded by sphere
!        mz -- elevation difference from sphere center to sphere border at 
!           current i,j location.
!        mz2 -- square of mz.
!        nodecount -- number of nodes in column contained in sphere
!        radsq -- square of search sphere radius
!        rnull -- real null value used throughout program 
!        xdem,ydem -- location of column in DEM coordinates
!        xcen,ycen,zcen -- center of search sphere relative to DEM origin
!        xdiffsq,ydiffsq -- used to calculate components of radius
!        zb(nx,ny) -- base of slip surface at each node (not column)
!        zdem(nx,ny) -- DEM elevations
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        USE CommonData
        USE GridData, ONLY: delxy,xcen,ycen,zcen,radsq,zdemnodes,nci1,nci2,nci3,nci4,nci5
        USE SetData, ONLY: insphere,zb,xdem,ydem,omini,omaxi,ominj,omaxj
        
        IMPLICIT NONE

        INTEGER, INTENT(in) :: ii,jj
        INTEGER, INTENT(inout) :: nodecount

        REAL(pr) :: xdiffsq,ydiffsq,mz2,mz,x,y
      
        x = xdem(ii,jj)
        y = ydem(ii,jj)
        xdiffsq = (x-xcen)*(x-xcen)
        ydiffsq = (y-ycen)*(y-ycen)
        mz2 = radsq-(xdiffsq+ydiffsq)  
!     Check whether DEM node is intersected by sphere. mz2 will
!     be less than zero if x or y of node is outside sphere.     
        !********************************************************************
        IF (mz2 .gt. 0.0_pr) THEN
          IF (jj.le.20) THEN
            mz = sqrt(mz2) -  nci1
          ELSE IF (jj.le.30) THEN
            mz = sqrt(mz2) -  nci2
          ELSE IF (jj.le.40) THEN
            mz = sqrt(mz2) -  nci3
          ELSE IF (jj.le.50) THEN
            mz = sqrt(mz2) -  nci4
          ELSE
            mz = sqrt(mz2) -  nci5
          END IF

!     If DEM elevation is above bottom of sphere at x,y location
          IF (zdemnodes(ii,jj).ge.(zcen-mz)) THEN
!     Mark case where elevation of the DEM surface is completely above the
!     sphere so that sets containing these columns can be discarded as  
!     truncated surfaces.
            IF (zdemnodes(ii,jj).ge.(zcen+mz)) THEN
              insphere(ii,jj) = -1
              zb(ii,jj) = rnull
            ELSE    
              insphere(ii,jj) = 1 
!     Find elevation of slip surface.
              zb(ii,jj) = zcen - mz 
            END IF
!     Keep track of min and max locations and total number 
!     of all intersected nodes.           
            omini = MIN(omini,ii)
            omaxi = MAX(omaxi,ii)
            omaxj = MAX(omaxj,jj)
            ominj = MIN(ominj,jj)
            nodecount = nodecount + 1 
          ELSE ! DEM surface is completely below sphere boundary.
            insphere(ii,jj) = 0
            zb(ii,jj) = rnull
          END IF
        ELSE ! DEM node is outside x,y bounds of sphere.
          insphere(ii,jj) = 0
          zb(ii,jj) = rnull
        END IF       
     
        END SUBROUTINE checkinsphere
