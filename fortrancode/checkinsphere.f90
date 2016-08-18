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
        USE GridData, ONLY: delxy,xcen,ycen,zcen,radsq,zdemnodes,nci1,nci2,nci3,nci4,nci5,nci6,nci7,nci8,nci9,nci10,xst
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
          mz = sqrt(mz2)
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
              mz = sqrt(mz2)
              ! zb(ii,jj) = zcen - mz
              ! IF (ii.le.85.5) THEN
              !   zb(ii,jj) = zcen - mz -  nci1
              ! ELSE IF (ii.ge.85.5 .and. ii.le.89.3) THEN
              !   zb(ii,jj) = zcen - mz -  nci2
              ! ELSE IF (ii.ge.89.3 .and. ii.le.93.1) THEN
              !   zb(ii,jj) = zcen - mz -  nci3
              ! ELSE IF (ii.ge.93.1 .and. ii.le.96.9) THEN
              !   zb(ii,jj) = zcen - mz -  nci4
              ! ELSE
              !   zb(ii,jj) = zcen - mz - nci5
              ! END IF
              IF (ii.le.xst) THEN
                xst = ii
              END IF

              IF (ii.le.(xst+nci10)) THEN
                zb(ii,jj) = zcen - mz
              ELSE IF (ii.ge.(xst+nci10) .and. ii.le.(xst+nci10*2.0)) THEN
                zb(ii,jj) = zcen - mz -  nci1
              ELSE IF (ii.ge.(xst+nci10*2.0) .and. ii.le.(xst+nci10*3.0)) THEN
                zb(ii,jj) = zcen - mz -  nci2
              ELSE IF (ii.ge.(xst+nci10*3.0) .and. ii.le.(xst+nci10*4.0)) THEN
                zb(ii,jj) = zcen - mz -  nci3              
              ELSE IF (ii.ge.(xst+nci10*4.0) .and. ii.le.(xst+nci10*5.0)) THEN
                zb(ii,jj) = zcen - mz -  nci4              
              ELSE IF (ii.ge.(xst+nci10*5.0) .and. ii.le.(xst+nci10*6.0)) THEN
                zb(ii,jj) = zcen - mz -  nci5              
              ELSE IF (ii.ge.(xst+nci10*6.0) .and. ii.le.(xst+nci10*7.0)) THEN
                zb(ii,jj) = zcen - mz -  nci6              
              ELSE IF (ii.ge.(xst+nci10*7.0) .and. ii.le.(xst+nci10*8.0)) THEN
                zb(ii,jj) = zcen - mz -  nci7              
              ELSE IF (ii.ge.(xst+nci10*8.0) .and. ii.le.(xst+nci10*9.0)) THEN
                zb(ii,jj) = zcen - mz -  nci8              
              ELSE
                zb(ii,jj) = zcen - mz -  nci9              
              END IF

              ! IF (ii.le.(xst+nci10)) THEN
              !   zb(ii,jj) = zcen - mz
              ! ELSE IF (ii.ge.(xst+nci10) .and. ii.le.(xst+nci10*2.0)) THEN
              !   zb(ii,jj) = nci1
              ! ELSE IF (ii.ge.(xst+nci10*2.0) .and. ii.le.(xst+nci10*3.0)) THEN
              !   zb(ii,jj) = nci2
              ! ELSE IF (ii.ge.(xst+nci10*3.0) .and. ii.le.(xst+nci10*4.0)) THEN
              !   zb(ii,jj) = nci3              
              ! ELSE IF (ii.ge.(xst+nci10*4.0) .and. ii.le.(xst+nci10*5.0)) THEN
              !   zb(ii,jj) = nci4              
              ! ELSE IF (ii.ge.(xst+nci10*5.0) .and. ii.le.(xst+nci10*6.0)) THEN
              !   zb(ii,jj) = nci5              
              ! ELSE IF (ii.ge.(xst+nci10*6.0) .and. ii.le.(xst+nci10*7.0)) THEN
              !   zb(ii,jj) = nci6              
              ! ELSE IF (ii.ge.(xst+nci10*7.0) .and. ii.le.(xst+nci10*8.0)) THEN
              !   zb(ii,jj) = nci7              
              ! ELSE IF (ii.ge.(xst+nci10*8.0) .and. ii.le.(xst+nci10*9.0)) THEN
              !   zb(ii,jj) = nci8              
              ! ELSE
              !   zb(ii,jj) = nci9              
              ! END IF
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
