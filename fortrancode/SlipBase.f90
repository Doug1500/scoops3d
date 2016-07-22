       Subroutine SlipBase(m,i1,i2,j1,j2)
     
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!      This subroutine calculates elevation of slip surface base at
!      column midpoint
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!      Called by Scoops3D
!
!      VARIABLES  
!      delxy -- DEM grid resolution (delta x, delta y)
!      dx1,dx2,dy1,dy2 -- x and y grid spacing of partial column 
!      in(nx,ny) -- number of nodes of each column bounded by sphere or
!         equals -1 if truncated surface.
!      insphere(nx,ny) -- indicates whether DEM column node is bounded by sphere
!      inst(4) -- local variable for insphere
!      failsurf(nx,ny) -- array of failure surface elevations for each DEM column
!      ifailsurf -- flag for whether failure surface file is used
!      nonode -- flag for when there is no intersection between search
!         sphere and slope line between outside and inside nodes.
!         If nonode = 1, new node can not be calculated.
!      outnode(2) -- array of partial column nodes which fall outside
!         slip surface
!      radsq -- square of search sphere radius
!      subset(nx,ny) -- indicates set membership of DEM columns
!      xcen,ycen,zcen -- location of center of search sphere
!      xdem,ydem -- x and y locations of column
!      xmid -- x location of middle of column
!      x(4),y(4) -- x and y locations of each column node.
!      ymid -- y location of middle of column
!      xrad,yrad,zrad -- x,y, and z coordinates of equation for search sphere radius
!      zb(nx,ny) -- base of slip surface
!      zbase(4) -- local coordinates of base of column
!      zdem(nx,ny) -- DEM elevations
!      zdemnodes(nx,ny) -- averaged DEM elevations at each column node
!      zmid(i,j) -- elevation of slip surface at column midpoint
!      znodes(4) -- elevation of DEM nodes in column
!      zz -- square of zrad 
!
!
!            dx2 
!          4----3
!      dy1 |    | dy2   numbering scheme for corners of column base
!          |    | 
!          1----2
!            dx1
!      
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
      USE CommonData
      USE GridData, ONLY: delxy,xcen,ycen,zcen,radsq,nci1,nci2,nci3,nci4,nci5,zdemnodes,zdem,halfdelxy
      USE SetData, ONLY: subset,zb,insphere,in,xdem,ydem,zmid,outnodes
      USE FailSurfData, ONLY: ifailsurf,failsurf      
              
      IMPLICIT NONE        
      
      INTEGER, INTENT(in) :: m,i1,i2,j1,j2
      INTEGER :: outnode(2),i,j,nonode,kk,inst(4),nn

      REAL(pr) :: xmid,ymid,zbase(4),znodes(4),x(4),y(4)
      REAL(pr) :: dx1,dx2,dy1,dy2,zz,zrad,xrad,yrad

      CHARACTER*70 :: problemtype
      CHARACTER*120 :: errmessage,solution
      
      errmessage = ' '
      solution = ' '
      problemtype = 'in location of potential failure surface base'  
  
!     Loop through all columns in current set.                                         
      DO j = j1,j2      
        DO i = i1,i2
!     Bypass if column not within slip surface.
          IF (subset(i,j).ne.m) CYCLE
          
!     If defined failure surface exists and this is a full column
!       spherical failure surface, zmid, has already been calculated in Subsets.  
!       Set zmid to lower of defined or spherical failure surface at column midpoint.
          IF (ifailsurf.eq.1.and.in(i,j).eq.4) THEN
            IF (zmid(i,j).lt.failsurf(i,j)) THEN
              zmid(i,j) = failsurf(i,j)
            END IF
          ELSE

!     Assign local coordinates.
            x(1) = xdem(i,j)
            y(1) = ydem(i,j)
            IF (in(i,j).lt.4) THEN
              x(2) = x(1)+delxy
              y(3) = y(1)+delxy
              zbase(1) = zb(i,j)
              zbase(2) = zb(i+1,j)
              zbase(3) = zb(i+1,j+1)
              zbase(4) = zb(i,j+1)
              znodes(1) = zdemnodes(i,j)
              znodes(2) = zdemnodes(i+1,j)
              znodes(3) = zdemnodes(i+1,j+1)
              znodes(4) = zdemnodes(i,j+1)
              inst(1) = insphere(i,j)
              inst(2) = insphere(i+1,j)
              inst(3) = insphere(i+1,j+1)
              inst(4) = insphere(i,j+1)
!    Initialize column variables.
              outnodes(i,j,1) = 0
              outnodes(i,j,2) = 0
              nn = 0
              dx1 = delxy
              dx2 = delxy
              dy1 = delxy
              dy2 = delxy
!    Determine which column nodes are not in surface for partial columns.         
              DO kk = 1,4
                IF (inst(kk).ne.1) THEN
                  nn = nn+1      
                  IF (nn.gt.2) THEN 
                    errmessage = 'included column has fewer than 2 nodes'
                    Call WriteError(1,errmessage,problemtype,'no','no ',0,' ')                   
                  END IF          
                  outnodes(i,j,nn) = kk
                  outnode(1) = outnodes(i,j,1)
                  outnode(2) = outnodes(i,j,2)
                END IF
              END DO
!     Find new coordinates for nodes in partially contained columns.              
              CALL newnodes (i,j,zbase,znodes,dx1,dx2,dy1,dy2,nonode)

!     If column actually only has one node, due to precision errors
!     in initial calculation, throw out column.
              IF (nonode.eq.1) THEN
                subset(i,j) = 0
                CYCLE
              END IF
            END IF

!     Find x,y center of column used to calculate dip and apparent dip of slip base.
            SELECT CASE (in(i,j))
              CASE (4)
                xmid = x(1) + halfdelxy
                ymid = y(1) + halfdelxy
              
              CASE (2)   
                IF ((outnode(1).eq.1.and.outnode(2).eq.2).or. &
                    (outnode(1).eq.1.and.outnode(2).eq.4)) THEN
                  xmid = x(2)-(dx1+dx2)/4.0_pr
                  ymid = y(3)-(dy1+dy2)/4.0_pr
                ELSE
                  xmid = x(1)+(dx1+dx2)/4.0_pr
                  ymid = y(1)+(dy1+dy2)/4.0_pr
                END IF

              CASE (3)
                SELECT CASE (outnode(1))
                 CASE (1)
                    xmid = x(2)-(dx1+dx2)/4.0_pr
                    ymid = y(3)-(dy1+dy2)/4.0_pr

                  CASE (2)
                    xmid = x(1)+(dx1+dx2)/4.0_pr
                    ymid = y(3)-(dy1+dy2)/4.0_pr

                  CASE (3)
                    xmid = x(1)+(dx1+dx2)/4.0_pr
                    ymid = y(1)+(dy1+dy2)/4.0_pr
  
                  CASE (4)
                    xmid = x(2)-(dx1+dx2)/4.0_pr
                    ymid = y(1)+(dy1+dy2)/4.0_pr
                END SELECT
            END SELECT

!    Find distance from column slip surface midpoint to search node.
            yrad = ymid-ycen
            xrad = xmid-xcen                      
            zz = radsq-xrad*xrad-yrad*yrad 
            IF (zz.lt.0.0_pr) THEN  
              errmessage = 'cannot calculate distance from column slip surface midpoint to search node'
              Call WriteError(1,errmessage,problemtype,'no','no ',0,' ')              
!              Print *,radsq,xrad,yrad,zz
            END IF                
            ! IF (i.le.85.5) THEN
            !   zrad = SQRT(zz) -  nci1
            ! ELSE IF (i.ge.85.5 .and. i.le.89.3) THEN
            !   zrad = SQRT(zz) -  nci2
            ! ELSE IF (i.ge.89.3 .and. i.le.93.1) THEN
            !   zrad = SQRT(zz) -  nci3
            ! ELSE IF (i.ge.93.1 .and. i.le.96.9) THEN
            !   zrad = SQRT(zz) -  nci4
            ! ELSE IF (i.ge.96.9 .and. i.le.100.7) THEN
            !   zrad = SQRT(zz) -  nci5
            ! END IF
            zrad = SQRT(zz)
            ! IF (i.ge.85.5) THEN
            !  zmid(i,j) = zcen - zrad - 5.0
            ! ELSE
             ! zmid(i,j) = zcen - zrad
            ! END IF
            IF (i.le.85.5) THEN
              zmid(i,j) = zcen - zrad -  nci1
            ELSE IF (i.ge.85.5 .and. i.le.89.3) THEN
              zmid(i,j) = zcen - zrad -  nci2
            ELSE IF (i.ge.89.3 .and. i.le.93.1) THEN
              zmid(i,j) = zcen - zrad -  nci3
            ELSE IF (i.ge.93.1 .and. i.le.96.9) THEN
              zmid(i,j) = zcen - zrad -  nci4
            ELSE
              zmid(i,j) = zcen - zrad -  nci5
            END IF
          END IF
          IF (zmid(i,j).gt.zdem(i,j)) zmid(i,j) = zdem(i,j)
        END DO
      END DO
      
      END SUBROUTINE SlipBase
