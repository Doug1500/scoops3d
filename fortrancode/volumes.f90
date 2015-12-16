      SUBROUTINE Volumes (nset,goodset)
        
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!     This subroutine calculates the volumes and areas of all  
!     subsets intersected by sphere. Calculations assume planar
!     slip surface defined by intersection of search sphere with 
!     column corners.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!     Called by Scoops3D
!
!     VARIABLES
!
!      area(nnset) -- surface area of subset
!      avetop(nx,ny) -- average elevation of 4 corner nodes of intersected 
!        column calculated using the average of the four surrounding cells
!        for full columns, and by finding the intersection of the slip surface
!        with the smoothed DEM surface for partial columns.
!      colxy(nx,ny) -- length of each column side (delxy for full columns)
!      delxy -- DEM grid resolution (delta x, delta y)
!      i,j -- DEM grid array location
!      ifailsurf -- flag for whether failure surface file is used
!      failsurf(nx,ny) -- array of failure surface elevations for each DEM column  
!      in(nx,ny) -- number of nodes of each column bounded by sphere or
!        equals -1 if truncated surface, -2 if contains boundary node.
!      insphere(nx,ny) -- indicates whether DEM column node is bounded by sphere;
!      m -- set number
!      mini(nnset),maxi(nnset) -- minimum and maximum i bounds of each subset
!      minj(nnset),maxj(nnset) -- minimum and maximum j bounds of each subset
!      ncol(nnset) -- number of columns in subset
!      nonode -- flag for when there is no intersection between search
!         sphere and slope line between outside and inside nodes or partial
!         column is truncated.
!         If nonode = 1, new node can not be calculated.
!         If nonode = 2, truncated column.
!      nset -- number of subsets
!      outnodes(nx,ny,2) -- corner location of column nodes which fall
!         outside the search sphere.
!      set(nnset) -- number of nodes in a contiguous set
!      subset(nx,ny) -- indicates set membership of DEM columns
!      volume(nnset) -- volume of subset
!      z(4) -- used to find average column length
!      zb(nx+1,ny+1) -- base of slip surface at column nodes
!      zcol -- length of column between slip surface and land surface
!      zdem(nx,ny) -- DEM elevations
!      zdemnodes(nx+1,ny+1) -- DEM node elevations
!      zmid(nx,ny) -- base of slip surface at midpoint of intersected column
!      znodes(4) -- Local variable for DEM node elevations
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        USE CommonData
        USE GridData, ONLY: zdem,delxy,nx,ny,zdemnodes,cellarea
        USE SetData, ONLY: volume,area,mini,minj,maxi,maxj,set,setflag,&
                           subset,zb,ncol,in,insphere,surfslope,vol,&
                           colxy,outnodes,avetop,zmid
        USE FosData, ONLY: single
        USE FailSurfData, ONLY: ifailsurf,failsurf
        
        IMPLICIT NONE

        INTEGER, INTENT(in) :: nset
        INTEGER, INTENT(out) :: goodset
        INTEGER :: m,i,j,nn,inst(4),ii,nonode
        
        REAL(pr) :: dx1,dx2,dy1,dy2,zbase(4),dz(4),parea,znodes(4)
        REAL(pr) :: side1,side2,dzin
        
        REAL(pr), EXTERNAL :: colvol,trivol 

        CHARACTER*70 :: problemtype
        CHARACTER*120 :: errmessage,solution
      
        errmessage = ' '
        solution = ' '
        problemtype = 'calculating volume'  

        DO m = 1,nset
                  
          IF (m.eq.1) goodset = 0
          area(m) = 0.0_pr
          volume(m) = 0.0_pr
          ncol(m) = 0
          IF (set(m).ne.0) THEN
            DO j = minj(m),maxj(m)
              DO i = mini(m),maxi(m)
                IF (i.gt.nx.or.j.gt.ny) CYCLE
                IF (subset(i,j).eq.m.and.in(i,j).ge.2) THEN
                  outnodes(i,j,1) = 0
                  outnodes(i,j,2) = 0
                  nn = 0
                  dz = 0.0_pr
                  dzin = 0.0_pr
                  dx1 = delxy
                  dx2 = delxy
                  dy1 = delxy
                  dy2 = delxy
 
!   If user-defined failure surface exists and is above sphere elevation at column center
!   set failure base elevation to defined failure surface. Otherwise use sphere base elevation.              
                  IF (ifailsurf.eq.1) THEN
                    IF (failsurf(i,j).eq.rnull) CYCLE  ! don't include colums with null failure surface.
                    IF (zmid(i,j).lt.failsurf(i,j).and.in(i,j).eq.4) THEN
                      zbase = failsurf(i,j)
                    ELSE
                      zbase(1) = zb(i,j)
                      zbase(2) = zb(i+1,j)
                      zbase(3) = zb(i+1,j+1)
                      zbase(4) = zb(i,j+1)
                    END IF
                  ELSE
                    zbase(1) = zb(i,j)
                    zbase(2) = zb(i+1,j)
                    zbase(3) = zb(i+1,j+1)
                    zbase(4) = zb(i,j+1)
                  END IF
                  
!    Set znodes to averaged node elevations for partial column calculations.
!    Set znodes to DEM elevation for full columns
                  IF (in(i,j).lt.4) THEN
                    znodes(1) = zdemnodes(i,j)
                    znodes(2) = zdemnodes(i+1,j)
                    znodes(3) = zdemnodes(i+1,j+1)
                    znodes(4) = zdemnodes(i,j+1)
                  ELSE
                    znodes = zdem(i,j)
                   avetop(i,j) = zdem(i,j)
                  END IF
                  
                  inst(1) = insphere(i,j)
                  inst(2) = insphere(i+1,j)
                  inst(3) = insphere(i+1,j+1)
                  inst(4) = insphere(i,j+1)
             
!    Determine which column nodes are not in surface for partial columns.
                  IF (in(i,j).lt.4) THEN
                    DO ii = 1,4
                      IF (inst(ii).ne.1) THEN
                        nn = nn+1      
                        IF (nn.gt.2) THEN 
                          errmessage = 'included column has fewer than 2 nodes'
                          Call WriteError(1,errmessage,problemtype,'no','no ',0,' ')                          
                        END IF          
                        outnodes(i,j,nn) = ii
                      END IF
                    END DO
!     Find new coordinates for nodes in partially contained columns.              
                    CALL newnodes (i,j,zbase,znodes,dx1,dx2,dy1,dy2,nonode)

!     If partial column is truncated, mark set or end run if single surface and only one set.
                    IF (nonode.eq.2) THEN
                      IF (single.eq.1.and.nset.eq.1) THEN
                        errmessage = 'no valid potential failure for single surface (surface is truncated)'
                        solution = 'check coordinates for center of sphere'
                        Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')                         
                      ELSE
                        IF (setflag(m).ne.2) setflag(m) = 1
                      END IF
                    END IF
!     If column actually only has one node, due to precision errors
!     in initial calculation, throw out column.
                    IF (nonode.eq.1) THEN
                      subset(i,j) = 0
                      CYCLE
                    END IF
                    avetop(i,j) = sum(znodes)/4.0_pr                    
                  END IF
              
                  ncol(m) = ncol(m) + 1   
!     Calculate volume of prismoid columns.
                  dz = znodes-zbase
                  WHERE (dz.lt.0.0_pr) dz = 0.0_pr                        
                  vol(i,j) = colvol(dx1,dx2,dy1,dy2,dz)
 
                  IF (in(i,j).eq.3) THEN
                    SELECT CASE (outnodes(i,j,1))
                     CASE (1)
                      IF (dx1.le.0.0_pr) THEN
                        side1 = dx1
                        side2 = dy1
                        dzin = dz(2)
                      ELSE
                        side1 = dy1
                        side2 = dx1
                        dzin = dz(4)
                      END IF
                     CASE (2)
                      IF (dx1.le.0.0_pr) THEN
                        side1 = dx1
                        side2 = dy2
                        dzin = dz(1)
                      ELSE
                        side1 = dy2
                        side2 = dx1
                        dzin = dz(3)
                      END IF
                     CASE (3)
                      IF (dx2.le.0.0_pr) THEN
                        side1 = dx2
                        side2 = dy2
                        dzin = dz(4)
                      ELSE
                        side1 = dy2
                        side2 = dx2
                        dzin = dz(2)
                      END IF
                     CASE (4)
                      IF (dx2.le.0.0_pr) THEN
                        side1 = dx2
                        side2 = dy1
                        dzin = dz(3)
                      ELSE
                        side1 = dy1
                        side2 = dx2
                        dzin = dz(1)
                      END IF
                    END SELECT
                      
                    vol(i,j) = vol(i,j) + trivol(side1,side2,dzin)
                   parea = cellarea - 0.5_pr*(delxy + side1)*(delxy - side2)
                  ELSE
!     Calculate area of slide on DEM surface.  
                    parea = .250_pr * (dx1+dx2)*(dy1+dy2)
                  END IF
                  volume(m) = volume(m) + vol(i,j)
                  area(m) = area(m)+ parea  ! for 2D area
                  colxy(i,j,1) = dx1
                  colxy(i,j,2) = dx2
                  colxy(i,j,3) = dy1
                  colxy(i,j,4) = dy2

                END IF  ! IF (subset(i,j).eq.m.and.in(i,j).ne.-1) 

              END DO  ! loop on i
            END DO  ! loop on j
          END IF  ! IF (set(m).ne.0)
          IF (volume(m).gt.0.0_pr) goodset = 1
        END DO  ! loop on m
        
        
        END SUBROUTINE volumes
        
        
             
      FUNCTION colvol(ddx1,ddx2,ddy1,ddy2,z)
      
        USE CommonData
        USE GridData, ONLY: delxy
        IMPLICIT NONE
        
        REAL(pr) :: colvol
        REAL(pr), INTENT(in) :: ddx1,ddx2,ddy1,ddy2,z(4)
        REAL(pr) :: s0,s1,s2
        
        IF ((ddx1.eq.delxy.or.ddx1.le.0.0_pr).and.(ddx2.eq.delxy.or.ddx2.le.0.0_pr)) THEN
          s0 = .5_pr*(z(1)+z(4))*ddy1
          s1 = .125_pr*(z(1)+z(2)+z(3)+z(4)) &
                            *(ddy1+ddy2)
          s2 = .5_pr*(z(2)+z(3))*ddy2            
        ELSE
          s0 = .5_pr*(z(1)+z(2))*ddx1
          s1 = .125_pr*(z(1)+z(2)+z(3)+z(4)) &
                             *(ddx1+ddx2)
          s2 = .5_pr*(z(3)+z(4))*ddx2            
        END IF
           
        colvol = (1._pr/6._pr)*(s0+4._pr*s1+s2)*delxy
        
      END FUNCTION colvol
      
      FUNCTION trivol(ddx1,ddx2,z)
      
        USE CommonData
        USE GridData, ONLY: delxy
        IMPLICIT NONE
        
        REAL(pr) :: trivol
        REAL(pr), INTENT(in) :: ddx1,ddx2,z
           
        trivol = (1._pr/6._pr)*(z*abs(ddx1))*(delxy-ddx2)
        
      END FUNCTION trivol
