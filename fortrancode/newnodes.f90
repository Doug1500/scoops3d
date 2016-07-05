        SUBROUTINE newnodes (k,l,zbase,znodes,dx1,dx2,dy1,dy2,nonode)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!       This subroutine calculates the location of the intersection of
!       the search sphere with the boundary of a column that is
!       partially contained in the slip surface.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!      Called by Volumes
!
!      VARIABLES
!
!      delxy -- DEM grid resolution (delta x, delta y)
!      dx1,dx2 -- x-grid spacing of new column relative to inside nodes
!      dxy -- local variable for column resolution used to solve for
!         new node location
!      dy1,dy2 -- y resolution of column using new nodes
!      in(k,l) -- # of nodes contained in slip surface
!      m -- slope of DEM grid between inside and outside nodes
!      nonode -- flag for when there is no intersection between search
!         sphere and slope line between outside and inside nodes or partial
!         column is truncated.
!         If nonode = 1, new node can not be calculated.
!         If nonode = 2, truncated column.
!      outnode(2) -- array of column node numbers outside search sphere
!      p1,p2 -- each part of quadratic solution; add or subtract to 
!         find correct solution
!      rad -- radius of trial sphere
!      radsq -- square of search sphere radius
!      xcen,ycen,zcen -- center of search sphere relative to DEM origin
!      x1 -- x location of column nodes 1 and 4
!      x2 -- x location of column nodes 2 and 3
!      y1 -- y location of column nodes 1 and 2
!      y2 -- y location of column nodes 3 and 4
!      z1,z2 -- elevation of DEM at new nodes which are at the 
!         intersection of the search sphere with the DEM surface
!      zbase(4) -- elevation of base of slip surface of each node in
!         column 
!      znodes(4) -- elevation of DEM nodes in column 
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
        USE CommonData
        USE GridData, ONLY: xcen,ycen,zcen,delxy,rad,radsq
        USE SetData, ONLY: in,xdem,ydem,outnodes
        
        IMPLICIT NONE

        INTEGER, INTENT(in) :: k,l
        INTEGER, INTENT(out) :: nonode

        REAL(pr), INTENT(inout) :: zbase(4),znodes(4),dx1,dx2,dy1,dy2
        REAL(pr) :: x1,x2,y1,y2,z1,z2
        INTEGER :: outnode1,outnode2
     

        x1 = xdem(k,l)
        y1 = ydem(k,l)
        x2 = x1 + delxy
        y2 = y1 + delxy
        outnode1 = outnodes(k,l,1)
        outnode2 = outnodes(k,l,2)
        
        nonode = 0
 
!     If 2 column nodes are not intersected by search sphere calculate
!     intersection of slip surface with DEM surface.
        IF (in(k,l).eq.2) THEN
          IF (((outnode1.eq.1).and.(outnode2.eq.4)).or.&
              ((outnode1.eq.2).and.(outnode2.eq.3))) THEN
            CALL getnode (x1,x2,y1,znodes(1),znodes(2),ycen,xcen,dx1,z1,nonode)
            IF (nonode.ne.0) RETURN
            CALL getnode (x1,x2,y2,znodes(4),znodes(3),ycen,xcen,dx2,z2,nonode)
            IF (nonode.ne.0) RETURN
            IF (outnode1.eq.1) THEN
!    Find distance from column node in surface to intersection with sphere.
              dx1 = delxy - dx1
              dx2 = delxy - dx2
!    Set slip base and top elevation of out of bounds node numbers equal to 
!    DEM elevation at new node.
              zbase(1) = z1
              zbase(4) = z2
              znodes(1) = z1
              znodes(4) = z2
            END IF
            IF (outnode1.eq.2) THEN
              zbase(2) = z1
              zbase(3) = z2
              znodes(2) = z1
              znodes(3) = z2
            END IF
          END IF

          IF (((outnode1.eq.1).and.(outnode2.eq.2)).or.&
              ((outnode1.eq.3).and.(outnode2.eq.4))) THEN
            CALL getnode (y1,y2,x1,znodes(1),znodes(4),xcen,ycen,dy1,z1,nonode)
            IF (nonode.ne.0) RETURN
            CALL getnode (y1,y2,x2,znodes(2),znodes(3),xcen,ycen,dy2,z2,nonode)
            IF (nonode.ne.0) RETURN
            IF (outnode1.eq.1) THEN
              dy1 = delxy - dy1
              dy2 = delxy - dy2
              zbase(1) = z1
              zbase(2) = z2
              znodes(1) = z1
              znodes(2) = z2
            END IF
            IF (outnode1.eq.3) THEN
              zbase(3) = z2
              zbase(4) = z1
              znodes(3) = z2
              znodes(4) = z1
            END IF
          END IF
!     If a new dx or dy equals zero, the node is not actually in the
!     sphere and this column has only one node in, so don't use this column.
          IF (dy1.le.0.0_pr.or.dy2.le.0.0_pr.or.&
              dx1.le.0.0_pr.or.dx2.le.0.0_pr) nonode = 1
        END IF 

!     If only one node is uncontained, need to determine new node that
!     creates larger quadrilateral.
        IF (in(k,l).eq.3) THEN
          IF (outnode1.eq.1) THEN
            CALL getnode (x1,x2,y1,znodes(1),znodes(2),ycen,xcen,dx1,z1,nonode)
            IF (nonode.eq.2) RETURN
            IF (nonode.eq.1) dx1 = delxy
            CALL getnode (y1,y2,x1,znodes(1),znodes(4),xcen,ycen,dy1,z2,nonode)
            IF (nonode.eq.2) RETURN
            IF (nonode.eq.1) dy1 = delxy
!    If no good new nodes found don't use this column.
            IF (dx1.eq.delxy.and.dy1.eq.delxy) THEN
              nonode = 1
              RETURN
            ELSE
              nonode = 0
            END IF
!    Find distance from column node in surface to intersection with sphere.
            dx1 = delxy - dx1
            dy1 = delxy - dy1
!    Determine which new node results in larger column and keep that one.
            IF (dx1.ge.dy1) THEN
              dy1 = -dy1
              znodes(1) = z1
              zbase(1) = z1
            ELSE
              dx1 = -dx1
              znodes(1) = z2
              zbase(1) = z2
            END IF                 
          END IF
          IF (outnode1.eq.2) THEN
            CALL getnode (x1,x2,y1,znodes(1),znodes(2),ycen,xcen,dx1,z1,nonode)
            IF (nonode.eq.2) RETURN
            IF (nonode.eq.1) dx1 = 0.0_pr
            CALL getnode (y1,y2,x2,znodes(2),znodes(3),xcen,ycen,dy2,z2,nonode)
            IF (nonode.eq.2) RETURN
            IF (nonode.eq.1) dy2 = delxy
!    If no good new nodes found don't use this column.
            IF (dx1.eq.0.0_pr.and.dy2.eq.delxy) THEN
              nonode = 1
              RETURN
            ELSE
              nonode = 0
            END IF
            dy2 = delxy - dy2
            IF (dx1.ge.dy2) THEN
              dy2 = -dy2
              znodes(2) = z1
              zbase(2) = z1
            ELSE
              dx1 = -dx1
              znodes(2) = z2
              zbase(2) = z2
            END IF                 
          END IF
          IF (outnode1.eq.3) THEN
            CALL getnode (y1,y2,x2,znodes(2),znodes(3),xcen,ycen,dy2,z1,nonode)
            IF (nonode.eq.2) RETURN
            IF (nonode.eq.1) dy2 = 0.0_pr
            CALL getnode (x1,x2,y2,znodes(4),znodes(3),ycen,xcen,dx2,z2,nonode)
            IF (nonode.eq.2) RETURN
            IF (nonode.eq.1) dx2 = 0.0_pr
!    If no good new nodes found don't use this column.
            IF (dy2.eq.0.0_pr.and.dx2.eq.0.0_pr) THEN
              nonode = 1
              RETURN
            ELSE
              nonode = 0
            END IF
            IF (dx2.ge.dy2) THEN
              dy2 = -dy2
              znodes(3) = z2
              zbase(3) = z2
            ELSE
              dx2 = -dx2
              znodes(3) = z1
              zbase(3) = z1
            END IF                 
          END IF
          IF (outnode1.eq.4) THEN
            CALL getnode (x1,x2,y2,znodes(4),znodes(3),ycen,xcen,dx2,z1,nonode)
            IF (nonode.eq.2) RETURN
            IF (nonode.eq.1) dx2 = delxy
            CALL getnode (y1,y2,x1,znodes(1),znodes(4),xcen,ycen,dy1,z2,nonode)
            IF (nonode.eq.2) RETURN
            IF (nonode.eq.1) dy1 = 0.0_pr
!    If no good new nodes found don't use this column.
            IF (dx2.eq.delxy.and.dy1.eq.0.0_pr) THEN
              nonode = 1
              RETURN
            ELSE
              nonode = 0
            END IF
            dx2 = delxy - dx2
            IF (dx2.ge.dy1) THEN
              dy1 = -dy1
              znodes(4) = z1
              zbase(4) = z1
            ELSE
              dx2 = -dx2
              znodes(4) = z2
              zbase(4) = z2
            END IF                 
          END IF
        END IF
        
        END SUBROUTINE newnodes

  
        SUBROUTINE getnode (xya,xyb,xy0,za,zb,xycen,yxcen,dxy,z0,nonode)
 
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!       This subroutine calculates the location of a new column node at 
!       intersection of slip base with DEM column boundary by equating the
!       sphere equation with the equation for the line between the excluded 
!       and included nodes and solving the resulting quadratic equation
!       for the x or y intersection of the column side, depending on which
!       side is being checked.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!      Called by Newnodes
!       
!      VARIABLES
!
!      a,b,c -- coefficients of quadratric used to solve for new node
!      del -- variable used to take into account precision errors
!      m -- DEM surface gradient along column border which is intersected
!         by sphere.
!      minus -- negative correction to location of new node
!      nonode -- flag for when there is no intersection between search
!         sphere and slope line between outside and inside nodes or partial
!         column is truncated.
!         If nonode = 1, new node can not be calculated.
!         If nonode = 2, truncated column.
!      plus -- positive correction to location of new node
!      radsq -- square of search sphere radius
!      xcen,ycen,zcen -- center of search sphere relative to DEM origin
!      xya,xyb -- local x or y variables used to solve for new node
!      xycen,yxcen -- local variables used to solve for new node (x or y 
!         of search sphere center)
!      xy0,yx0 -- x or y location of new node
!      z0 -- z location of new node
!      za,zb -- local z variables used to solve for new node
!      zsphere -- top z elevation of sphere at new x,y node location
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            USE CommonData
            USE GridData, ONLY: radsq,xcen,ycen,zcen
            USE SetData, ONLY: setflag
            
            IMPLICIT NONE
            
            INTEGER, INTENT(out) :: nonode            

            REAL(pr), INTENT(in) :: xya,xyb,xy0,za,zb,xycen,yxcen
            REAL(pr), INTENT(out) :: dxy,z0
            REAL(pr) :: plus,minus,yx0,p1,p2,del
            REAL(pr) :: m,c,a,b,mtimesxya,bsqmin4atimesc

            nonode = 0
            m = (zb-za)/(xyb-xya)
            mtimesxya = m*xya                     
            c = ((za-zcen)-m*xya)*((za-zcen)-m*xya)+(xycen-xy0)*(xycen-xy0)&
                +yxcen*yxcen-radsq            
            a = 1.0_pr+m*m
            b = 2.0_pr*(m*(za-mtimesxya-zcen) - yxcen)
            p1 = -b/(2.0_pr*a)
            bsqmin4atimesc = b*b-4.0_pr*a*c
            IF ((bsqmin4atimesc).lt.0.0_pr) THEN
              nonode = 1
              RETURN
            END IF
            p2 = (1.0_pr/(2.0_pr*a)) * sqrt(bsqmin4atimesc)
            plus = p1+p2
            minus = p1-p2
!     Calculate del to use for taking into account precision errors. 
!     Count very close points as coinciding with the node point.
            del = .0001_pr*(xyb-xya)
            IF (plus.ge.xya-del.and.plus.le.xyb+del) THEN
              yx0 = plus
            ELSE
              IF (minus.ge.xya-del.and.minus.le.xyb+del) THEN
                yx0 = minus
              ELSE
                nonode = 1
                RETURN
              END IF
            END IF
            IF (yx0.lt.xya) yx0 = xya
            IF (yx0.gt.xyb) yx0 = xyb
            dxy = yx0 - xya
            z0 = m*dxy + za

!    If elevation of sphere center is below the elevation of the new node, mark this as a
!      truncated surface so it will not be calculated.                  
            IF (z0.gt.zcen) then
              nonode = 2  ! Set contains truncated column.
            END IF
            
            END SUBROUTINE getnode
  
