        SUBROUTINE Line2d (i,j,x1,x2,x3,x4,y1,y2,y3,y4,i2dcen,j2dcen,x2dcen,&
                           y2dcen,cosang,sinang,mline,outnode1,outnode2,&
                           dx1,dx2,dy1,dy2,line,zrad2d,tandip2d,&
                           x2d,y2d)
     
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!     This subroutine determines whether the 2-D slip line intersects the
!     column and calculates the length of the slip line.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!     Called by Ordinary and Bishop
!
!     VARIABLES
!
!      ang --  angle of slip converted to radians
!      angle -- assumed angle of slip (azimuth)
!      cosang -- cosine of slip direction angle
!      deltax -- width of column used to calculate length of intersecting
!        line using slope from slip angle.
!      delxy -- DEM grid resolution (delta x, delta y)
!      dx1,dx2,dy1,dy2 -- x and y grid spacing of new column relative to inside nodes
!      i2dcen,j2dcen -- i,j location of center of search sphere intersection 
!         with DEM surface
!      in(nx,ny) -- number of nodes of each column bounded by sphere
!      line -- length of 2-D slip line in current column
!      mline -- slope of slip direction line relative to x and y grid
!      outnode(2) -- array of partial column nodes which fall outside
!         slip surface
!      rnull -- real null value used throughout program 
!      sinang -- sin of slip direction angle
!      x(4) -- x location of each column corner node
!      x2dcen,y2dcen -- x,y location of center of search sphere intersection 
!         with DEM surface
!      xl,yl -- intersection of slip line with grid line for given j or i.
!      xl1,yl1 -- x,y intersection of slip line with one side of column
!      xl2,yl2 -- x,y intersection of slip line with another side of column
!      xlow,xhi -- min and max x of the two column intersections 
!      xy1,xy2 -- x-intersection of slip line with horizontal lines at j and j+1
!      y(4) -- y location of each column corner node
!      ylow,yhi -- min and max y of the two column intersections 
!
!
!              4---3
!              |   |    numbering scheme for corners of column
!              |   |
!              1---2
!
!                3
!               ---
!              |   |    numbering scheme for sides of column
!             4|   |2
!               ---
!                1
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      USE CommonData
      USE GridData, ONLY: halfdelxy,delxy,angle,xcen,ycen,zcen,radsq
      USE SetData, ONLY: in
            
      IMPLICIT NONE                 

      INTEGER, INTENT(in) :: i,j,i2dcen,j2dcen
      INTEGER, INTENT(in) :: outnode1,outnode2

      REAL(pr), INTENT(in) :: cosang,sinang,mline,x2dcen,y2dcen,dx1,dx2,dy1,dy2 
      REAL(pr), INTENT(inout) :: x1,x2,x3,x4,y1,y2,y3,y4
      REAL(pr), INTENT(out) :: line,zrad2d,tandip2d,x2d,y2d
      REAL(pr) :: deltax,xl1,xl2,yl1,yl2,xn1,xn2,xn3,xn4,yn1,yn2,yn3,yn4
      REAL(pr) :: xlow,xhi,ylow,yhi,zz,xrad2d,yrad2d
      REAL(pr) :: xy1,xy2,xl,yl,xdiff,ydiff,del

      CHARACTER*70 :: problemtype
      CHARACTER*120 :: errmessage,solution
      
      errmessage = ' '
      solution = ' '
      problemtype = 'calculating line for 2D factor of safety'   
            
      line = 0.0_pr
      xl1 = rnull
      xl2 = rnull
      yl1 = rnull
      yl2 = rnull
      xn1 = x1
      xn2 = x2
      xn3 = x3
      xn4 = x4
      yn1 = y1
      yn2 = y2
      yn3 = y3
      yn4 = y4
      del = .0001*delxy
   
!     If slip direction perpendicular to DEM grid and not in center column, return.
      IF ((sinang.gt.-0.00001_pr.and.sinang.lt.0.00001_pr.and.j.ne.j2dcen).or.&
          (cosang.gt.-0.00001_pr.and.cosang.lt.0.00001_pr.and.i.ne.i2dcen)) RETURN
            
!     If all four column nodes included   
      IF (in(i,j).eq.4) THEN
!     If slip direction is 0 or 180 degrees, calculate length of line.
        IF (sinang.gt.-0.00001_pr.and.sinang.lt.0.00001_pr) THEN
          line = delxy
          xl1 = xn1
          xl2 = xn2
          yl1 = y2dcen
          yl2 = y2dcen
        ELSE
!     If slip direction is 90 or 270 degrees
          IF (cosang.gt.-0.00001_pr.and.cosang.lt.0.00001_pr) THEN
            line = delxy
            xl1 = x2dcen
            xl2 = x2dcen
            yl1 = yn1
            yl2 = yn3
          ELSE  ! For all other slip directions
          
!     Find x-intersection of slip direction line and horizontal lines through j and j+1
            xy1 = (yn1-y2dcen)/mline + x2dcen
            xy2 = (yn3-y2dcen)/mline + x2dcen
            
!     Determine whether current column is intersected by the slip direction line, by
!     checking if any of the four sides the of column are intersected.
            IF ((xy1.lt.xn1).and.((xy2.gt.xn1).and.(xy2.lt.xn2))) THEN
!     Line intersects sides 3 and 4.
              deltax = xy2 - xn1
              line = deltax/cosang
              xl1 = xn1
              xl2 = xy2
              yl1 = (xn1-x2dcen)*mline + y2dcen
              yl2 = yn3 
            END IF
                 
            IF (((xy1.gt.xn1).and.(xy1.lt.xn2)).and.(xy2.gt.xn2)) THEN
!     Line intersects sides 1 and 2.
              deltax = xn2 - xy1
              line = deltax/cosang
              xl1 = xy1
              xl2 = xn2
              yl1 = yn1
              yl2 = (xn2-x2dcen)*mline + y2dcen
            END IF
             
            IF (((xy1.gt.xn1).and.(xy1.lt.xn2)).and.(xy2.lt.xn1)) THEN
!     Line intersects sides 1 and 4.
              deltax = xy1 - xn1
              line = deltax/cosang 
              xl1 = xn1
              xl2 = xy1
              yl1 = (xn1-x2dcen)*mline + y2dcen
              yl2 = yn1
            END IF
             
            IF ((xy1.gt.xn2).and.((xy2.gt.xn1).and.(xy2.lt.xn2))) THEN   
!     Line intersects sides 2 and 3.             
              deltax = xn2 - xy2
              line = deltax/cosang 
              xl1 = xy2
              xl2 = xn2
              yl1 = yn3
              yl2 = (xn2-x2dcen)*mline + y2dcen 
            END IF
            
            IF (((xy1.le.xn1).and.(xy2.ge.xn2)).or. &
               ((xy1.ge.xn2).and.(xy2.le.xn1))) THEN  
!     Line intersects sides 2 and 4.   
              line = delxy/cosang
              xl1 = xn1
              xl2 = xn2
              yl1 = (xn1-x2dcen)*mline + y2dcen
              yl2 = (xn2-x2dcen)*mline + y2dcen  
            END IF
             
            IF (((xy1.ge.xn1).and.(xy1.le.xn2)).and. &
               ((xy2.ge.xn1).and.(xy2.le.xn2))) THEN
!     Line intersects sides 1 and 3.  
              line = delxy/sinang
              xl1 = xy1
              xl2 = xy2
              yl1 = yn1
              yl2 = yn3   
            END IF
          END IF   
        END IF    

      ELSE  !  If partial column, check intersection with each side of quadrilateral.
!     If 2 nodes outside sphere assign x and y locations of new nodes.
        IF (in(i,j).eq.2) THEN
          IF (outnode1.eq.1.and.outnode2.eq.2) THEN
            yn1 = yn4-dy1
            yn2 = yn3-dy2
          END IF
          IF (outnode1.eq.2.and.outnode2.eq.3) THEN
            xn2 = xn1+dx1
            xn3 = xn4+dx2
          END IF
          IF (outnode1.eq.3.and.outnode2.eq.4) THEN
            yn3 = yn2+dy2
            yn4 = yn1+dy1
          END IF
          IF (outnode1.eq.1.and.outnode2.eq.4) THEN
            xn1 = xn2-dx1
            xn4 = xn3-dx2
          END IF
        END IF
!     If 1 node outside sphere assign x and y location of new node.
        IF (in(i,j).eq.3) THEN
          IF (outnode1.eq.1) THEN
            xn1 = xn2-dx1
            yn1 = yn4-dy1
          END IF
          IF (outnode1.eq.2) THEN
            xn2 = xn1+dx1
            yn2 = yn3-dy2
          END IF
          IF (outnode1.eq.3) THEN
            xn3 = xn4+dx2
            yn3 = yn2+dy2
          END IF
          IF (outnode1.eq.4) THEN
            xn4 = xn3-dx2
            yn4 = yn1+dy1
          END IF
        END IF   

!     If slip direction is not 0 or 180 degrees, check intersection with side 1.
        IF (sinang.le.-0.00001_pr.or.sinang.ge.0.00001_pr) THEN        
!     Find slip line intersection location with line defined by nodes 1 and 2.
          CALL Intersection (xn1,yn1,xn2,yn2,x2dcen,y2dcen,mline,cosang,xl,yl)        
          xlow = xn1
          xhi = xn2
!     Determine which node has bigger y coordinate.
          IF (yn1.gt.yn2) THEN
            ylow = yn2
            yhi = yn1
          ELSE
            ylow = yn1
            yhi = yn2
          END IF
!     Determine if slip line intersection falls within column.
          IF (((xl.ge.xlow-del).and.(xl.le.xhi+del)).and.((yl.ge.ylow-del).and.(yl.le.yhi+del))) THEN
!     Define first column side intersection as being in side 1.
            xl1 = xl
            yl1 = yl
          END IF
        END IF          

!    If slip direction is not 90 or 270 degrees check intersection with side 2.
        IF (cosang.le.-0.00001_pr.or.cosang.ge.0.00001_pr) THEN          
!    Find slip line intersection location with line defined by nodes 2 and 3.             
          CALL Intersection (xn2,yn2,xn3,yn3,x2dcen,y2dcen,mline,cosang,xl,yl)
          ylow = yn2
          yhi = yn3
!    Determine which node has bigger x coordinate.
          IF (xn2.gt.xn3) THEN
            xlow = xn3
            xhi = xn2
          ELSE
            xlow = xn2
            xhi = xn3
          END IF   
!     Determine if slip line intersection falls within column.     
          IF (((xl.ge.xlow-del).and.(xl.le.xhi+del)).and.((yl.ge.ylow-del).and.(yl.le.yhi+del))) THEN
!     If no column side intersections yet found, define first one in side 2.
            IF (xl1.lt.0.0_pr) THEN
              xl1 = xl
              yl1 = yl
            ELSE
!     If one intersection already found, second one is in side 2.
              xl2 = xl
              yl2 = yl
            END IF
          END IF
        END IF

!    If slip direction is not 0 or 180 degrees, check intersection with side 3.
        IF (xl2.lt.0.0_pr.and.sinang.le.-0.00001_pr.or.sinang.ge.0.00001_pr) THEN  
!    Find slip line intersection location with line defined by nodes 3 and 4.        
          CALL Intersection (xn3,yn3,xn4,yn4,x2dcen,y2dcen,mline,cosang,xl,yl)       
          xlow = xn4
          xhi = xn3
!    Determine which node has bigger y coordinate.
          IF (yn3.gt.yn4) THEN
            ylow = yn4
            yhi = yn3
          ELSE
            ylow = yn3
            yhi = yn4
          END IF        
          IF (((xl.ge.xlow-del).and.(xl.le.xhi+del)).and.((yl.ge.ylow-del).and.(yl.le.yhi+del))) THEN
!     If no column side intersections yet found, define first one in side 3.
            IF (xl1.lt.0.0_pr) THEN
              xl1 = xl
              yl1 = yl
            ELSE
!     If one intersection already found, second one is in side 3.
              xl2 = xl
              yl2 = yl
            END IF
          END IF
        END IF   

!    If slip direction is not 90 or 270 degrees check intersection with side 4.
        IF (xl2.lt.0.0_pr.and.cosang.le.-0.00001_pr.or.cosang.ge.0.00001_pr) THEN     
!    Find slip line intersection location with line defined by nodes 4 and 1.         
          CALL Intersection (xn4,yn4,xn1,yn1,x2dcen,y2dcen,mline,cosang,xl,yl)  
          ylow = yn1
          yhi = yn4
!    Determine which node has bigger x coordinate.
          IF (xn1.gt.xn4) THEN
            xlow = xn4
            xhi = xn1
          ELSE
            xlow = xn1
            xhi = xn4
          END IF        
          IF (((xl.ge.xlow-del).and.(xl.le.xhi+del)).and.((yl.ge.ylow-del).and.(yl.le.yhi+del))) THEN
!     If no column side intersections yet found, then this column is not intersected.
            IF (xl1.lt.0.0_pr) THEN
              RETURN
            ELSE
!     Second column intersection is in side 4.
              xl2 = xl
              yl2 = yl
            END IF
          END IF 
        END IF
      END IF
        
      IF (xl2.ge.0.0_pr.or.yl2.ge.0.0_pr) THEN
!    Determine x and y location of center of 2d column.
          xdiff = xl2-xl1
          ydiff = yl2-yl1
          x2d = xl1 + xdiff * 0.5_pr
          y2d = yl1 + ydiff * 0.5_pr  
!    Determine elevation of slip base at center of 2d column.        
          xrad2d = x2d-xcen
          yrad2d = y2d-ycen 
          zz = radsq - xrad2d*xrad2d - yrad2d*yrad2d
          IF (zz.lt.0.0_pr) THEN
            errmessage = 'cannot calculate column midpoint'
            CLOSE (33)
            Call WriteError(1,errmessage,problemtype,'no','no ',0,' ')               
!            PRINT *,radsq,x2d,y2d,zz,xl1,xl2,yl1,yl2
          END IF
          zrad2d = sqrt(zz)          
!     Use partial derivative of equation for the 2D search sphere 
!     at column center to compute true dip of base.
          tandip2d = -(xrad2d*cosang+yrad2d*sinang)/zrad2d
  
!     Find length of line using pythagorean theorem.
          IF (abs(line).le.0.00001_pr) line = sqrt(xdiff*xdiff + ydiff*ydiff)                      
      END IF     

      END SUBROUTINE Line2d        


      SUBROUTINE Intersection (xx1,yy1,xx2,yy2,x2dcen,y2dcen,mline,cosang,xi,yi)
     
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!      This subroutine calculates the intersection of the slip direction
!      line with the line segment defined by the column corner nodes xx1,yy1 .
!      and xx2,yy2. The intersection is (x,y).
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!      Called by Line2d
!
!      VARIABLES
!
!      det -- determinant
!      mline -- slope of slip line relative to positive x of DEM grid
!      x2dcen,y2dcen -- location of center of search sphere intersection 
!            with DEM surface
!      xx1,yy1,xx2,yy2 -- coordinates of grid line along one side of column
!      x,y -- intersection of slip line with grid line
!
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        USE CommonData
        
        IMPLICIT NONE

        REAL(pr), INTENT(in) :: xx1,yy1,xx2,yy2,mline,x2dcen,y2dcen,cosang
        REAL(pr), INTENT(out) :: xi,yi
        REAL(pr) :: a(2,2),b(2,2),a1,a2,b1,b2,c1,c2,det

        CHARACTER*70 :: problemtype
        CHARACTER*100 :: errmessage,solution
      
        errmessage = ' '
        solution = ' '
        problemtype = 'calculating line for 2D factor of safety calculation'        
                
!     If slip direction is 90 or 270 degrees
        IF (cosang.gt.-0.00001_pr.and.cosang.lt.0.00001_pr) THEN
          b1 = yy1-yy2
          a1 = xx2-xx1
          c1 = (xx1*yy2) - (xx2*yy1)
          a2 = 0.0_pr
          b2 = 1.0_pr
          c2 = -x2dcen
        ELSE
          a1 = yy2-yy1
          b1 = xx1-xx2
          c1 = (xx2*yy1) - (xx1*yy2)
          a2 = -mline
          b2 = 1.0_pr
          c2 = mline*x2dcen - y2dcen
        END IF

        a(1,1) = a1
        a(1,2) = b1
        a(2,1) = a2
        a(2,2) = b2        

!  Compute the determinant.

        det = a(1,1) * a(2,2) - a(1,2) * a(2,1)

        IF (det.eq.0.0_pr) THEN
          errmessage = 'error in determining line intersections'
          CLOSE (33)
          Call WriteError(1,errmessage,problemtype,'no','no ',0,' ')         
        END IF

!  Compute the entries of the inverse matrix using an explicit formula.

        b(1,1) = + a(2,2) / det
        b(1,2) = - a(1,2) / det
        b(2,1) = - a(2,1) / det
        b(2,2) = + a(1,1) / det
 
!  If the inverse exists, then multiply the inverse times -C to get the intersection point.
        IF (cosang.gt.-0.00001_pr.and.cosang.lt.0.00001_pr) THEN
          yi = - b(1,1) * c1 - b(1,2) * c2
          xi = - b(2,1) * c1 - b(2,2) * c2
        ELSE
          xi = - b(1,1) * c1 - b(1,2) * c2
          yi = - b(2,1) * c1 - b(2,2) * c2
        END IF
                        
      END SUBROUTINE Intersection
        
        
        
        

