        SUBROUTINE Radius(zd,newrad,va1max,radiniprev,rmax,ierr,goodradcnt,badradcnt)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!       This subroutine computes the radius of the slip surface for 
!       the first sphere at the search grid point, adjusting the 
!       radius on subsequent calls to meet the volume or area constraint.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!        Called by Scoops3D
!
!        VARIABLES
!
!        a0,a1,a2 -- parameters used to calculate radius
!        a2sq,a2cube -- parameters used to calculate radius
!        armin,armax -- minimum and maximum area bounds for valid slip surface
!        ddxcensq,ddycensq,ddzcensq -- square of center of DEM column to 
!            search center
!        delxy -- DEM grid resolution (delta x, delta y)
!        diff -- parameter used to calculate radius
!        dr -- differential of r found from derivative of of
!           volume with respect to r.
!        dt -- difference between vmin or armin and set volume
!        dtold -- previous dt
!        dxcen -- x center of DEM grid
!        dycen -- y center of DEM grid
!        goal -- vmin or armin + half of tol
!        ia -- flag for area control of scoop size
!        ierr -- error flag indicating initial radius problem (2) 
!           or too many subsets (1)
!        iv -- flag for volume control of scoop size
!        ncount -- total number of changes in radius estimates
!        newrad -- indicates whether valid initial range subset has been found.
!           (1=no, 0=yes, 2=no valid sets after 10 radius adjustments.)
!        nx -- number of DEM cells in x direction
!        ny -- number of DEM cells in y direction
!        oldzcen -- last zcen, used to calculate new radius when searching new k at same i,j
!        r -- parameter used to calculate radius
!        rnull -- default value for DEM grid point outside boundary
!        rad -- radius of search sphere
!        rad1 -- complex solution for radius
!        radini - initial radius
!        radiniprev -- used as first radius estimate at new search lattice elevation.
!           Equals previous radius plus elevation change, or = -1 if radius not yet
!           calculated at current search lattice x-y location.
!        radold -- last radius guess with non-zero volume
!        rnull -- real null value used throughout program 
!        q -- parameter used to calculate radius
!        s1,s2 -- complex parameters used to calculate radius
!        t1 -- complex parameter used to calculate radius
!        tol -- tolerance on primary control criterion for calculating initial 
!           radius.  (i.e. initial volume must fall between volume and volume+tol)
!        vacriterion -- indicates primary control for intial radius (v or a
!           for volume or area)
!        va1max -- max. volume or area (primary criterion) of a subset
!        vmin,vmax -- minimum and maximum volume bounds for valid slip surface
!        xcen,ycen,zcen -- center of search sphere relative to DEM origin
!        z -- distance used to calculate new radius.
!        zavgmax -- average distance of surface of columns in subset to sphere center
!        zd -- elevation of DEM grid at same x,y as sphere center
!        zmin -- minimum DEM elevation
!        zzsrchres -- resolution of search grid on z-axis
!       
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        USE CommonData
        USE GridData, ONLY: xcen,ycen,zcen,rad,delxy,zmin,nx,ny,pi
        USE SearchData, ONLY: zavgmax,vmin,vmax,armin,armax,tol,&
                              vacriterion,goal
        
        IMPLICIT NONE
  
        INTEGER, INTENT(inout) :: newrad,ierr
        INTEGER :: ncount,intersect
        INTEGER :: goodradcnt,badradcnt
        
        REAL(pr), INTENT(inout) :: radiniprev
        REAL(pr), INTENT(in) :: va1max,zd,rmax
        REAL(pr) :: dt,dr,radini,ddzcensq,zdiff
        REAL(pr) :: radold,oldzcen
        REAL(pr) :: dxcen,dycen,z,ddxcensq,ddycensq
        REAL(pr) :: a0,a1,a2,q,r,diff,a2sq,a2cube

        COMPLEX(pr) :: s1,s2,t1,rad1
                
        SAVE radold,ncount
        SAVE radini,z,oldzcen,intersect
        dt = 0.0_pr
     
!     If this is first radius estimate at current search lattice node (i,j,k)
        IF (rad.lt.0.0_pr) THEN        
!     If radius has never been estimated at this search grid i,j location       
          IF (radiniprev.lt.0.0_pr) THEN          
!     If search node is outside DEM, estimate the initial radius as the 
!     distance to the center of the DEM grid with minimum DEM elevation, zmin.
            IF (zd.eq.-888.0_pr) THEN
              dxcen = REAL(((nx+1)/2)- 1,pr) * delxy
              dycen = REAL(((ny+1)/2)-1,pr) * delxy
              ddxcensq = (dxcen-xcen)*(dxcen-xcen)
              ddycensq = (dycen-ycen)*(dycen-ycen)
              ddzcensq = (zmin-zcen)*(zmin-zcen)
              rad = SQRT(ddxcensq + ddycensq + ddzcensq)
            ELSE
!     If search node is above null DEM value, estimate z as distance
!     from search lattice elevation to minimum DEM elevation. Otherwise use 
!     distance to DEM elevation directly below search center.
              IF (zd.eq.rnull) THEN
                z = zcen - zmin
              ELSE
                z = zcen - zd
              END IF            
!     If search grid point is within DEM boundary and volume
!     is primary criterion calculate initial guess for radius assuming a flat
!     DEM surface at distance z from the search node and the volume equation
!     for a truncated sphere.
              IF (vacriterion.eq.'v') THEN         
                a0 = .50_pr * (z*z*z) - (3.0_pr*vmin)/(2.0_pr*pi)
                a1 = 0.0_pr
                a2 = (-1.50_pr)*z
                a2sq = a2*a2
                a2cube = a2*a2*a2
                q = (1.0_pr/3.0_pr)*a1 - (1.0_pr/9.0_pr)*a2sq
                r = (1.0_pr/6.0_pr)*(a1*a2 - 3.0_pr*a0) - (1.0_pr/27.0_pr)*a2cube
                t1 = (q*q*q + r*r)
                t1 = SQRT(t1)
                s1 = (r + t1)**(1.0_pr/3.0_pr)
                diff = r-t1
                IF (ABS(diff).gt.1.E-2) THEN
                  s2 = (r - t1)**(1.0_pr/3.0_pr)
                ELSE
                  s2 = 0.00_pr
                END IF
                rad1 = (s1+s2) - (a2/3.0_pr)
                rad = REAL(rad1)   
              ELSE   ! if area is primary criterion
                rad = SQRT(armin/pi + z*z)
              END IF
            END IF  ! if (zd.eq.-888.)
!     If radius was previously calculated at this i,j (but not this search node elevation), 
!     add change in search node elevation for new radius initial estimate.            
          ELSE 
            zdiff = zcen - oldzcen
            IF (vacriterion.eq.'v') THEN          
              dr = 0.5_pr * zdiff *(radiniprev + z) / radiniprev
            ELSE
              dr = z * zdiff/ radiniprev
            END IF
            rad = radiniprev + dr
          END IF  ! if (radiniprev,eq,-1)
            
          radini = rad
          radold = rad
          ncount = 0
          intersect = 0
        
!     If not first estimate of radius adjust value to meet criterion.
        ELSE
          ncount = ncount + 1
!     If grid not intersected, increment radius.
          IF (va1max.eq.0.0_pr) THEN
            IF (intersect.eq.0) THEN  ! If never intersected
              radold = rad
              rad = radini + 1.50_pr**(ncount+1)*delxy
              IF (rad.gt.rmax) rad = rmax - delxy 
            ELSE  ! If intersected at least once, but volume now zero
              rad = (rad + radold)/2.0_pr
            END IF
          ELSE
            intersect = 1
!     Calculate difference from target volume or area depending on primary criterion.
!     (goal=vmin or armin + half of tol)
            dt = va1max - goal  

!     Otherwise find new radius using dv/dr or da/dr, depending on primary criterion.              
            IF (vacriterion.eq.'v') THEN
              dr = -dt/(2.0_pr* pi * rad * (rad-zavgmax))
            ELSE
              dr = -dt/(2.0_pr* pi * rad)
            END IF 
            z = zavgmax                        
            radold = rad
            rad = rad + dr
!     If the new radius would be negative, use half of last radius.
            IF (rad.le.0.0_pr) rad = 0.5_pr * radold    
          END IF      
          
!     If can't hit goal criterion range and volume is non-zero
!     choose largest value of radius and continue by setting newrad
!     to 2, so won't get sent back into Radius.f to recalculate.          
          IF (ncount.gt.10) THEN
            badradcnt = badradcnt + 1
            IF (dt.gt.0.0_pr) THEN
              rad = MIN(rad,radold)
            ELSE
              rad = MAX(rad,radold)
            END IF
            newrad = 2
            ierr=1
          ELSE
            goodradcnt = goodradcnt + ncount 
          END IF
        END IF  !  if (rad.eq.-1)
              
        oldzcen = zcen 
        IF (ierr.ne.1) THEN
          radiniprev = rad
        ELSE
          radiniprev = -1
        END IF
            

        END SUBROUTINE radius
            
