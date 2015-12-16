        SUBROUTINE Fos (nset,iretro,fgrid)      

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!       This subroutine calls subroutine Ordinary or Bishop, depending on
!       solution method chosen by user, to calculate the factor of safety
!       for each scoop and tracks various minimum factors of safety depending
!       on options chosen.  
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!       VARIABLES
!
!        arclength -- arc length of current 2d slice
!        area(nset) -- surface area of slide
!        area2d -- cross-sectional area of current 2d slice
!        angle -- assumed angle of slip (azimuth)
!        colfile -- flag for opening error file reporting cols < limcol
!        deginc -- degree increment to search for minimum FOS
!        degmax -- degrees to search on either side of estimated
!           fall line (describes half width of search segment)
!        delr -- radius increment for search
!        delxy -- DEM grid resolution (delta x, delta y)
!        delz -- z resolution of search grid
!        failsurf(nx,ny) -- array of failure surface elevations for each DEM column
!        fall -- estimated fall line of scoop based on DEM slopes
!        fallfailsurf -- estimated fall line of scoop based on failure surface slope
!        fangle -- angle of minimum factor of safety for a subset
!        fangle2d -- angle of minimum FOS 2d slice for a subset
!        farclength -- arc length of minimum FOS 2d slice for a subset
!        farea2d -- area of minimum FOS 2d slice for a subset
!        fgrid -- minimum factor of safety at search grid point
!        filtcount(nx,ny) -- counter of number of times each column was filtered by
!           minma
!        fminma -- min m-alpha at current minimum fos slip angle
!        foscut -- factor of safety cutoff for removing slip surface from DEM
!        fosmin -- minimum factor of safety for a subset at DEM column
!        fosmin2d -- 2-d minimum factor of safety for a subset at DEM column
!        fosmin2d3d -- 3-d minimum factor of safety associated with min FOS 2-d slice
!           for a subset at DEM column
!        fosmin3d2d -- 2-d minimum factor of safety associated with min FOS 3-d scoop
!           for a subset at DEM column
!        foso2d -- 2-d factor of safety at trial angle
!        fsangle(i,j) -- slip angle of surface associated with min FOS
!        fsangle2d(nx,ny) -- array of fall angles of minimum 2d slice at each DEM column
!        fsarclength(nx,ny) -- array of arclengths of minimum 2d slice at each DEM column
!        fsarea2d(nx,ny) -- array of areas of minimum 2d slice at each DEM column
!        fsmin(i,j) -- minimum factor of safety at each DEM column
!        fsmin2d(i,j) -- 2-d minimum factor of safety at each DEM column
!        fsmin2d3d(nx,ny) -- 3-d minimum factor of safety associated with min 2d slice
!           at each DEM column
!        fsmin3d2d(nx,ny) -- 2-d minimum factor of safety associated with min 3d scoop
!           at each DEM column
!        fsrad(i,j) -- radius of sphere associated with min FOS at each DEM column
!        fsrad2d(nx,ny) -- array of radii of minimum 2d FOS slice at each DEM column
!        fsvol(i,j) -- volume of minimum fos scoop for DEM column i,j
!        fswidth2d(nx,ny) -- array of widths of minimum 2d slice at each DEM column
!        fsx(nx,ny),fsy(nx,ny),fsz(nx,ny) -- search sphere center of minimum FOS 
!           scoop at each DEM column
!        fsx2d(nx,ny),fsy2d(nx,ny),fsz2d(nx,ny) -- search sphere center of minimum FOS 
!           2d slice at each DEM column
!        fwidth2d -- width of minimum FOS 2d slice for a subset
!        i1,i2,j1,j1 -- Array bounds for current set
!        icrit(nx,ny) -- array of volume or area flags indicating proximity to limits
!            of criteria range for each DEM column
!        icritprox -- volume or area boundary proximity flag (i.e. 1 if vol<vmin+tol, 2 if
!           vol>vmax-tol)
!        ifailsurf -- flag for whether failure surface file is used
!        ifallfail -- flag for whether searchable fall line was found based on failure surface dip
!        ifos2d -- flag for calculating 2-d factors of safety
!        ilattice -- flag for creating file foslattice.out for minimum FOS at each search
!           grid node, 1=create file in ijk format, 2=create file in xyz&
!           format
!        i2dcolcen,j2dcolcen -- i,j location of center of search sphere intersection 
!           with DEM surface
!        isubsurf -- flag for creating 3-d file of FOS below DEM surface&
!          1=create file in ijk format, 2=create file in xyz format
!        in(nx,ny) -- number of nodes of each column bounded by sphere or
!          equals -1 if truncated surface.
!        insphere(nx,ny,4) -- indicates whether DEM column node is bounded by sphere;
!          nodes are numbered counterclockwise from lower left corner
!        iretro -- counter for number of retrogressions
!        iseg -- integer associated with current search angle
!        limcol -- optimal number of columns per slip surface.  Fewer columns
!           will generate error in problems.out file
!        linein(nx,ny) -- array of flags for whether 2-D slide directon line intersects
!           column
!        m -- set number
!        mcol(nx,ny) -- number of columns in minimum FOS scoop at each DEM column
!        mcol2d(nx,ny) -- number of columns in minimum 2d FOS slice at each DEM column
!        method -- B for Bishops Simplified, O for Ordinary (Fellenius) method for
!           calculation of factor of safety.
!        mflag -- flag indicating whether filter bounds were exceeded
!        minma -- min absolute value of m-alpha of all iterations for a surface
!        mini(nnset),maxi(nnset) -- minimum and maximum i bounds of each subset
!        minj(nnset),maxj(nnset) -- minimum and maximum j bounds of each subset
!        ncol(nset) -- number of columns in current scoop
!        ncol2d -- number of columns in current 2-D slice
!        nrange(nnset) -- array of flags for whether set fits volume or area criteria 
!           range
!        nretro -- number of failure retrogressions to compute
!        numdir -- number of search slip directions for each failure
!        nset -- number of subsets found
!        ntry -- total number of slip surfaces for which FOS is calculated
!        nx -- number of DEM cells in x direction
!        ny -- number of DEM cells in y direction
!        nz -- number of nodes used for 3-D FOS values when isubsurf=1or2
!        oangle -- slip angle associated with overall minimum factor of safety
!        oangle2d -- slip angle associated with overall minimum 2d factor of safety
!        oarclength -- arc length of overall minimum 2d factor of safety
!        oarclength3d2d -- arc length of 2d surface associated with overall min 3d scoop
!        oarea -- area of overall minimum scoop
!        oarea2d -- area of overall minimum 2d slice
!        oarea2d3d -- area of 3d surface associated with overall min 2d slice
!        oarea3d2d -- area of 2d surface associated with overall min 3d scoop
!        ocol -- number of columns in minimum scoop
!        ocol2d -- number of columns in overall minimum 2-D scoop
!        ocol2d3d -- number of columns in 3d surface associated with min 2d slice
!        ocol3d2d -- number of columns in 2d surface associated with min 3d slice
!        ofos -- overall minimum factor of safety
!        ofos2d -- 2-d FOS associated with overall minimum 3-d factor of safety
!        ofos2d3d -- FOS of 3d slice associated with overall minimum 2d FOS
!        ofos3d2d -- FOS of 2d slice associated with overall minimum 3d FOS
!        ofsx,ofsy,ofsz -- search sphere center of minimum FOS scoop
!        ofsx2d,ofsy2d,ofsz2d -- search sphere center of minimum 2d FOS scoop
!        omset -- set number of overall minimum scoop
!        orad -- radius of sphere associated with overall minimum factor of safety
!        orad2d -- radius of arc associated with overall minimum 2d factor of safety
!        osliparea -- slip surface area of overall min 3d scoop
!        osliparea2d3d -- slip surface area of 3d slice associated with overall minimum 2d FOS
!        ovangle -- slip angle of largest volume scoop with FOS < cutoff
!        ovarea -- area of largest volume scoop with FOS < cutoff
!        ovfos -- factor of safety of largest volume scoop with FOS < cutoff
!        ovmset -- set number of largest volume scoop with FOS < cutoff
!        ovol, oarea - volume and area of minimum FOS scoop
!        ovol2d3d -- volume of 3d surface associated with overall min 2d slice
!        ovrad -- radius of largest volume scoop with FOS < cutoff
!        ovfsx,ovfsy,ovfsz -- search sphere center of largest volume scoop with 
!           FOS < cutoff
!        ovsliparea -- slip surface area of largest volume scoop with 
!           FOS < cutoff
!        ovvol,ovarea -- volume and area of largest volume scoop with FOS < cutoff
!        ovwt -- weight of largest volume scoop with FOS < cutoff
!        owidth2d -- wisth of overall min 2d FOS
!        owidth3d2d -- width of 2d surface associated with overall min 3d scoop
!        owt -- weight of overall min 3d scoop
!        owt2d3d -- weight of 3d slice associated with overall minimum 2d FOS
!        ozb(nx,ny) -- array of DEM after scoops < FOS cut off removed.
!        remove -- A= all scoops < FOS cutoff will be removed from new
!           DEM.  L= only largest volume scoop < FOS cutoff will be 
!           removed. Other= the scoop with the minimum FOS will be removed.
!        rad -- radius of search sphere
!        rnull -- real null value used throughout program 
!        sdeg -- starting angle
!        setflag(nnset) -- indicates valid subset (0), or invalid by containing
!           truncated node (1), or adjacent to DEM boundary (2), or empty (-1)
!        single -- flag for calculating single slip surface
!        sliparea -- area of slip surface
!        subset(nx,ny) -- indicates set membership of DEM columns
!        tol -- tolerance on primary control criterion for calculating initial 
!           radius.  (i.e. initial volume must fall between volume and volume+tol)
!        vacriterion -- indicates primary control for intial radius (v or a
!           for volume or area)
!        vaflag -- flag=1 if exceed maximum volume/area allowed,-1 if under volume/area
!        vmin,vmax -- minimum and maximum volume bounds for valid slip surface
!        volume(nset) -- volume of scoop
!        weight -- weight of slide
!        width2d -- width of current 2d slice
!        xcolcount -- number of critical surfaces with total columns < limcol
!        xcount --  number of surfaces eliminated due to m-alpha limits or nonconvergence
!        xcen,ycen,zcen -- center of search sphere relative to DEM origin
!        x2dcen,y2dcen -- perpendicular projection of x,y location of center of search sphere  
!           intersection with DEM surface onto axis of rotation
!        x2dcolcen,y2dcolcen -- x,y location of center of search sphere intersection 
!           with DEM surface
!        xfailslope(nx,ny),yfailslope(nx,ny) -- x and y slope vectors of input failure surface
!        xslope(nx,ny),yslope(nx,ny) -- x and y slope vectors of DEM surface
!        xsum,ysum -- sum of DEM slope vectors of full columns included in slip surface
!        xfailsum,yfailsum -- sum of failure surface slope vectors of full columns
!        zmid(i,j) -- elevation of slip surface at center of intersected column
!        zbot -- elevation of bottom of 3-d fos grid array
!        zdem(i,j) -- DEM elevations
!        zfos(nx,ny,nz) -- 3-d array of minimum factor of safety below DEM surface
!        ztop -- z elev 1 node above DEM surface
!
!        OUTPUT FILES USED
!        unit filename
!          32 'fosretro_out.txt' -- file of retrogression data. Opened in Readin
!               and written in Scoops, Readin and Fos if Remove=A and nretro>0.
!          33 'filter_out.txt' -- file of data describing surfaces which fell outside
!               of user-defined m-alpha and nonconverging surfaces.  
!               Opened and written in FOS if filter = 1.
!          34 'ncolerr_out.txt' -- file of information about critical surfaces with
!               fewer columns than user-defined limit, limcol. Opened and
!               written in Fos if colfile = 1.
!          40 'spheresltcut_out.txt' --  file of sphere center coordinates, radius, 
!                fall angle, # of columns in scoop, vol, and factor of safety of the
!                 scoops with F<fostol. Opened in readin and written in fos, only if remove='A'.
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        USE GridData, ONLY: xcen,ycen,zcen,zdem,angle,delz,rad,radsq,nx,ny,nz,delxy,pi,&
                           lengthunits,xll,yll,zcenrot,xcenrot,ycenrot,halfdelxy
        USE SetData, ONLY: in,xslope,yslope,area,ncol,nrange,mini,minj,maxi,maxj,&
                           subset,setflag,colfile,xcolcount,zmid,&
                           sliparea,weight,volume,insphere
        USE SearchData, ONLY : vmin,tol,vmax,armin,armax,vacriterion,irotcen
        USE FOSData
        USE FailSurfData, ONLY: failsurf,xfailslope,yfailslope,ifailsurf
                    
        IMPLICIT NONE

        INTEGER, INTENT(in) :: nset,iretro
        INTEGER :: i,j,m,k,i1,i2,j1,j2,i2dcolcen,j2dcolcen
        INTEGER :: iseg,icritprox,mflag,ifallfail
     
        REAL(pr), INTENT(inout) :: fgrid
        REAL(pr) :: fall,fosmin,ztop,z,fminma,x2dcen,y2dcen,xdiff,ydiff
        REAL(pr) :: xsum,ysum,fangle,sdeg,farea2d,farclength2d,fncol2d
        REAL(pr) :: fosmin2d,fosmin2d3d,fosmin3d2d,fangle2d,fwidth2d
        REAL(pr) :: fwidth3d2d,farea3d2d,farclength3d2d,fncol3d2d,rad2d
        REAL(pr) :: fallfailsurf,xfailsum,yfailsum,x2dcolcen,y2dcolcen,yprime      
        REAL(pr) :: felfosmin,cosang,sinang,ang
       
        CHARACTER*10 :: cnum
        CHARACTER*4 :: fosmeth
        CHARACTER*1 :: charetro
        
!     Initialize variables
        x2dcen = 0.0_pr
        y2dcen = 0.0_pr 
        i2dcolcen = 0
        j2dcolcen = 0
        rad2d =  0.0_pr
        mflag = 0
        cnum = '0123456789'
        charetro = cnum(iretro+1:iretro+1)
        sdeg = 365.0_pr
        IF (irotcen.ne.1) THEN
          xcenrot = xcen
          ycenrot = ycen
          zcenrot = zcen
        END IF
        
!!!  Fellenius comparison addition
        felfoso = 0.0_pr 
!!!
        
!     Loop through all sets in current slip surface
        DO m = 1,nset
          fosmin = nullhi
          foso = nullhi 
!     If not valid subset, go on to next one.
          IF (nrange(m).eq.0.and.single.ne.1) CYCLE              
          xsum = 0.0_pr
          ysum = 0.0_pr  
          IF (ifailsurf.eq.1) THEN        
            xfailsum = 0.0_pr
            yfailsum = 0.0_pr 
            ifallfail = 0      
          END IF 
          
!    Set loop range as minimum and maximum locations of columns
!    in set m.              
          j1 = minj(m)
          j2 = maxj(m)
          i1 = mini(m)
          i2 = maxi(m)               
          
!     If single run determine whether slip direction is specified.
!     If angle is greater than 360 then Scoops will calculate fall line.
          IF (single.eq.1) THEN
            IF (setflag(m).ne.0) CYCLE
            numdir = 1
            IF (abs(angle).le.360.0_pr) THEN 
!     Use input angle as slip direction
              sdeg = angle
            ELSE ! use fall line as slip direction
              degmax = 0.0_pr
              deginc = 0.0_pr
            END IF
          END IF

!      Calculate fall line if not given as part of single search.
          IF (sdeg.gt.360.0_pr) THEN         
            DO j = j1,j2
              DO i = i1,i2           
!     Sum DEM surface slope vectors to find resultant slope 
!     vector using set columns that are fully intersected.
                IF (subset(i,j).eq.m) THEN
                  IF (in(i,j).eq.4) THEN
                    xsum = xsum + xslope(i,j)
                    ysum = ysum + yslope(i,j)
                    IF (ifailsurf.eq.1) THEN
!     If using defined failure surface, sum failure surface slope vectors to 
!      find resultant slope vector using set columns that are fully intersected.
                      IF (zmid(i,j).lt.failsurf(i,j)) THEN
                        xfailsum = xfailsum + xfailslope(i,j) 
                        yfailsum = yfailsum + yfailslope(i,j)                       
                      END IF
                    END IF
                  END IF
                END IF  
              END DO
            END DO
!    If surface is flat, move to next set.
            IF (ysum.eq.0.0_pr.and.xsum.eq.0.0_pr) CYCLE

!     Estimate fall line (convert to degrees) and angle search range.
            fall = ATAN2(ysum,xsum)
            fall = fall*(180.0_pr/pi)
            IF (fall.lt.0.0_pr) fall = fall + 360.0_pr
            sdeg = fall - degmax
          END IF
	        
!     Estimate fall line of defined failure surface.
          IF (ifailsurf.eq.1) THEN 
            IF (yfailsum.ne.0.0_pr.or.xfailsum.ne.0.0_pr) THEN
              fallfailsurf = ATAN2(yfailsum,xfailsum)
              fallfailsurf = fallfailsurf*(180.0_pr/pi)
              IF (fallfailsurf.lt.0.0_pr) fallfailsurf = fallfailsurf + 360.0_pr
!      Determine if failure surface fall line is within +/- 90 deg. of DEM fall line
!      and use as second search direction if so.
              IF (fallfailsurf.gt.fall-90.or.fallfailsurf.lt.fall+90) THEN
                IF (fallfailsurf-fall.gt.0.1_pr) THEN
                  ifallfail = 1
                  numdir = numdir + 1
                END IF
              END IF
            END IF
          END IF

!     To calculate 2D center of rotation find approximate center of all intersected columns using
!     max and min column locations. Add 1 to find center of min and max nodes.	        
          IF (ifos2d.eq.1) THEN 
            fosmin2d = nullhi 
            i2dcolcen = INT(REAL(i1+i2+1)*0.5_pr+0.001_pr*delxy)
            j2dcolcen = INT(REAL(j1+j2+1)*0.5_pr+0.001_pr*delxy) 
            x2dcolcen = (REAL(i1+i2-1)*0.5_pr) * delxy 
            y2dcolcen = (REAL(j1+j2-1)*0.5_pr) * delxy 
            xdiff = x2dcolcen - xcenrot
            ydiff = y2dcolcen - ycenrot
          END IF
	        
!     Find minimum factors of safety for each slip direction.
          DO iseg = 1,numdir
            IF (iseg.eq.numdir.and.ifailsurf.eq.1.and.ifallfail.eq.1) THEN
              angle = fallfailsurf
            ELSE
              angle = sdeg + REAL(iseg-1,pr)*deginc
            END IF
            ang = angle*pi/180.0_pr  
            cosang = COS(ang)
            sinang = SIN(ang)
            IF (ifos2d.eq.1) THEN 
!     Project center of intersected columns onto axis of rotation
!     to estimate 2D center of rotation. This center will be used
!     for 2-D FOS solutions such that all 2-D fall lines will
!     pass through this point.
              yprime = xdiff*(-sinang)+ydiff*cosang
              x2dcen = -yprime*sinang + xcenrot
              y2dcen = yprime*cosang + ycenrot
              rad2d = sqrt(radsq - yprime*yprime)
            END IF
 
            SELECT CASE (method)
              CASE ('F','f','O','o')
                CALL Ordinary(m,i1,i2,j1,j2,ang,cosang,sinang,&
                               i2dcolcen,j2dcolcen,x2dcen,y2dcen,rad2d)
                
              CASE ('B','b') 
                CALL Bishop(m,i1,i2,j1,j2,ang,cosang,sinang,&
                            i2dcolcen,j2dcolcen,x2dcen,y2dcen,rad2d,mflag)                
            END SELECT

!     Count how many FOS trials.            
            ntry = ntry + 1
            
!     If surface fell outside of m-alpha or FOS limits or did not converge 
!     in Bishop write scoop parameters to filter.out file.
            IF (mflag.ne.0) THEN
!     If filter file not yet opened, open it and write header lines. Initialize 
!     count of filtered search nodes.
              IF (mfile.eq.0) THEN
                IF (iretro.gt.0) THEN 
                  OPEN (33,STATUS='replace',file = outputdir(1:LEN_TRIM(outputdir))//&                 
                     filin(bfile:nfile)//'_filter'//charetro//'_out.txt')
                ELSE
                  OPEN (33,STATUS='replace',file = outputdir(1:LEN_TRIM(outputdir))//&
                     filin(bfile:nfile)//'_filter_out.txt') 
                END IF
!                IF (method.eq.'b') THEN
                  WRITE (33,1500)
                  WRITE (33,1510)
                  IF (absminma.ne.0) WRITE (33,1520) absminma 
                  WRITE (33,1540)
                  IF (LEN_TRIM(lengthunits).eq.1) THEN
                    WRITE (33,1551) lengthunits,lengthunits,lengthunits,lengthunits
                  ELSE
                    IF (LEN_TRIM(lengthunits).eq.2) THEN
                       WRITE (33,1552) lengthunits,lengthunits,lengthunits,lengthunits
                    ELSE 
                      WRITE (33,1550)
                    END IF  
                 END IF                   
                mfile = 1
                xcount = 0
              END IF
!    Write filtered search center data and calculated FOS to filter file.
              IF (method.eq.'b') THEN
                WRITE (33,1560) xcen+xll,ycen+yll,zcen,rad,angle,minma,foso,felfoso
              ELSE
                WRITE (33,1620) xcen+xll,ycen+yll,zcen,rad,angle,foso
              END IF
!    Set calculated FOS to null value so they will not be compared to
!    other factors of safety.
              IF (foso.ne.100.0_pr) foso = 111.0_pr
              IF (foso2d.ne.100.0_pr) foso2d = 111.0_pr
              xcount = xcount + 1
            END IF  !  IF (mflag.ne.0)

!     If removing all surfaces with FOS<foscut, save bottom elevation to ozb.
            IF (remove.eq.'A') THEN
              IF (foso.lt.foscut) THEN
               WRITE (40,1400) xcen+xll, ycen+yll,zcen,rad,angle,volume(m),area(m),foso
                IF (nretro.gt.0) WRITE (32,1300) foso,xcen+xll,&
                     ycen+yll,zcen,rad,angle,volume(m),area(m)
                DO i=i1,i2
                  DO j=j1,j2
                    IF (zmid(i,j).lt.ozb(i,j).and.subset(i,j).eq.m) THEN
                      ozb(i,j) = zmid(i,j)
                    END IF
                  END DO
                END DO
              END IF
            END IF  !  IF (remove.eq.'A')

!     Save minimum 3-D factors of safety and associated 2-D values if necessary
!     for all angles tried at this search sphere.
            IF (foso.lt.fosmin) THEN
              fosmin = foso
!!!  Fellenius comparison addition
              IF (method.eq.'b') felfosmin = felfoso
!!!
              fangle = angle
              fminma = minma    
              IF (ifos2d.eq.1) THEN ! Most 2d values change with slip direction
                fosmin3d2d = foso2d  
                farea3d2d = area2d
                farclength3d2d = arclength
                fwidth3d2d = width2d
                fncol3d2d = ncol2d
              END IF  
            END IF  ! IF (foso.lt.fosmin)
            
!     If calculating 2-D factors of safety, save minimum 2-D factor of safety
!     and associated 3-D factor of safety for all angles.
            IF (ifos2d.eq.1) THEN
              IF (foso2d.lt.fosmin2d) THEN
                fosmin2d = foso2d
                fosmin2d3d = foso
                fangle2d = angle 
                fwidth2d = width2d
                farea2d = area2d
                farclength2d = arclength
                fncol2d = ncol2d
              END IF
            END IF   ! IF (ifos2d.eq.1)         
                        
          END DO  ! loop on iseg
                            
          IF (single.ne.1) THEN 
!     If # columns less than limcol, write search location to ncolerr.out.            
            IF (ncol(m).lt.limcol) THEN
              IF (colfile.eq.0) THEN  
                IF (iretro.gt.0) THEN              
                  OPEN (34,STATUS='replace',file = outputdir(1:LEN_TRIM(outputdir))//&
                   filin(bfile:nfile)//'_ncolerr'//charetro//'_out.txt')
                ELSE
                  OPEN (34,STATUS='replace',file = outputdir(1:LEN_TRIM(outputdir))//&
                        filin(bfile:nfile)//'_ncolerr_out.txt')
                END IF
                WRITE (34,1000)
                IF (method.eq.'b') THEN
                     fosmeth = 'Bish'
                ELSE
                    fosmeth = 'Ord'
                END IF                             
                IF (LEN_TRIM(lengthunits).eq.1) THEN
                    WRITE (34,1110) lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,fosmeth
                ELSE
                   IF (LEN_TRIM(lengthunits).eq.2) THEN
                      WRITE (34,1120) lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,fosmeth
                   ELSE
                      WRITE (34,1100) fosmeth
                  END IF  
               END IF                             
                colfile = 1
                xcolcount = 0
              END IF
              xcolcount = xcolcount + 1
              WRITE (34,1200) ncol(m),xcen+xll,ycen+yll,zcen,rad,volume(m),area(m),foso
            END IF  ! IF (ncol(m).lt.limcol)
          END IF  ! IF (single.ne.1)

!     Keep track of minimum factor of safety at each search lattice
!     node and write to file if ilattice = 1.
          IF (single.ne.1.and.(ilattice.eq.1.or.icritlattice.eq.1)) THEN             
            IF (fosmin.le.fgrid) THEN
              fgrid = fosmin   
            END IF
          END IF
            
          IF (fosmin.eq.nullhi) CYCLE

!     Check if 3-D factor of safety is less than previous minimum 
!     for each DEM column.
          DO j= j1,j2
            DO i = i1,i2
              IF (subset(i,j).eq.m) THEN
                IF (fosmin.lt.fsmin(i,j)) THEN
                  fsmin(i,j) = fosmin
!!!  Fellenius comparison addition
                  IF (method.eq.'b') felfsmin(i,j) = felfosmin
!!!
                  IF (ifos2d.eq.1) fsmin3d2d(i,j) = fosmin3d2d
                  fsrad(i,j) = rad
                  fsangle(i,j) = fangle
                  fsx(i,j) = xcen
                  fsy(i,j) = ycen
                  fsz(i,j) = zcen
                  mcol(i,j) = ncol(m)
                  fsvol(i,j) = volume(m)
                  fsarea(i,j) = area(m)
                  icritprox = 0
                  IF (single.ne.1) THEN
!    Mark surfaces that were within one tolerance value of volume and
!    area bounds.
                    IF (vacriterion.eq.'v') THEN
                      IF (volume(m).lt.vmin+tol) icritprox = 1
                      IF (volume(m).gt.vmax-tol) icritprox = 2                        
                    ELSE
                      IF (area(m).lt.armin+tol) icritprox = 1
                      IF (area(m).gt.armax-tol) icritprox = 2
                    END IF
                  END IF
                  icrit(i,j) = icritprox
                END IF  ! IF (fosmin.lt.fsmin(i,j))
                
!     Track minimum 2-D FOS as well.
                IF (ifos2d.eq.1) THEN
                  IF ((linein(i,j).eq.1).and.(fosmin2d.lt.fsmin2d(i,j))) THEN
                    fsmin2d(i,j) = fosmin2d
                    fsmin2d3d(i,j) = fosmin2d3d
                    fsrad2d(i,j) = rad
                    fsangle2d(i,j) = fangle2d
                    fsx2d(i,j) = xcen
                    fsy2d(i,j) = ycen
                    fsz2d(i,j) = zcen
                    mcol2d(i,j) = fncol2d
                    fsarea2d(i,j) = farea2d
                    fswidth2d(i,j) = fwidth2d
                    fsarclength(i,j) = farclength2d
                  END IF
                END IF  ! IF (ifos2d.eq.1)
!     Keep running sum of number of times each column is filtered.               
                IF (mflag.ne.0) filtcount(i,j) = filtcount(i,j) + 1                 

                IF (isubsurf.eq.1.or.isubsurf.eq.2.or.isubsurf.eq.3) THEN    
!     Find z stopping point at 5 nodes above actual DEM.
                  ztop = zdem(i,j) + 5.0_pr*delz
!     Now check z values above base of trial surface.
                  DO k = 1,nz
                    z = zbot + REAL(k-1,pr)*delz
!     Compare FOS, replace if z >= zmid and  z <= stop elev
!     and FOS less.
                    IF (z.ge.zmid(i,j).and.z.le.ztop) THEN
                      IF (fosmin.lt.zfos(i,j,k)) zfos(i,j,k)=fosmin
                    END IF
                  END DO 
                END IF  ! IF (isubsurf.eq.1.or.isubsurf.eq.2.or.isubsurf.eq.3)
              END IF   ! IF (subset(i,j).eq.m) 
            END DO  ! loop on i
          END DO  !  loop on j

          if (fosmin.lt.0.0) then
            print *,'WARNING - fosmin<0!!',fosmin,xcen+xll,ycen+xll,zcen,rad
            WRITE (39,*) 'WARNING - fosmin<0!!',fosmin,xcen+xll,ycen+xll,zcen,rad
          end if

!     Check if overall factor of safety is minimum.
          IF (fosmin.lt.ofos) THEN
            ofos = fosmin
!!!  Fellenius comparison addition
            IF (method.eq.'b') felminfos = felfosmin
            orad = rad
            oangle = fangle
            ofsx = xcen
            ofsy = ycen
            ofsz = zcen
            ovol = volume(m)
            oarea = area(m)
            ocol = ncol(m)
            omset = m
            ominma = fminma
            owt = weight
            osliparea = sliparea
            IF (ifos2d.eq.1) THEN
              ofos3d2d = fosmin3d2d
              oarea3d2d = farea3d2d
              oarclength3d2d = farclength3d2d
              owidth3d2d = fwidth3d2d
              ocol3d2d = fncol3d2d
            END IF         
          END IF
          IF (ifos2d.eq.1) THEN
            IF (fosmin2d.lt.ofos2d) THEN
              ofos2d = fosmin2d
              orad2d = rad
              oangle2d = fangle2d
              ofsx2d = xcen
              ofsy2d = ycen
              ofsz2d = zcen
              oarea2d = farea2d
              oarclength = farclength2d
              owidth2d = fwidth2d
              ocol2d = fncol2d
              ofos2d3d = fosmin2d3d
              oarea2d3d = area(m)
              ovol2d3d = volume(m)
              ocol2d3d = ncol(m)
              owt2d3d = weight
              osliparea2d3d = sliparea
            END IF
          END IF

!     Keep track of largest vol scoop with FOS < cut off (option remove = L).
          IF (remove.eq.'L') THEN
            IF (fosmin.lt.foscut.and.volume(m).gt.ovvol) THEN
              ovfos = fosmin
              ovrad = rad
              ovangle = fangle
              ovfsx = xcen
              ovfsy = ycen
              ovfsz = zcen
              ovvol = volume(m)
              ovarea = area(m)
              ovmset = m
              ovwt = weight
              ovsliparea = sliparea
            END IF
          END IF      
         sdeg = 365.0_pr
        END DO  ! loop on subsets m
  
1000   FORMAT ('Trial surfaces with columns < limcol')
1100   FORMAT ('# cols       xcen          ycen         zcen      radius    volume      area       F_',A4)
1110   FORMAT ('# cols     xcen_',A1,'          ycen_',A1,'       zcen_',A1,'   radius_',A1,' volume_',A1,'^3  area_',A1,&
                            '^2  ','   F_',A4)
1120   FORMAT ('# cols     xcen_',A2,'         ycen_',A2,'      zcen_',A2,'  radius_',A2,' volume_',A2,'^3 area_',A2,&
                          '^2 ','  F_',A4)
1200   FORMAT (i5,2f16.4,2f10.3,2es11.3,2x,f8.4)      
1290   FORMAT (3(' ',f10.2),' ',f11.4,' ',f8.4,' ',2es11.3,f7.4)
1300   FORMAT (f7.4,2f16.4,' ',f10.2,' ',f11.4,' ',f8.4,' ',2es11.3)

1400   FORMAT (2f16.4,' ',f10.2,' ',f11.4,' ',f8.4,' ',2es11.3,f7.4)

1500   FORMAT ('This file contains data for trial slip surfaces filtered by&
                & absminma')
1510   FORMAT ( ' and nonconverging surfaces.')
1520   FORMAT ('Minimum absolute value m-alpha: ',f10.5) 
1540   FORMAT ('Converging F are written here, but not used in search for minimum F.')
1550   FORMAT ('           xcen           ycen      zcen    radius     angle     minma      F_Bish    F_Ord')     
1551   FORMAT ('         xcen_',A1,'         ycen_',A1,'    zcen_',A1,'   radius_',A1,'   angle     minma     F_Bish      F_Ord')     
1552   FORMAT ('        xcen_',A2,'         ycen_',A2,'   zcen_',A2,'   radius_',A2,'   angle     minma     F_Bish    F_Ord')     
1560   FORMAT (2f16.4,3f10.3,es11.3,2f10.3) 
1620   FORMAT (2f16.4,3f10.3,f10.3)
       
       END SUBROUTINE Fos
