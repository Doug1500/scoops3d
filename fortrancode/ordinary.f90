        SUBROUTINE Ordinary (m,i1,i2,j1,j2,ang,cosang,sinang,i2dcen,j2dcen,&
                             x2dcen,y2dcen,rad2d)
     
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!       This subroutine calculates local and overall factors of safety
!       using the Ordinary (Fellenius) method.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!       Called by Fos
!
!       VARIABLES
!
!      ang --  angle of slip converted to radians
!      angle -- assumed angle of slip (azimuth)
!      arclength -- arc length of current 2d slice
!      area2d -- cross-sectional area of current 2d slice
!      avetop(nx,ny) -- average elevation of 4 corner nodes of intersected 
!        column calculated using the average of the four surrounding cells
!        for full columns, and by finding the intersection of the slip surface
!        with the smoothed DEM surface for partial columns.
!      cee(nmat) -- cohesion of each layer
!      cflag -- indicates whether 3-d cohesion data are provided (0=yes,1=no)
!      coh -- cohesion at base of column
!      cohes(nx,ny,nzstr) -- array of 3-d cohesion data
!      colarea2d -- area of 2d column
!      colxy(nx,ny) -- length of each column side (delxy for full columns)
!      cosang -- cosine of slip direction angle
!      cosdip -- cosine of apparent dip of slip surface in column
!      costruedip -- cosine of true dip of slip surface in a column
!      dcohdz(nx,ny,strmaxk) -- array of cohesion gradients between vertical nodes
!      delxy -- DEM grid resolution (delta x, delta y)
!      delz -- difference between elevation of slip surface base
!         and nearest pressure data node.
!      delzp -- pressure head data vertical resolution
!      dflag -- indicates whether 3-d unit weight data are provided (0=yes,1=no)
!      dfricdz(nx,ny,strmaxk) -- array of friction gradients between vertical nodes
!      df -- local driving force
!      dip -- slope of column base in direction of slip angle (apparent dip)
!      dlay -- thickness of layer at each column
!      dlayd,dlayw -- wet and dry thickness of a layer in each column 
!      dpdz(nx,ny,nz) -- array of pressure gradients
!      drive -- overall driving force
!      duwtdz(nx,ny,strmaxk) -- array of unit weight gradients between vertical nodes
!      duwtlaydz(nx,ny,nmat-1) -- gradient of depth-weighted unit weights between layers
!      dx1,dx2,dy1,dy2 -- x and y grid spacing of partial column relative to inside nodes
!      dzdx,dzdy -- used to calculate dip
!      dzstr -- z spacing of 3-d strength data
!      earm -- moment arm of earthquake moment, taken from center of
!         column to elevation of search center
!      eq -- earthquake acceleration coefficient
!      failsurf(nx,ny) -- array of failure surface elevations for each DEM column            
!      failsurfslope(nx,ny) -- slope of input failure surface
!      fflag -- indicates whether 3-d friction data are provided (0=yes,1=no)
!      filter -- indicates whether solution filters are applied (malpha and fos)
!      foso(numdir) -- factor of safety of current scoop at each slip angle
!      foso2d(numdir) -- 2-d factor of safety of current scoop at each slip angle
!      fxa(nmat) -- curve-fitting parameter, a from the equation for the soil-water characteristic 
!            curve defined by Fredlund and Xing (1994) 
!      fxn(nmat) -- curve-fitting parameter, n from the equation for the soil-water characteristic 
!            curve defined by Fredlund and Xing (1994) 
!      fxm(nmat) -- curve-fitting parameter, m from the equation for the soil-water characteristic 
!            curve defined by Fredlund and Xing (1994) 
!      fxr(nmat) -- curve-fitting parameter, psi_r (soil suction associated with residual moisture 
!            content) from the equation for the soil-water characteristic curve defined by 
!            Fredlund and Xing (1994)
!      fsurf -- moment arm of hydrostatic pressure applied to submerged columns
!      gamr(nmat,3) -- array of total, partially saturated, and saturated
!         unit weights for each layer
!      gamw -- water unit weight in units of problem
!      gamsurf -- unit weight of surface load in units of problem   
!      gsameall -- flag indicating whether all densities are equivalent
!      gsameeach -- flag indicating whether densities in each layer are equiv
!      hsforce -- hydrostatic force on a partial column with surface water layer
!      i1,i2,j1,j1 -- Array bounds for current set
!      i2dcen,j2dcen -- i,j location of center of search sphere intersection 
!         with DEM surface
!      ifailbase -- flag for whether slip base is from failure surface or sphere
!      ifailsurf -- flag for whether failure surface file is used
!      ifos2d -- flag for calculating 2-d factors of safety
!      in(nx,ny) -- number of nodes of each column bounded by sphere
!      isurfwat -- flag for analysis with surface water layer
!      isurfwater -- flag for whether column has surface water layer above the DEM.
!      irot -- local flag for whether irotcen is implemented
!      irotcen -- flag for specfied rotational center for single surface
!      iwater -- flag for method of specifying 3D pressure heads or piezometric
!         surface, where 0 = no water pressures,1 =  ru approximation, 2 =  piezometric surface, 
!         3=input 3-d pressure file, 4=3D variably saturated file containing
!         pressure head and moisture content, 5=3D variably saturated file,
!         moisture content from vanGenuchten SWCC, 6=3D variably saturated file,
!         moisture content from Fredlund and Xing SWCC
!      kklo -- 3D pressure array location just below slip surface
!      l -- layer number
!      last -- previous layer found in column
!      layer(nmat,nx,ny) -- array of bottom elevations for material layers
!      line -- length of 2-D slip line in current column
!      lineflag -- flag for determining whether 2-D line/column intersection
!         already found for current j. Skip calculation if so (=-1).
!      linein(nx,ny) -- array of flags for whether 2-D slide directon line intersects
!         column
!      linterp -- linear interpolation flag; 1=use linear interpolation, other=use
!         nearest node below slip surface base for strength data.
!      m -- set number
!      maxpk(nx,ny) -- array of highest k value of pressure data at each DEM column
!      maxstrk(nx,ny) -- array of highest k value of strength data at each DEM column
!      mline -- slope of slip direction line relative to x and y grid
!      moist -- determines which element of soil unit weight array to use: 1-total weight
!         2-partially saturated, 3-fully saturated.
!      ncol(nnset) -- number of columns in subset
!      ncol2d -- number of columns in current 2-D slice
!      nmat -- number of material layers
!      numdir -- number of search slip directions for a single failure
!      outnode(2) -- array of partial column nodes which fall outside
!         slip surface
!      outnodes(nx,ny,2) -- corner location of column nodes which fall
!         outside the search sphere.
!      parea -- area of column projected on DEM grid.  Parea equals delxy**2 
!         except for partial columns, which have a quadrilateral shape
!      pgrid -- flag for pressure head grid, 1=regular, 2=irregular grid
!      ph -- pressure head on slip surface base in column
!      piezo(nx,ny) -- array of piezometric elevation at each DEM column
!      pmaxk -- maximum number of 3-d pressure values at each DEM column
!      pressh(nx,ny,pmaxk) -- array of pressure head data
!      presshpz(nx,ny,pmaxk) -- array of z elevations of pressure head data
!      pvol -- volume of intersected column below piezometric surface
!      pzsurf -- flag for whether piezometric surface was found,1=yes,0=no
!      r -- moment arm for resisting moment
!      rr -- square of moment arm
!      radsq -- square of search sphere radius
!      resist -- overall resisting moment
!      rfc(nx,ny) -- cohesional component of resisting force
!      rff(nx,ny) -- frictional component of resisting force
!      rnull -- real null value used throughout program 
!      ru(nmat) -- pore pressure ratio approximation
!      ru2 -- pore pressure at slip surface base
!      sinang -- sin of slip direction angle
!      sindip -- sin of apparent dip of slip base in slide direction
!      single -- flag for calculating single slip surface
!      sliparea -- area of slip surface
!      subset(nx,ny) -- indicates set membership of DEM columns
!      str3d -- flag for using 3-d strength file
!      strgrid -- flag for strength data grid; 1=regular, 2=irregular grid
!      strmaxk -- number of 3-d strength values at each DEM column
!      strz(nx,ny,strmaxk) -- array of z locations of irregular strength data
!      surfwat(nx,ny) -- array of surface water elevations at each DEM cell
!      tanfric -- tangent of friction angle at base of column
!      tanphi(nnmat) -- tangent of friction angle for each layer
!      tfric(nx,ny,strmaxk) -- array of 3-d friction data
!      vga(nmat) -- curve-fitting parameter, alpha from the equation for the 
!           soil-water characteristic curve defined by vanGenuchten (1980) 
!      vgn(nmat) -- curve-fitting parameter, n, from the equation for the soil-water 
!            characteristic curve defined by vanGenuchten (1980) 
!      vol(nx,ny) -- volume of current column in slip surface
!      volfrac2d -- ratio of 2-d area to volume for column used to calculate 2-d FOS
!      watht -- depth of surface water in column
!      watwt -- weight of surface water in column
!      weight -- weight of slide
!      width2d -- width of current 2d slice
!      wt -- weight of column
!      x2dcen,y2dcen -- x,y location of center of search sphere intersection 
!         with DEM surface
!      xcen,ycen,zcen -- location of center of search sphere
!      xdem,ydem -- x and y locations of column
!      xfailslope(nx,ny),yfailslope(nx,ny) -- x and y slope vectors of input failure surface
!      xmid -- x location of middle of column
!      xprime -- x component of column midpoint projected on axes transformed by slip angle
!      x(4),y(4) -- x and y locations of each column node.
!      ymid -- y location of middle of column
!      xrad,yrad,zrad -- x,y, and z coordinates of equation for search sphere radius
!      zdiff -- difference between strength or pressure elevation and elevation 
!         of base of slip surface
!      zlay(nmat) -- local variable for layer elevation in a column
!      zmid(i,j) -- elevation of slip surface at column midpoint
!      zmidbase -- local variable for elevation of slip surface at column midpoint
!      zpiez -- local column variable for elevation of piezometric
!         surface at each column 
!      zpmin -- minimum elevation of pressure head data
!      zstrmin -- minimum elevation of strength data
!      ztop -- average of node elevations for column
!      zz -- square of zrad
!
!            dx2 
!          4----3
!      dy1 |    | dy2   numbering scheme for corners of column base
!          |    | 
!          1----2
!            dx1
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      USE CommonData
      USE GridData, ONLY: delxy,xcen,ycen,zcen,angle,radsq,zdem,pi,halfdelxy,&
            zdemnodes,xcenrot,ycenrot,zcenrot,cellarea
      USE MaterialData, ONLY: nmat,gsameall,gsameeach,gamw,eq,layer,cee,tanphi,&
            gamr,ru,uwtlay,duwtlaydz,gamsurf,thetasat,thetares,&
            fxa,fxn,fxm,fxr,vga,vgn
      USE SetData, ONLY: subset,insphere,in,ncol,xdem,ydem,sliparea,weight,zmid,&
            vol,colxy,outnodes,avetop,xslope,yslope
      USE WaterData  
      USE StrengthData, ONLY: strmaxk,str3d,cflag,fflag,dflag,strgrid,linterp,&
            maxstrk,zstrmin,dzstr,cohes,tfric,strz,dcohdz,dfricdz,duwtdz,uwt3d,uwt
      USE FOSData, ONLY: numdir,foso,foso2d,mfile,diter,ifos2d,arclength,&
            area2d,ncol2d,width2d,single,filter,foslocal,linein
      USE FailSurfData, ONLY: ifailsurf,failsurfslope,failsurf,xfailslope,yfailslope     
      USE SearchData, ONLY:  irotcen           
      
      IMPLICIT NONE            
      
      CHARACTER*200 :: filout
      
      INTEGER, INTENT(in) :: m,i1,i2,j1,j2,i2dcen,j2dcen
      INTEGER :: l,kklo,last,irot
      INTEGER :: i,j,kk,lineflag,isurfwater 
      INTEGER :: outnode1,outnode2,ifailbase 

      REAL(pr), INTENT(in) :: x2dcen,y2dcen,rad2d,cosang,sinang,ang
      REAL(pr) :: xmid,ymid,coh,fosl,wt,xprime,gamma,f,pore
      REAL(pr) :: zdiff,zpiez,zmidbase,xprime2d
      REAL(pr) :: sindip,cosdip
      REAL(pr) :: ztop,pvol,ru2,ph,colarea2d
      REAL(pr) :: dx1,dx2,dy1,dy2,parea,rff,rfc,dlayd,dlayw 
      REAL(pr) :: xrad,yrad,tanfric,earm,rr,r,rm
      REAL(pr) :: zz,zrad,dzdx,dzdy,dip,zlay(nmat),truearea,zrad2d
      REAL(pr) :: rff2d,rfc2d,mline,line,costruedip,volfrac2d      
      REAL(pr) :: resist,drive,rf,df,resist2d,drive2d,rf2d,df2d
      REAL(pr) :: z2d,cosdip2d,x1,x2,x3,x4,y1,y2,y3,y4 
      REAL(pr) :: x2d, y2d, xrad2d, yrad2d,surfdip,zsurf,normf,dip2dbase
      REAL(pr) :: watht,watwt,hsforce,fsurf,rrsurf,rmsurf,gammasurf
      REAL(pr) :: consta,tanalpha,tanomega,fdip,tandip2d,cossurfdip
      REAL(pr) :: a,n,mswcc,rswcc,cfx,tsat,tres,theta,securrent 
      
      CHARACTER*70 :: problemtype
      CHARACTER*120 :: errmessage,solution
      
      errmessage = ' '
      solution = ' '
      problemtype = 'in factor of safety calculation'        
      sliparea = 0.0_pr
      weight = 0.0_pr
      IF (cosang.ne.0.0_pr) THEN
        mline = TAN(ang)
      ELSE
        mline = rnull
      END IF                                                          
      drive = 0.0_pr 
      resist = 0.0_pr
      lineflag = 0
      resist2d = 0.0_pr
      drive2d = 0.0_pr  
      area2d = 0.0_pr
      arclength = 0.0_pr
      width2d = 0.0_pr
      ncol2d = 0
      earm = 0.0_pr
      ru2 = 0.0_pr
      IF (irotcen.eq.1) THEN
        irot = 1
      ELSE
        irot = 0
      END IF

      IF (foslocal.eq.1) THEN
        filout = outputdir(1:LEN_TRIM(outputdir))//filin(bfile:nfile)//'_foslocal_out.txt'
        OPEN (16,file = filout)
        WRITE (16,1000)
      END if
                                         
!     Loop through all columns in current set.                                          
      DO j = j1,j2
        lineflag = 0
        DO i = i1,i2

!     Bypass if column not within slip surface.
          IF (subset(i,j) .ne. m) CYCLE
          
!     Bypass if failure surface is defined, but value is rnull.
          IF (ifailsurf.eq.1) THEN
            IF (failsurf(i,j).eq.rnull) CYCLE
          END IF

!     Assign local variables.
          x1 = xdem(i,j)
          x2 = x1+delxy
          x3 = x2
          x4 = x1
          y1 = ydem(i,j)
          y2 = y1
          y3 = y1+delxy
          y4 = y3
          dx1 = colxy(i,j,1)
          dx2 = colxy(i,j,2)
          dy1 = colxy(i,j,3)
          dy2 = colxy(i,j,4)
          outnode1 = outnodes(i,j,1)
          outnode2 = outnodes(i,j,2)
!      avetop is defined in Volumes subroutine. It equals zdem for full columns
!      and the average of column corners for partial nodes.
          ztop = avetop(i,j)
          IF (pzsurf.eq.1) THEN
            zpiez = piezo(i,j)
          END IF
          IF (nmat.gt.1) THEN
            DO l=1,nmat-1
              zlay(l) = layer(l,i,j)
            END DO
          END IF
            
!    Initialize column variables.
          ifailbase = 0
          isurfwater = 0
          tanfric = 0.0_pr
          coh = 0.0_pr
          wt = 0.0_pr
          watwt = 0.0_pr
          watht = 0.0_pr
          pore = 0.0_pr
          IF (ifos2d.eq.1)THEN
            line = 0.0_pr
            linein(i,j) = 0
          END IF

!     If there is a defined failure surface determine whether it defines the column base or
!     the sphere does (whichever is higher for full columns. Only use sphere for partial columns.)
          IF (ifailsurf.eq.1) THEN
            IF (zmid(i,j).le.failsurf(i,j).and.in(i,j).eq.4) THEN
              ifailbase = 1
              zmid(i,j) = failsurf(i,j)
              zmidbase = zmid(i,j)
            END IF
          END IF
          
!     If there is a surface water layer determine whether above this column
          IF (isurfwat.eq.1) THEN
            IF (surfwat(i,j).gt.ztop) isurfwater = 1
          END IF
    
!     Find x,y center of column.
          SELECT CASE (in(i,j))
            CASE (4)  ! Full column
              xmid = x1 + halfdelxy
              ymid = y1 + halfdelxy
              
            CASE (2)   !  Partial column with 2 nodes in
              IF ((outnode1.eq.1.and.outnode2.eq.2).or. &
                  (outnode1.eq.1.and.outnode2.eq.4)) THEN
                xmid = x2-(dx1+dx2)*0.25_pr
                ymid = y3-(dy1+dy2)*0.25_pr
              ELSE
                xmid = x1+(dx1+dx2)*0.25_pr
                ymid = y1+(dy1+dy2)*0.25_pr
              END IF

            CASE (3)  !  Partial column with 3 nodes in
!     Shorter intersected side was set to the negative of the segment length to identify. 
!     Both segments are used for column center approximation. Excluded triangular segment
!     is deleted from full column area to find partial column projected area.            
              SELECT CASE (outnode1)
                CASE (1)
                  xmid = x2-(abs(dx1)+dx2)*0.25_pr
                  ymid = y3-(abs(dy1)+dy2)*0.25_pr
                  parea = cellarea-0.5_pr*(delxy-abs(dx1))*(delxy-abs(dy1))
                  IF (dx1.le.0.0_pr) THEN
                    dx1 = delxy
                  ELSE
                    dy1 = delxy
                  END IF
                CASE (2)
                  xmid = x1+(abs(dx1)+dx2)*0.25_pr
                  ymid = y3-(dy1+abs(dy2))*0.25_pr
                  parea = cellarea-0.5_pr*(delxy-abs(dx1))*(delxy-abs(dy2))
                  IF (dx1.le.0.0_pr) THEN
                    dx1 = delxy
                  ELSE
                    dy2 = delxy
                  END IF
                CASE (3)
                  xmid = x1+(dx1+abs(dx2))*0.25_pr
                  ymid = y1+(dy1+abs(dy2))*0.25_pr
                  parea = cellarea-0.5_pr*(delxy-abs(dx2))*(delxy-abs(dy2))
                  IF (dx2.le.0.0_pr) THEN
                    dx2 = delxy
                  ELSE
                    dy2 = delxy
                 END IF   
                CASE (4)
                  xmid = x2-(dx1+abs(dx2))*0.25_pr
                  ymid = y1+(abs(dy1)+dy2)*0.25_pr
                  parea = cellarea-0.5_pr*(delxy-abs(dx2))*(delxy-abs(dy1))
                  IF (dx2.le.0.0_pr) THEN
                    dx2 = delxy
                  ELSE
                    dy1 = delxy
                  END IF
              END SELECT
          END SELECT     

!    Find distance from column slip surface midpoint to sphere center.
          yrad = ymid-ycen
          xrad = xmid-xcen 
	               
            zz = radsq-xrad*xrad-yrad*yrad
            IF (zz.lt.0.0_pr) THEN
              errmessage = 'cannot calculate column midpoint'
              CLOSE (33)
              Call WriteError(1,errmessage,problemtype,'no','no ',0,' ')                  
            END IF            
            zrad = SQRT(zz)
            zmid(i,j) = zcen - zrad
            zmidbase = zmid(i,j)
          
            IF (zmidbase.ge.ztop) CYCLE
	                      
!     Use partial derivative of equation for the search sphere 
!     at column center to compute cosine of true dip of base (costruedip) 
!     and slope of base in slide direction (apparent dip).
!     (Downslope apparent dip is positive in the direction of the slide.)
            dzdx = -xrad/zrad
            dzdy = -yrad/zrad             

          costruedip=1.0_pr/SQRT(1.0_pr+ dzdx*dzdx + dzdy*dzdy)              
          tanalpha = dzdx*cosang + dzdy*sinang             
          dip = ATAN(tanalpha) 
          sindip = SIN(dip)
          cosdip = cos(dip) 
           
!     Find projected area for full columns (dx*dy) and partial columns with 2 nodes in.
!     Projected area for columns with 3 nodes was already calculated.
          IF (in(i,j).ne.3) THEN                                                                           
            parea = .250_pr * (dx1+dx2)*(dy1+dy2)
          END IF      
!     Accumulate slip surface area 
          truearea = parea/costruedip       
          sliparea = sliparea + truearea           
            
!     Find dip and apparent dip in slide direction of DEM surface if water on column.
          IF (isurfwater.eq.1) THEN  
!     Find slope of DEM in direction of slide
            dzdx = xslope(i,j)/(8.0_pr*delxy)
            dzdy = yslope(i,j)/(8.0_pr*delxy)
!     Adjust sign of slope so downslope apparent dip is positive in the direction of the slide	        
            IF (xslope(i,j).ne.0.0_pr) THEN
              IF (angle.lt.90.0_pr.or.angle.gt.270.0_pr) THEN
                IF (xslope(i,j).lt.0.0_pr) THEN
                  dzdx = -dzdx
                  dzdy = -dzdy
                END IF
              ELSE
                IF (xslope(i,j).gt.0.0_pr) THEN
                 dzdx = -dzdx
                 dzdy = -dzdy
                END IF
              END IF 
            ELSE
              IF (angle.gt.0.0_pr.and.angle.lt.180.0_pr) THEN
                IF (yslope(i,j).lt.0.0_pr) THEN
                  dzdx = -dzdx
                  dzdy = -dzdy
                END IF
              ELSE
                IF (yslope(i,j).gt.0.0_pr) THEN
                 dzdx = -dzdx
                 dzdy = -dzdy
                END IF
              END IF 
            END IF           
            tanomega = dzdx*cosang + dzdy*sinang   
            surfdip = ATAN(tanomega)
            cossurfdip = COS(surfdip)
          END IF

!      Determine whether slip direction line intersects column for 2-D
!      solution and calculate length of line and slip base in the column.
!      Only check if line/column intersection not yet found for this j
!      or found in adjacent column.
          IF (ifos2d.eq.1.and.lineflag.ge.0) THEN 
            CALL line2d(i,j,x1,x2,x3,x4,y1,y2,y3,y4,i2dcen,j2dcen,x2dcen,&
                        y2dcen,cosang,sinang,mline,outnode1,outnode2,dx1,&
                        dx2,dy1,dy2,line,zrad2d,tandip2d,x2d,y2d)
            line = abs(line)
            IF (line.gt.0.0_pr) THEN
              IF (ifailbase.eq.1) THEN
                dip2dbase = dip
                z2d = zmidbase
              ELSE
                dip2dbase = ATAN(tandip2d)
                z2d = zcen - zrad2d
              END IF            
              width2d = width2d + line
              cosdip2d = COS(dip2dbase) 
              arclength = arclength + abs(line/cosdip2d)
              ncol2d = ncol2d + 1
              IF (ztop.gt.z2d) THEN
                colarea2d = line*(ztop-z2d)                      
                area2d = area2d + colarea2d
!     Find ratio of 2D area to volume of columns intersected by slip direction.              
                volfrac2d = colarea2d/vol(i,j)
                lineflag = 1
                linein(i,j) = 1
              ELSE
                colarea2d = 0.0_pr
                volfrac2d = 0.0_pr
                lineflag = 1
              END IF
            ELSE
              IF (lineflag.eq.1) lineflag = -1
            END IF
          END IF            

!     If 3D strength file was read...
          IF (str3d.eq.1) THEN
!     If slip base is below bottom data point exit program.
            IF (zmidbase.lt.strz(i,j,1)) THEN
              errmessage = 'slip surface below defined strength data'
              solution = 'increase depth of strength data or decrease volume of potential failures'              
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                                                
              PRINT *,'    DEM cell location, i= ',i,'j= ',j
              PRINT *,'    elevation of slip surface = ',zmidbase
              PRINT *,'    minimum elevation of strength data = ',strz(i,j,1)
              WRITE (39,*) '     DEM cell location,    i =  ',i,'j = ',j
              WRITE (39,*) '     elevation of slip surface = ',zmidbase
              WRITE (39,*) '     minimum elevation of strength data = ',strz(i,j,1)             
              CLOSE (20)
              CLOSE (33)
              CLOSE (39)
              STOP            
            END IF
 
!     If uniform grid ...
            IF (strgrid.eq.1) THEN                   
!     Find first strength node below column center base elevation                
              kklo = INT((zmidbase-zstrmin)/dzstr) + 1 
!     If non-uniform grid
            ELSE 
              kklo = maxstrk(i,j)                                          
              DO kk = 2,maxstrk(i,j)
!     Find first strength node below column center base elevation
                IF (strz(i,j,kk).gt.zmidbase) THEN 
                  kklo = kk - 1
                  EXIT
                END IF 
              END DO  
            END IF    !  IF (strgrid.eq.1)  
            zdiff = zmidbase - strz(i,j,kklo)
!     If not using linear interpolation use strength values from first node below base.                   
            IF (linterp.ne.1) THEN          
              IF (fflag.eq.0) tanfric = tfric(i,j,kklo)
              IF (cflag.eq.0) coh = cohes(i,j,kklo)                  
            ELSE  
!     If interpolating, use strength values at first node below base and strength gradient.
              IF (fflag.eq.0) tanfric = tfric(i,j,kklo) + dfricdz(i,j,kklo)*zdiff
              IF (cflag.eq.0) coh = cohes(i,j,kklo) + dcohdz(i,j,kklo)*zdiff
            END IF 
!     Use gradient from bottom node to find depth-averaged unit weight at slip base.
            IF (dflag.eq.0) wt = vol(i,j) * (duwtdz(i,j,kklo)*zdiff + uwt3d(i,j,kklo))                                 
          END IF  ! IF (str3d.eq.1) 
  
!     If one layer, assign column cohesion, friction and ru values.
          IF (nmat.eq.1) THEN
            IF (cflag.eq.1) coh = cee(1)
            IF (fflag.eq.1) tanfric = tanphi(1)
            IF (iwater.eq.1) ru2 = ru(1) 
            IF (iwater.eq.5.or.iwater.eq.6) THEN
                tsat = thetasat(1)
                tres = thetares(1)
                SELECT CASE (iwater)
                    CASE(5)  ! vanGunuchten curve fit  
                       a = vga(1)
                       n = vgn(1)
                    CASE(6)  ! Fredlund and Xing curve fit 
                       a = fxa(1)
                       n = fxn(1)
                       mswcc = fxm(1)
                       rswcc = fxr(1)                          
!                    CASE(7)  ! Gardner curve fit
!                       a = gra(1)
!                       n = grn(1) 
                END SELECT                                                                
            END IF    
          END IF
                     
!     If only 1 material layer or all densities for all layers are same,
!     calculate weight of column.
          IF (iwater.lt.4.and.dflag.eq.1.and.((nmat.eq.1).or.(gsameall.eq.1))) THEN
            IF (pzsurf.eq.0.or.gsameall.eq.1) THEN
              wt = vol(i,j)*gamr(1,moist)
            ELSE
              IF (zpiez.le.zmidbase) THEN
                wt = vol(i,j)*gamr(1,2)
              ELSE
                IF (zpiez.ge.ztop) THEN
                  wt = vol(i,j)*gamr(1,3)
                ELSE
                  dlayd = ztop-zpiez
                  dlayw = zpiez-zmidbase
                  wt = vol(i,j) * (dlayd*gamr(1,2) + dlayw*gamr(1,3))/(ztop-zmidbase)
                END IF
              END IF
            END IF
          END IF

!     If multiple material layers, compute weight of each column based on 
!     depth-weighted average unit weights.
          IF (nmat.gt.1) THEN 
            zdiff = 0.0_pr
            last = nmat
            zlay(nmat) = zmidbase
            DO l = 1,nmat 
!     If current layer is intersected by column, check if base of material is 
!     above slip surface.  
              IF (zlay(l).ne.rnull.and.zlay(l).lt.ztop) THEN
                zdiff = zmidbase - zlay(l)
              ELSE
                CYCLE
              END IF
!     If layer bottom is above slip surface, go to next layer
              IF (zdiff.lt.0.0_pr) THEN
                last = l
!     If layer bottom is below slip surface, calculate column weight
              ELSE 
!     Skip if already calculated weight or using unsaturated weights
                IF (iwater.lt.4.and.gsameall.ne.1) THEN
                  IF (l.ne.nmat.and.last.ne.nmat) THEN
                    wt = vol(i,j) * (uwtlay(i,j,l) + duwtlaydz(i,j,l)*zdiff) 
                  ELSE 
                    IF ((pzsurf.eq.0).or.(gsameeach.eq.1)) THEN
                      IF (last.ne.nmat) THEN
                        zdiff = zlay(last)-zmidbase
                        wt = vol(i,j) * (zdiff*gamr(l,moist) + &
                           uwtlay(i,j,last)*(ztop-zlay(last)))/(ztop-zmidbase)
                      ELSE
                        wt = vol(i,j) * gamr(l,moist)
                      END IF
                    ELSE
                      IF (pzsurf.eq.1) THEN                   
                        IF (zlay(last).le.zpiez)THEN
                          IF (last.ne.nmat) THEN
                            zdiff = zlay(last)-zmidbase
                            wt = vol(i,j) * (zdiff*gamr(l,3) + &
                               uwtlay(i,j,last)*(ztop-zlay(last)))/(ztop-zmidbase)
                          ELSE
                            wt = vol(i,j) * gamr(l,3)
                          END IF
                        ELSE 
                          IF (zpiez.le.zmidbase) THEN
                            IF (last.ne.nmat) THEN
                              zdiff = zlay(last)-zmidbase
                              wt = vol(i,j) * (zdiff*gamr(l,2) + &
                                 uwtlay(i,j,last)*(ztop-zlay(last)))/(ztop-zmidbase)
                            ELSE
                              wt = vol(i,j) * gamr(l,2)
                            END IF
                          ELSE
                            IF (last.ne.nmat) THEN
                              dlayd = zlay(last)-zpiez
                              dlayw = zpiez-zmidbase
                              wt = vol(i,j) * (dlayd*gamr(l,2) + dlayw*gamr(l,3) +&
                                 uwtlay(i,j,last)*(ztop-zlay(last)))/(ztop-zmidbase)
                            ELSE
                              IF (zpiez.ge.ztop) THEN
                                wt = vol(i,j) * gamr(l,3)
                              ELSE
                                dlayd = ztop-zpiez
                                dlayw = zpiez-zmidbase
                                wt = vol(i,j) * (dlayd*gamr(l,2) + dlayw*gamr(l,3))/(ztop-zmidbase)
                              END IF
                            END IF
                          END IF
                        END IF
                      END IF
                    END IF
                  END IF
                END IF  
!     Determine friction and cohesion values of slip surface       
                IF (fflag.eq.1) tanfric = tanphi(l)
                IF (cflag.eq.1) coh = cee(l)
                IF (iwater.eq.1) ru2 = ru(l)
                IF (iwater.eq.5.or.iwater.eq.6) THEN
                  tsat = thetasat(l)
                  tres = thetares(l)                  
                  SELECT CASE (iwater)
                   CASE(5)  ! vanGunuchten curve fit  
                       a = vga(l)
                       n = vgn(l)
                   CASE(6)  ! Fredlund and Xing curve fit 
                       a = fxa(l)
                       n = fxn(l)
                       mswcc = fxm(l)
                       rswcc = fxr(l)                          
!                    CASE(7)  ! Gardner curve fit
!                       a = gra(l)
!                       n = grn(l)
                  END SELECT                                                                                             
                END IF                                
                EXIT
              END IF
            END DO  !  layer loop 
          END IF  !  If (nmat.gt.1)        

!     Calculate resisting force.
!     If using piezometric surface file find volume of saturated portion
!     of column.
          IF (iwater.eq.2) THEN
            IF (zpiez.ne.rnull) THEN
!     If piezometric surface is above DEM surface, assume entire
!     column is saturated.  Add weight of water column if there is surface 
!     water layer.
              IF (zpiez.ge.ztop) THEN
                pvol = vol(i,j)
              ELSE
                pvol = vol(i,j) - (ztop-zpiez)*parea
              END IF
            END IF
            IF (pvol.lt.0.0_pr) THEN
              pvol = 0.0_pr
              pore = 0.0_pr
            ELSE
                pore = pvol*gamw
            END IF
          ELSE
            
!     If using 3d pressure head data file and data is regularly spaced in z,
!     interpolate data to estimate pressure heads at slip surface.   
            IF (iwater.gt.2) THEN
              IF (zmidbase.lt.presshpz(i,j,1)) THEN              
                errmessage = 'slip surface below defined pressure head data'
                solution = 'increase depth of pressure head data or decrease volume of potential failures'
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                     
                PRINT *,'    DEM cell location, i = ',i,'j = ',j  
                PRINT *,'    elevation of slip surface = ',zmidbase
                PRINT *,'    minimum elevation of pressure head data = ',presshpz(i,j,1)          
                WRITE (39,*) '     DEM cell location,     i =',i,'j = ',j
                WRITE (39,*) '     elevation of slip surface = ',zmidbase
                WRITE (39,*) '     minimum elevation of pressure head data = ',presshpz(i,j,1) 
                CLOSE (20)
                CLOSE (33) 
                CLOSE (39) 
                STOP                
              END IF            
              IF (pgrid.eq.1) THEN  ! uniform pressure grid in z              
!     Find pressure head node just below average slip surface
                kklo = INT((zmidbase-zpmin)/delzp)+1
              ELSE              
!     If irregularly spaced pressure head data (pgrid=2)              
!     find pressure node just below base.
                kklo = maxpk(i,j)
                DO kk = 2,maxpk(i,j)
                  IF (presshpz(i,j,kk).gt.zmidbase) THEN
                    kklo = kk - 1
                    EXIT
                  END IF
                END DO                    
              END IF
              zdiff = zmidbase - presshpz(i,j,kklo)
              ph = zdiff*dpdz(i,j,kklo) + pressh(i,j,kklo)
!     If variably saturated use gradient from bottom node to find depth-averaged unit weight at slip base.
              IF (iwater.eq.4.or.iwater.eq.5.or.iwater.eq.6) THEN
                wt = vol(i,j) * (duwtdz(i,j,kklo)*zdiff + uwt3d(i,j,kklo))              
!     If variably saturated and negative pressure, multiply by effective saturation
                IF (ph.lt.0.0_pr) THEN
! vanGunuchten curve fit  
                  IF (iwater.eq.5) theta = tres + (tsat-tres)/(1+(a*abs(ph*gamw))**n)**(1-1/n) 
! Fredlund and Xing curve fit 
                  IF (iwater.eq.6) THEN
                       cfx = 1 - log(1+abs(ph*gamw)/rswcc)/log(1+(1000000/rswcc))
      		       theta = cfx * tsat / (log(exp(1.) + (abs(ph*gamw)/a)**n)**mswcc ) 
                  END IF                            
!              ! Gardner curve fit
!                       theta = tres + (tsat-tres)/(1+(a*ph)^n))                                                 
                  IF (iwater.eq.4) THEN
                       ph = ph * (zdiff * dsedz(i,j,kklo) + se(i,j,kklo)) 
                  ELSE
                       securrent = (theta - tres)/(tsat - tres)
                       IF (iwater.eq.6.and.abs(ph)*gamw.gt.rswcc) THEN
                         ph = 0.0_pr                         
                       ELSE    
                         ph = ph * securrent
                       END IF
                  END IF
                END IF  
              END IF
!      If pressure head on base is negative and not using the variable saturation option, 
!      set pressure head to zero for force calculation.
              IF (iwater.eq.3.and.ph.lt.0.0_pr) ph = 0.0_pr
              IF (iwater.eq.6.and.abs(ph).gt.rswcc) ph = 0.0_pr              
              pore = ph*gamw*parea                                   
            ELSE
!      If calculating resisting force using hydrostatic piezometric
!      surface approximation (iwater=1) or no water data (iwater=0,ru=0).
              pore = wt*ru2        
            END IF  !  if (iwater.eq.3)
          END IF  !  if (iwater.eq.2)        

!     Calculate moment arm for earthquake force from rotational
!     center to center of column.
          IF (eq.gt.0.0_pr) earm = zcenrot - (ztop + zmidbase)*0.5_pr        
                                         
!     Calculate other moment arms using rotational center
          xrad =  xmid - xcenrot
          yrad =  ymid - ycenrot
          zrad =  zcenrot - zmidbase
          xprime = -(xrad*cosang+yrad*sinang)
          rr = zrad*zrad + xprime*xprime
          IF (rr.le.0.0_pr) THEN
            errmessage = 'cannot calculate moment arm'
            CLOSE (33)
            Call WriteError(1,errmessage,problemtype,'no','no ',0,' ')             
          END IF
          rm = SQRT(rr) 
          r = rm
        
          consta = cosdip*cosdip/costruedip
          normf = consta*(wt*(1-eq*tanalpha) + watwt*(1+tanalpha*tanomega))
          rff =  (normf-pore/costruedip)*tanfric
          df = wt*(xprime + eq*earm)          

          IF (rff.lt.0.0_pr) rff = 0.0_pr
          rfc = coh*truearea
          rf = (rff + rfc)*r

!     Accumulate total driving force
          drive = drive + df                                                     

!     Accumulate total resisting force                                    
          resist = resist + rf 
        
!!     Accumulate forces for 2-D factors of safety
          IF (ifos2d.eq.1.and.line.gt.0.0_pr) THEN
            consta = cosdip2d*cosdip2d/costruedip
            normf = consta*(wt*(1-eq*tandip2d) + watwt*(1+tandip2d*tanomega))  
            rff2d = volfrac2d*(normf-pore/cosdip2d)*tanfric           
            IF (eq.gt.0.0_pr) earm = zcenrot - (ztop + z2d)*0.5_pr
            zrad2d = zcenrot - z2d
            xrad2d = x2d - xcenrot
            yrad2d = y2d - ycenrot
            xprime2d = -(xrad2d*cosang+yrad2d*sinang)
            rfc2d = coh*line/cosdip2d
            df2d = volfrac2d*wt*(xprime2d + eq*earm)
            IF (rff2d.lt.0.0_pr) rff2d = 0.0_pr
            rf2d = rff2d+rfc2d
            rr = zrad2d*zrad2d + xprime2d*xprime2d
            IF (rr.le.0.0_pr) THEN
              errmessage = 'cannot calculate 2D moment arm'
              CLOSE (33)
              Call WriteError(1,errmessage,problemtype,'no','no ',0,' ')              
            END IF
            rm = SQRT(rr)
            r = rm
 
!     Accumulate total driving force
            drive2d = drive2d + df2d                                                
!     Accumulate total resisting force                                    
            resist2d = resist2d + rf2d * r
          END IF
          
!     Accumulate total weight
          weight = weight + wt          
          
!     Calculate local factors of safety if required.
          IF (foslocal.eq.1) THEN    
            IF (df.eq.0.0_pr) THEN                                        
              fosl = 111.0_pr
            ELSE                                                       
              fosl = rf/df
            END IF
            WRITE (16,2000) i,j,truearea,rff,rfc,rf,df,fosl
2000    FORMAT (2i4,6es12.4)
          END IF
  
        END DO   !  loop on i                                                     
      END DO     !  loop on j                                                     
                                   

      IF (resist.le.0.0_pr) resist = 0.0_pr                                                       
      IF (drive.lt.0.010_pr*resist.or.drive.le.0.0_pr) THEN
        foso = 100.0_pr
      ELSE
        foso = resist/drive
      END IF
      
      IF (ifos2d.eq.1) THEN
        IF (resist2d.le.0.0_pr) resist2d = 0.0_pr
        IF (drive2d.lt.0.010_pr*resist2d.or.drive2d.le.0.0_pr) THEN
          foso2d = 100.0_pr
        ELSE
          foso2d = resist2d/drive2d
        END IF
      END IF
     

      CLOSE (16)
    
1000  FORMAT ('   i   j  sliparea    res. frict.  res. coh.  total res.  total driv. res./driv.')
                                   
      END SUBROUTINE Ordinary                                                              
