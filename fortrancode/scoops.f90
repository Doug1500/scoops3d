       PROGRAM Scoops3D
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!      This U.S. Geological Survey program, Scoops3D, analyzes 3D slope stability 
!      throughout a digital landscape.  It determines the minimum factor of safety 
!      associated with each cell in a digital elevation model (DEM).  
!
!      More information on the theoretical basis of the program, program operation, 
!      practical considerations, and testing of Scoop3D can be found in:  
!         Reid, M.E., Christian, S.B., Brien, D.L., and Henderson, S.T., 2014, 
!         Scoops3D-Software to Analyze Three-Dimensional Slope Stability Throughout
!         a Digital Landscape: U.S. Geological Survey Techniques and Methods, book 14, 
!         chap. A1  [http://pubs.usgs.gov/tm/14/a01].
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!      Program coding: Sarah Christian, Dianne Brien, Mark Reid
! 
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!      LIST OF VARIABLES
!
!      area(nnset) -- surface area of subset
!      armin,armax -- minimum and maximum area bounds for valid slip surface
!      colfile -- flag for opening error file reporting cols < limcol
!      dr -- radius increment for search
!      delxy -- DEM grid resolution
!      demflag -- indicates whether file newDEM will be generated
!      dzdx,dzdy -- surface gradients in x and y direction of DEM
!      error -- file allocation flag
!      failsurf(nx,ny) -- array of failure surface elevations for each DEM column            
!      failsurfslope(nx,ny) -- slope of input failure surface
!      felfsmin(nx,ny) -- array of minimum Ordinary (Fellenius) FOS at each DEM column
!      fgrid -- minimum factor of safety at search grid point
!      filtcount(nx,ny) -- counter of number of times each DEM cell was filtered
!         by minma.
!      filter -- indicates whether solution filters are applied (malpha and fos)
!      foscut -- factor of safety cutoff for removing slip surface from DEM
!      foso -- factor of safety of current scoop at each slip angle
!      foso2d -- 2D factor of safety of current scoop at each slip angle
!      fs1 through fs8 -- value of discontinuity surface used for slope calculation, 
!         use center cell value for any null cells in surrounding 8 cells.
!      fsangle(nx,ny) -- array of fall angles of minimum scoop at each DEM column
!      fsangle2d(nx,ny) -- array of fall angles of minimum 2d slice at each DEM column
!      fsarclength(nx,ny) -- array of arclengths of minimum 2d slice at each DEM column
!      fsarea2d(nx,ny) -- array of areas of minimum 2d slice at each DEM column
!      fsmin(nx,ny) -- array of minimum FOS at each DEM column
!      fsminold(i,j) -- last minimum factor of safety at each DEM column used to
!         compare change in FOS during seed iterations to see if criterion is met
!      fsmin2d(nx,ny) -- 2D array of minimum FOS at each DEM column
!      fsmin2d3d(nx,ny) -- 3D minimum factor of safety associated with min 2D slice
!         at each DEM cell
!      fsmin3d2d(nx,ny) -- 2D minimum factor of safety associated with min 3D scoop
!         at each DEM column
!      fsrad(nx,ny) -- array of radius of minimum scoop at each DEM column
!      fsrad2d(nx,ny) -- array of radii of minimum 2D FOS slice at each DEM column
!      fsvol(nx,ny) -- volume of minimum FOS scoop at each DEM column
!      fswidth2d(nx,ny) -- array of widths of minimum 2D slice at each DEM column
!      fsx(nx,ny),fsy(nx,ny),fsz(nx,ny) -- search sphere center of minimum FOS 
!         scoop at each DEM cell
!      fsx2d(nx,ny),fsy2d(nx,ny),fsz2d(nx,ny) -- search sphere center of minimum FOS 
!         2D slice at each DEM column
!      halfdelxy -- 1/2 DEM grid spacing
!      icen,jcen,kcen -- center of search sphere in i,j coordinates
!      icrit(nx,ny) -- array of volume or area flags indicating proximity to limits
!          of criteria range for each DEM column
!      ierr -- error flag indicating initial radius problem (2) 
!         or too many subsets (1)
!      icritlattice -- flag for creating file critfoslattice_out.3D for minimum FOS of 
!          critical surfaces at each search  grid node
!      ilatt,jlatt -- indices for search lattice files
!      ilattice -- flag for creating file foslattice_out.3D for minimum FOS at each 
!         search grid node
!      ifailsurf -- flag for whether failure surface file is used
!      in(nx,ny) -- number of nodes of each column bounded by sphere or
!        equals -1 if truncated surface.
!      insphere(nx+1,ny+1) -- indicates whether DEM column node is bounded by sphere
!      inull -- integer null value used throughout program
!      irefine -- flag for implementing coarse to fine search iterations
!      iretro -- counter for number of retrogressions
!      iseed -- grid refinement iteration count also used to fill srchlatt
!         lattice. iseed=1 represents initial coarse grid locations
!      ismin,ismax -- minimum and maximum boundaries of search grid on x-axis,
!         specified relative to DEM cells i=1 through nx
!      isqout -- flag for creating search quality output files
!      isrchmin,isrchmax -- array bounds of searchout output file
!      isrchnum,jsrchnum -- the number of search lattice nodes in x and y directions
!      isubsurf -- flag for creating 3D file of FOS below DEM surface,
!          1=create file in ijk format, 2=create file in xyz format
!          3=create .3D format
!       iwater -- flag indicating method for modeling water pressure
!          0 = no water pressures, 1 = ru approximation, 2 = piezometric surface, 
!          3 = 3d pressure head file, 4 = 3d pressure head file with relative moisture contents
!          for variable saturation calculation
!      jsmin,jsmax -- minimum and maximum boundaries of search grid on y-axis,
!       specified relative to DEM cells j=1 through ny
!      jsrchmin,jsrchmax -- array bounds of searchout output file
!      kseed -- array number representing search lattice elevation at current i,j
!      kset -- number of valid subsets
!      layer(nmat,nx,ny) -- array of bottom elevations for material layers
!      linein(nx,ny) -- array of flags for whether 2D slide directon line intersects
!         column
!      mcol(nx,ny) -- number of columns in minimum FOS scoop at each DEM column
!      mfile -- flag indicating whether filter file was opened
!      mini(nnset),maxi(nnset) -- minimum and maximum i bounds of each subset
!      minj(nnset),maxj(nnset) -- minimum and maximum j bounds of each subset
!      multres -- resolution multiplier for coarse search lattice
!      nbdy(nx+1,ny+1) -- array of flags for identifying nodes adjacent to 
!         DEM boundary cells
!      newrad -- indicates whether valid initial range subset has been found
!         (1=no, 0=yes, 2=no valid sets after 10 radius adjustments.)
!      nkseed(nxsrch,nysrch) -- number of k values for current seed at each i,j
!         search node
!      nmat -- number of material layers
!      nnset -- maximum number of subsets
!      nrange(nnset) -- array of flags for whether set fits volume or area criteria 
!         range
!      nres -- search lattice loop multiplier
!      nretro -- number of failure retrogressions to compute
!      nseed -- number of new search nodes generated from current seeds.
!      nseg -- number of search slip directions for each search sphere
!      nset -- number of subsets found
!      nsetmax -- flag for exceeding number of subsets allowed nnset
!      nsrchpt -- total number of search points
!      nsrchres -- resolution of finest search grid relative to DEM resolution.
!         Acts as a multiplier of DEM resolution delxy
!      nstate -- flag used by subroutine Criteria to indicate whether valid sets
!         were found at a particular search node.
!         100 = recalculate radius to find first good set at a search node.
!         300 = increase radius by dr and continue. This flag occurs if good
!               sets were found, or only sets not grown out of range are too
!               small to be in range.
!         500 = go to next search grid point because no good sets possible.
!      ntry -- total number of slip surfaces for which FOS is calculated
!      nullhi -- positive real null value used for initializations and set in CommonData.
!      nxsrchout,nysrchout -- number of search nodes in x and y directions in output
!         search grid file, includes intersection of search range with DEM range.
!      nx -- number of DEM cells in x direction
!      ny -- number of DEM cells in y direction
!      nz -- number of nodes used for 3D FOS values when isubsurf=1 or 2 or 3
!      oangle -- slip angle associated with overall minimum factor of safety
!      ocol -- number of columns in minimum scoop
!      ofos -- overall minimum factor of safety
!      ofsx,ofsy,ofsz -- search sphere center of minimum FOS scoop
!      omset -- set number of overall minimum scoop
!      orad -- angle of sphere associated with overall minimum factor of safety
!      outputdir -- optional directory path to place output files into
!      ovangle -- slip angle of largest volume scoop with FOS < cutoff
!      ovmset -- set number of largest volume scoop with FOS < cutoff
!      ovol, oarea - volume and area of minimum  FOS scoop
!      ovfos -- factor of safety of largest volume scoop with FOS < cutoff
!      ovrad -- radius of largest volume scoop with FOS < cutoff
!      ovfsx,ovfsy,ovfsz -- search sphere center of largest volume scoop with 
!         FOS < cutoff
!      ovvol,ovarea -- volume and area of largest volume scoop with FOS < cutoff
!      ozb(nx,ny) -- array containing elevation of DEM after scoops are removed 
!         depending on value of variable 'remove' (A, L, or M)
!      piezo(nx,ny) -- array of piezometric elevation at each DEM column
!      r2 -- variable used for calculating maximum radius
!      rad -- radius of search sphere
!      radiniprev -- used as first radius estimate at new search lattice elevation.
!         Equals previous radius plus elevation change, or = -1 if radius not yet
!         calculated at current search lattice x-y location.
!      radsq -- square of search sphere radius
!      remove -- A= all scoops < FOS cutoff will be removed from new
!         DEM.  L= only largest volume scoop < FOS cutoff will be 
!         removed. M= the scoop with the minimum FOS will be removed,
!         N= no scoops removed.
!      rmax -- maximum radius allowed
!      rnull -- real negative null value used throughout program and set in CommonData.
!      srchlatt(isrchnum,jsrchnum,ksmax) -- search lattice array with seed
!         iteration number at each searched node. 
!      searchgrid - array to specify search grid points. Values less than 0 are
!         not searched.
!      searchout -- output array of area covered by search grid and indicating
!         which grid points were searched and whether they fell inside or outside
!         the DEM.
!      single -- flag for calculating single slip surface
!      srchct -- running sum of number of search grid points completed
!      srchfile -- flag indicating whether using search file to specify valid
!         search grid nodes
!      srchpct -- percentage of search file completed
!      srchprint -- flag to report percent of completed search
!      srchpt -- flag for valid search grid location, 1=valid,0=not
!      str3d -- flag for using 3D strength file
!      subset(nx,ny) -- indicates set membership of DEM columns
!      surfslope(nx,ny) -- estimated slope of DEM surface at each DEM cell
!      vacriterion -- indicates primary control for intial radius (v or a
!         for volume or area)
!      va1max -- max. volume or area (primary criterion) of a subset
!      version -- version of Scoops
!      vmin,vmax -- minimum and maximum volume bounds for valid slip surface
!      volume(nnset) -- volume of subset
!      xcen,ycen,zcen -- center of search sphere relative to DEM origin
!      xcolcount -- number of critical surfaces with total # of columns < limcol
!      xcount -- number of surfaces eliminated due to m-alpha limits and/or nonconvergence
!      xdiff,ydiff -- variables used to calculate maximum radius
!      xlast,ylast,zlast -- last search grid node tried
!      xll,yll -- x and y origin of DEM grid, read in header lines
!      xllcorner,yllcorner -- character version of these header lines
!      xdem(nx,ny),ydem(nx,ny) -- x and y location of DEM column bottom left node
!      xslope(nx,ny),yslope(nx,ny) -- x and y slope vectors of DEM surface
!      xfailslope(nx,ny),yfailslope(nx,ny) -- x and y slope vectors of input failure surface
!      zavgmax -- average distance of surface of columns in subset to sphere center
!      zb(nx+1,ny+1) -- base of slip surface at each node
!      zmid(nx,ny) -- base of slip surface at center of intersected columns
!      zd -- elevation of DEM grid at same x,y as sphere center
!      zdem(nx,ny) -- DEM elevations
!      zfos(nx,ny,nz) -- 3D array of minimum factor of safety below DEM surface
!      zmin -- minimum DEM elevation
!      zsmin,zsmax -- minimum and maximum elevations of search grid
!      zsrchlatt(isrchnum,jsrchnum,ksmax) -- search lattice array with 
!         elevations of new search centers generated from seed.
!      zsrchres -- resolution of fine search grid on z-axis
!      zzsrchres -- resolution of current search grid on z-axis
!  
!      INPUT FILES
!       unit filename
!         10 search filename -- file that defines search grid. Opened and
!              closed in Readsearch.
!         10 include area filename -- file that defines area to include in all
!              failure surfaces. Opened and closed in Readinclude.
!         12 input filename -- file of basic input information. Opened in
!              Readin.
!         13 DEM filename -- file of DEM input information. Opened and
!              closed in Readdem.          
!         13 3D pressure filename -- file containing 3D pressure input.
!              Opened and closed in Readpressh.
!         13 3D strength filename -- file of 3D strength for each DEM column.
!              Opened and closed in Readstrength.
!         14 piezometric filename -- file of piezometric elevations at each
!              DEM column.  Opened and closed in Readpiezo.
!         14 strat layers filenames -- files of layer boundary elevations.
!              Each layer should have own file with extension giving the 
!              layer number, i.e. 'layer.1'. Opened and closed in Readstrat.
!
!      OUTPUT FILES
!       unit filename
!         15 'foslattice_out.3D' -- file of minimum factor of safety   
!               at each search grid point. Opened and written in Scoops 
!              if ilattice = 1.
!         16 'foslocal_out.txt' -- file of local resisting force, components
!              of the resisting force, driving force and factors
!              of safety for each column in minimum slip surface.
!              Opened and written in Ordinary if foslocal = 1.
!         17 'filtergrid_out.txt' -- file of number of times each DEM cell was
!              filtered due to absmina.  Opened and written in
!              Writeout if filter = 1 and mflag = 1.
!         19 'fos2d_out.asc' -- file of minimum 2D factors of safety at each
!              DEM column. Opened and written in Writeout if ifos2d=1.
!         20 'inputfilename_out.txt' -- echo of input and results containing
!              overall minimum slip surface data, written in subroutines
!              Readin, Readpiezo, Readpressh, Readsearch, Readstrat,
!              Readstrength, and Writeout.
!         21 'spheres_out.okc' -- file of sphere center coordinates, radius, 
!              fall angle, # of columns in scoop, vol, and factor of safety 
!              of the least stable scoop at each DEM column. Opened and written
!              in Writeout.
!         22 'critcheck_out.asc' -- file of volume (or area) checks.  If volume of minimum
!              FOS scoop is < vmin+tol shows 1, if > vmax-tol shows 2,
!              if not close to volume boundary shows 0. Opened and written
!              in Writeout if isqout = 1.
!	        23 'subsurffos_out.txt' or '.vtk' -- file of minimum FOS in 3D array, written
!              if isubsurf = 1 or 2 or 3. Opened and written in subroutine writeout.
!         24 'newDEM_out.asc' -- file of replace DEM with certain scoops removed 
!              depending on value of flag remove. Opened and written in Writeout.
!         26 'fos3d_out.asc' -- file of minimum factors of safety at each
!              DEM column. Opened and written in Writeout.
!         27 'fosvol_out.asc' or 'fosarea_out.asc' -- file of volumes or areas (depending on 
!              primary criteria) for min FOS at each DEM column. Opened and written in 
!              Writeout.
!         28 'searchgrid_out.asc' -- search grid file.  Record of which nodes were searched
!              and whether they fell within the DEM bounds or not. Opened
!              and written in Scoops if isqout = 1.
!         29 'numcols_out.asc' -- file of # of columns for min FOS at each 
!              DEM column. Opened and written in Writeout if isqout = 1.
!         32 'fosretro_out.txt' -- file of retrogression data. Opened in Readin
!              and written in Scoops, Readin and Fos if Remove=A and nretro>0.
!         33 'filter_out.txt' -- file of data describing surfaces which fell outside
!              of user-defined m-alpha and/or FOS limits and nonconverging surfaces.  
!              Opened and written in FOS if filter = 1.
!         34 'ncolerr_out.txt' -- file of information about critical surfaces with
!              fewer columns than user-defined limit, limcol. Opened and
!              written in Fos if colfile = 1.
!         35 'arcs2d_out.okc' -- file of sphere center coordinates, radius, fall angle,
!              and 2D factor of safety of the least stable 2D solution at each DEM column. 
!              Opened and written in Writeout if ifos2d=1.
!         36 'slope_out.asc' -- slope of DEM. Opened and written in Writeout.
!         37 'fos3drel_out.asc' -- file of relative minimum factors of safety at each
!              DEM column. Relative factor of safety is defined as factor of safety divided
!              by the global minimum factor of safety (ofos).  Opened and written in Writeout 
!              if irelfos=1.
!         38 'fos2drel_out.txt' -- file of relative minimum 2D factors of safety at each
!              DEM column. 2D relative factor of safety is defined as factor of safety divided
!              by the global minimum 2D factor of safety (ofos2d).Opened and written in Writeout 
!              if ifos2d=1 and irelfos=1.
!         39 'errors.txt' -- file containing error messages
!         40 'spheresltcut_out.txt' --  file of sphere center coordinates, radius, 
!                fall angle, # of columns in scoop, vol, and factor of safety of the
!                 scoops with F<fostol. Opened in readin and written in fos, only if remove='A'.
!         41 'failsurfslope_out.asc' -- file of defined failure surface slopes
!         42 'failsurfdepth_out.asc' -- file of defined failure surface depths
!         43 'ordfos3d_out.asc' -- file of minimum factors of safety computed by the
!              Ordinary (Fellenius) method at each DEM column. Opened and written in Writeout.
!              if method = 'b'
!         44 'boundcheck_out.asc' -- file of check for boundary limitations in
!              the search lattice. Opened and written in Writeout if isqout = 1
!         45 'critfoslattice_out.3D' -- file of minimum factor of safety   
!              for critical surfaces at each search grid point. Opened and written in Scoops 
!              if icritlattice = 1.
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        USE CommonData
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE GridData, ONLY: nx,ny,nz,delxy,xcen,ycen,zcen,nci1,nci2,nci3,nci4,nci5,zmin,rad,xll,yll,demflag,&
                            radsq,zdem,xllcorner,yllcorner,lengthunits,&
                            halfdelxy,zdemnodes,nbdy
        USE MaterialData, ONLY: nmat,layer
        USE WaterData, ONLY: iwater,piezo,isurfwat,surfwat
        USE StrengthData, ONLY: str3d
        USE SearchData, ONLY: srchfile,ismin,ismax,jsmin,jsmax,nsrchres,irefine,nkseed,&
                            zsmin,zsmax,zsrchres,zavgmax,vmin,vmax,armin,armax,ksmax,&
                            isqout,vacriterion,searchgrid,nsrchpt,zsrchlatt,srchlatt,multres,&
                            isrchmin,jsrchmin,isrchmax,jsrchmax,nxsrchout,nysrchout,searchout,&
                            iincflag,includegrid,irotcen
        USE FOSData, ONLY: remove,ntry,nretro,isubsurf,single,ofos,orad,oangle,omset,&
                           ofsx,ofsy,ofsz,ovol,oarea,ovfos,ovrad,ovangle,ocol,ovmset,&
                           ovfsx,ovfsy,ovfsz,ovvol,ovarea,foscut,fsmin,fsrad,dr,xcount,&
                           fsangle,foso,fsx,fsy,fsz,fsvol,zfos,ozb,mcol,icrit,numdir,mfile,&
                           fsmin2d,foso2d,ifos2d,ofos2d3d,ofos3d2d,fsmin3d2d,fsmin2d3d,&
                           fsx2d,fsy2d,fsz2d,fsrad2d,fsangle2d,fsarea2d,mcol2d,ofos2d,&
                           fsarea,fswidth2d,fsarclength,filtcount,filter,icritlattice,ilattice,linein,&
                           fsminold,irelfos,method,&
!!!     Fellenius comparison addition
                           felfsmin
!!!
        USE SetData, ONLY: xslope,yslope,surfslope,xdem,ydem,nnset,subset,nsetmax,&
                           nrange,in,insphere,zb,area,volume,mini,maxi,minj,maxj,&
                           colfile,xcolcount,zmid,outnodes,colxy,vol,avetop
        USE BishopArrays
        USE FailSurfData, ONLY: failsurf,xfailslope,yfailslope,failsurfslope,ifailsurf,&
                                dArray1,dArray2,dArray1_2d,dArray2_2d
 
        IMPLICIT NONE

        INTEGER :: goodradcnt,badradcnt
        INTEGER :: nstate,newrad,error,nseed,nres,kseed
        INTEGER :: iretro,ierr,icen,jcen,kcen,zcount,goodset
        INTEGER :: i,j,k,nn,m,l,mm,ll,kset,nset,nzres,skip
        INTEGER :: srchct,srchpcnt,srchprint,isrchnum,jsrchnum,ilatt,jlatt
        INTEGER(int2) :: iseed

        REAL(pr) :: rmax,zd,radiniprev,va1max,r2,zdn(4)
        REAL(pr) :: xdiff,ydiff,zavg(nnset),fs1,fs2,fs3,fs4,fs5,fs6,fs7,fs8
        REAL(pr) :: dzdx,dzdy,fgrid,zzsrchres,zlast
        REAL(pr), EXTERNAL :: Newzcen
        
        CHARACTER*8 :: date
        CHARACTER*10 :: time        
        CHARACTER*50 :: version
        CHARACTER*10 :: cnum
        CHARACTER*1 :: charetro 
        CHARACTER*4 :: fosmeth    
        CHARACTER*70 :: problemtype
        CHARACTER*120 :: errmessage,solution
        
        errmessage = ' '
        solution = ' '
        problemtype = 'Scoops3D main program'  

        version = '1.0'
 
!     Read basic input data.
        CALL Readin (version)

!     Allocate arrays and initialize variables and arrays.
        ALLOCATE (xslope(nx,ny),yslope(nx,ny),fsmin(nx,ny),fsangle(nx,ny),STAT=error)
        ALLOCATE (fsrad(nx,ny),fsx(nx,ny),fsy(nx,ny),fsz(nx,ny),fsvol(nx,ny),STAT=error)
        ALLOCATE (icrit(nx,ny),mcol(nx,ny),xdem(nx+1,ny+1),ydem(nx+1,ny+1),STAT=error)
        ALLOCATE (surfslope(nx,ny),fsarea(nx,ny),filtcount(nx,ny),zmid(nx,ny),STAT=error)
        ALLOCATE (subset(nx,ny),in(nx,ny),insphere(nx+1,ny+1),zb(nx+1,ny+1),STAT=error)
        ALLOCATE (outnodes(nx,ny,2),colxy(nx,ny,4),vol(nx,ny),avetop(nx,ny),STAT=error)  
        IF (ifos2d.eq.1) THEN
          ALLOCATE (fsrad2d(nx,ny),fsx2d(nx,ny),fsy2d(nx,ny),fsz2d(nx,ny),STAT=error)
          ALLOCATE (fsarea2d(nx,ny),mcol2d(nx,ny),fsangle2d(nx,ny),STAT=error)
          ALLOCATE (fsmin2d(nx,ny),fsmin3d2d(nx,ny),fsmin2d3d(nx,ny),STAT=error)
          ALLOCATE (fswidth2d(nx,ny),fsarclength(nx,ny),linein(nx,ny),STAT=error)    
        END IF
        IF (irefine.eq.1) THEN
          isrchnum = INT((ismax-ismin)/nsrchres)+1
          jsrchnum = INT((jsmax-jsmin)/nsrchres)+1 
          ALLOCATE (fsminold(nx,ny),srchlatt(isrchnum,jsrchnum,ksmax),STAT=error)
          ALLOCATE (zsrchlatt(isrchnum,jsrchnum,ksmax),STAT=error)
          ALLOCATE (nkseed(isrchnum,jsrchnum),STAT=error)
        END IF
        IF (method.eq.'b'.or.method.eq.'B') THEN
          ALLOCATE (rf(nx,ny),sindip(nx,ny),costruedip(nx,ny),STAT=error)
          ALLOCATE (tanfric(nx,ny),STAT=error)
!!!    Fellenius comparison change
          ALLOCATE (felfsmin(nx,ny),STAT=error)
!!!
         
          IF (ifos2d.eq.1) THEN
            ALLOCATE (dip2d(nx,ny),rf2d(nx,ny),STAT=error)
          END IF         
          IF (ifailsurf.eq.1.or.irotcen.eq.1) THEN
            ALLOCATE (dArray1(nx,ny),dArray2(nx,ny),STAT=error)
            IF (ifos2d.eq.1) ALLOCATE (dArray1_2d(nx,ny),dArray2_2d(nx,ny),STAT=error)
          END IF
        END IF
        IF (ifailsurf.eq.1) THEN
          ALLOCATE (xfailslope(nx,ny),yfailslope(nx,ny),failsurfslope(nx,ny),STAT=error)
        END IF

        IF (error.ne.0) THEN
          errmessage = 'arrays not allocated successfully' 
          solution = 'reduce memory requirements. See "Practical Considerations" chapter of Scoops3D manual'
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')         
        END IF          
        goodradcnt = 0
        badradcnt = 0
        iretro = 0
        halfdelxy = .5_pr*delxy

        isrchmin = 1
        jsrchmin = 1
        isrchmax = nx
        jsrchmax = ny  
!     Allocate 3D array of FOS.
        IF (isubsurf.eq.1.or.isubsurf.eq.2.or.isubsurf.eq.3) THEN
          ALLOCATE (zfos(nx,ny,nz), STAT=error)
          IF (error.ne.0) THEN
            errmessage = 'subsurface factor of safety array not allocated successfully' 
            solution = 'reduce memory requirements - select "isubsurf"=0. See "Practical Considerations" chapter of Scoops3D manual'
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')                   
          END IF    
        END IF
       
!     Allocate search array for searchgrid.out
        IF (isqout.eq.1) THEN
          IF (single.ne.1) THEN
            IF (ismin.lt.1) isrchmin = ismin
            IF (jsmin.lt.1) jsrchmin = jsmin
            IF (ismax.gt.nx) isrchmax = ismax
            IF (jsmax.gt.ny) jsrchmax = jsmax
          END IF            
          nxsrchout=(isrchmax-isrchmin)+1
          nysrchout=(jsrchmax-jsrchmin)+1
          ALLOCATE (searchout(isrchmin:isrchmax,jsrchmin:jsrchmax),STAT=error)
          IF (error.ne.0) THEN
            errmessage = 'search array not allocated successfully' 
            solution = 'reduce memory requirements. See "Practical Considerations" chapter of Scoops3D manual'
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')                   
          END IF 
!     Initialize output search array to -1 if inside DEM boundaries, or -11 if outside.
          DO j = jsrchmin,jsrchmax
            DO i = isrchmin,isrchmax
              IF ((i.ge.1).and.(j.ge.1).and.(i.le.nx).and.(j.le.ny)) THEN
                IF (zdem(i,j).ne.rnull) THEN
                  searchout(i,j) = -1_int2
                ELSE
                  searchout(i,j) = -11_int2
                END IF
              ELSE
                searchout(i,j) = -11_int2
              END IF   
            END DO
          END DO      
          IF (single.eq.1) THEN
           IF (xcen.gt.0) THEN
            icen = INT(xcen/delxy) + 1
           ELSE
            icen = INT(xcen/delxy)
           END IF
           IF (ycen.gt.0) THEN 
            jcen = INT(ycen/delxy) + 1
           ELSE
            jcen = INT(ycen/delxy)
           END IF  
           IF (icen.ge.isrchmin.and.icen.le.isrchmax.and.&
                jcen.le.jsrchmin.and.jcen.ge.jsrchmax)&           
            searchout(icen,jcen) = 2_int2
          END IF
        END IF     

!     Calculate x and y location of nodes and components of slope and other useful column values.       
25      DO j = 1,ny
          DO i = 1,nx
            xdem(i,j) = REAL(i-1,pr)*delxy
            ydem(i,j) = REAL(j-1,pr)*delxy 
            IF (i.eq.nx) THEN
              xdem(nx+1,j) = xdem(nx,j)+delxy
              ydem(nx+1,j) = ydem(nx,j)+delxy
            END IF
            IF (j.eq.ny) THEN
              xdem(i,ny+1) = xdem(i,ny)+delxy
              ydem(i,ny+1) = ydem(i,ny)+delxy
            END IF
            IF (i.gt.1.and.j.gt.1.and.i.lt.nx.and.j.lt.ny) THEN
              IF(nbdy(i,j).ne.1.and.nbdy(i+1,j).ne.1.and.&
                nbdy(i+1,j+1).ne.1.and.nbdy(i,j+1).ne.1) THEN
                xslope(i,j)= (zdem(i-1,j+1)+2.0_pr*zdem(i-1,j)+zdem(i-1,j-1))&
                          -(zdem(i+1,j+1)+2.0_pr*zdem(i+1,j)+zdem(i+1,j-1))
                yslope(i,j)= (zdem(i-1,j-1)+2.0_pr*zdem(i,j-1)+zdem(i+1,j-1))&
                          -(zdem(i-1,j+1)+2.0_pr*zdem(i,j+1)+zdem(i+1,j+1))                           
                dzdx = xslope(i,j)/(8.0_pr*delxy)
                dzdy = yslope(i,j)/(8.0_pr*delxy)
                surfslope(i,j) = ATAN(sqrt(dzdx*dzdx+dzdy*dzdy)) ! in radians
                IF (ifailsurf.eq.1) THEN
                  IF (failsurf(i,j).ne.rnull) THEN
!     If any surrounding discontinuity surface cells are null, set equal to center cell for slope calculation.                   
                    IF (failsurf(i-1,j+1).eq.rnull) THEN
                      fs1 = failsurf(i,j)
                    ELSE
                      fs1 = failsurf(i-1,j+1)
                    END IF
                    IF (failsurf(i-1,j).eq.rnull) THEN
                      fs2 = failsurf(i,j)
                    ELSE
                      fs2 = failsurf(i-1,j)
                    END IF
                    IF (failsurf(i-1,j-1).eq.rnull) THEN
                      fs3 = failsurf(i,j)
                    ELSE
                      fs3 = failsurf(i-1,j-1)
                    END IF                      
                    IF (failsurf(i+1,j+1).eq.rnull) THEN
                      fs4 = failsurf(i,j)
                    ELSE
                      fs4 = failsurf(i+1,j+1)
                    END IF                      
                    IF (failsurf(i+1,j).eq.rnull) THEN
                      fs5 = failsurf(i,j)
                    ELSE
                      fs5 = failsurf(i+1,j)
                    END IF                      
                    IF (failsurf(i+1,j-1).eq.rnull) THEN
                      fs6 = failsurf(i,j)
                    ELSE
                      fs6 = failsurf(i+1,j-1)
                    END IF                      
                    IF (failsurf(i,j-1).eq.rnull) THEN
                      fs7 = failsurf(i,j)
                    ELSE
                      fs7 = failsurf(i,j-1)
                    END IF                      
                    IF (failsurf(i,j+1).eq.rnull) THEN
                      fs8 = failsurf(i,j)
                    ELSE
                      fs8 = failsurf(i,j+1)
                    END IF                      
                    xfailslope(i,j)= (fs1 + 2.0_pr*fs2 + fs3)&
                              - (fs4 + 2.0_pr*fs5 + fs6)
                    yfailslope(i,j)= (fs3 +2.0_pr*fs7 + fs6)&
                              - (fs1 + 2.0_pr*fs8 + fs4)                           
                    dzdx = xfailslope(i,j)/(8.0_pr*delxy)
                    dzdy = yfailslope(i,j)/(8.0_pr*delxy)
                    failsurfslope(i,j) = ATAN(sqrt(dzdx*dzdx+dzdy*dzdy))  ! in radians
                  ELSE
                    xfailslope(i,j) = 0.0_pr
                    yfailslope(i,j) = 0.0_pr
                    failsurfslope(i,j) = 0.0_pr
                  END IF
                END IF
              ELSE
                xslope(i,j)= 0.0_pr
                yslope(i,j)= 0.0_pr
                surfslope(i,j) = rnull
                IF (ifailsurf.eq.1) THEN
                  xfailslope(i,j)= 0.0_pr
                  yfailslope(i,j)= 0.0_pr
                  failsurfslope(i,j) = rnull
                END IF
              END IF
            ELSE
              xslope(i,j)= 0.0_pr
              yslope(i,j)= 0.0_pr
              surfslope(i,j) = rnull
              IF (ifailsurf.eq.1) THEN
                xfailslope(i,j)= 0.0_pr
                yfailslope(i,j)= 0.0_pr
                failsurfslope(i,j) = rnull
              END IF
            END IF     
          END DO
        END DO
        
!     If printing mimimum FOS at each search lattice node, open output file
!     and write header information.                    
        IF (ilattice.eq.1.or.icritlattice.eq.1) THEN
          IF (iretro.gt.0) THEN    
            IF (ilattice.eq.0) THEN         
               OPEN (15,STATUS='scratch')     
            ELSE
              OPEN (15,STATUS='replace',file = outputdir(1:LEN_TRIM(outputdir))//&
                   filin(bfile:nfile)//'_foslattice'//charetro//'_out.3D')             
           END IF    
          ELSE     
            IF (ilattice.eq.0) THEN              
               OPEN (15,STATUS='scratch')           
            ELSE    
                OPEN (15,STATUS='replace',file = outputdir(1:LEN_TRIM(outputdir))//&
                    filin(bfile:nfile)//'_foslattice_out.3D')           
            END IF    
          END IF
          IF (method.eq.'b') then 
             fosmeth = 'Bish'
          ELSE
             fosmeth = 'Ord '
         END IF   
          IF (LEN_TRIM(lengthunits).eq.1) THEN
             WRITE (15,1210) lengthunits,lengthunits,lengthunits,fosmeth
          ELSE
            IF (LEN_TRIM(lengthunits).eq.2) THEN
             WRITE (15,1215) lengthunits,lengthunits,lengthunits,fosmeth
            ELSE
             WRITE (15,1220) fosmeth
            END IF  
         END IF  
         
        END IF
       
!     If using retrogression, write number to 'fosretro.out' if remove=A.
!     Set variables for retrogression file names
        IF (nretro.gt.0.AND.single.ne.1) THEN
          IF (remove.eq.'A') WRITE (32,1000) iretro
          cnum = '0123456789'
          charetro = cnum(iretro+1:iretro+1)              
        END IF 

!     Initialize variables needed in each retrogression
        fsmin = nullhi
!!!   Fellenius comparison change
        IF (method.eq.'b'.or.method.eq.'B') felfsmin = nullhi
!!!

        fsarea = rnull
        fsangle = rnull
        fsrad = rnull
        fsx = rnull
        fsy = rnull
        fsz = rnull
        ofos = nullhi
        ofos2d = nullhi
        orad = 0.0_pr
        oangle = 0.0_pr
        ofsx = 0.0_pr
        ofsy = 0.0_pr
        ofsz = 0.0_pr
        omset = 0
        ovol = 0.0_pr
        oarea = 0.0_pr
        ocol = 0
        ovfos = nullhi
        ovrad = 0.0_pr
        ovangle = 0.0_pr
        ovmset = 0
        ovfsx = 0.0_pr
        ovfsy = 0.0_pr
        ovfsz = 0.0_pr
        ovvol = vmin
        ovarea = 0.0_pr
        nrange = 0
        srchct = 0
        xcount = 0
        xcolcount = 0 
        filtcount = 0
        IF (ifos2d.eq.1) THEN
          fsmin2d = nullhi
          fsmin2d3d = nullhi
          fsmin3d2d = nullhi
          fsangle2d = rnull
          fsrad2d = rnull
          fsx2d = rnull
          fsy2d = rnull
          fsz2d = rnull
          fsarea2d = rnull
          mcol2d = inull
          fswidth2d = rnull
          fsarclength = rnull
          linein = 0.0_pr
        END IF
        fsvol = rnull
        icrit = -1
        mcol = inull
        mfile = 0
        colfile = 0
        nsetmax = 0
        ntry = 0
        IF (isubsurf.eq.1.or.isubsurf.eq.2.or.isubsurf.eq.3) zfos = nullhi
        in = 0
        subset = 0
        insphere = 0
        iseed = 0_int2
        nseed = 0            

        IF (irefine.eq.1) THEN  ! Initial coarse search  
           nres = multres * nsrchres
           zzsrchres = REAL(multres,pr) * zsrchres
           DO i=1,isrchnum
            DO j=1,jsrchnum
             DO k=1,ksmax
               srchlatt(i,j,k) = 0_int2
             END DO
            END DO
           END DO         
           fsminold = rnull
           nzres = multres
        ELSE
           zzsrchres = zsrchres  
           nres = nsrchres   
           nzres = 1 
        END IF

!     If not doing single run for factors of safety
        IF (single.ne.1) THEN 
          PRINT *,filin(bfile:LEN_TRIM(filin))," - Starting search using Scoops3D"                  
! Do until restol is reached at every DEM cell or no new search seeds are produced
          DO      
            srchct = 0 
            srchprint = 0

!     If doing coarse to fine iterative search          
            IF (irefine.eq.1.and.iseed.ge.1_int2) THEN    
              multres = max(int(multres/2),1)
              nres = multres * nsrchres
              nzres = multres
              zzsrchres = REAL(nzres,pr)*zsrchres
              CALL SeedSearch(iseed,nseed,nres,nzres,zzsrchres,isrchnum,jsrchnum)
              IF (nseed.eq.0) EXIT    
            END IF
                      
            IF (irefine.eq.1.and.iseed.ne.0_int2) nsrchpt = nseed  
            IF (iseed.eq.0_int2) THEN
              skip = nres 
            ELSE
              skip = nsrchres
            END IF              
!     Loop through search grid searching for minimum FOS.
            DO jcen = jsmin,jsmax,skip
              DO icen = ismin,ismax,skip             
!     Check if using search file and then whether icen,jcen is a valid location.  
!     If using search lattice refinement method, search file only applies to
!     coarse search.        
                IF (srchfile.eq.1.and.iseed.eq.0_int2) THEN
                  IF (searchgrid(icen,jcen).lt.1) CYCLE
                END IF
!     If doing lattice refinement, calculate indices of search lattice arrays.                
                IF (irefine.eq.1) THEN
                  ilatt = INT((icen-ismin)/nsrchres)+1
                  jlatt = INT((jcen-jsmin)/nsrchres)+1
!     Check whether at a valid i,j search location.                  
                  IF (iseed.gt.1_int2) THEN
                    IF (nkseed(ilatt,jlatt).lt.1) CYCLE  
                  END IF
                END IF

!     Calculate percent of search completed.              
                IF (srchct.gt.0) THEN
                  srchpcnt = srchct*100/nsrchpt
                  srchpcnt = INT(srchpcnt/10)*10
                    
                  IF (srchpcnt.gt.srchprint.or.(srchpcnt.eq.10.and.srchprint.eq.0)) THEN
!     In cases with a search grid and search refinement, the value for
!       srchct is approximate and can result in a value of srchpcnt > 100
                    IF (srchpcnt.gt.100) srchpcnt = 100                   
                    srchprint = MAX(srchprint + 10,srchpcnt)                                    
                    IF (irefine.eq.0) THEN
                      PRINT &
                        "(A20,' - Search node:',i6,',',i6,'; search grid ,',i4,' % completed,')",&
                          filin(bfile:LEN_TRIM(filin)),icen,jcen,srchpcnt
                      PRINT "('      ',i10,' trial surfaces analyzed')",ntry     
                    ELSE
                      IF (iseed.eq.0_int2) THEN
                        PRINT &
                           "(A20,' - Search node:',i6,',',i6,'; coarse search ,',i4,' % completed,')",&
                           filin(bfile:LEN_TRIM(filin)),icen,jcen,srchpcnt
                        PRINT "('      ',i10,' trial surfaces analyzed')",ntry        
                      ELSE
                        PRINT &
                           "(A20,' - Search node:',i6,',',i6,'; fine search # ',i3,' ,',i4,' % completed,')",&
                           filin(bfile:LEN_TRIM(filin)),icen,jcen,iseed-1_int2,srchpcnt
                        PRINT "('      ',i10,' trial surfaces analyzed')",ntry        
                      END IF
                    END IF
                  END IF 
                END IF
                IF (irefine.eq.1.and.iseed.gt.1_int2) THEN
                  srchct = srchct + nkseed(ilatt,jlatt)
                  kseed = 1
                ELSE      
                  srchct = srchct + 1
                END IF
                           
!     Find starting elevation of search node.
                IF (iseed.gt.1_int2) THEN
                  zcen = REAL(zsrchlatt(ilatt,jlatt,kseed),pr)
                ELSE
                  zcen = zsmin
                END IF   
                   
!     If above DEM and ground elevation is higher than zsmin, 
!     find starting search elevation at first search lattice node above DEM.
                IF (icen.gt.0.and.jcen.gt.0.and.icen.le.nx.and.jcen.le.ny) THEN
                  IF (zdem(icen,jcen).ne.rnull) THEN

!     If creating output search file, set array value equal to 2 if within
!     DEM boundaries, or 22 if outside boundaries or zdem = rnull.
                    IF (isqout.eq.1) searchout(icen,jcen) = 2_int2 
                    IF  (zsmin.lt.zdem(icen,jcen).and.iseed.eq.0_int2) THEN
                      nn = INT((zdem(icen,jcen)-zsmin)/zzsrchres) + 1
                      zcen = zsmin + (REAL(nn,pr)*zzsrchres)                     
                    END IF
                  ELSE
                    IF (isqout.eq.1) searchout(icen,jcen) = 22_int2
                  END IF  ! IF (zdem(icen,jcen).ne.rnull)
                ELSE
                  IF (isqout.eq.1) searchout(icen,jcen) = 22_int2
                END IF  ! IF (icen.gt.0.and.jcen.gt.0.and.icen.le.nx.and.jcen.le.ny)
                
!     Set flag for no prior radius estimate at current search i,j location.                
                radiniprev = -1      
!     Find x and y coordinates of search node based on DEM column (1,1) lower left  
!     node at x=0., y=0. Assume search node is at center of DEM cell.        
                xcen = REAL(icen-1,pr)*delxy + halfdelxy
                ycen = REAL(jcen-1,pr)*delxy + halfdelxy 
                zlast = zcen 
                fgrid = nullhi
!    Loop over range of zcen (search lattice elevations).                
                DO       
                  kcen = NINT((zcen-zsmin)/zsrchres) + 1                                                     
                  
!    Keep track of minimum factor of safety at each search lattice
!    node and write to file if ilattice = 1
                  IF (ilattice.eq.1.or.icritlattice.eq.1) THEN
!    If at a new search lattice point, write minimum fos data for last
!    grid point to file.
                    IF (zcen.ne.zlast) THEN 
                      WRITE (15,1110) xcen+xll,ycen+yll,zlast,fgrid
                      IF (fgrid.ge.nullhi-1.0_pr) radiniprev = -1                        
                      fgrid = nullhi 
                    END IF      
                  END IF
!    If higher than maximum search elevation move to next point.                
                  IF (zcen.gt.zsmax) EXIT 
                  IF (irefine.eq.1.and.iseed.eq.0_int2) srchlatt(ilatt,jlatt,kcen) = -1_int2                  
                  zlast = zcen           
                  rmax = 0.0_pr
!     Calculate maximum radius to furthest border.
                  DO m = 1,2
                    DO l = 1,2
                      mm = (m-1)*(nx-1) + 1
                      ll = (l-1)*(ny-1) + 1
                      xdiff = xcen-(REAL(mm-1,pr)*delxy+halfdelxy)
                      ydiff = ycen-(REAL(ll-1,pr)*delxy+halfdelxy)                  
!     If DEM has null value at corners, approximate elevation value as zmin
!     (the minimum DEM elevation) to calculate a maximum radius.
                      IF (zdem(mm,ll).eq.rnull) THEN
                        r2 = xdiff*xdiff + ydiff*ydiff + &
                            (zcen-zmin)*(zcen-zmin)
                      ELSE
                        r2 = xdiff*xdiff + ydiff*ydiff + & 
                            (zcen-zdem(mm,ll))*(zcen-zdem(mm,ll))
                      END IF                
                      rmax = MAX(rmax,r2)
                    END DO
                  END DO

                  rmax = SQRT(rmax)                

!    Initialize variables used to calculate initial radius guess.
                  ierr = 0
                  newrad = 1
                  rad = -1.0_pr
                  nset = 1      
                  IF ((icen.lt.1.or.icen.gt.nx).or.(jcen.lt.1.or.jcen.gt.ny)) THEN
                    zd = -888.0_pr
                  ELSE
                    zd = zdem(icen,jcen)
                  END IF
                  zavgmax = rnull
                  va1max = 0.0_pr
            
!     Calculate radius of trial sphere.                 
100               CALL Radius (zd,newrad,va1max,radiniprev,rmax,ierr,goodradcnt,badradcnt)
                  radsq = rad*rad  
                  
                        
!     If radius is too large and hit max number of iterations, move to next point.
                  IF (ierr.eq.1.and.rad.gt.rmax) THEN
                    IF (iseed.gt.1_int2) kseed = kseed + 1
                    zcen = Newzcen(ilatt,jlatt,zcen,zzsrchres,iseed,kseed)
                    CYCLE  !  zcen DO loop
                  END IF

!     Loop over range of radii                  
                  DO  
!     Find all sets of contiguous DEM columns that are intersected
!     by the trial sphere.
                    CALL Subsets (newrad,nset,kset,zavg,ierr)
                    
!     If more subsets than parameter nnset skip search node.
                    IF (ierr.eq.1) THEN
                      IF (iseed.gt.1_int2) kseed = kseed + 1
                      zcen = Newzcen(ilatt,jlatt,zcen,zzsrchres,iseed,kseed)
                      EXIT  ! leave radius DO loop   
                    END IF

!     If no intersected sets recalculate radius.
                    IF (kset.eq.0) THEN
                      va1max = 0.0_pr
                      go to 100
                    END IF
            
!     Calculate volumes and areas of sets, and check whether they meet
!     criteria.  
                    CALL Volumes(nset,goodset)

!     If all intersected sets thrown out in volumes recalculate radius.
                    IF (goodset.eq.0) THEN
                      va1max = 0.0_pr
                      go to 100
                    END IF               

                    IF ((vacriterion.eq.'V').or.(vacriterion.eq.'v')) THEN
                      CALL Criteria (nset,newrad,zavg,volume,area,vmin,vmax,armin,armax,&
                                     va1max,nstate)            
                    ELSE           
                      CALL Criteria (nset,newrad,zavg,area,volume,armin,armax,vmin,vmax,&
                                     va1max,nstate)            
                    END IF
  
                    IF (nstate.eq.0) THEN                 
!     Calculate factor of safety for all sections within correct
!     volume interval and keep minimum factor of safety stored.        
                       CALL Fos (nset,iretro,fgrid)
                    ELSE  ! if need to recalculate radius and have not hit max # tries.
                      IF (nstate.eq.100.and.newrad.ne.0) go to 100
                    END IF
                    rad = rad + dr
                    radsq = rad*rad
                      
!     If no good sets possible with increasing radius based on currently intersected set
!     details or radius exceeds rmax, move to next search node. 
!!!   If it is desired to continue searching to rmax, simply delete the (nstate.eq.500)                    
                    IF ((nstate.eq.500).or.(rad.gt.rmax)) THEN
                      IF (iseed.gt.1_int2) kseed = kseed + 1
                      zcen = Newzcen(ilatt,jlatt,zcen,zzsrchres,iseed,kseed)
                      EXIT  ! leave radius DO loop
                    END IF
                  END DO !  loop over radii            
                END DO !  (zcen.lt.zsmax)                        
              END DO  ! loop on icen of search grid
            END DO  ! loop on jcen of search grid
           
            IF (irefine.eq.0) THEN
            PRINT *, filin(bfile:LEN_TRIM(filin)),' - search , ',100.,'% complete,',ntry,' trial surfaces'
              EXIT
            ELSE
              IF (iseed.eq.0_int2) THEN
                PRINT *, filin(bfile:LEN_TRIM(filin)),' - coarse search , ',100.,'% complete,',ntry,' trial surfaces'
                iseed = 1_int2
              ELSE
                PRINT "(A,' - fine search # ',i3,', ',f4.0,'% complete,',i10,' trial surfaces')",&
                     filin(bfile:LEN_TRIM(filin)),iseed-1_int2,100.,ntry
              END IF
            END IF
           

          END DO  ! iteration loop on seeds   
 
        ELSE   
!     If option chosen to calculate FOS for single failure option  
             
          radsq = rad*rad
          CALL Subsets (newrad,nset,kset,zavg,ierr)
          CALL Volumes (nset,goodset) 
          CALL Fos (nset,iretro,fgrid)          
          IF ((ilattice.eq.1.or.icritlattice.eq.1)) WRITE (15,1110) xcen+xll,ycen+yll,zcen,ofos
        END IF

!     Create a new DEM with scoop removed depending on options chosen.
!     If all scoops below FOS cutoff are to be removed, array ozb is assigned
!     in suboutine FOS as scoops are found to be below FOS cutoff. 
        IF (ofos.ne.nullhi) THEN 
          IF (remove.ne.'A'.and.remove.ne.'N') THEN
            subset = 0
            IF (remove.eq.'L') THEN
!     Assign search center variables values for largest scoop under FOS cutoff.
              xcen = ovfsx
              ycen = ovfsy
              zcen = ovfsz
              rad = ovrad
              radsq = rad*rad
              m = ovmset
            ELSE
              IF (remove.eq.'M') THEN
!     Assign search center variables for overall minimum FOS scoop.
                xcen = ofsx
                ycen = ofsy
                zcen = ofsz
                rad = orad
                radsq = rad*rad
                m = omset
              END IF
            END IF
!     Assign removed scoop base elevation values to array zmid using chosen search center.
            CALL Subsets (newrad,nset,kset,zavg,ierr)
            CALL SlipBase(m,mini(m),maxi(m),minj(m),maxj(m))
            DO j = minj(m),maxj(m)
              DO i = mini(m),maxi(m)
                IF (subset(i,j).eq.m) ozb(i,j) = zmid(i,j)
              END DO
            END DO
          END IF 
        END IF
!     Place new DEM in original DEM array zdem.        
        IF (demflag.eq.1) zdem = ozb        
 
        CALL Writeout (iretro)              

!     Check to see if more retrogressions are needed and
!     if min FOS is less than foscut and continue if so.
        IF (nretro.gt.0.and.iretro.lt.nretro.and.ofos.lt.foscut.and.single.ne.1) THEN
          iretro = iretro + 1
!     Average 4 surrounding cells to find node elevations for continuous 
!     DEM surface approximation      
          DO j = ny+1,1,-1
            DO i = 1,nx+1
              zcount = 0
              zdn = 0
              IF (i.gt.1) THEN
                IF (j.gt.1) THEN
                  IF (zdem(i-1,j-1).ne.rnull) THEN
                    zdn(1) = zdem(i-1,j-1)
                    zcount = zcount + 1
                  END IF
                END IF
                IF (j.lt.ny+1) THEN
                  IF (zdem(i-1,j).ne.rnull) THEN
                    zdn(4) = zdem(i-1,j)
                    zcount = zcount + 1
                  END IF
                END IF
              END IF
              IF (i.lt.nx+1) THEN
                IF (j.gt.1) THEN
                  IF (zdem(i,j-1).ne.rnull) THEN
                    zdn(2) = zdem(i,j-1)
                    zcount = zcount + 1
                  END IF
                END IF
                IF (j.lt.ny+1) THEN
                  IF (zdem(i,j).ne.rnull) THEN
                    zdn(3) = zdem(i,j)
                    zcount = zcount + 1
                  END IF
                END IF
              END IF
              IF (zcount.gt.0) zdemnodes(i,j) = SUM(zdn)/REAL(zcount,pr)         
            END DO
          END DO
          GOTO 25
        END IF                

        CLOSE(15)
        CLOSE (16)
        CLOSE (20)
        CLOSE (32)
        CLOSE (39)
        CLOSE(45)

        CALL DATE_AND_TIME(date,time) 
             
        PRINT *, filin(bfile:LEN_TRIM(filin)),' - Successful execution of Scoops3D version number ',version
        PRINT "('Date and Time: ',A2,'/',A2,'/',A4,'  ',A2,':',A2,':',A2)",&
                 date(5:6),date(7:8),date(1:4),time(1:2),time(3:4),time(5:6)
                 
!       PRINT *, 'Number of bad initial volumes: ',badradcnt
!        PRINT *, 'Number of radius changes: ',goodradcnt
!        PRINT *, 'Number of trial surfaces: ',ntry

1000    FORMAT ('Retrogression # = ',i5,/3x,'FOS','        ',&
                 &'x','          ','y','           ','z','       ',&
                 &'rad','        ','angle','     ','vol','     ','area') 
1110    FORMAT (3f30.4,f11.4)  
1210    FORMAT ('x_',A1,' y_',A1,' z_',A1,' F_',A4)                
1215    FORMAT ('x_',A2,' y_',A2,' z_',A2,' F_',A4)                
1220    FORMAT ('x y z F_',A4)     

                 	
        END PROGRAM Scoops3D
        
        
        
             
        FUNCTION Newzcen(ilatt,jlatt,zcen,zzsrchres,iseed,kseed)
      
          USE CommonData
          USE SearchData, ONLY: nkseed,zsmax,zsrchlatt
          IMPLICIT NONE
        
          INTEGER, INTENT(in) :: kseed,ilatt,jlatt
          INTEGER(int2), INTENT(in) :: iseed
          REAL(pr), INTENT(in) :: zcen,zzsrchres
          REAL(pr) :: Newzcen
        
          IF (iseed.gt.1_int2) THEN
            IF (kseed.le.nkseed(ilatt,jlatt)) THEN
              Newzcen = REAL(zsrchlatt(ilatt,jlatt,kseed),pr)
            ELSE
              Newzcen = zsmax + zzsrchres
            END IF
          ELSE          
            Newzcen = zcen + zzsrchres
          END IF
        
        END FUNCTION Newzcen
