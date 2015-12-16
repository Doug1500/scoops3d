MODULE CommonData
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Used in all subroutines to control precision of all real variables.
!
!        VARIABLES
!
!        bfile - location in character string of first character of 
!           input file name prefix, defined by searching from the
!           end of the string for the occurence of a slash character 
!        double -- represents double precision
!        filin -- name of primary input file of parameters
!        int2 -- 2 byte integer
!        inull -- integer null value used throughout program
!        nfile -- location in character string of last character of 
!           input file name prefix, defined by searching from the 
!           end of the string for the occurence of a period
!        nullhi -- positive high null value used throughout program
!        outputdir -- optional directory path to place output files into
!        pr -- precision of all declared real variables in program
!        rnull --  real null value used throughout program
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  SAVE
  INTEGER, PARAMETER :: double = SELECTED_REAL_KIND(15)
  INTEGER, PARAMETER :: pr = double
  INTEGER, PARAMETER :: int2 = SELECTED_INT_KIND(4)
  INTEGER :: nfile,bfile
  INTEGER, PARAMETER :: inull=-9999
  REAL(pr), PARAMETER :: rnull=-9999.0_pr,nullhi=9999.0_pr
  CHARACTER*200 :: filin,outputdir

END MODULE CommonData


MODULE GridData

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Used in subroutines readin, readdem, readpiezo, readpressh, readstrat,
!  radius, subsets, volumes, fos, bishop, fellenius, checkinsphere, 
!  checksets, newnodes, readstrength, and writeout.
!
!        VARIABLES
!
!        angle -- assumed angle of slip (azimuth)
!        cellarea -- area of DEM grid cell (delxy*delxy)
!        delxy -- DEM grid resolution (delta x, delta y)
!        delz -- z resolution of search grid
!        demflag -- indicates whether file replacedem will be generated.
!        halfdelxy -- 1/2 DEM grid spacing
!        lengthunits -- units of length, used for labels in output files 
!        nbdy(nx+1,ny+1) -- array of flags for identifying nodes adjacent to 
!           DEM boundary cells
!        nx -- number of DEM cells in x direction
!        ny -- number of DEM cells in y direction
!        nz -- number of nodes used for 3-D FOS values when isubsurf=1or2
!        rad -- radius of search sphere
!        radsq -- square of search sphere radius
!        xcen,ycen,zcen -- center of search sphere relative to DEM origin
!        xll,yll -- x and y origin of DEM grid, read in header lines
!        xllcorner,yllcorner -- character version of these header lines
!        zdem(nx,ny) -- DEM elevations
!        zdemnodes(nx+1,ny+1) -- DEM node elevations averaged from srrounding cells
!        zmin -- minimum DEM elevation
!        zmax -- maximum DEM elevation
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  USE CommonData
  IMPLICIT NONE
  SAVE
  
  INTEGER :: nx,ny,nz,demflag
  
  REAL(pr) :: delxy,halfdelxy,xcen,ycen,zcen,zmin,delz
  REAL(pr) :: rad,radsq,angle,zmax,xll,yll
  REAL(pr) :: xcenrot,ycenrot,zcenrot,cellarea
  
  CHARACTER*2 :: lengthunits
  CHARACTER*60 :: xllcorner,yllcorner
  
  REAL(pr), PARAMETER :: pi=3.141592653589793_pr
  
  INTEGER, ALLOCATABLE :: nbdy(:,:)
  REAL(pr), ALLOCATABLE :: zdem(:,:),zdemnodes(:,:)
END MODULE GridData



MODULE MaterialData

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Used in subroutines readin, readstrat, fos, bishop, and fellenius. 
!
!        VARIABLES
!
!        activelay(nmat-1) -- number of active (nonnull) cells in each layer
!        cee(nmat) -- cohesion of each layer
!        duwtlaydz(nx,ny,nmat-1) -- gradient of depth-weighted unit weights between layers
!        eq -- earthquake acceleration coefficient
!        failsurf(nx,ny) -- array of user-defined failure surface elevations at each DEM cell
!        fsurf -- flag for whether failure surface was found,1=yes,0=no
!        fxa(nmat) -- curve-fitting parameter, a from the equation for the soil-water characteristic 
!            curve defined by Fredlund and Xing (1994) 
!        fxn(nmat) -- curve-fitting parameter, n from the equation for the soil-water characteristic 
!            curve defined by Fredlund and Xing (1994) 
!        fxm(nmat) -- curve-fitting parameter, m from the equation for the soil-water characteristic 
!            curve defined by Fredlund and Xing (1994) 
!        fxr(nmat) -- curve-fitting parameter, psi_r (soil suction associated with residual moisture 
!            content) from the equation for the soil-water characteristic curve defined by 
!            Fredlund and Xing (1994)
!        gamr(nmat,3) -- array of total, partially saturated, and saturated
!           unit weights for each layer
!        gamw -- water unit weight in units of problem 
!        gamsurf -- unit weight of surface load in units of problem
!        gsameall -- flag indicating whether all unit weights are equivalent
!        gsameeach -- flag indicating whether unit weights in each layer are equivalent
!           to each other
!        layer(nmat,nx,ny) -- array of bottom elevations for material layers
!        maxlayer(nmat-1) --maximum value of each stratigraphic layer file
!        minlayer(nmat-1) -- min value of each stratigraphic layer file
!        nmat -- number of material layers
!        ru(nmat) -- pore pressure ratio approximation
!        tanphi(nmat) -- tangent of friction angle for each layer
!        thetares(nmat) -- residual water content for each layer
!        thetasat(nmat) -- saturated water content for each layer
!        uwtlay(nx,ny,nmat-1) -- depth-weighted average unit weight at each layer
!           bottom with correct moist for water conditions
!        vga(nmat) -- curve-fitting parameter, alpha from the equation for the 
!           soil-water characteristic curve defined by vanGenuchten (1980) 
!        vgn(nmat) -- curve-fitting parameter, n, from the equation for the soil-water 
!            characteristic curve defined by vanGenuchten (1980) 
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  USE CommonData
  IMPLICIT NONE
  SAVE
  
  INTEGER :: gsameall,gsameeach
  INTEGER :: nmat,fsurf
  INTEGER, ALLOCATABLE :: activelay(:)
  
  REAL(pr) :: gamw,eq,gamsurf
  
  REAL(pr), ALLOCATABLE :: layer(:,:,:),uwtlay(:,:,:),cee(:),tanphi(:)
  REAL(pr), ALLOCATABLE :: gamr(:,:),ru(:),duwtlaydz(:,:,:),thetasat(:)
  REAL(pr), ALLOCATABLE :: minlayer(:),maxlayer(:),failsurf(:,:),thetares(:)
  REAL(pr), ALLOCATABLE :: fxa(:),fxn(:),fxm(:),fxr(:),vga(:),vgn(:)
END MODULE MaterialData

  
  
MODULE WaterData

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Used in subroutines readin, readpiezo, readpressh, readstrat, fos, 
!  bishop, and fellenius.
!
!        VARIABLES
!
!        delzp -- vertical resolution of pressure head data
!        dpdz(nx,ny,pmaxk) -- array of pressure head gradients
!        dsedz(nx,ny,pmaxk) -- array of effective saturation gradients
!        dthetazdz(nx,ny,pmaxk) -- array of volumetric water content gradients
!        iwater -- flag indicating method for modeling water pressure.
!           0=no water pressures, 1=ru approximation, 2= piezometric surface, 
!           3=input 3-d pressure file, 4=3D variably saturated file containing
!           pressure head and moisture content, 5=3D variably saturated file,
!           moisture content from vanGunuchten SWCC, 6=3D variably saturated file,
!           moisture content from Fredlund and Xing SWCC
!        isurfwat -- flag for analysis with surface water layer
!        maxpk(nx,ny) -- array of highest k value of pressure data at each DEM cell
!        moist -- determines which element of soil unit weight array to use: 1-total weight
!           2-partially saturated, 3-fully saturated.
!        pgrid -- flag for pressure head grid, 1=regular, 2=irregular grid
!        piezo(nx,ny) -- array of piezometric elevation at each DEM cell
!        pmaxk -- maximum number of 3-d pressure values at each cell
!        pressh(nx,ny,pmaxk) -- array of pressure head data
!        presshpz(nx,ny,pmaxk) -- array of z elevations of pressure head data
!        pzsurf -- flag for whether piezometric surface was found,1=yes,0=no
!        se(nx,ny,pmaxk) -- array of effective saturation 
!        surfwat(nx,ny) -- array of surface water elevations at each DEM cell
!        thetaz(nx,ny,pmaxk) -- array of volumetric water content
!        zpmin -- minimum elevation of pressure head values
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  USE CommonData
  IMPLICIT NONE
  SAVE
  
  INTEGER :: iwater,pzsurf,pgrid,isurfwat
  INTEGER :: pmaxk,moist
  INTEGER, ALLOCATABLE :: maxpk(:,:)
  
  REAL(pr):: delzp,zpmin
  
  REAL(pr), ALLOCATABLE :: pressh(:,:,:),thetaz(:,:,:),surfwat(:,:)
  REAL(pr), ALLOCATABLE :: piezo(:,:),dpdz(:,:,:),presshpz(:,:,:),dthetazdz(:,:,:)
  REAL(pr), ALLOCATABLE :: se(:,:,:),dsedz(:,:,:)
END MODULE WaterData

  

MODULE StrengthData

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Used in subroutines readin, readstrength, bishop, and fellenius.
!
!        VARIABLES
!
!        cflag -- indicates whether 3-d cohesion data are provided; 0=no,1=yes
!        ceect -- total number of cohesion values read
!        ceeunits -- units for cohesion, used for labels in output files
!        cohes(nx,ny,strmaxk) -- array of 3-d cohesion data
!        dcohdz(nx,ny,strmaxk) -- array of cohesion gradients between vertical nodes
!        dflag -- indicates whether 3-d unit weight data are provided; 0=no,1=yes
!        dfricdz(nx,ny,strmaxk) -- array of friction gradients between vertical nodes
!        duwtdz(nx,ny,strmaxk) -- array of unit weight gradients between vertical nodes 
!        dzstr -- z spacing of 3-d strength data
!        fflag -- indicates whether 3-d friction data are provided; 0=no,1=yes
!        fricct -- total number of friction angles read
!        linterp -- linear interpolation flag; 1=use linear interpolation, other=use
!          nearest node below slip surface base for strength data.
!        maxstrk(nx,ny) -- array of highest k value of strength data at each DEM cell
!        mincee,maxcee -- overall min and max cohesion
!        minfric,maxfric -- overall min and max friction angles
!        minuwt(3),maxuwt(3) -- overall min and max unit weights for each saturation
!        str3d -- flag for using 3-d strength file, 1=3D strengths
!        strgrid -- flag for strength data grid; 1=regular, 2=irregular grid
!        strmaxk -- number of 3-d strength values at each cell
!        strz(nx,ny,strmaxk) -- array of z locations of irregular strength data
!        tfric(nx,ny,strmaxk) -- array of 3-d friction data
!        uwt(nx,ny,strmaxk,3) -- array of unit weights
!        uwtct -- total number of unit weights read
!        uwtunits -- units for unit weight, used for labels in output files
!        uwt3d(nx,ny,strmaxk) -- unit weight at each cell (depth weighted average)
!          with correct moist for water conditions
!        zstrmin -- minimum elevation of strength data
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  USE CommonData
  IMPLICIT NONE
  SAVE
  
  INTEGER :: str3d,cflag,fflag,dflag,strgrid,linterp
  INTEGER :: ceect,fricct,uwtct,strmaxk
  INTEGER, ALLOCATABLE :: maxstrk(:,:)
  
  REAL(pr) :: zstrmin,dzstr
  REAL :: mincee,maxcee,minfric,maxfric,minuwt(3),maxuwt(3)
  
  REAL(pr), ALLOCATABLE :: cohes(:,:,:),tfric(:,:,:)
  REAL(pr), ALLOCATABLE :: strz(:,:,:),dcohdz(:,:,:),dfricdz(:,:,:),duwtdz(:,:,:)
  REAL(pr), ALLOCATABLE :: uwt3d(:,:,:),uwt(:,:,:,:)
  
  CHARACTER*8 ceeunits, uwtunits
  
END MODULE StrengthData
 


MODULE SearchData

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Used in subroutines readin, radius, criteria, subsets, and fos.
!
!        VARIABLES
!
!        allin -- flag to indicate if include columns are contained by slip surface
!        armin,armax -- minimum and maximum area bounds for valid slip surface
!        fostol -- maximum percent difference in FOS between seed iterations to
!           not generate a new seed
!        goal -- vmin or armin + half of tol
!        iincflag -- flag to require an area that must be included in all potential failure
!           surfaces, defined by an inclusion area file, 1 = inclusion area file provided
!        includegrid -- Defines columns of the DEM that must be included in all
!           potential failure surfaces (1 if included, 0 if not)
!        iincludemin,iincludemax -- minimum and maximum i range of include columns
!        irefine -- flag for implementing coarse to fine search iterations
!        irotcen -- flag indicating rotational center was specified for single surface
!        ismin,ismax -- minimum and maximum boundaries of search grid on x-axis,
!           specified relative to DEM columns i=1 through nx
!        isqout -- flag for creating search quality output files, 1=create files
!        isrchmin,isrchmax -- array bounds of searchout output file.
!        jincludemin,jincludemax -- minimum and maximum j range of include columns
!        jsmin,jsmax -- minimum and maximum boundaries of search grid on y-axis,
!           specified relative to DEM rows j=1 through ny
!        jsrchmin,jsrchmax -- array bounds of searchout output file
!        ksmax -- maximum boundary of search grid on z-axis calculated from 
!           user-specified zsmax and z resolution
!        minclude -- set number that contains the inclusion area
!        multres -- resolution multiplier for coarse search lattice
!        nincpt -- number of active cells in include file
!        nkseed(nxsrch,nysrch) -- number of k values for current seed at each i,j
!           search node
!        nsrchpt -- total number of search points
!        nsrchres -- resolution of search grid relative to DEM resolution.  Acts
!           as a multiplier of DEM resolution delxy
!        nxsrchout,nysrchout -- number of search nodes in x and y directions in output
!           search grid file, includes intersection of search range with DEM range.
!           only, 0 = primary and secondary criteria
!        searchgrid - array to specify search grid points
!           which grid points were searched and whether they fell inside or outside
!           the DEM.
!        searchout -- output array of area covered by search grid and indicating
!           which grid points were searched and whether they fell inside or outside
!           the DEM.
!        srchfile -- flag indicating whether using search file to specify valid
!           search grid nodes, 1=search file provided
!        srchlatt(nxsrch,nysrch,nzsrch+1) -- search lattice array with seed
!           iteration number at each searched node. 
!        tol -- tolerance on primary control criterion for calculating initial 
!           radius.  (i.e. initial volume must fall between volume and volume+tol)
!        vacriterion -- indicates primary control for intial radius (v or a
!           for volume or area)
!        vmin,vmax -- minimum and maximum volume bounds for valid slip surface
!        zavgmax -- approximation used to help find initial radius estimate
!        zsmin,zsmax -- minimum and maximum elevations of search grid
!        zsrchlatt(nxsrch,nysrch,nzsrch) -- search lattice array with 
!           elevations of new search centers generated from seed.
!        zsrchres -- resolution of search grid on z-axis
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  USE CommonData
  IMPLICIT NONE
  SAVE
  
  CHARACTER*1 :: vacriterion
  
  INTEGER :: ismin,ismax,jsmin,jsmax,nsrchpt
  INTEGER :: ksmax,isrchmin,jsrchmin,isrchmax,jsrchmax
  INTEGER :: nxsrchout,nysrchout,nsrchres,multres
  INTEGER :: srchfile,isqout,irefine,iincflag,allin,minclude,nincpt,irotcen
  INTEGER :: iincludemin,iincludemax, jincludemin,jincludemax
  INTEGER, ALLOCATABLE :: nkseed(:,:),searchgrid(:,:),includegrid(:,:)
  INTEGER(int2), ALLOCATABLE :: srchlatt(:,:,:),searchout(:,:)
  
  REAL(pr) :: armin,armax,vmin,vmax,tol,zavgmax,zsmin,zsmax,zsrchres,goal,fostol
  REAL, ALLOCATABLE :: zsrchlatt(:,:,:)


END MODULE SearchData


MODULE FOSData

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Used in subroutines readin, fos, and writeout.
!
!        VARIABLES
!        absminma -- abs. value of minimum m-alpha allowed in Bishop FOS calculation
!        arclength -- arc length of current 2d slice
!        area2d -- cross-sectional area of current 2d slice
!        deginc -- degree increment for additional slip directions to search for minimum FOS
!        degmax -- degrees to search on either side of estimated
!           fall line (describes half width of search segment)
!        dr -- radius increment for search
!        diter -- convergence requirement in terms of percent difference between
!           last two iterations for factor of safety using Bishop's Method 
!        filtcount(nx,ny) -- counter of number of times each cell was filtered by
!           minma.
!        filter -- indicates whether solution filters are applied (malpha and fos)
!        ifos2d -- flag for calculating 2-d factors of safety
!        foscut -- factor of safety cutoff for removing slip surface from DEM
!        foslocal -- flag for writing local factors of safety to file, only
!           available for single failure option
!        foso(numdir) -- factor of safety of current scoop at each slip angle
!        foso2d(numdir) -- 2-d factor of safety of current scoop at each slip angle
!        fsangle(nx,ny) -- array of fall angles of minimum scoop at each DEM cell
!        fsangle2d(nx,ny) -- array of fall angles of minimum 2d slice at each DEM cell
!        fsarclength(nx,ny) -- array of arclengths of minimum 2d slice at each DEM cell
!        fsarea2d(nx,ny) -- array of areas of minimum 2d slice at each DEM cell
!        fsmin(nx,ny) -- array of minimum FOS at each DEM cell
!        fsminold(i,j) -- last minimum factor of safety at each DEM cell used to
!           compare change in FOS during seed iterations to see if criterion is met
!        fsmin2d(nx,ny) -- 2-d minimum factor of safety at each DEM cell
!        fsmin2d3d(nx,ny) -- 3-d minimum factor of safety associated with min 2d slice
!           at each DEM cell
!        fsmin3d2d(nx,ny) -- 2-d minimum factor of safety associated with min 3d scoop
!           at each DEM cell
!        fsrad(nx,ny) -- array of radius of minimum scoop at each DEM cell
!        fsrad2d(nx,ny) -- array of radius of minimum 2d slice at each DEM cell
!        fsvol(nx,ny) -- volume of minimum FOS scoop at each DEM cell
!        fswidth2d(nx,ny) -- array of widths of minimum 2d slice at each DEM cell
!        fsx(nx,ny),fsy(nx,ny),fsz(nx,ny) -- search sphere center of minimum FOS 
!           scoop at each DEM cell
!        fsx2d(nx,ny),fsy2d(nx,ny),fsz2d(nx,ny) -- search sphere center of minimum FOS 
!           2d slice at each DEM cell
!        icritlattice -- flag for creating file critfoslattice_out.3D for minimum FOS of 
!          critical surfaces at each search  grid node
!        ilattice -- flag for creating file foslattice_out.3D for minimum FOS at each search
!           grid node
!        irelfos -- flag for creating relative factor of safety grids,1=create files
!        isubsurf -- flag for creating 3-d file of FOS below DEM surface, 1=create file&
!          1=create file in ijk format, 2=create file in xyz format
!        icrit(nx,ny) -- array of volume or area flags indicating proximity to limits
!            of criteria range for each DEM cell
!        limcol -- optimal number of columns per slip surface.  Fewer columns
!           will generate error in ncolerr.out file
!        linein(nx,ny) -- array of flags for whether 2-D slide directon line intersects
!           column
!        mcol(nx,ny) -- number of columns in minimum FOS scoop at each DEM cell
!        mcol2d(nx,ny) -- number of columns in minimum 2d FOS slice at each DEM cell
!        mfile -- flag indicating whether filter file was opened
!        method -- B for Bishops Simplified, F for Fellenius method for
!           calculation of factor of safety.
!        mindip -- minimum apparent dip of slip surface column base
!        minma -- min absolute value of m-alpha of all iterations for a surface
!        ncol2d -- number of columns in current 2-D slice
!        nretro -- number of failure retrogressions to compute
!        numdir -- number of search slip directions for a single failure
!        ntry -- total number of slip surfaces for which FOS is calculated
!        oangle -- slip angle associated with overall minimum factor of safety
!        oangle2d -- slip angle associated with overall minimum 2d factor of safety
!        oarclength -- arc length of overall minimum 2d factor of safety
!        oarclength3d2d -- arc length of 2d surface associated with overall min 3d scoop
!        oarea -- area of overall minimum scoop
!        oarea2d -- area of overall minimum 2d slice
!        oarea2d3d -- area of 3d surface associated with overall min 2d slice
!        oarea3d2d -- area of 2d surface associated with overall min 3d scoop
!        ocol -- number of columns in overall minimum scoop
!        ocol2d -- number of columns in overall minimum 2-D scoop
!        ocol2d3d -- number of columns in 3d surface associated with min 2d slice
!        ocol3d2d -- number of columns in 2d surface associated with min 3d scoop
!        ofos -- overall minimum factor of safety
!        ofos2d -- 2-d FOS associated with overall minimum 3-d factor of safety
!        ofos2d3d -- FOS of 3d slice associated with overall minimum 2d FOS
!        ofos3d2d -- FOS of 2d slice associated with overall minimum 3d FOS
!        ofsx,ofsy,ofsz -- search sphere center of minimum FOS scoop
!        ofsx2d,ofsy2d,ofsz2d -- search sphere center of minimum 2d FOS scoop
!        ominma -- overall min absolute value of m-alpha of all iterations for a surface
!        omset -- set number of overall minimum scoop
!        orad -- radius of sphere associated with overall minimum factor of safety
!        orad2d -- radius of arc associated with overall minimum 2d factor of safety
!        osliparea -- slip surface area of overall min 3d scoop
!        osliparea2d3d -- slip surface area of 3d slice associated with overall minimum 2d FOS
!        ovangle -- slip angle of largest volume scoop with FOS < cutoff
!        ovarea -- area of largest volume scoop with FOS < cutoff
!        ovmset -- set number of largest volume scoop with FOS < cutoff
!        ovfos -- factor of safety of largest volume scoop with FOS < cutoff
!        ovrad -- radius of largest volume scoop with FOS < cutoff
!        ovfsx,ovfsy,ovfsz -- search sphere center of largest volume scoop with 
!           FOS < cutoff
!        ovol, oarea - volume and area of minimum  FOS scoop
!        ovol2d3d -- volume of 3d surface associated with overall min 2d slice
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
!           removed. M= the scoop with the minimum FOS will be removed,
!           N= no scoops removed.
!        single -- flag for calculating single slip surface
!        width2d -- width of current 2d slice
!        xcount -- number of surfaces eliminated due to m-alpha limits and/or nonconvergence
!        zbot -- elevation of bottom of 3-d fos grid array
!        zfos(nx,ny,nz) -- 3-d array of minimum factor of safety below DEM surface
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  USE CommonData
  IMPLICIT NONE
  SAVE
  
  CHARACTER*1 :: remove,method
  INTEGER :: icritlattice,ilattice,single,foslocal
  INTEGER :: filter,isubsurf,irelfos,mfile,ifos2d
  INTEGER :: nretro,omset,ovmset
  INTEGER :: numdir,limcol,xcount,ocol2d,ncol2d,ocol3d2d,ocol2d3d,ocol,ntry
  INTEGER, ALLOCATABLE :: linein(:,:),icrit(:,:)
  INTEGER, ALLOCATABLE :: mcol(:,:),mcol2d(:,:),filtcount(:,:)
  
  REAL(pr) :: dr,degmax,deginc,ofos,foso,ovvol,foscut,zbot,diter,ofos2d
  REAL(pr) :: absminma,minma,mindip,arclength,area2d,width2d,foso2d
  REAL(pr) :: orad,ovrad,ovfsx,ovfsy,ovfsz,ofsx,ofsy,ofsz
  REAL :: ominma,ovwt,ovsliparea,osliparea2d3d,owt2d3d
  REAL :: ofos2d3d,ofos3d2d,oarea3d2d,oarclength3d2d
  REAL :: owidth3d2d,oarea2d3d,ovol2d3d,owidth2d,osliparea,owt
  REAL :: oangle,ovarea,ovol,oarea,ovangle
  REAL :: ovfos,orad2d,oangle2d,ofsx2d,ofsy2d,ofsz2d,oarea2d,oarclength
  
  REAL(pr), ALLOCATABLE :: fsmin(:,:),fsminold(:,:),fsmin2d(:,:),zfos(:,:,:),ozb(:,:)
  REAL, ALLOCATABLE :: fsx(:,:),fsy(:,:),fsz(:,:),fsvol(:,:),fsarclength(:,:)
  REAL, ALLOCATABLE :: fsrad2d(:,:),fsangle2d(:,:),fswidth2d(:,:)
  REAL, ALLOCATABLE :: fsx2d(:,:),fsy2d(:,:),fsz2d(:,:),fsmin2d3d(:,:),fsmin3d2d(:,:)
  REAL, ALLOCATABLE :: fsarea2d(:,:)
  REAL, ALLOCATABLE :: fsrad(:,:),fsangle(:,:),fsarea(:,:)
!!! Added Fellenius variables  
  REAL, ALLOCATABLE :: felfsmin(:,:)
  REAL :: felfoso,felminfos
END MODULE FOSData


MODULE SetData

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Used in subroutines subsets, volumes, criteria, fos, bishop, fellenius
! checkinsphere, checksets, and newnodes.
!
!        VARIABLES
!
!        area(nnset) -- surface area of subset
!        avetop(nx,ny) -- average elevation of 4 corner nodes of intersected 
!          column calculated using the average of the four surrounding cells
!          for full columns, and by finding the intersection of the slip surface
!          with the smoothed DEM surface for partial columns.
!        colfile -- flag for opening error file reporting cols < limcol
!        colxy(nx,ny) -- length of each column side (delxy for full columns)
!        in(nx,ny) -- number of nodes of each column bounded by sphere
!        insphere(nx,ny) -- indicates which DEM column nodes are bounded by sphere
!        mini(nnset),maxi(nnset) -- minimum and maximum i bounds of each subset
!        minj(nnset),maxj(nnset) -- minimum and maximum j bounds of each subset
!        ncol(nnset) -- number of columns in subset
!        nnset -- maximum number of subsets
!        nrange(nnset) -- array of flags for whether set fits volume or area criteria 
!           range, 1=fits criteria range
!        nsetmax -- flag for exceeding number of subsets allowed nnset
!        omini,omaxi,ominj,omaxj -- minimum and maximum bounds of all sets included
!           in search sphere
!        outnodes(nx,ny,2) -- corner location of column nodes which fall
!           outside the search sphere.
!        set(nnset) -- number of columns in a contiguous set
!        setflag(nnset) -- indicates valid subset (0), or invalid by containing
!           truncated node (1), or adjacent to DEM boundary (2), or empty (-1)
!        sliparea -- area of slip surface
!        subset(nx,ny) -- indicates set membership of DEM columns
!        surfslope(nx,ny) -- slope of dem surface
!        vol(nx,ny) -- volume of column in subset
!        volume(nnset) -- volume of subset
!        weight -- weight of slide
!        xcolcount -- number of critical surfaces with total # of columns < limcol
!        xdem(nx+1,ny+1),ydem(nx+1,ny+1) -- x and y locations of DEM column nodes
!        xslope(nx,ny),yslope(nx,ny) -- x and y slope vectors of DEM surface
!        zb(nx,ny) -- base of slip surface at each column node
!        zmid(nx,ny) -- base of slip surface at midpoint of intersected column
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  USE CommonData
  IMPLICIT NONE
  SAVE
  
  INTEGER, PARAMETER :: nnset=10
  
  INTEGER :: omini,omaxi,ominj,omaxj,xcolcount
  INTEGER :: nsetmax
  INTEGER :: mini(nnset),minj(nnset),maxi(nnset),maxj(nnset),set(nnset),ncol(nnset)
  INTEGER :: nrange(nnset),setflag(nnset),colfile
  INTEGER, ALLOCATABLE :: outnodes(:,:,:),subset(:,:),insphere(:,:),in(:,:)
  
  REAL(pr),ALLOCATABLE :: xslope(:,:),yslope(:,:),surfslope(:,:),zb(:,:),avetop(:,:)
  REAL(pr),ALLOCATABLE :: xdem(:,:),ydem(:,:),zmid(:,:),vol(:,:),colxy(:,:,:)
  REAL(pr) :: area(nnset),volume(nnset),sliparea,weight
  
END MODULE SetData

MODULE FailSurfData

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Used in subroutines readfailsurf, subsets, volumes, criteria, fos, bishop, fellenius
! checkinsphere, checksets, and newnodes.
!
!        VARIABLES
!
!        failsurf(nx,ny) -- array of failure surface elevations for each DEM column    
!        failsurfdepth(nx,ny) -- array of failure surface depths for each DEM column        
!        failsurfslope(nx,ny) -- slope of input failure surface
!        ifailsurf -- flag for whether failure surface file is used
!        xfailslope(nx,ny),yfailslope(nx,ny) -- x and y slope vectors of input failure surface surface
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  USE CommonData
  IMPLICIT NONE
  SAVE
  
  INTEGER :: ifailsurf
  
  REAL(pr),ALLOCATABLE :: xfailslope(:,:),yfailslope(:,:),failsurfslope(:,:),failsurf(:,:)
  REAL(pr),ALLOCATABLE :: dArray1(:,:),dArray2(:,:),failsurfdepth(:,:)
  REAL(pr),ALLOCATABLE :: dArray1_2d(:,:),dArray2_2d(:,:) 
END MODULE FailSurfData

MODULE ResultsStats

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Used in subroutine writeout
!
! Minimums and maximums of results arrays
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  USE CommonData
  IMPLICIT NONE
  SAVE
  
  INTEGER :: numdatapts,minmcol,maxmcol,minfiltcount,maxfiltcount,minmcol2d,&
             maxmcol2d

  REAL :: minslope,maxslope,minfsmin,maxfsmin,minfsrelmin,maxfsrelmin,&
              minfsvol,maxfsvol,minfsarea,maxfsarea,minzdem,maxzdem,&
              minfsmin2d,maxfsmin2d,minfsrelmin2d,maxfsrelmin2d,minfsmin3d2d,&
              maxfsmin3d2d,minfsmin2d3d,maxfsmin2d3d, minfsz2d,maxfsz2d,&
              minfsx2d,maxfsx2d,minfsy2d,maxfsy2d,minfsrad2d,maxfsrad2d,&
              minfsangle2d,maxfsangle2d,minfsarclength,maxfsarclength,&
              minfswidth2d,maxfswidth2d,minfsz,maxfsz,minfsx,maxfsx,minfsy,&
              maxfsy,minfsrad,maxfsrad,minfsangle,maxfsangle,minfsarea2d,&
              maxfsarea2d,mincol2d,maxcol2d,minfelfsmin,maxfelfsmin,&
              minfailslope,maxfailslope,minfaildepth,maxfaildepth
  
END MODULE ResultsStats


MODULE BishopArrays

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Used in subroutine Bishop
!
!        costruedip(nx,ny) -- cosine of the true dip of the sliding base at the
!           center of each column.
!        dip2d -- Dip in the slide direction at the center of each 2d column.
!        rf(nx,ny) -- resisting force of each column
!        rf2d(nx,ny) -- 2D resisting force of each column
!        sindip(nx,ny) -- sine of the dip of the sliding base in the direction of 
!           the slide at the center of each column.
!        tanfric(nx,ny) -- friction angle at the base of each column. 
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  USE CommonData
  IMPLICIT NONE
  SAVE
  
  REAL(pr),ALLOCATABLE :: rf(:,:),sindip(:,:),costruedip(:,:),rf2d(:,:)
  REAL(pr),ALLOCATABLE :: tanfric(:,:),dip2d(:,:)

END MODULE BishopArrays
