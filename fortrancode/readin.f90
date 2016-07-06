        SUBROUTINE readin (version)
	
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!	This subroutine reads primary input file for Scoops and calls all 
!       other data input subroutines. Input data is echoed to output file.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!       Called by Scoops3D
!
!	      VARIABLES
!
!       absminma -- abs. value of minimum m-alpha allowed in FOS calculation
!       activelay(nmat-1) -- number of active (nonnull) cells in each layer
!       angle -- assumed angle of slip (azimuth)
!       armin,armax -- minimum and maximum area bounds for valid slip surface&
!       cee(nmat) -- cohesion of each layer
!       ceect -- total number of cohesion values read in 3-D file
!       ceeunits -- units for cohesion, used for labels in output files
!       cflag -- indicates whether 3-d cohesion data are provided (0=yes,1=no)
!       ddeg -- degree interval to search for minimum fos
!       deldeg -- degrees to search on either side of estimated
!       fall line (describes half width of search segment)
!       dr -- radius increment for search
!       delxy -- DEM grid resolution (delta x, delta y)
!       delz -- z resolution of search grid
!       demfile -- dem file name
!       demflag -- indicates whether file newDEM will be generated.
!       dflag -- indicates whether 3-d unit weight data are provided (0=yes,1=no)
!       diter -- convergence tolerance between last two iterations for factor of safety 
!           using Bishop's Method 
!       duwtlaydz(nx,ny,nmat-1) -- gradient of depth-weighted unit weights between layers
!       eq -- earthquake acceleration coefficient. 
!       fflag -- indicates whether 3-d friction data are provided (0=yes,1=no)
!       filter -- intdicates whether solution filters are applied (malpha)
!       foscut -- factor of safety cutoff for removing slip surface from DEM
!       foslocal -- flag for writing local factors of safety to file, only
!           available for single failure option
!       fostol -- maximum percent difference in FOS between seed iterations to
!           not generate a new seed
!       fricct -- total number of friction angles read in 3-D file
!       fxa(nmat) -- curve-fitting parameter, a from the equation for the soil-water characteristic 
!            curve defined by Fredlund and Xing (1994) 
!       fxn(nmat) -- curve-fitting parameter, n from the equation for the soil-water characteristic 
!            curve defined by Fredlund and Xing (1994) 
!       fxm(nmat) -- curve-fitting parameter, m from the equation for the soil-water characteristic 
!            curve defined by Fredlund and Xing (1994) 
!       fxr(nmat) -- curve-fitting parameter, psi_r (soil suction associated with residual water 
!            content) from the equation for the soil-water characteristic curve defined by 
!            Fredlund and Xing (1994)
!       gamr(nmat,3) -- array of total, partially saturated, and saturated
!           unit weights for each layer
!       gamw -- water unit weight in units of problem 
!       gamsurf -- unit weight of surface load in units of problem  
!       goal -- vmin or armin + half of tol
!       gsameall -- flag indicating whether all unit weights are equivalent
!       gsameeach -- flag indicating whether unit weights in each layer are equivalent
!       icritlattice -- flag for creating file critfoslattice_out.3D for minimum FOS of 
!           critical surfaces at each search  grid node
!       ifos2d -- flag for calculating 2-d factors of safety
!       isupported -- function to indicate feature is supported in current version of SCOOPS
!           isupported = 0 indicates that feature is not supported,1 indicates that feature is supported
!       iincflag -- flag that requires an area be included in all potential failure
!           surfaces, defined by an inclusion area file, 1 = inclusion area file provided
!       ilattice -- flag for creating file foslattice_out.3D for minimum FOS at each search
!           grid node
!       inull -- integer null value used throughout program
!       irefine -- flag for implementing coarse to fine search iterations
!       irotcen -- flag indicating rotational center was specified for single surface
!       isflag -- output generation flag for whether isqout line was given
!       ismin,ismax -- minimum and maximum boundaries of search grid on x-axis,
!           specified relative to DEM columns i=1 through nx
!       isqout -- flag for creating search quality output files
!       isubsurf -- flag for creating 3-d file of FOS below DEM surface
!           1=create file in ijk format, 2=create file in xyz format,
!           3=create file in .vtk format
!       isurfwat -- flag for surface water layer
!       iv -- flag for volume control of scoop size, = 1 if 
!           volume is primary criterion, 0 if not
!       iwater -- flag indicating method for modeling water pressure.
!           0 = no water pressures, 1 = ru approximation, 2 = piezometric surface, 
!           3=input 3-d pressure file, 4=3D variably saturated file containing
!           pressure head and water content, 5=3D variably saturated file,
!           water content from vanGenuchten SWCC, 6=3D variably saturated file,
!           water content from Fredlund and Xing SWCC
!       jsmin,jsmax -- minimum and maximum boundaries of search grid on y-axis,
!           specified relative to DEM rows j=1 through ny
!       ksmax -- maximum number of search grid nodes on z-axis
!       layer(nmat,nx,ny) -- array of bottom elevations for material layers
!       lnum(nmat) -- layer identification number
!       lengthunits -- units of length, used for labels in output files
!       limcol -- optimal number of columns per slip surface.  Fewer columns
!           will generate error in ncolerr.out file
!       linterp -- linear interpolation flag; 1=use linear interpolation, other-use
!           nearest node for strength data.
!       maxlayer(nmat-1) --maximum value of each stratigraphic layer file
!       maxuwt(3) -- maximum unit weights for 3-D strengths file
!       method -- B for Bishops Simplified, O for Ordinary (Fellenius) method for
!           calculation of factor of safety.
!       mincee,maxcee -- overall min and max cohesion
!       minfric,maxfric -- overall min and max friction angles
!       minlayer(nmat-1) -- min value of each stratigraphic layer file
!       minuwt(3) -- minimum unit weights for 3-D strengths file
!       multres -- resolution multiplier for coarse search lattice
!       ndir -- number of search slip directions for each failure, input by user
!       nmat -- number of material layers
!       nretro -- number of failure retrogressions to compute
!       numdir -- number of search slip directions for each failure, calculated by SCOOPS
!       nsrchpt -- total number of search points
!       nsrchres -- resolution of search grid relative to DEM resolution.  Acts
!           as a multiplier of DEM resolution delxy
!       nxsrch,nysrch -- number of search nodes in x and y directions
!       nx -- number of DEM columns
!       ny -- number of DEM rows
!       nz -- number of nodes used for 3-D FOS values when isubsurf=1, 2 or 3&
!           only, 0 = primary and secondary criteria
!       outputdir -- optional directory path to place output files into
!       phi(nnmat) -- friction angle of each material in stratigraphy
!       pcount -- number of DEM cells with piezometric surface found
!       pnum -- number of horizontal cells with pressure values
!       pzsurf -- flag for whether piezometric surface was found,1=yes,0=no
!       rad -- radius of search sphere
!       remove -- A= all scoops < FOS cutoff will be removed from new
!           DEM.  L= only largest volume scoop < FOS cutoff will be 
!           removed. M= the scoop with the minimum FOS will be removed,
!           N= no scoops removed.
!       rnull -- negative real null value used for initializations and set in CommonData.
!       ru(nmat) -- pore pressure ratio approximation
!       single -- flag for calculating single slip surface
!       srch-- 'fi'=search file,'si'=single set of parameters,'bo'=search box specified
!       str3d -- flag for using 3-d strength file
!       strnum -- number of horizontal cells with strength data given
!       surfwfile -- surface water file name
!       tanphi(nmat) -- tangent of friction angle for each layer
!       thetares(nmat) -- residual  water content for each layer
!       thetasat(nmat) -- saturated  water content for each layer
!       title -- title of primary input file. Should be first line of input file.
!       tol -- tolerance on primary control criterion for calculating initial 
!           radius.  (i.e. initial volume must fall between volume and volume+tol)
!       uwtct -- total number of unit weights read in 3-D file
!       uwtlay(nx,ny,nmat-1) -- depth-weighted average unit weight at each layer
!           bottom with correct moist for water conditions
!       uwtunits -- units for unit weight, used for labels in output files
!       vacriterion -- indicates primary control for intial radius (v or a
!           for volume or area)
!       version - scoops version number and date of last modification
!       vga(nmat) -- curve-fitting parameter, alpha from the equation for the 
!           soil-water characteristic curve defined by vanGenuchten (1980) 
!       vgn(nmat) -- curve-fitting parameter, n, from the equation for the soil-water 
!            characteristic curve defined by vanGenuchten (1980) 
!       vmin,vmax -- minimum and maximum volume bounds for valid slip surface
!       water --  pressure head method. 'NO'=none,'RU','PZ','3D', or 'VS'
!       xcen,ycen,zcen -- center of search sphere relative to DEM origin
!       xcenrot,ycenrot,zcenrot -- rotational center of slip surface relative to DEM origin
!       xll,yll -- x and y origin of DEM grid, read in header lines
!       zbot -- elevation of bottom of 3-d fos grid array
!       zsmin,zsmax -- minimum and maximum elevations of search grid
!       zsrchres -- resolution of search grid on z-axis
!       zfrac -- fraction of delxy to make delz for 3D FOS output
!       zmin -- minimum DEM elevation
!       zmax -- maximum DEM elevation
!
!
!      INPUT FILES
!       unit filename
!         12 'inputfilename.#' -- file of basic input information
!
!      OUTPUT FILES
!       unit filename
!         20 'inputfilename_out.txt' -- echo of input and results containing
!              overall minimum slip surface data, written in subroutines
!              Readin, Readpiezo, Readpressh, Readsearch, Readstrat,
!              Readstrength, and Writeout.
!         32 'fosretro_out.txt' -- file of retrogression data. Opened in Readin
!              and written in Scoops, Readin and Fos if Remove=A and nretro>0.
!         39 'errors_out.txt' -- file containing error messages
!         40 'spheresltcut_out.txt' --  file of sphere center coordinates, radius, 
!              fall angle, # of columns in scoop, vol, and factor of safety of the
!              scoops with F<fostol. Opened in readin and written in fos, only if remove='A'.
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        USE CommonData
        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE GridData, ONLY: nx,ny,nz,delxy,zmin,delz,xll,yll,xcen,ycen,zcen,nci1,nci2,nci3,nci4,nci5,&
                           rad,angle,pi,zmax,demflag,zdem,lengthunits,xcenrot,ycenrot,zcenrot
        !!!!!!!!!!!!!!!!!!!!!!!!!!
        USE MaterialData, ONLY: nmat,gsameall,gsameeach,gamw,eq,cee,tanphi,gamr,gamsurf,ru,&
                           thetares,thetasat,layer,minlayer,maxlayer,activelay,uwtlay,duwtlaydz,&
                           fxa,fxn,fxm,fxr,vga,vgn
        USE WaterData, ONLY: iwater,pzsurf,piezo,zpmin,delzp,pmaxk,pgrid,pressh,moist,&
                           isurfwat,surfwat,thetaz
        USE StrengthData, ONLY: cflag,ceeunits,fflag,dflag,str3d,zstrmin,dzstr,strmaxk,&
                           linterp,strgrid,cohes,minuwt,maxuwt,mincee,maxcee,minfric,&
                           maxfric,uwtct,uwtunits,ceect,fricct
        USE SearchData, ONLY: vacriterion,srchfile,ismin,ismax,jsmin,jsmax,nsrchres,goal,&
                           isqout,armin,armax,vmin,vmax,tol,zsmin,zsmax,zsrchres,nsrchpt,&
                           ksmax,irefine,multres,fostol,iincflag,includegrid,irotcen
        USE FOSData, ONLY: remove,method,diter,numdir,limcol,irelfos,icritlattice,ilattice,&
                           isubsurf,nretro,single,dr,degmax,deginc,foscut,zbot,foslocal,absminma,&
                           ifos2d,filter
        USE FailSurfData, ONLY: ifailsurf,failsurf,failsurfdepth
 
        IMPLICIT NONE

        INTEGER :: n,ios,ierr,ios2,error,nxsrch,nysrch,isflag,pcount
        INTEGER :: strnum,pnum,nsrch,nsrchptini,ndir
        INTEGER, ALLOCATABLE :: gsame(:),lnum(:)
        LOGICAL  :: isupported
        
        REAL(pr) :: zfrac,ditertemp
        REAL(pr), ALLOCATABLE :: phi(:)

        CHARACTER*50, INTENT(in) :: version
        CHARACTER*8 :: date
        CHARACTER*6 :: srch
        CHARACTER*60 :: heading,units
        CHARACTER*220 :: demfile,filout,prfile,pzfile,sfile,stratfile,surfwfile
        CHARACTER*220 :: strgthfile,junk,errfilout,failsurffile,includefile
        CHARACTER*10 :: time
        CHARACTER*2 :: water
        CHARACTER*120 :: title
        CHARACTER*3 :: strcoords,prcoords
        CHARACTER*4 :: fosmeth        
        CHARACTER*70 :: problemtype
        CHARACTER*120 :: errmessage,solution
        
        errmessage = ' '
        solution = 'check Scoops3D manual for description of required parameters'
        problemtype = 'reading input file'  
        n = 1   
        nretro = 0
        ierr = 0
        nsrchres = inull
        zsrchres = rnull
        ismin = inull
        ismax = inull
        jsmin = inull
        jsmax = inull
        zsmin = rnull
        zsmax = rnull        
        dr = rnull
        ndir = inull
        numdir = inull
        deginc = rnull
        degmax = rnull
        zfrac = rnull
        armin = rnull
        armax = rnull
        vmin = rnull
        vmax = rnull
        vacriterion = 'X'
        foslocal = 0
        icritlattice = 0
        ilattice = 0 
        isubsurf = 0
        irelfos = 0
        tol = rnull
        limcol = inull
        linterp = 0
        nmat = 0
        method = 'X'
        iwater = 0
        gamw = 0.0_pr
        gamsurf = 0.0_pr
        eq = rnull
        diter = 0.0001
        ditertemp = rnull
        single = 0
        xcen = rnull
        ycen = rnull
        zcen = rnull
        rad = rnull
        !!!!!!!!!!!!!!!!!!!!!!!!
        nci1 = rnull
        nci2 = rnull
        nci3 = rnull
        nci4 = rnull
        nci5 = rnull
        !!!!!!!!!!!!!!!!!!!!!!!!
        angle = rnull
        isqout = 0
        str3d = 0
        pzsurf = 0
        ios = 1
        ios2 = 1
        cflag = 1
        fflag = 1
        dflag = 1 
        absminma = 0.0_pr
        ifos2d = 0
        demflag = 0
        filter = 0
        isflag = 0
        irefine = 0
        multres = 1
        fostol = rnull
        remove = 'N'
        demfile = 'X'
        filout = 'X'
        errfilout = 'X'
        prfile = 'X'
        pzfile = 'X'
        surfwfile = 'X'
        sfile = 'X'
        stratfile = 'X'
        strgthfile = 'X'
        failsurffile = 'X'
        includefile = 'X'
        outputdir = ''
        title = 'X'
        srch = 'X'
        water = 'X'
        lengthunits = '  '
        ceeunits = '        '
        uwtunits = '        '
        ifailsurf = 0
        irotcen = 0
        iincflag = 0
        isurfwat = 0
        
        PRINT *
        PRINT *, 'Executing Scoops3D'
        PRINT *

        DO
          PRINT *, 'Input file name?'
          
          ! READ (*,1000,IOSTAT=ios) filin
          ! OPEN (12,STATUS = 'old',FILE = filin,IOSTAT=ios2)

          filin = '/home/yewintun/mygo/src/scoops3d/test16.scp'        
          OPEN (12,STATUS = 'old',FILE = filin,IOSTAT=ios2)

          ios = 0

          IF (ios.eq.0.and.ios2.eq.0) EXIT
          n=n+1
          IF (n.gt.5) THEN
            PRINT *,'*** ERROR - Scoops3D  execution terminated'          
            PRINT *, '*** error opening input file ',filin
            STOP
          END IF
        END DO

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!   Read main input file.    
        PRINT *, 'Reading input file: ',filin        

!  Search input file for ouput directory path
        DO
          READ (12,1000,IOSTAT=ios) heading
          IF (ios.ne.0) EXIT
          IF (heading(1:2).eq.'ou'.or.heading(1:2).eq.'OU') THEN
           if (heading(1:6).eq.'output'.or.heading(1:6).eq.'OUTPUT') THEN
            READ(12,1000,IOSTAT=ios2) outputdir         
            IF (LEN_TRIM(outputdir).gt.220) THEN
              errmessage = 'Scoops3D allows a maximum of 220 characters for path names'
              Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')              
            END IF         
            IF (ios2.ne.0) THEN
              errmessage = 'error in output directory name'
              solution = 'check that output directory is a valid path name'
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')             
              ierr = 1 
              CYCLE
            END IF   
           END IF   
          END IF                           
         END DO
         
!  Count the number of characters in the input file name before the file extension
        nfile = SCAN(filin,'.',BACK=.true.)
        IF (nfile.eq.0) THEN
          nfile = LEN_TRIM(filin)
        ELSE 
          nfile=nfile-1
        END IF

!  Count the number of characters in the input file name before the file path name
        bfile = SCAN(filin,'/',BACK=.true.)
        IF (bfile.eq.0) bfile = SCAN(filin,'\',BACK=.true.)
        bfile=bfile+1

!  Concatenate output directory path to output file name and open output file
        outputdir = ADJUSTL (outputdir) 
        IF (outputdir.ne.''.and.&
            outputdir(LEN_TRIM(outputdir):LEN_TRIM(outputdir)).ne.'/'.and.&
            outputdir(LEN_TRIM(outputdir):LEN_TRIM(outputdir)).ne.'\') THEN
          ierr = 1
          errmessage = 'error in output directory name'
          solution = 'check that output path name ends with a slash'
          Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                     
        END IF      
        filout = outputdir(1:LEN_TRIM(outputdir))//filin(bfile:nfile)//'_out.txt'
        errfilout = outputdir(1:LEN_TRIM(outputdir))//filin(bfile:nfile)//'_errors_out.txt'
        OPEN (20,STATUS='replace',FILE = filout,IOSTAT=ios)
        OPEN (39,STATUS='replace',FILE = errfilout,IOSTAT=ios)
        IF (ios.ne.0) THEN
          errmessage = 'error opening output files'
          solution = 'check for existence of output directory'
          Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')          
          PRINT *,'sumary output file name',filout
          PRINT *,'error output file name',errfilout
          WRITE (20,*) 'sumary output file name',filout
          WRITE (20,*) 'error output file name',filout,errfilout        
          STOP
        END IF           
       
        CALL DATE_AND_TIME(date,time)

!  Rewind input file to obtain other variables
        REWIND(12) 

!  Read in other variables in input file
        DO 
          READ (12,1000,IOSTAT=ios) heading
          solution = 'check Scoops3D manual for description of required parameters'                  
          IF (ios.ne.0) EXIT
          SELECT CASE (heading(1:2))
          
          CASE ('ti','TI')
            if (heading(1:5).eq.'title'.or.heading(1:5).eq.'TITLE') READ (12,1000,IOSTAT=ios2) title
        
!   Search grid and size control parameters
!       First line should be definition of search option
          CASE ('sr','SR')
            IF (heading(1:4).eq.'srch'.or.heading(1:4).eq.'SRCH') THEN
              READ (12,*,IOSTAT=ios2) srch
              IF (ios2.ne.0) THEN            
                errmessage = 'invalid or missing "srch" parameter'
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                                   
                ierr = 1
              END IF 
            END IF
            
          IF (srch(1:6).eq.'single'.or.srch(1:6).eq.'SINGLE') THEN
!       If single parameters                            
              READ (12,*)
              READ(12,*,IOSTAT=ios2) xcen,ycen,zcen,rad,angle,nci1,nci2,nci3,nci4,nci5
              srch = 'single'
              single = 1
              ifos2d = 1  ! Always generate 2-D data with single option  
!  generate local fos for Fellenius option paired with single option.              
              IF (method.eq.'o') foslocal = 1   
              IF (rad.le.0) THEN
                errmessage = 'invalid value for "rad"; "rad" must be > 0'
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')              
                ierr = 1
              END IF                     
!       If search file to be specified           
          ELSE IF (srch(1:4).eq.'file'.or.srch(1:4).eq.'FILE') THEN
              srchfile = 1
              srch = 'file'
              READ (12,*)
              READ (12,*,IOSTAT=ios2) ismin,jsmin,nsrchres       
!       If search box to be specified           
          ELSE IF (srch(1:3).eq.'box'.or.srch(1:3).eq.'BOX') THEN
              srch = 'box'
              READ (12,*)
              READ (12,*,IOSTAT=ios2) ismin,jsmin,ismax,jsmax,nsrchres
!       If invalid characters read
          ELSE
              errmessage = 'invalid "srch" parameter'
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                
              ierr = 1
          END IF
                     
          IF ((srch.eq.'box'.or.srch.eq.'file').and.nsrchres.le.0) THEN
            errmessage = 'invalid value for "nsrchres"; "nsrchres" must be > 0'
            Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')              
            ierr = 1
          END IF                            

          IF (ios2.ne.0) THEN
            errmessage = 'invalid or missing value following "srch" parameter-line ID'
            Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')              
            ierr = 1
          END IF        
        
          CASE ('ro','RO')   
            IF (heading(1:3).eq.'rot'.or.heading(1:3).eq.'ROT') THEN          
             IF (isupported(heading(1:2))) THEN                
              READ (12,*,IOSTAT=ios2) xcenrot,ycenrot,zcenrot
              irotcen = 1                                            
              IF (single.ne.1) THEN
                errmessage = '"rotcen" line will be ignored except with single surface'
                Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')                 
              END IF
              IF (ios2.ne.0) THEN
                errmessage = 'invalid or missing value following "rotcen" parameter-line ID'
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')              
                ierr = 1
              END IF
             ELSE  ! write error message and return to readin to stop execution
              errmessage = 'an unsupported feature was selected in main parameter input file'    
              Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')
              WRITE (*,*) 'unsupported feature = ', heading
              WRITE (39,*)'unsupported feature = ',heading   
              STOP                                
             END IF 
            END IF  

        
          CASE ('zs','ZS')
           IF (heading(1:5).eq.'zsmin'.or.heading(1:5).eq.'ZSMIN') THEN
            READ (12,*,IOSTAT=ios2) zsmin,zsmax,zsrchres
            IF (single.eq.1) THEN
              errmessage = '"zsmin" line will be ignored with single surface'
              Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')                          
            ELSE IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "zsmin" parameter-line ID'
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')             
              ierr = 1
            END IF  
            IF (zsmin.gt.zsmax) THEN
              errmessage = 'minimum elevation > maximum elevation for search lattice' 
              solution = 'assign "zsmin" < "zsmax"'            
              Call WriteError(0,errmessage,problemtype,solution,'no ',0,' ')
              solution = 'check Scoops3D manual for description of required parameters'              
            END IF 
!      Calculate max number of search lattice elevations
            IF (zsrchres.gt.0.) THEN
              ksmax = NINT((zsmax-zsmin)/zsrchres)+1
            ELSE
              errmessage = '"zsrchres" must be > 0'
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')              
              ierr = 1
            END IF
           END IF 
            
          CASE ('ir','IR')
            SELECT CASE (heading(1:4))
              CASE ('iref','IREF')
                READ (12,*,IOSTAT=ios2) irefine,multres,fostol
                IF (ios2.ne.0) THEN
                  errmessage = 'invalid or missing value following "irefine" parameter-line ID'
                  Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                                 
                  ierr = 1
                  CYCLE
                END IF 
                IF (single.eq.1.and.irefine.eq.1) THEN
                  errmessage = '"irefine" line will be ignored with single surface'
                  Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')                 
                  irefine = 0
                END IF
                IF (multres.le.1.and.irefine.eq.1) THEN
                  irefine = 0 
                  errmessage = '"multres" should be > 1; search will continue at "nsrchres", "zsrchres" resolution'
                  Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')                  
                ELSE
                  IF (fostol.le.0) THEN
                    errmessage = '"fostol" must be > 0' 
                    Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')     
                    ierr = 1                  
                  END IF
                  fostol = fostol/100.0_pr  ! convert percent to fraction
                END IF  
                        
              CASE ('irel','IREL')
                READ(12,*,IOSTAT=ios2) irelfos
                IF (ios2.ne.0) THEN
                  errmessage = 'invalid or missing value following "irelfos" parameter-line ID'
                  Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                   
                  ierr = 1
                  CYCLE
                END IF  
            END SELECT                             
                       
          CASE ('dr','DR')
            READ (12,*,IOSTAT=ios2) dr,deginc,degmax,ndir
            IF (ios2.ne.0) THEN 
              IF (dr.eq.rnull.OR.degmax.eq.rnull) THEN              
                errmessage = 'invalid or missing value following "dr" parameter-line ID'
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')              
                ierr = 1
                CYCLE                
              ELSE
                IF (ndir.eq.inull) THEN
                  ios2 = 0
                  BACKSPACE(12) 
                ELSE
                  errmessage = 'invalid or missing value following "dr" parameter-line ID'              
                  Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                 
                  ierr = 1
                  CYCLE                                
                END IF 
              END IF 
            END IF            
            IF (dr.le.0) THEN  
              errmessage = 'invalid value for "dr"; "dr" must be > 0'
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')     
              ierr = 1
            END IF  
            IF (deginc.lt.0) THEN  
              errmessage = 'invalid value for "deginc"; "deginc" must be > 0'
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')     
              ierr = 1
            END IF  
            IF (degmax.lt.0) THEN  
              errmessage = 'invalid value for "degmax"; "degmax" must be > 0' 
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')     
              ierr = 1
            END IF              
            IF (ndir.eq.0) THEN
              BACKSPACE(12)
              READ (12,*,IOSTAT=ios2) dr,deginc,degmax 
            END IF                            
            IF (single.eq.1) THEN
              errmessage = '"dr" line will be ignored with single surface'
              Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')               
            END IF            
!       Calculate number of slip directions.            
            IF ((deginc.le.degmax).and.(deginc.ne.0.)) THEN
              numdir = NINT((2.0_pr*degmax)/deginc) + 1
            ELSE
              numdir = 1
              degmax = 0.0_pr
            END IF  
            IF (degmax.gt.0.and.MOD(degmax,deginc).ne.0) THEN
              errmessage = 'error following "dr" parameter; "deginc" must be a multiple of "degmax"'              
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                     
              ierr=1             
            END IF        
!       Compare number of slip direction with number of slip directions specified in input file                     
            IF (ierr.eq.0.and.ndir.ne.inull) THEN
              IF (ndir.ne.inull.and.ndir.ne.numdir) THEN
                errmessage = '"numdir" is incorrect; Scoops3d will use calculated number of directions and continue'
                Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')                               
                PRINT *,'user specified # of directions =',ndir,'calculated # of directions = ',numdir 
                WRITE (39,*) 'user specified # of directions =',ndir,'calculated # of directions = ',numdir 
              END IF
            END IF              
                               
          CASE ('va','VA')
           IF (heading(1:6).eq.'vacrit'.or.heading(1:6).eq.'VACRIT') THEN
            READ (12,*,IOSTAT=ios2) vacriterion,armin,armax,vmin,vmax,tol,limcol            
            IF ((vacriterion.eq.'V'.or.vacriterion.eq.'v').and.vmin.lt.0) THEN  
              errmessage = 'invalid value for "vmin"; "vmin" must be > 0'              
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')     
              ierr = 1
            END IF  
            IF ((vacriterion.eq.'A'.or.vacriterion.eq.'a').and.armin.lt.0) THEN  
              errmessage = 'invalid value for "armin"; "armin" must be > 0'              
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')     
              ierr = 1
            END IF 
            IF (tol.le.0) THEN  
              errmessage = 'invalid value for "tol"; "tol" must be > 0'              
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')     
              ierr = 1
            END IF                                 
            IF (limcol.le.0) THEN  
              errmessage = 'invalid value for "limcol"; "limcol" must be > 0'              
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')     
              ierr = 1
            END IF                                                              
            IF (single.eq.1) THEN
              errmessage = '"vacriterion" line will be ignored with single surface'
              Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')              
            END IF
            IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "vacriterion" parameter-line ID'              
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')              
              ierr = 1
            END IF
           END IF 
            
!   Failure surface flag (not required)     
          CASE ('fl','FL')
           IF (heading(1:6).eq.'flsurf'.or.heading(1:6).eq.'FLSURF') THEN          
            IF (isupported(heading(1:2))) THEN                              
              READ (12,*,IOSTAT=ios2) ifailsurf  
              IF (ios2.ne.0) THEN
                errmessage = 'invalid or missing value following "flsurf" parameter-line ID'              
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                
                ierr = 1
                CYCLE
              END IF
              IF (ifailsurf.ne.0) ifailsurf = 1  
            ELSE  ! write error message and return to readin to stop execution
              errmessage = 'an unsupported feature was selected in main parameter input file'    
              Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')
              WRITE (*,*) 'unsupported feature = ', heading
              WRITE (39,*)'unsupported feature = ',heading   
              STOP                               
            END IF
           END IF   
        
!   inclusion area flag (not required)     
          CASE ('ii','II')
           IF (heading(1:8).eq.'iincflag'.or.heading(1:8).eq.'IINCFLAG') THEN          
            IF (isupported(heading(1:2))) THEN                                      
              READ (12,*,IOSTAT=ios2) iincflag                     
              IF (ios2.ne.0) THEN
                errmessage = 'invalid or missing value following "iincflag" parameter-line ID'              
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                
                ierr = 1
                CYCLE
              END IF                               
              IF (iincflag.ne.0) iincflag = 1
            ELSE  ! write error message and return to readin to stop execution
              errmessage = 'an unsupported feature was selected in main parameter input file'    
              Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')
              WRITE (*,*) 'unsupported feature = ', heading
              WRITE (39,*)'unsupported feature = ',heading   
              STOP                                    
            END IF    
           END IF 
            
!   Stability method parameters       
          CASE ('me','ME')
           IF (heading(1:6).eq.'method'.or.heading(1:6).eq.'METHOD') THEN         
            READ (12,*,IOSTAT=ios2) method  
            IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "method" parameter-line ID'              
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')             
              ierr = 1
              CYCLE
            END IF       
            IF (method.eq.'B'.or.method.eq.'b') THEN
              method = 'b'
!             diter = diter/100.0_pr  ! convert percent to fraction
            ELSE
              IF (method.eq.'F'.or.method.eq.'f'.or.method.eq.'O'.or.method.eq.'o') THEN
                method = 'o'
                IF (single.eq.1) foslocal = 1 ! Always generate local fos for
!     Fellenius option paired with single option.
              ELSE
                errmessage = 'invalid "method" parameter; must choose Bishop or ordinary method'              
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')               
                ierr = 1
              END IF           
            END IF
           END IF 

!   user-specified value for diter can supercede SCOOPS default       
          CASE ('di','DI')
           IF (heading(1:5).eq.'diter'.or.heading(1:5).eq.'DITER') THEN           
            READ (12,*,IOSTAT=ios2) ditertemp  
            IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "diter" parameter-line ID'              
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')             
              ierr = 1
              CYCLE
            END IF       
            IF (method.eq.'B'.or.method.eq.'b') THEN
! if the user has specified a valid value for diter, this supercedes the default value
              IF (ditertemp.ne.rnull.and.ditertemp.gt.0.and.ditertemp.lt.0.01) THEN
               diter = ditertemp
              ELSE          
                IF (ditertemp.le.0) THEN
                  errmessage = 'invalid "diter" parameter; "diter" must be greater than 0 and < 0.01'              
                  Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                 
                  ierr = 1
                END IF  
!              diter = diter/100.0_pr  ! convert percent to fraction 
              END IF
            ELSE            
              errmessage = '"diter" line will be ignored when "method" is not Bishop'
              Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')                                    
            END IF 
           END IF       

!   Solution filters (not required)
          CASE ('ab','AB')
           IF (heading(1:8).eq.'absminma'.or.heading(1:8).eq.'ABSMINMA') THEN           
            READ (12,*,IOSTAT=ios2) absminma
            IF (absminma.lt.0) THEN  
              errmessage = 'invalid value for "absminma"; "absminma" must be > 0'              
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')     
              ierr = 1
            END IF                            
            IF (method.ne.'b') THEN
              errmessage = '"absminma" parameter will be ignored when "method" is not Bishop'
              Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')                  
            END IF
            IF (single.ne.1.and.method.eq.'b') filter = 1
            IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "absminma" parameter-line ID'              
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')             
              ierr = 1
            END IF
           END IF        

!   2-D Option (not required)  
          CASE ('IF','if') 
           IF (heading(1:6).eq.'ifos2d'.or.heading(1:6).eq.'IFOS2D') THEN      
            IF (isupported(heading(1:2))) THEN                   
              READ(12,*,IOSTAT=ios2) ifos2d
              IF (single.eq.1) ifos2d = 1
              IF (ios2.ne.0) THEN
                errmessage = 'invalid or missing value following "ifos2d" parameter-line ID'              
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')               
                ierr = 1
              END IF
            ELSE  ! write error message and return to readin to stop execution
              errmessage = 'an unsupported feature was selected in main parameter input file'    
              Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')
              WRITE (*,*) 'unsupported feature = ', heading
              WRITE (39,*)'unsupported feature = ',heading     
              STOP
            END IF    
           END IF                 

!   Groundwater configuration
          CASE ('wa','WA')
           IF (heading(1:5).eq.'water'.or.heading(1:5).eq.'WATER') THEN           
            READ (12,*,IOSTAT=ios2) water
            errmessage = 'invalid or missing value following "water" parameter-line ID' 
            IF (ios2.ne.0) THEN             
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')             
              ierr = 1  
              CYCLE
            END IF                          

            SELECT CASE (water(1:2))

             CASE ('no','NO')
               iwater = 0

             CASE ('ru','RU')
               iwater = 1
 
             CASE ('pz','PZ')            
              iwater = 2
              pzsurf = 1
              BACKSPACE (12)
              READ (12,*,IOSTAT=ios2) water,gamw
              IF (ios2.ne.0) THEN              
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')               
                ierr = 1 
                CYCLE 
              END IF              

             CASE ('3d','3D')
              iwater = 3
              BACKSPACE (12)
              READ (12,*,IOSTAT=ios2) water,gamw
              IF (ios2.ne.0) THEN              
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')               
                ierr = 1  
              END IF              

             CASE ('vs','VS')
              iwater = 4
              BACKSPACE (12)
              READ (12,*,IOSTAT=ios2) water,gamw
              IF (ios2.ne.0) THEN             
                 Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                 
                 ierr = 1  
              END IF
 
             CASE ('vg','VG')
              iwater = 5
              BACKSPACE (12)
              READ (12,*,IOSTAT=ios2) water,gamw
              IF (ios2.ne.0) THEN              
                 Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                 
                 ierr = 1  
              END IF

             CASE ('fx','FX')
              iwater = 6
              BACKSPACE (12)
              READ (12,*,IOSTAT=ios2) water,gamw
              IF (ios2.ne.0) THEN             
                   Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                 
                   ierr = 1  
              END IF    

             CASE DEFAULT
                 errmessage = 'invalid "water" parameter'              
                 Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')               
                 ierr = 1 
                   
            END SELECT 
            IF (iwater.ge.2.and.iwater.le.6) THEN
              IF (gamw.lt.0.0_pr) THEN  
                   errmessage = 'invalid "gamw" parameter; "gamw" must be >= 0'
                   Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                 
                   ierr = 1                    
              END IF
            END IF 
          END IF
        
          CASE ('st','ST')
           IF (heading(1:5).eq.'str3d'.or.heading(1:5).eq.'STR3D') THEN            
            READ (12,*,IOSTAT=ios2) str3d,linterp
            IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "str3D" parameter-line ID'              
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')             
              ierr=1 
            END IF  
           END IF
           
!  Units for length, unit weight, and cohesion
         CASE ('le','LE')   
          IF (heading(1:11).eq.'lengthunits'.or.heading(1:11).eq.'LENGTHUNITS') THEN         
           READ (12,1000,IOSTAT=ios) units 
           IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "lengthunits" parameter-line ID'              
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ') 
              ierr=1 
           ELSE                       
              units = ADJUSTL (units)
              n = SCAN(units,' ',BACK=.false.)          
              lengthunits = units(1:n-1)
              units = units(n:LEN_TRIM(units))
              units = ADJUSTL (units)          
              n = SCAN(units,' ',BACK=.false.)           
              ceeunits = units(1:n-1)
              units = units(n:LEN_TRIM(units))
              units = ADJUSTL (units)                    
              n = SCAN(units,' ',BACK=.false.)
              uwtunits = units(1:n-1)
           END IF
          END IF
           	                 
!   Material properties (either 3-D or layers, not both)
          CASE ('nm','NM')
           IF (heading(1:4).eq.'nmat'.or.heading(1:4).eq.'NMAT') THEN            
            READ (12,*,IOSTAT=ios2) nmat
            IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "nmat" parameter-line ID'              
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ') 
              ierr = 1
            END IF
            IF (nmat.lt.1) THEN  
                errmessage = 'invalid value for "nmat"; "nmat" must be > 0'              
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')     
                ierr = 1
                CYCLE
            END IF             
!            IF (nmat.lt.1) CYCLE
            READ (12,1000,IOSTAT=ios) heading  
            IF (heading(1:4).ne.'lnum'.and.heading(1:4).ne.'LNUM') THEN
              errmessage = 'invalid or missing "lnum" parameter-line ID'              
              CLOSE(12)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ') 
            END IF         
            IF (ALLOCATED (cee)) DEALLOCATE (cee)
            IF (ALLOCATED (phi)) DEALLOCATE (phi)
            IF (ALLOCATED (tanphi)) DEALLOCATE (tanphi)
            IF (ALLOCATED (ru)) DEALLOCATE (ru)
            IF (ALLOCATED (thetares)) DEALLOCATE (thetares)
            IF (ALLOCATED (thetasat)) DEALLOCATE (thetasat)
            IF (ALLOCATED (fxa)) DEALLOCATE (fxa)
            IF (ALLOCATED (fxn)) DEALLOCATE (fxn) 
            IF (ALLOCATED (fxm)) DEALLOCATE (fxm) 
            IF (ALLOCATED (fxr)) DEALLOCATE (fxr) 
            IF (ALLOCATED (vga)) DEALLOCATE (vga)   
            IF (ALLOCATED (vgn)) DEALLOCATE (vgn)       
            IF (ALLOCATED (gamr)) DEALLOCATE (gamr)
            IF (ALLOCATED (lnum)) DEALLOCATE (lnum)
            IF (ALLOCATED (gsame)) DEALLOCATE (gsame)
            ALLOCATE (cee(nmat),phi(nmat),tanphi(nmat),ru(nmat),&
                     gamr(nmat,3),lnum(nmat),gsame(nmat),STAT=error)
            IF (iwater.gt.3) THEN
               ALLOCATE (thetares(nmat),thetasat(nmat),STAT=error)     
               IF (iwater.eq.5) ALLOCATE (vga(nmat),vgn(nmat),STAT=error) 
               IF (iwater.eq.6) ALLOCATE (fxa(nmat),fxn(nmat),fxm(n),fxr(n),STAT=error)
            END IF   
            IF (error.ne.0) THEN
              errmessage = 'strength arrays not allocated successfully' 
              solution = 'reduce memory requirements. See "Practical Considerations" chapter of Scoops3D manual'
              CLOSE(12)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ') 
            END IF    
            cee = -1.0_pr
            phi = -1.0_pr        
            gamr = -1.0_pr
            lnum = inull
            ru = 0.0_pr
            IF (iwater.gt.3) THEN
              thetares = -1.0_pr
              thetasat = -1.0_pr                 
            END IF  
            IF (iwater.eq.5) THEN
              vga = rnull
              vgn = rnull
            END IF 
            IF (iwater.eq.6) THEN
              fxa = rnull
              fxn = rnull
              fxm = rnull
              fxr = rnull
            END IF  
            gsame = 0
            gsameeach = 1
            gsameall = 1
            IF (str3d.eq.1.and.nmat.gt.1) THEN
              errmessage = 'layer files cannot be used in combination with "str3d" = 1' 
              solution = 'select "nmat" = 1 when using 3D material properties file'
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ') 
              solution = 'check Scoops3D manual for description of required parameters'                   
              ierr=1
              nmat = 1
            END IF

            DO n=1,nmat                              
              SELECT CASE (iwater)
                CASE(0)       ! no water pressure
                  READ (12,*,IOSTAT=ios2) lnum(n),cee(n),phi(n),gamr(n,1)                         
                CASE(1)      ! ru approximation
                  READ (12,*,IOSTAT=ios2) lnum(n),cee(n),phi(n),gamr(n,1),ru(n)             
                  IF (ru(n).lt.0) THEN
                    errmessage = 'invalid or missing value following "lnum" parameter-line ID, material # ' 
                    solution = 'assign valid "ru" value for each layer'
                    Call WriteError(2,errmessage,problemtype,solution,'int',lnum(n),' ')                   
                    ierr = 1
                    CYCLE
                  END IF                                                            
                CASE(2,3)  ! peizometric surface or 3d pressure file
                  READ (12,*,IOSTAT=ios2) lnum(n),cee(n),phi(n),&
                     gamr(n,2),gamr(n,3)
                CASE(4)    ! 3d pressure file with relative water contents and curve fit parameters                   
                        READ (12,*,IOSTAT=ios2) lnum(n),cee(n),phi(n),&
                          gamr(n,3),thetares(n),thetasat(n)
                CASE(5)  ! 3d pressure file with vanGenuchten SWCC
                        READ (12,*,IOSTAT=ios2) lnum(n),cee(n),&
                          phi(n),gamr(n,3),thetares(n),thetasat(n),vga(n),vgn(n) 
                          IF (vga(n).eq.rnull.or.vgn(n).eq.rnull.or.vga(n).le.0.or.vgn(n).le.0) THEN
                           errmessage = 'invalid or missing value following "lnum" parameter-line ID, material # ' 
                           solution = '"nmat" line must contain positive numeric values for vga, vgn '
                           Call WriteError(2,errmessage,problemtype,solution,'int',lnum(n),' ')                           
                           ierr = 1
                           CYCLE
                         END IF                                                
                CASE(6)    ! 3d pressure file with Fredlund and Xing SWCC
                        READ (12,*,IOSTAT=ios2) lnum(n),cee(n),&
                          phi(n),gamr(n,3),thetares(n),thetasat(n),fxa(n),fxn(n),fxm(n),fxr(n)                                  
                          IF (fxa(n).eq.rnull.or.fxn(n).eq.rnull.or.fxm(n).eq.rnull.or.fxr(n).eq.rnull.or.&
                                fxa(n).le.0.or.fxn(n).le.0.or.fxm(n).le.0.or.fxr(n).le.0) THEN
                           errmessage = 'error or missing value following "lnum" parameter-line ID, material # ' 
                           solution = '"nmat" line must contain positive numeric values for fxa, fxn, fxm, fxr'
                           Call WriteError(2,errmessage,problemtype,solution,'int',lnum(n),' ')                        
                           ierr = 1
                           CYCLE
                          END IF                                                                                                                
              END SELECT  
 
              solution = 'check Scoops3D manual for description of required parameters'
              IF (lnum(n).ne.n) THEN
                    errmessage = 'layer identification numbers "lnum" should be listed in ascending order, lnum = ' 
                    Call WriteError(2,errmessage,problemtype,solution,'int',lnum(n),' ')                   
                    ierr = 1
                    CYCLE
              END IF                  
              IF (cee(n).lt.0) THEN
               IF (str3d.ne.1.or.cee(n).ne.-1) THEN  
                errmessage = 'invalid value for "cee"; "cee" must be >= 0'              
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')     
                ierr = 1
               END IF 
              END IF 
              IF (phi(n).lt.0) THEN
               IF (str3d.ne.1.or.phi(n).ne.-1) THEN    
                errmessage = 'invalid value for "phi"; "phi" must be >= 0'              
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')     
                ierr = 1
               END IF 
              END IF
 
              IF (iwater.eq.0.or.iwater.eq.1) THEN                        
                    IF (gamr(n,1).lt.0) THEN 
                     IF (str3d.ne.1.or.gamr(1,1).ne.-1) THEN 
                      errmessage = 'invalid value for "gamt"; "gamt" must be >= 0'              
                      Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')     
                      ierr = 1
                     END IF 
                    END IF 
              ELSE IF (iwater.eq.2.or.iwater.eq.3) THEN   ! peizometric surface or 3d pressure file             
                    IF (gamr(n,2).lt.0.or.gamr(n,3).lt.0) THEN 
                     IF (str3d.ne.1.or.gamr(1,1).ne.-1) THEN 
                      errmessage = 'invalid value for "gamps" or "gams"; "gamps" and "gams" must be >= 0'              
                      Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')     
                      ierr = 1
                     END IF 
                    END IF             
              ELSE IF (iwater.eq.4.or.iwater.eq.5.or.iwater.eq.6) THEN
                 IF (gamr(n,3).lt.0) THEN  
                    errmessage = 'invalid value for "gams"; "gams" must be >= 0'              
                    Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')     
                    ierr = 1
                 END IF                       
                 IF (thetares(n).lt.0.or.thetasat(n).lt.0) THEN
                    errmessage = 'invalid value for "thetares" or "thetasat"; values must be >= 0, material # ' 
                    Call WriteError(2,errmessage,problemtype,solution,'int',lnum(n),' ')                  
                    ierr = 1
                    CYCLE
                 END IF                      
              END IF
              
              IF (ios2.ne.0) THEN
                errmessage = 'invalid or missing value following "lnum" parameter-line ID, material #' 
                Call WriteError(2,errmessage,problemtype,solution,'int',n,' ')               
                ierr = 1
                CYCLE
              END IF
!       Set flags for constant parameters when using 3D strength input.
!       (cee(1) and phi(1) equal -1 if in 3D strength file)                              
              IF (cee(1).eq.-1.0_pr) cflag = 0
              IF (phi(1).eq.-1.0_pr) THEN
                fflag = 0
              ELSE
                tanphi(n) = phi(n)*pi/180.0_pr
                tanphi(n) = TAN(tanphi(n))
              END IF
              
!     Determine whether unit weights are identical for unsat, partial sat and sat 
!     conditions for each layer and for all layers.
              IF (gamr(1,2).lt.0.0_pr.and.str3d.eq.1) THEN
                dflag = 0
                gsameeach = 0
                gsameall = 0
              ELSE            
                SELECT CASE (iwater)
                  CASE(0,1,4,5,6)       ! no water pressure, ru, partially saturated
                    gsame(n) = 1                     
                  CASE(2,3)       ! peizometric surface or 3d pressure file
                    IF  (gamr(n,2).eq.gamr(n,3)) gsame(n) = 1                                      
                END SELECT
              END IF
            END DO  !  Loop on material layers
            
            IF (dflag.eq.1) THEN 
              DO n=1,nmat
                IF (gsame(n).eq.0) THEN
                  gsameall = 0
                  gsameeach = 0
                  EXIT
                ELSE
                  IF (n.gt.1) THEN                
                    SELECT CASE (iwater)
                      CASE(0,1)       ! no water pressure, ru, partially saturated
                        IF (gamr(n,1).ne.gamr(n-1,1)) gsameall=0                   
                      CASE(2,3)       ! peizometric surface or 3d pressure file
                        IF (gamr(n,2).ne.gamr(n-1,2)) gsameall=0                                      
                      CASE(4,5,6)       ! peizometric surface or 3d pressure file
                        IF (gamr(n,3).ne.gamr(n-1,3)) gsameall=0                                   
                    END SELECT                
                  END IF
                END IF 
              END DO 
            END IF
           END IF 
                                     
           
!   Earthquake loading            
          CASE ('eq','EQ')
            READ(12,*,IOSTAT=ios2) eq
            IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "eq" parameter-line ID' 
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                
              ierr = 1
            END IF
            IF (eq.lt.0) THEN 
                errmessage = 'invalid value for "eq"; "eq" must be >= 0'              
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')     
                ierr = 1
            END IF             

!   Output control            
          CASE ('is','IS')           
            SELECT CASE (heading(1:4))
              CASE ('isqo','ISQO')
                READ(12,*,IOSTAT=ios2) isqout
                IF (ios2.ne.0) THEN
                  errmessage = 'invalid or missing value following "isqout" parameter-line ID' 
                  Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                  
                  ierr = 1
                  CYCLE
                END IF
              CASE ('isub','ISUB')
                READ(12,*,IOSTAT=ios2) isubsurf,zfrac 
                IF (ios2.ne.0) THEN
                  errmessage = 'invalid or missing value following "isubsurf" parameter-line ID' 
                  Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                   
                  ierr = 1
                  CYCLE
                END IF
                IF (zfrac.le.0) THEN  
                    errmessage = '"zfrac" must be > 0; no subsurface file will be created'              
                    Call WriteError(0,errmessage,problemtype,solution,'no ',0,' ')     
                END IF                               
                
              CASE ('isur','ISUR')
                READ(12,*,IOSTAT=ios2) isurfwat,gamsurf 
                IF (ios2.ne.0) THEN
                  errmessage = 'invalid or missing value following "isurfwat" parameter-line ID' 
                  Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                   
                  ierr = 1
                  CYCLE
                 END IF                 
             END SELECT   
    
          CASE ('ic','IC')
           IF (heading(1:5).eq.'icrit'.or.heading(1:5).eq.'ICRIT') THEN           
            READ(12,*,IOSTAT=ios2) icritlattice
            IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "icritlattice" parameter-line ID' 
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')               
              ierr = 1
              CYCLE
            END IF                   
           END IF               

          CASE ('il','IL')
           IF (heading(1:8).eq.'ilattice'.or.heading(1:8).eq.'ILATTICE') THEN           
            READ(12,*,IOSTAT=ios2) ilattice
            IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "ilattice" parameter-line ID' 
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')               
              ierr = 1
              CYCLE
            END IF   
           END IF                       

          CASE ('re','RE')
           IF (heading(1:6).eq.'remove'.or.heading(1:6).eq.'REMOVE') THEN           
            READ (12,*,IOSTAT=ios2) remove,foscut
            IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "remove" parameter-line ID' 
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')               
              ierr = 1
              CYCLE
            END IF
            IF (foscut.lt.0) THEN  
              errmessage = 'invalid value for "foscut"; "foscut" should be >= 0'              
              Call WriteError(0,errmessage,problemtype,solution,'no ',0,' ')     
              ierr = 1
            END IF              
           END IF              

          CASE ('nr','NR')
           IF (heading(1:6).eq.'nretro'.or.heading(1:6).eq.'NRETRO') THEN          
            IF (isupported(heading(1:2))) THEN                               
              READ(12,*,IOSTAT=ios2) nretro
              IF (ios2.ne.0) THEN
                errmessage = 'invalid or missing value following "nretro" parameter-line ID' 
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                 
                ierr = 1
                CYCLE
              END IF      
            ELSE  ! write error message and return to readin to stop execution
              errmessage = 'an unsupported feature was selected in main parameter input file'    
              Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')
              WRITE (*,*) 'unsupported feature = ', heading
              WRITE (39,*)'unsupported feature = ',heading          
              STOP
            END IF
           END IF 
 
!   Input File names

          CASE ('de','DE')
           IF (heading(1:3).eq.'dem'.or.heading(1:3).eq.'DEM') THEN          
            READ(12,1000,IOSTAT=ios2) demfile
            IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "DEM file" parameter-line ID' 
              solution = 'specify a valid DEM file name'               
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')             
              ierr = 1
              CYCLE
            END IF
           END IF
            
          CASE ('se','SE')
           IF (heading(1:6).eq.'search'.or.heading(1:6).eq.'SEARCH') THEN          
            READ(12,1000,IOSTAT=ios2) sfile
            IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "search file" parameter-line ID' 
              solution = 'specify a valid search file name' 
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')               
              ierr = 1
              CYCLE
            END IF    
           END IF
           
          CASE ('pi','PI')
           IF (heading(1:5).eq.'piezo'.or.heading(1:5).eq.'PIEZO') THEN
            READ(12,1000,IOSTAT=ios2) pzfile
            IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "piezometric file" parameter-line ID' 
              solution = 'specify a valid piezometric-surface file name' 
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                
              ierr = 1
              CYCLE
            END IF
           END IF
           
          CASE ('pr','PR')
           IF (heading(1:8).eq.'pressure'.or.heading(1:8).eq.'PRESSURE') THEN
            READ(12,1000,IOSTAT=ios2) prfile
            IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "pressure head file" parameter-line ID' 
              solution = 'specify a valid 3D pressure-head file name'
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')               
              ierr = 1
              CYCLE
            END IF   
           END IF
            
          CASE ('la','LA')
           IF (heading(1:5).eq.'layer'.or.heading(1:5).eq.'LAYER') THEN
            READ(12,1000,IOSTAT=ios2) stratfile
            IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "layer file" parameter-line ID' 
              solution = 'specify a valid layer file name'
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')              
              ierr = 1
              CYCLE
            END IF         
           END IF
            
          CASE ('ma','MA')
           IF (heading(1:8).eq.'material'.or.heading(1:8).eq.'MATERIAL') THEN          
            READ(12,1000,IOSTAT=ios2) strgthfile
            IF (ios2.ne.0) THEN
              errmessage = 'invalid or missing value following "material properties file" parameter-line ID' 
              solution = 'specify a valid 3D material properties file name'
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                 
              ierr = 1
              CYCLE
            END IF   
           END IF
            
          CASE ('fa','FA')
           IF (heading(1:7).eq.'failure'.or.heading(1:7).eq.'FAILURE') THEN          
            IF (isupported(heading(1:2))) THEN                                  
              READ(12,1000,IOSTAT=ios2) failsurffile            
              IF (ios2.ne.0) THEN
                errmessage = 'invalid or missing value following "failsurffile" parameter-line ID' 
                solution = 'specify a valid failure surface file name'
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                 
                ierr = 1
                CYCLE
              END IF 
            ELSE  ! write error message and return to readin to stop execution
              errmessage = 'an unsupported feature was selected in main parameter input file'    
              Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')
              WRITE (*,*) 'unsupported feature = ', heading
              WRITE (39,*)'unsupported feature = ',heading     
              STOP
            END IF    
           END IF
 
         CASE ('in','IN')
           IF (heading(1:7).eq.'include'.or.heading(1:7).eq.'INCLUDE') THEN         
            IF (isupported(heading(1:2))) THEN                                           
              READ(12,1000,IOSTAT=ios2) includefile
              IF (ios2.ne.0) THEN
                errmessage = 'invalid or missing value following "includefile" parameter-line ID' 
                solution = 'specify a valid inclusion area file name'
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                 
                ierr = 1
                CYCLE
              END IF 
            ELSE  ! write error message and return to readin to stop execution
              errmessage = 'an unsupported feature was selected in main parameter input file'    
              Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')
              WRITE (*,*) 'unsupported feature = ', heading
              WRITE (39,*)'unsupported feature = ',heading                                   
              STOP     
            END IF    
           END IF
           
         CASE ('su','SU')
           IF (heading(1:7).eq.'surface'.or.heading(1:7).eq.'SURFACE') THEN         
            IF (isupported(heading(1:2))) THEN                                           
              READ(12,1000,IOSTAT=ios2) surfwfile
              IF (ios2.ne.0) THEN
                errmessage = 'invalid or missing value following "surfwfile" parameter-line ID' 
                solution = 'specify a valid surface-water file name'
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                    
                ierr = 1
                CYCLE
              END IF 
            ELSE  ! write error message and return to readin to stop execution
               errmessage = 'an unsupported feature was selected in main parameter input file'    
               Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')
               WRITE (*,*) 'unsupported feature = ', heading
               WRITE (39,*)'unsupported feature = ',heading             
               STOP     
            END IF
           END IF
           
! output directory was read in earlier, just skip over this line
          CASE ('ou','OU')
            READ(12,1000,IOSTAT=ios2) junk
            
          END SELECT

        END DO
       
! End read main input file.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!   Read additional input files and allocate arrays as required.

        IF (str3d.eq.1.and.nmat.lt.1) THEN
          nmat = 1
          ALLOCATE (cee(nmat),phi(nmat),tanphi(nmat),ru(nmat),&
                     gamr(nmat,3),lnum(nmat),gsame(nmat),STAT=error)
          IF (error.ne.0) THEN
            errmessage = 'material property arrays not allocated successfully' 
            solution = 'reduce memory requirements; see "Practical Considerations" chapter of Scoops3D manual'
            CLOSE(12)
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')           
          END IF    
          cee = -1.0_pr
          phi = -1.0_pr        
          gamr = -1.0_pr
          ru = 0.0_pr
          gsame = 0
          gsameeach = 1
          gsameall = 1                
          IF (iwater.ne.1) THEN
            cflag = 0
            dflag = 0
            fflag = 0
          ELSE
            errmessage = 'missing "ru" values ' 
            solution = '"nmat" parameter line must contain values for "ru"'
            Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                 
            ierr = 1
          END IF 
        END IF        
!     Read DEM data.
        PRINT *, 'Opening DEM file: ',demfile
        IF (LEN_TRIM(demfile).gt.220) THEN
          errmessage = 'Scoops3D allows a maximum of 220 characters for file names'
          Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')            
        END IF                    
        CALL Readdem (demfile)

!     Set initial search resolution        
        IF (irefine.eq.1.and.multres.gt.1) THEN
          nsrch = multres * nsrchres
        ELSE
          nsrch = nsrchres
        END IF
!     If reading search file.
        IF (srchfile.eq.1) THEN       
          PRINT *, 'Opening search file: ',sfile
          IF (LEN_TRIM(sfile).gt.220) THEN
            errmessage = 'Scoops3D allows a maximum of 220 characters for file names'
            Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')             
          END IF                    
          CALL Readsearch (sfile,nxsrch,nysrch,nsrchpt) 
!     Find approximation of number of search points if reading search grid
!     but have nsrch > 1.
          IF (nsrch.gt.1) THEN 
            nsrchptini = nsrchpt
            nsrchpt = nsrchpt/(nsrch*nsrch)
          END IF       
          IF (nsrchpt.lt.1) nsrchpt = 1
!     Calculate number of search grid points.       
        ELSE
          nxsrch=(ismax-ismin) + 1
          nysrch=(jsmax-jsmin) + 1
          IF (nsrch.gt.1) THEN
            IF (nxsrch.gt.nsrch) THEN
              nxsrch = INT((nxsrch-1)/nsrch) + 1
            ELSE
              nxsrch = 1
            END IF
            IF (nysrch.gt.nsrch) THEN
              nysrch = INT((nysrch-1)/nsrch) + 1
            ELSE
              nysrch = 1
            END IF            
          END IF
          nsrchpt = nxsrch*nysrch
        END IF
        
!     Allocate, initialize, and read failure surface data.
        IF (ifailsurf.eq.1) THEN
          PRINT *, 'Opening failure surface file: ',failsurffile
          IF (LEN_TRIM(failsurffile).gt.220) THEN
             errmessage = 'Scoops3D allows a maximum of 220 characters for file names'
             Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')             
          END IF                   
          ALLOCATE (failsurf(nx,ny),failsurfdepth(nx,ny),STAT=error) 
          IF (error.ne.0) THEN          
            errmessage = 'failure surface arrays not allocated successfully' 
            solution = 'reduce memory requirements; see "Practical Considerations" chapter of Scoops3D manual'
            CLOSE(12)
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
          END IF            
          CALL Readfailsurf(failsurffile)
        END IF

!     If reading inclusion area file.
        IF (iincflag.eq.1) THEN        
          PRINT *, 'Opening include area file: ',includefile
          IF (LEN_TRIM(includefile).gt.220) THEN
            errmessage = 'Scoops3D allows a maximum of 220 characters for file names'
            Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')             
          END IF                   
          ALLOCATE (includegrid(nx,ny),STAT=error) 
          CALL Readinclude (includefile) 
        END IF          

!     Allocate, initialize, and read piezometric data.
        IF (iwater.eq.2) THEN
          PRINT *, 'Opening piezometric file: ',pzfile
          IF (LEN_TRIM(pzfile).gt.220) THEN
            errmessage = 'Scoops3D allows a maximum of 220 characters for file names'
            Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')            
          END IF                         
          ALLOCATE (piezo(nx,ny), STAT=error)
          IF (error.ne.0) THEN          
            errmessage = 'piezometric data arrays not allocated successfully' 
            solution = 'reduce memory requirements; see "Practical Considerations" chapter of Scoops3D manual'
            CLOSE(12)
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')    
          END IF            
          piezo = rnull
          CALL Readpiezo(pzfile)
        END IF

!     Allocate, initialize, and read pressure head data file.
        IF (iwater.gt.2) THEN        
          PRINT *, 'Opening pressure head file: ',prfile
           IF (LEN_TRIM(prfile).gt.220) THEN
            errmessage = 'Scoops3D allows a maximum of 220 characters for file names'
            Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')            
          END IF                   
          CALL Readpressh(prfile,prcoords,pcount,pnum)
        END IF

!     If no piezometric surface
        IF (pzsurf.eq.0) THEN
!     If no water data or ru approximation use total unit weights      
          IF (iwater.eq.0.or.iwater.eq.1) moist = 1
!     If using 3D pressure data and no piezometric surface detected
!     use partially saturated conditions                   
          IF (iwater.eq.3) moist = 2
          IF (iwater.eq.4.or.iwater.eq.5.or.iwater.eq.6) moist = 3
        ELSE
          moist = 2
        END IF 

!     Allocate, initialize, and read stratigraphic layer data, if more than one material.
        IF (nmat.gt.1) THEN   
          PRINT *, '0pening stratigraphic layers file: ',stratfile  
          IF (LEN_TRIM(stratfile).gt.220) THEN
            errmessage = 'Scoops3D allows a maximum of 220 characters for file names'
            Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')           
          END IF          
          ALLOCATE (layer(nmat-1,nx,ny),uwtlay(nx,ny,nmat-1),STAT=error)
          ALLOCATE (duwtlaydz(nx,ny,nmat-1),STAT=error)
          ALLOCATE (minlayer(nmat-1),maxlayer(nmat-1),activelay(nmat-1), STAT=error)
          IF (error.ne.0) THEN
            errmessage = 'layer arrays not allocated successfully' 
            solution = 'reduce memory requirements; see "Practical Considerations" chapter of Scoops3D manual'
            CLOSE(17)
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')            
          END IF
          layer = rnull 
          uwtlay = rnull
          minlayer = -rnull
          maxlayer = rnull
          activelay = 0
          duwtlaydz = 0.0_pr
          CALL Readstrat(stratfile)
        END IF

!     Read 3D strength data.
        IF (str3d.eq.1) THEN   
          PRINT *, 'Opening 3D strength parameters file: ',strgthfile 
          IF (LEN_TRIM(strgthfile).gt.220) THEN
            errmessage = 'Scoops3D allows a maximum of 220 characters for file names'
            Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')            
          END IF                         
          CALL Readstrength(strgthfile,strcoords,strnum)
        END IF

!     Process water content data to precalculate depth-weighted unit weights
        IF (iwater.eq.4.or.iwater.eq.5.or.iwater.eq.6) Call iwater4wt 

!     Read surface water layer data.
        IF (isurfwat.eq.1) THEN 
          IF (iwater.eq.0.or.iwater.eq.1) THEN
            errmessage = 'can not combine submergence with ru or no water pressures' 
            solution = 'check Scoops3D manual for description of valid groundwater configurations'
            Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')              
            ierr = 1
          END IF   
          PRINT *, 'Opening surface water elevation file: ',surfwfile 
          IF (LEN_TRIM(surfwfile).gt.220) THEN
            errmessage = 'Scoops3D allows a maximum of 220 characters for file names'
            Call WriteError(0,errmessage,problemtype,'no','no ',0,' ')             
          END IF
          ALLOCATE (surfwat(nx,ny), STAT=error)
          IF (error.ne.0) THEN
            errmessage = 'surface water array not allocated successfully' 
            solution = 'reduce memory requirements; see "Practical Considerations" chapter of Scoops3D manual'
            CLOSE(12)
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')            
          END IF            
          surfwat = rnull                         
          CALL Readsurfwater(surfwfile)
        END IF        
        
! End read additional input files.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  

!     Check that essential input was read in and other error checking:

        solution = 'check Scoops3D manual for description of required parameters'
        IF (srch.eq.'file'.or.srch.eq.'box') THEN
          IF (ismin.eq.inull) THEN
            errmessage = 'invalid or missing "ismin" parameter -- required with file or box search' 
            Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')             
            ierr = 1
          END IF
        ELSE IF (srch.eq.'single') THEN
          IF (xcen.eq.rnull) THEN
              errmessage = 'invalid or missing "xcen" parameter -- required with single surface' 
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                
              ierr = 1
          END IF
        END IF

        IF (srch.eq.'file'.or.srch.eq.'box') THEN
          IF (zsmin.eq.rnull) THEN
            errmessage = 'invalid or missing "zsmin" parameter -- required with file or box search' 
            Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')              
            ierr = 1
          END IF
          IF (dr.eq.rnull) THEN
            errmessage = 'invalid or missing "dr" parameter -- required with file or box search' 
            Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')              
            ierr = 1
          END IF
          IF (vacriterion.eq.'X') THEN
            errmessage = 'invalid or missing "vacriterion" parameter  -- required with file or box search' 
            Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')              
            ierr = 1
          END IF
        END IF
        IF (method.eq.'X') THEN
          errmessage = 'invalid or missing "method" parameter' 
          Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')            
          ierr = 1
        END IF
        IF (water.eq.'X') THEN
          errmessage = 'invalid or missing "water" parameter' 
          Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')            
          ierr = 1
        END IF
        IF (str3d.eq.0.and.nmat.lt.1) THEN
          errmessage = 'invalid or missing material strength data' 
          solution = 'check Scoops3D manual for description of required parameters -- specify layers or 3D strength file'
          Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')        
          ierr = 1
        END IF

!     Echo input to output file
        WRITE (20,1100)
        WRITE (20,1200)
        WRITE (20,1300)
        WRITE (20,1400)
        WRITE (20,1500) version
        WRITE (20,1600)
        WRITE (20,1900) filout
        WRITE (20,1700) date(5:6),date(7:8),date(1:4),time(1:2),time(3:4),time(5:6)
        IF (title.ne.'X') WRITE (20,1800) title
        
        WRITE (20,2000) 
        WRITE (20,2100) 
        WRITE (20,2200) demfile
        IF (srchfile.eq.1) WRITE (20,2300) sfile
        IF (iwater.eq.2) WRITE (20,2400) pzfile
        IF (iwater.eq.3) WRITE (20,2500) prfile
        IF (iwater.eq.4.or.iwater.eq.5.or.iwater.eq.6) WRITE (20,2510) prfile        
        IF (nmat.gt.1) WRITE (20,2600) stratfile
        IF (str3d.eq.1) WRITE (20,2700) strgthfile
        IF (ifailsurf.eq.1) WRITE (20,2750) failsurffile
        IF (iincflag.eq.1) WRITE (20,2760) includefile
        IF (isurfwat.eq.1) WRITE (20,2765) surfwfile
        WRITE (20,2800) filin
        
        WRITE (20,2000)
        WRITE (20,2900)
        WRITE (20,3000)
        WRITE (20,3100)
        WRITE (20,3200) demfile
        WRITE (20,3300) nx,ny
        WRITE (20,3350) nx*ny
        WRITE (20,3360) COUNT(zdem.ne.rnull)
        IF (lengthunits.eq.'  ') THEN
          WRITE (20,3405) delxy
          WRITE (20,3505) zmin
          WRITE (20,3605) zmax 
          WRITE (20,3705) xll,yll 
        ELSE
          WRITE (20,3400) lengthunits,delxy
          WRITE (20,3500) lengthunits,zmin
          WRITE (20,3600) lengthunits,zmax 
          WRITE (20,3700) lengthunits,xll,yll 
        END IF    

        IF (lengthunits.ne.'  ') THEN
          WRITE (20,3000)
          WRITE (20,8580)
          WRITE (20,8583) 
          WRITE (20,8586) lengthunits,ceeunits,uwtunits
        END IF  
 
        WRITE (20,3000)
        WRITE (20,8600)
        IF (str3d.eq.1) THEN
          WRITE (20,8750)
          IF (nmat.eq.1) THEN
            WRITE (20,8800) nmat            
            IF (ceeunits.ne.'        ') THEN    
              SELECT CASE (iwater)
                CASE(0)       ! no water pressure
                  WRITE (20,8811)                
                  WRITE (20,8901) 
                  WRITE (20,8911) ceeunits,uwtunits
                  WRITE (20,9001) lnum(1),cee(1),phi(1),gamr(1,1)                         
                CASE(1)      ! ru approximation
                   WRITE (20,8811)                              
                   WRITE (20,8902)  
                   WRITE (20,8911) ceeunits,uwtunits
                   WRITE (20,9002) lnum(1),cee(1),phi(1),gamr(1,1),ru(1)                                 
                CASE(2,3)  ! peizometric surface or 3d pressure file
                   WRITE (20,8812)                
                   WRITE (20,8813)                              
                   WRITE (20,8903)
                   WRITE (20,8912) ceeunits,uwtunits,uwtunits
                   WRITE (20,9002) lnum(1),cee(1),phi(1), gamr(1,2),gamr(1,3)
                CASE(4) ! 3d pressure file with relative water contents
                   WRITE (20,8814)                
                   WRITE (20,8815)                   
                   WRITE (20,8904)
                   WRITE (20,8911) ceeunits,uwtunits
                   WRITE (20,9003) lnum(1),cee(1),phi(1),gamr(1,3),thetares(1),thetasat(1)                       
                CASE(5)  ! 3d pressure file with vanGenuchten SWCC
                   WRITE (20,8814)                
                   WRITE (20,8815)                 
                   WRITE (20,8907)
                   WRITE (20,8911) ceeunits,uwtunits
                   WRITE (20,9005) lnum(1),cee(1),phi(1),gamr(1,3),thetares(1),thetasat(1),&
                          vga(1),vgn(1)              
                CASE(6)  ! 3d pressure file with Fredlund and Xing SWCC   
                   WRITE (20,8814)                
                   WRITE (20,8815)             
                   WRITE (20,8905)
                   WRITE (20,8911) ceeunits,uwtunits
                   WRITE (20,9004) lnum(1),cee(1),phi(1),gamr(1,3),thetares(1),thetasat(1),&
                          fxa(1),fxn(1),fxm(1),fxr(1)                                                               
              END SELECT 
            ELSE           
              SELECT CASE (iwater)
                CASE(0)       ! no water pressure
                  WRITE (20,8811)                      
                  WRITE (20,8901) 
                  WRITE (20,9001) lnum(1),cee(1),phi(1),gamr(1,1)                         
                CASE(1)      ! ru approximation
                  WRITE (20,8811)                    
                  WRITE (20,8902)  
                  WRITE (20,9002) lnum(1),cee(1),phi(1),gamr(1,1),ru(1)                                 
                CASE(2,3)  ! peizometric surface or 3d pressure file
                  WRITE (20,8812)                
                  WRITE (20,8813)                        
                  WRITE (20,8903)
                  WRITE (20,9002) lnum(1),cee(1),phi(1), gamr(1,2),gamr(1,3)
                CASE(4)    ! 3d pressure file with relative water contents
                  WRITE (20,8814)                
                  WRITE (20,8815)                   
                  WRITE (20,8904)
                  WRITE (20,9003) lnum(1),cee(1),phi(1),gamr(1,3),thetares(1),thetasat(1) 
                CASE(5) 
                  WRITE (20,8814)                
                  WRITE (20,8815)                   
                  WRITE (20,8907)
                  WRITE (20,9005) lnum(1),cee(1),phi(1),gamr(1,3),thetares(1),thetasat(1),&
                          vga(1),vgn(1)                         
                CASE(6) 
                  WRITE (20,8814)                
                  WRITE (20,8815)                   
                  WRITE (20,8905)
                  WRITE (20,9004) lnum(1),cee(1),phi(1),gamr(1,3),thetares(1),thetasat(1),&
                          fxa(1),fxn(1),fxm(1),fxr(1)                                                                                                    
              END SELECT 
            END IF                       

            IF (cflag.ne.0.and.fflag.ne.0.and.dflag.ne.0) THEN
              errmessage = 'no parameters specified for 3D material properties file' 
              solution = 'check Scoops3D manual for description of 3D material properties file'
              Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                         
              ierr = 1
            END IF
            WRITE (20,9100)
          END IF
          WRITE (20,9300) strgthfile
          IF (strgrid.eq.1) THEN
            WRITE (20,9400)
            IF (lengthunits.eq.'  ') THEN
              WRITE (20,9505) dzstr
              WRITE (20,9605) zstrmin            
            ELSE 
              WRITE (20,9500) lengthunits,dzstr
              WRITE (20,9600) lengthunits,zstrmin
            END IF
            WRITE (20,9650) strnum
            WRITE (20,9700) strmaxk
            IF (linterp.eq.1) THEN
              WRITE (20,9900) 'yes'
            ELSE
              WRITE (20,9900) 'no'
            END IF
            WRITE (20,8300) strcoords
          ELSE  !  irregular grid
            WRITE (20,9450) 
            WRITE (20,9650) strnum
            WRITE (20,9700) strmaxk
            IF (linterp.eq.1) THEN
              WRITE (20,9900) 'yes'
            ELSE
              WRITE (20,9900) 'no'
            END IF
            WRITE (20,8300) strcoords
          END IF  
          IF (cflag.eq.0) THEN
            WRITE (20,9800) ceect
            IF (ceeunits.eq.'        ') THEN
              WRITE (20,9815) mincee
              WRITE (20,9825) maxcee
            ELSE              
              WRITE (20,9810) ceeunits,mincee
              WRITE (20,9820) ceeunits,maxcee
            END IF  
          END IF
          IF (fflag.eq.0) THEN
            WRITE (20,9830) fricct
            WRITE (20,9840) minfric
            WRITE (20,9850) maxfric
          END IF
          IF (dflag.eq.0) THEN
            WRITE (20,9860) uwtct/3
            IF (uwtunits.eq.'        ') THEN 
               SELECT CASE (iwater)
                  CASE(0,1)       ! no water pressure or ru
                    WRITE (20,9875)                        
                  CASE(2,3)       ! peizometric surface or 3d pressure file
                    WRITE (20,9876)                                  
                  CASE(4,5,6)  ! 3d pressure file - we currently do not allow this case
                    WRITE (20,9877)                  
                END SELECT            
            ELSE
                SELECT CASE (iwater)
                  CASE(0,1)       ! no water pressure or ru
                    WRITE (20,9870)  uwtunits                  
                  CASE(2,3)       ! peizometric surface or 3d pressure file
                    WRITE (20,9871)   uwtunits                               
                  CASE(4,5,6)  ! 3d pressure file - we currently do not allow this case
                    WRITE (20,9872) uwtunits                 
                END SELECT                       
            END IF  
            SELECT CASE (iwater)
                  CASE(0,1)       ! no water pressure or ru
                     WRITE (20,9880) minuwt(1)                 
                  CASE(2,3)       ! peizometric surface or 3d pressure file
                     WRITE (20,9881) minuwt(2:3)                              
                  CASE(4,5,6)  ! 3d pressure file - we currently do not allow this case
                     WRITE (20,9880) minuwt(3:3)                 
            END SELECT                       

            IF (uwtunits.eq.'        ') THEN
               SELECT CASE (iwater)
                  CASE(0,1)       ! no water pressure or ru
                    WRITE (20,9895)                        
                  CASE(2,3)       ! peizometric surface or 3d pressure file
                    WRITE (20,9896)                                  
                  CASE(4,5,6)  ! 3d pressure file - we currently do not allow this case w/3d strengths
                    WRITE (20,9897)                  
                END SELECT                        
            ELSE
                SELECT CASE (iwater)
                  CASE(0,1)       ! no water pressure or ru
                    WRITE (20,9890)  uwtunits                      
                  CASE(2,3)       ! peizometric surface or 3d pressure file
                    WRITE (20,9891)   uwtunits                               
                  CASE(4,5,6)  ! 3d pressure file - we currently do not allow this case w/3d strengths
                    WRITE (20,9892) uwtunits                 
                END SELECT                                   
            END IF   
         
            SELECT CASE (iwater)
                  CASE(0,1)       ! no water pressure or ru
                     WRITE (20,9880) maxuwt(1)                 
                  CASE(2,3)       ! peizometric surface or 3d pressure file
                     WRITE (20,9881) maxuwt(2:3)                              
                  CASE(4,5,6)  ! 3d pressure file - we currently do not allow this case w/3d strengths
                     WRITE (20,9880) maxuwt(3:3)                 
            END SELECT                                 
          END IF
        ELSE  !  if specifying stratigraphic layer files
          WRITE (20,8700)
          WRITE (20,8800) nmat          
          IF (ceeunits.ne.'        ') THEN    
            SELECT CASE (iwater)
              CASE(0)       ! no water pressure
                WRITE (20,8811)                    
                WRITE (20,8901) 
                WRITE (20,8911) ceeunits,uwtunits                     
              CASE(1)      ! ru approximation
                 WRITE (20,8811)                    
                 WRITE (20,8902)  
                 WRITE (20,8911) ceeunits,uwtunits                               
              CASE(2,3)  ! peizometric surface or 3d pressure file
                 WRITE (20,8812)                
                 WRITE (20,8813)                    
                 WRITE (20,8903)
                 WRITE (20,8912) ceeunits,uwtunits,uwtunits
              CASE(4)    ! 3d pressure file with relative water contents
                WRITE (20,8814)                
                WRITE (20,8815)                   
                WRITE (20,8904)
                WRITE (20,8911) ceeunits,uwtunits                      
              CASE(5) ! 3d pressure file with vanGenuchten SWCC  
                WRITE (20,8814)                
                WRITE (20,8815)                
                WRITE (20,8907)
                WRITE (20,8911) ceeunits,uwtunits      
              CASE(6)  ! 3d pressure file with Fredlund and Xing SWCC
                WRITE (20,8814)                
                WRITE (20,8815)                  
                WRITE (20,8905)
                WRITE (20,8911) ceeunits,uwtunits                     
            END SELECT 
          ELSE            
            SELECT CASE (iwater)
              CASE(0)       ! no water pressure
                WRITE (20,8811)                 
                WRITE (20,8901)                         
              CASE(1)      ! ru approximation
                WRITE (20,8811)                    
                WRITE (20,8902)                                 
              CASE(2,3)  ! peizometric surface or 3d pressure file
                WRITE (20,8812)                
                WRITE (20,8813)                    
                WRITE (20,8903)
              CASE(4)    ! 3d pressure file with relative water contents
                WRITE (20,8814)                
                WRITE (20,8815)                  
                WRITE (20,8904)                      
              CASE(5) ! 3d pressure file with vanGenuchten SWCC  
                WRITE (20,8814)                
                WRITE (20,8815)                  
                WRITE (20,8907)    
              CASE(6)  ! 3d pressure file with Fredlund and Xing SWCC  
                WRITE (20,8814)                
                WRITE (20,8815)                 
                WRITE (20,8905)                                                    
              END SELECT                                     
 
          END IF                                       

          DO n = 1,nmat 
              SELECT CASE (iwater)
                CASE(0)       ! no water pressure
                  WRITE (20,9001) lnum(n),cee(n),phi(n),gamr(n,1)                         
                CASE(1)      ! ru approximation
                  WRITE (20,9002) lnum(n),cee(n),phi(n),gamr(n,1),ru(n)                                 
                CASE(2,3)  ! peizometric surface or 3d pressure file
                  WRITE (20,9002) lnum(n),cee(n),phi(n), gamr(n,2),gamr(n,3)
                CASE(4)    ! 3d pressure file with relative water contents                 
                     WRITE (20,9003) lnum(n),cee(n),phi(n),gamr(n,3),thetares(n),thetasat(n) 
                CASE(5) ! 3d pressure file with vanGenuchten SWCC                    
                     WRITE (20,9005) lnum(n),cee(n),phi(n),gamr(n,3),thetares(n),thetasat(n),&
                        vga(n),vgn(n)                        
                CASE(6) ! 3d pressure file with Fredlund and Xing SWCC                  
                     WRITE (20,9004) lnum(n),cee(n),phi(n),gamr(n,3),thetares(n),thetasat(n),&
                        fxa(n),fxn(n),fxm(n),fxr(n)                                                                                                                       
               END SELECT                                                 
          END DO          

          IF (nmat.gt.1) THEN
            WRITE (20,9200) stratfile
            DO n = 1,nmat-1
              WRITE (20,9210) n
              WRITE (20,9220) activelay(n)
              IF (lengthunits.eq.'  ') THEN
                WRITE (20,9235) minlayer(n)
                WRITE (20,9245) maxlayer(n)
              ELSE   
                WRITE (20,9230) lengthunits,minlayer(n)
                WRITE (20,9240) lengthunits,maxlayer(n)
              END IF  
            END DO
          END IF
        END IF

        IF (ifailsurf.eq.1) THEN
          WRITE (20,3000)
          WRITE (20,9910) failsurffile
          WRITE (20,9920) COUNT(failsurf.ne.rnull)
          IF (lengthunits.eq.'  ') THEN
            WRITE (20,9935) MINVAL(failsurf,MASK=failsurf.ne.rnull)
          ELSE
            WRITE (20,9930) lengthunits,MINVAL(failsurf,MASK=failsurf.ne.rnull)
          END IF  
          IF (lengthunits.eq.'  ') THEN
            WRITE (20,9945) MAXVAL(failsurf,MASK=failsurf.ne.rnull)
          ELSE
            WRITE (20,9940) lengthunits,MAXVAL(failsurf,MASK=failsurf.ne.rnull)
          END IF
        END IF
        WRITE (20,3000)
        WRITE (20,7200)
        SELECT CASE (iwater)
         CASE (0)
           WRITE (20,7300)
         CASE (1)
           WRITE (20,7310)
           WRITE (20,7315)
         CASE (2)
           WRITE (20,7320)
           IF (uwtunits.eq.'        ') THEN
             WRITE (20,7405) gamw
           ELSE
             WRITE (20,7400) uwtunits,gamw
           END IF  
           WRITE (20,7500) pzfile
           WRITE (20,7700) COUNT(piezo.ne.rnull)
           IF (lengthunits.eq.'  ') THEN
             WRITE (20,7755) MINVAL(piezo,MASK=piezo.ne.rnull)
           ELSE
             WRITE (20,7750) lengthunits,MINVAL(piezo,MASK=piezo.ne.rnull)
           END IF  
           IF (lengthunits.eq.'  ') THEN
             WRITE (20,7765) MAXVAL(piezo,MASK=piezo.ne.rnull)
           ELSE
             WRITE (20,7760) lengthunits,MAXVAL(piezo,MASK=piezo.ne.rnull)
           END IF
         CASE (3)
           WRITE (20,7330)
           IF (uwtunits.eq.'        ') THEN                      
             WRITE (20,7405) gamw
           ELSE
             WRITE (20,7400) uwtunits,gamw
           END IF    
           WRITE (20,7800) prfile
           IF (pgrid.eq.1) THEN  ! regularly spaced data
             WRITE (20,7900)
             IF (lengthunits.ne.'  ') THEN 
               WRITE (20,8005) delzp
             ELSE
               WRITE (20,8000) lengthunits,delzp
             END IF 
             IF (lengthunits.eq.'  ') THEN   
               WRITE (20,8105) zpmin 
             ELSE
               WRITE (20,8100) lengthunits,zpmin
             END IF  
             WRITE (20,8150) pnum
             WRITE (20,8160) pcount
             WRITE (20,8200) pmaxk
             IF (lengthunits.eq.'  ') THEN
               WRITE (20,8555) MINVAL(pressh)
               WRITE (20,8565) MAXVAL(pressh)
             ELSE
               WRITE (20,8550) lengthunits,MINVAL(pressh)
               WRITE (20,8560) lengthunits,MAXVAL(pressh)
             END IF
             WRITE (20,8300) prcoords
           ELSE  ! irregularly spaced data
             WRITE (20,8400)
             WRITE (20,8150) pnum
             WRITE (20,8160) pcount
             WRITE (20,8200) pmaxk
             IF (lengthunits.eq.'  ') THEN
               WRITE (20,8555) MINVAL(pressh)
               WRITE (20,8565) MAXVAL(pressh)
             ELSE             
               WRITE (20,8550) lengthunits,MINVAL(pressh)
               WRITE (20,8560) lengthunits,MAXVAL(pressh)
             END IF               
             WRITE (20,8300) prcoords
           END IF           
         CASE (4,5,6)
           WRITE (20,7335)
           IF (uwtunits.eq.'        ') THEN                      
             WRITE (20,7405) gamw
           ELSE
             WRITE (20,7400) uwtunits,gamw
           END IF    
           SELECT CASE (iwater)
              CASE(4)                   
                WRITE (20,7414)                       
              CASE(5)                   
                WRITE (20,7417)  
              CASE(6)                   
                WRITE (20,7415)                                                        
           END SELECT                
           WRITE (20,7810) prfile              
           IF (pgrid.eq.1) THEN  ! regularly spaced data
             WRITE (20,7900)
             IF (lengthunits.ne.'  ') THEN 
               WRITE (20,8005) delzp
             ELSE
               WRITE (20,8000) lengthunits,delzp
             END IF 
             IF (lengthunits.eq.'  ') THEN   
               WRITE (20,8105) zpmin 
             ELSE
               WRITE (20,8100) lengthunits,zpmin
             END IF  
             WRITE (20,8150) pnum
             WRITE (20,8160) pcount
             WRITE (20,8200) pmaxk
             IF (lengthunits.eq.'  ') THEN
               WRITE (20,8555) MINVAL(pressh)
               WRITE (20,8565) MAXVAL(pressh)
             ELSE
               WRITE (20,8550) lengthunits,MINVAL(pressh)
               WRITE (20,8560) lengthunits,MAXVAL(pressh)
             END IF 
             IF (iwater.eq.4) THEN         
               WRITE (20,8570) MINVAL(thetaz)
               WRITE (20,8575) MAXVAL(thetaz)             
               WRITE (20,8300) prcoords  
             ELSE
               WRITE (20,8300) prcoords  
               WRITE(20,8568)            
               WRITE (20,8570) MINVAL(thetaz)
               WRITE (20,8575) MAXVAL(thetaz)            
             END IF                         
           ELSE  ! irregularly spaced data
             WRITE (20,8400)
             WRITE (20,8150) pnum
             WRITE (20,8160) pcount
             WRITE (20,8200) pmaxk
             IF (lengthunits.eq.'  ') THEN
               WRITE (20,8555) MINVAL(pressh)
               WRITE (20,8565) MAXVAL(pressh)
             ELSE             
               WRITE (20,8550) lengthunits,MINVAL(pressh)
               WRITE (20,8560) lengthunits,MAXVAL(pressh)
             END IF  
             IF (iwater.eq.4) THEN         
               WRITE (20,8570) MINVAL(thetaz)
               WRITE (20,8575) MAXVAL(thetaz)             
               WRITE (20,8300) prcoords  
             ELSE
               WRITE (20,8300) prcoords  
               WRITE(20,8568)            
               WRITE (20,8570) MINVAL(thetaz)
               WRITE (20,8575) MAXVAL(thetaz)            
             END IF                         
           END IF           
                                 
        END SELECT
          
        IF (eq.ne.rnull) THEN
          WRITE (20,3000)
          WRITE (20,9950)
          WRITE (20,10000) eq
        ELSE ! set eq to 0 for calculations
          eq = 0.0_pr
        END IF

        WRITE (20,3000)
        WRITE (20,3750)
        IF (method.eq.'o') THEN
          WRITE (20,3760)
        ELSE
          WRITE (20,3770)
          IF (ditertemp.ne.rnull)  WRITE (20,3780) diter
        END IF          

        IF (filter.eq.1) THEN
          WRITE (20,3000)
          WRITE (20,10100)
          WRITE (20,10200)
          WRITE (20,10250) absminma
        END IF
        WRITE (20,3000)  

        IF (srch.eq.'box'.or.srch.eq.'file') THEN         
          IF (srch.eq.'box') THEN
            WRITE (20,3810) 
          ELSE  ! srch = file
            WRITE (20,3800)
          END IF
          WRITE (20,5600)
          IF ((vacriterion.eq.'V').or.(vacriterion.eq.'v')) THEN
            vacriterion='v'
            WRITE (20,5700)
            goal = vmin + .50_pr*tol
          ELSE
            IF ((vacriterion.eq.'A').or.(vacriterion.eq.'a')) THEN
              vacriterion='a'             
              WRITE (20,5750)
              goal = armin + .50_pr*tol
            ELSE
              errmessage = 'invalid value for "vacriterion"  -- area or volume must be selected' 
              solution = 'check Scoops3D manual for description of required parameters'
              IF (ierr.eq.0) Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')              
              ierr = 1
            END IF
          END IF
          IF (lengthunits.eq.'  ') THEN
            IF (armin.eq.0.0_pr.and.armax.eq.0.0_pr) THEN
              WRITE (20,5910)
            ELSE
             WRITE (20,5805) armin,armax
            END IF   
            IF (vmin.eq.0.0_pr.and.vmax.eq.0.0_pr) THEN
              WRITE (20,5915)
            ELSE           
              WRITE (20,5905) vmin,vmax
            END IF  
          ELSE  
            IF (armin.eq.0.0_pr.and.armax.eq.0.0_pr) THEN
              WRITE (20,5910)
            ELSE   
              WRITE (20,5800) lengthunits         
              WRITE (20,5802) armin,armax
            END IF   
            IF (vmin.eq.0.0_pr.and.vmax.eq.0.0_pr) THEN
              WRITE (20,5915)
            ELSE     
              WRITE (20,5900) lengthunits,vmin,vmax
            END IF     
          END IF
          WRITE (20,6000) tol
          WRITE (20,6100) 
          WRITE (20,6200) limcol  
          WRITE (20,5000)
          WRITE (20,5200)
          WRITE (20,5300) deginc
          WRITE (20,5400) degmax
          WRITE (20,5500) numdir
          WRITE (20,3900)
          WRITE (20,3910)          
          IF (lengthunits.eq.'  ') THEN
            WRITE (20,4205) zsmin
            WRITE (20,4305) zsmax            
          ELSE
            WRITE (20,4200) lengthunits,zsmin
            WRITE (20,4300) lengthunits,zsmax
          END IF
          IF (lengthunits.eq.'  ') THEN
            WRITE (20,4585) zsrchres
          ELSE
            WRITE (20,4580) lengthunits,zsrchres
          END IF
          IF (lengthunits.eq.'  ') THEN
            WRITE (20,5105) dr
          ELSE
            WRITE (20,5100) lengthunits,dr
          END IF
          WRITE (20,4310)
          WRITE (20,4312) ismin,jsmin             
          IF (srch(1:2).eq.'bo') THEN 
            WRITE (20,4314) ismax,jsmax
          END IF  
          WRITE (20,4510) nsrchres             
          IF (irefine.eq.1) THEN
             WRITE (20,4400)
             WRITE (20,4410) multres
             WRITE (20,4420) fostol*100.0
          END IF
          IF (srch(1:2).eq.'fi') THEN        
            WRITE (20,4520) sfile
            WRITE (20,4530) nxsrch,nysrch
            WRITE (20,4540) nxsrch*nysrch
            WRITE (20,4550) nsrchptini
            WRITE (20,4560) ismax,jsmax
            WRITE (20,4570) lengthunits,delxy*nsrchres        
          END IF 
        ELSE  ! srch = single
          WRITE (20,3820)
          IF (lengthunits.eq.'  ') THEN
            WRITE (20,6305) 
            WRITE (20,6400) xcen,ycen,zcen,rad
            IF (irotcen.eq.1) THEN
              WRITE (20,6306) 
              WRITE (20,6401) xcenrot,ycenrot,zcenrot
            END IF
          ELSE
            WRITE (20,6300) lengthunits,lengthunits,lengthunits,lengthunits
            WRITE (20,6400) xcen,ycen,zcen,rad
            IF (irotcen.eq.1) THEN
              WRITE (20,6301) 
              WRITE (20,6401) xcenrot,ycenrot,zcenrot
            END IF
          END IF
!          WRITE (20,6410)
          WRITE (20,*)
! convert xcen, ycen to coordinates relative to the lower left corner of the DEM          
          xcen = xcen - xll
          ycen = ycen - yll
          IF (irotcen.eq.1) THEN
            xcenrot = xcenrot - xll
            ycenrot = ycenrot - yll
          END IF                              
          IF (abs(angle).le.360.0_pr) THEN
            WRITE (20,6450) angle
          ELSE
            WRITE (20,6500) 
            WRITE (20,6600) angle
          END IF
          IF (method.eq.'f') THEN
            WRITE (20,6700)
            WRITE (20,6800) foslocal
          END IF
        END IF
                   
        WRITE (20,3000)
        WRITE (20,10400)
        IF (isqout.eq.1) WRITE (20,10500) isqout
        IF (irelfos.eq.1) WRITE (20,10510) irelfos
        IF (icritlattice.eq.1) WRITE (20,10515) icritlattice        
        IF (ilattice.eq.1) WRITE (20,10520) ilattice
        IF (isubsurf.eq.1.or.isubsurf.eq.2.or.isubsurf.eq.3) THEN
          IF (single.eq.1) THEN
             WRITE (20,10555)
          ELSE
            WRITE (20,10530) isubsurf
            WRITE (20,10540) zfrac   
            delz = zfrac*delxy
!       Compute limits for 3-D array of FOS.
!       Bottom of 3-D array is arbitrarily set to minimum DEM elevation - 0.2*relief.
            zbot = zmin - ((zmax-zmin)*0.20_pr)
!       Make sure to go above DEM z by at least 1 extra node.
            nz = int((zmax-zbot)/delz) + 2
          END IF
        END IF
        SELECT CASE (remove)
         CASE ('A','a')
           remove = 'A'
           demflag = 1           
         CASE ('L','l') 
           remove = 'L'
           demflag = 1         
         CASE ('M','m') 
           remove = 'M'
           demflag = 1
         CASE ('N','n')
           remove = 'N'
           foscut = -9999.
         CASE DEFAULT
           errmessage = 'invalid value for "remove"  -- no potential failure masses will be removed from new DEM' 
           solution = 'check Scoops3D manual for description of required parameters'
           Call WriteError(0,errmessage,problemtype,solution,'no ',0,' ')
           remove = 'N'
           foscut = -9999.              
        END SELECT 
        WRITE (20,10410) remove
        IF (remove.eq.'N') WRITE (20,10420) 
        IF (remove.eq.'A') WRITE (20,10430)
        IF (remove.eq.'M') WRITE (20,10440) 
        IF (remove.eq.'L') WRITE (20,10450)
        IF (remove.ne.'N') WRITE (20,10460) foscut       
        IF (nretro.gt.0) WRITE (20,10900) nretro
        	
        IF (ALLOCATED (phi)) DEALLOCATE (phi)
        IF (ALLOCATED (lnum)) DEALLOCATE (lnum)
        IF (ALLOCATED (gsame)) DEALLOCATE (gsame) 
        IF (ALLOCATED (minlayer)) DEALLOCATE (minlayer) 
        IF (ALLOCATED (maxlayer)) DEALLOCATE (maxlayer) 
        IF (ALLOCATED (activelay)) DEALLOCATE (activelay)     
        	
        IF (ierr.eq.1) THEN
          CLOSE (20)
          CLOSE (39)
          CLOSE (12)
          STOP
        END IF
        IF (remove.eq.'A') THEN
          filout = outputdir(1:LEN_TRIM(outputdir))//filin(bfile:nfile)//'_spheresltcut_out.txt'
          OPEN (40,STATUS='replace',FILE = filout)
          WRITE (40,13000) foscut
          IF (method.eq.'b') THEN
            fosmeth = 'Bish'
          ELSE
            fosmeth = 'Ord '
          END IF                                         
          IF (LEN_TRIM(lengthunits).eq.1) THEN
            WRITE (40,13110) lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,fosmeth
          ELSE
            IF (LEN_TRIM(lengthunits).eq.2) THEN
              WRITE (40,13120) lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,fosmeth
            ELSE 
              WRITE (40,13100) fosmeth
            END IF 
          END IF 
            
          IF (nretro.gt.0) THEN
            filout = outputdir(1:LEN_TRIM(outputdir))//filin(bfile:nfile)//'_fosretro_out.txt'
            OPEN (32,STATUS='replace',FILE = filout)
            WRITE (32,13000) foscut
          END IF  
        END IF

        CLOSE (12)       

1000    FORMAT (A)
1100    FORMAT ('    *************************************************************')
1200    FORMAT ('                               Scoops3D')
1300    FORMAT ('          3D Slope Stability Throughout a Digital Landscape')
1400    FORMAT ('                        U.S. Geological Survey')
1500    FORMAT ('                             Version: ',A50)
1600    FORMAT ('    *************************************************************',//)
1700    FORMAT ('Start date and time: ',A2,'/',A2,'/',A4,2x,A2,':',A2,':',A2)
1800    FORMAT ('Description: ',A)
1900    FORMAT ('This file: ',A,/)

2000    FORMAT (/,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
2100    FORMAT ('I. INPUT FILES:')
2200    FORMAT ('  DEM file: ',A)
2300    FORMAT ('  Search-lattice file: ',A)
2400    FORMAT ('  Piezometric surface file: ',A)
2500    FORMAT ('  3D pressure-head file: ',A)
2510    FORMAT ('  3D variably saturated file: ',A)
2600    FORMAT ('  Material layer file: ',A)
2700    FORMAT ('  3D material properties file: ',A)
2750    FORMAT ('  Failure surface file: ',A)
2760    FORMAT ('  Include area file: ',A)
2765    FORMAT ('  Surface water file: ',A)
2800    FORMAT ('  Main parameter input file: ',A)

2900    FORMAT ('II. SIMULATION PARAMETERS:')
3000    FORMAT ('------------')
3100    FORMAT ('DEM')
3200    FORMAT ('Input file for topography: ',A)
3300    FORMAT ('  Dimensions of DEM grid (x,y):                                   ',2i6)
3350    FORMAT ('  Number of cells in DEM grid:                                        ',i8)
3360    FORMAT ('  Number of non-null cells in DEM grid:                               ',i8)
3400    FORMAT ('  Horizontal resolution of DEM grid (',A2,'):                             ',f8.2)
3405    FORMAT ('  Horizontal resolution of DEM grid:                                  ',f8.2)
3500    FORMAT ('  Minimum elevation of DEM (',A2,'):                                    ',f10.3)
3505    FORMAT ('  Minimum elevation of DEM:                                         ',f10.3)
3600    FORMAT ('  Maximum elevation of DEM (',A2,'):                                    ',f10.3)
3605    FORMAT ('  Maximum elevation of DEM:                                         ',f10.3)
3700    FORMAT ('  xllcorner and yllcorner (',A2,'):                 ',2f15.3)
3705    FORMAT ('  xllcorner and yllcorner:                      ',2f15.3)
3750    FORMAT ('LIMIT-EQUILIBRIUM METHOD')
3760    FORMAT ('Analysis method (method):                                  Ordinary (Fellenius)')
3770    FORMAT ('Analysis method (method):                                               Bishop')
3780    FORMAT ('Iteration tolerance for Bishop Method,   (diter):                   ',es10.3)

3800    FORMAT ('SEARCH METHOD (srch)                                                      file')
3810    FORMAT ('SEARCH METHOD (srch)                                                       box')
3820    FORMAT ('SEARCH METHOD (srch)                                                    single')
3900    FORMAT ('SEARCH-LATTICE EXTENT AND RESOLUTION')
3910    FORMAT ('VERTICAL EXTENT AND RESOLUTION')
4200    FORMAT ('Minimum elevation of search-lattice nodes (',A2,') (zsmin):             ',f10.3)
4205    FORMAT ('Minimum elevation of search-lattice nodes (zsmin):                  ',f10.3)
4300    FORMAT ('Maximum elevation of search-lattice nodes (',A2,') (zsmax):             ',f10.3)
4305    FORMAT ('Maximum elevation of search-lattice nodes (zsmax):                  ',f10.3)

4310    FORMAT ('HORIZONTAL EXTENT AND RESOLUTION')
4312    FORMAT ('Starting search-lattice horizontal node (ismin,jsmin):              ',2i5)
4314    FORMAT ('Ending search-lattice horizontal node (ismax,jsmax):                ',2i5)

4400    FORMAT ('COARSE-TO-FINE SEARCH PARAMETERS')
4410    FORMAT ('Horizontal and vertical multiplier for initial coarse search (multres):   ',i4)
4420    FORMAT ('Search iteration tolerance - percent change F (fostol):             ',f10.4)

4510    FORMAT ('Horizontal spacing - multiple of DEM resolution (nsrchres):               ',i4)
4520    FORMAT (/,'Input file name for search lattice (horizontal extent):      ',A)
4530    FORMAT ('    Dimensions of search lattice:                                   ',2i5)
4540    FORMAT ('    Number of horizontal nodes in search lattice:                     ',i8)
4550    FORMAT ('    Number of non-null horizontal nodes in search lattice:            ',i8)
4560    FORMAT ('    Ending search-lattice node:                                     ',2i5)
!4570    FORMAT ('    Resolution of  search lattice:                                      ',f8.2,/)
4570    FORMAT ('    Search-lattice horizontal spacing (',A2,'):                         ',f10.3,/)
4580    FORMAT ('Search-lattice vertical spacing (',A2,') (zsrchres):                    ',f10.3)
4585    FORMAT ('Search-lattice vertical resolution (zsrchres):                     ',f10.3)


5100    FORMAT ('Increment amount for potential failure surface sphere radius(',A2,')(dr):',f9.3)
5105    FORMAT ('Increment amount for potential failure surface sphere radius (dr):   ',f9.3)
5000    FORMAT ('SLIP DIRECTIONS')
5200    FORMAT ('Interval to search slip directions on each side of overall ')
5300    FORMAT (' fall direction of potential failure, in degrees (degmax):          ',f10.3)
5400    FORMAT ('Increment amount for slip direction, in degrees (deginc):           ',f10.3)
5500    FORMAT ('Calculated number of slip directions tested for each lattice node:        ',i4)

5600    FORMAT ('POTENTIAL FAILURE SIZE CONTROLS')
5700    FORMAT ('Primary constraint, volume or area (vacriterion):                       Volume')
5750    FORMAT ('Primary constraint, volume or area (vacriterion):                         Area')
5800    FORMAT ('Surface area range of potential failures (',A2,'^2) (armin, armax):')
5802    FORMAT ('                                                          ',2es10.3)
5805    FORMAT ('Surface area range of potential failures (armin, armax):  ',2es10.3)
5900    FORMAT ('Volume range of potential failures (',A2,'^3)(vmin, vmax):    ',2es10.3)
5905    FORMAT ('Volume range of potential failures (vmin, vmax):          ',2es10.3)
5910    FORMAT ('Surface area is not a criterion for size restriction.')
5915    FORMAT ('Volume is not a criteria for size restriction.') 
6000    FORMAT ('Tolerance amount for initial potential failure (tol):               ',es10.3)
6100    FORMAT ('Minimum number of active columns in potential failure required,')
6200    FORMAT ('      otherwise error message generated (limcol):                        ',i5)

6300    FORMAT ('         xcen(',A2,')        ycen(',A2,')   zcen(',A2,')  radius(',A2,')')
6301    FORMAT (' xcenrot(',A2,') ycenrot(',A2,') zcenrot(',A2,')')
6305    FORMAT ('          xcen             ycen      zcen      radius')
6306    FORMAT ('       xcenrot       ycenrot     zcenrot')
6400    FORMAT (2f16.4,f10.4,f12.5)
6401    FORMAT (2f16.4,f10.4)
!6410    FORMAT ('coordinates for xcen and ycen are given as distance from lower left corner')
6450    FORMAT ('Slip direction, relative to lattice,                                ',f10.4)
6500    FORMAT ('Slip direction, relative to lattice,')
6600    FORMAT ('       if >|360| then uses estimated fall line (angle):             ',f10.4) 
6700    FORMAT ('Create output file of local F for each DEM column, ')
6800    FORMAT ('   Ordinary method only. (foslocal, yes=1, no=0):                            ',i1)

7200    FORMAT ('GROUNDWATER CONFIGURATION')
7300    FORMAT ('Groundwater method (water):                                               None')
7310    FORMAT ('Groundwater method (water):                                                 Ru')
7315    FORMAT ('  See material properties for Ru values in each material')
7320    FORMAT ('Groundwater method (water):                           Piezometric surface file')
7330    FORMAT ('Groundwater method (water):                              3D pressure-head file')
7335    FORMAT ('Groundwater method (water):                         3D variably saturated file')
7400    FORMAT ('Unit weight of fluid (',A7,') (gamw):                                ',f8.3,/)
7405    FORMAT ('Unit weight of fluid (gamw):                                          ',f8.3,/)
!7410    FORMAT ('Curve-fitting method for soil-water characteristic curve (SWCC):')
7414    FORMAT ('Water content values specified in 3D variably saturated file')
7415    FORMAT ('Method for soil-water characteristic curve (SWCC):    Fredlund and Xing (1994)')
7417    FORMAT ('Method for soil-water characteristic curve (SWCC):        van Genuchten (1980)')
7500    FORMAT (/,'Input file name for piezometric surface: ',A)
7700    FORMAT ('    Number of non-null cells in piezometric grid:                     ',i8)
7750    FORMAT ('    Minimum piezometric elevation (',A2,'):                             ',es10.3)
7755    FORMAT ('    Minimum piezometric elevation:                                  ',es10.3)
7760    FORMAT ('    Maximum piezometric elevation (',A2,'):                             ',es10.3,/)
7765    FORMAT ('    Maximum piezometric elevation:                                  ',es10.3,/)
7800    FORMAT (/,'Input file name for 3D pressure head: ',A)
7810    FORMAT (/,'Input file name for 3D variably saturated pressure head: ',/,'      ',A)
7900    FORMAT ('    Pressure heads regularly spaced with depth')
8000    FORMAT ('    Vertical spacing (',A2,'):                                           ',f8.2)
8005    FORMAT ('    Vertical spacing:                                                 ',f8.2)
8100    FORMAT ('    Lowest elevation with pressure head value (',A2,'):                ',es10.3)
8105    FORMAT ('    Lowest elevation with pressure head value:                      ',es10.3)
8150    FORMAT ('    Number of columns with 3D pressure heads:                         ',i8)
8160    FORMAT ('    Number of columns with a piezometric surface:                     ',i8)
8200    FORMAT ('    Maximum number of pressure heads in a column:                       ',i6)
8300    FORMAT ('    Data coordinate system:                                                ',A3)
8400    FORMAT ('    Pressure heads irregularly spaced with depth')
!8500    FORMAT ('    Maximum number of pressure heads below a DEM cell:                  ',i6)
8550    FORMAT ('    Minimum pressure head (',A2,'):                                     ',es10.3)
8555    FORMAT ('    Minimum pressure head:                                          ',es10.3)
8560    FORMAT ('    Maximum pressure head (',A2,'):                                     ',es10.3)
8565    FORMAT ('    Maximum pressure head:                                          ',es10.3)
8568    FORMAT ('Values calculated from SWCC: ')
8570    FORMAT ('    Minimum volumetric water content:                               ',es10.3)
8575    FORMAT ('    Maximum volumetric water content:                               ',es10.3)
8580    FORMAT ('UNIT DESCRIPTORS (used for labels in output files)')
8583    FORMAT ('lengthunits   ceeunits    gammaunits')
8586    FORMAT (5x,A2,5x,2A12)

8600    FORMAT ('MATERIAL PROPERTIES')
8700    FORMAT ('Property method:                                                         layer')
8750    FORMAT ('Property method:                                                       3D file')
8800    FORMAT ('Number of layers (nmat):                                                   ',i3)

8811    FORMAT ('                      total unit wt.  ')
8812    FORMAT ('                           unit weight  ')
8813    FORMAT ('                        part.sat.   sat.   ')
8814    FORMAT ('                       unit wt    water content ')
8815    FORMAT ('                      saturated residual saturated ')

8901    FORMAT ('lnum    cee     phi      gamt  ')
8902    FORMAT ('lnum    cee     phi      gamt      ru')
8903    FORMAT ('lnum    cee     phi      gamps   gams ')
8904    FORMAT ('lnum    cee     phi     gams    thetares thetasat ')
8905    FORMAT ('lnum    cee     phi     gams    thetares thetasat    fxa      fxn      fxm      fxr')
!8906    FORMAT ('lnum   cee        phi   gams        thetares   thetasat  gra     grn')
8907    FORMAT ('lnum    cee     phi     gams    thetares   thetasat  vga      vgn ')

8911    FORMAT ('        ',A8,'        ',A8)
8912    FORMAT ('        ',A8,'        ',2(A8,1x))

9001    FORMAT (i3,' ',f9.2,2(f8.3,' '))
9002    FORMAT (i3,' ',f9.2,3(f8.3,' '))
9003    FORMAT (i3,' ',f9.2,4(f8.3,' '))
9004    FORMAT (i3,' ',f9.2,8(f8.3,' '))
9005    FORMAT (i3,' ',f9.2,6(f8.3,' '))

!8900    FORMAT ('lnum   cee        phi    gamr(1)  gamr(2)  gamr(3)     ru')
!8910    FORMAT ('       ',A8,'         ',3(A8,1x))
!9000    FORMAT (i3,' ',f9.2,5(f8.3,' '))

9100    FORMAT ('Note: values with -1 indicate data are contained in a 3D file')
9200    FORMAT (/,'Input file prefix for material layers: ',A)
9210    FORMAT ('  Layer # ',i4)
9220    FORMAT ('    Number of non-null cells in layer:                                  ',i6)
9230    FORMAT ('    Minimum bottom elevation (',A2,'):                                  ',es10.3)
9235    FORMAT ('    Minimum bottom elevation:                                       ',es10.3)
9240    FORMAT ('    Maximum bottom elevation (',A2,'):                                  ',es10.3)
9245    FORMAT ('    Maximum bottom elevation:                                       ',es10.3)
9300    FORMAT (/,'Input file name for 3D material properties: ',A)
9400    FORMAT ('    Regularly spaced values')
9450    FORMAT ('    Irregularly spaced values')
9500    FORMAT ('    Vertical spacing (',A2,'):                                        ',es10.3)
9505    FORMAT ('    Vertical spacing:                                               ',es10.3)
9600    FORMAT ('    Lowest elevation with strength data (',A2,'):                     ',es10.3)
9605    FORMAT ('    Lowest elevation with strength data:                            ',es10.3)
9650    FORMAT ('    Number of columns with strength values:                           ',i8)
9700    FORMAT ('    Maximum number of values in a column:                               ',i6)
9800    FORMAT ('    Number of cee values read                                         ',i8)
9810    FORMAT ('    Minimum cee value (',A8,'):                                    ',f9.2)
9815    FORMAT ('    Minimum cee value:                                               ',f9.2)
9820    FORMAT ('    Maximum cee value (',A8,'):                                    ',f9.2)
9825    FORMAT ('    Maximum cee value:                                               ',f9.2)
9830    FORMAT ('    Number of phi values read                                         ',i8)
9840    FORMAT ('    Minimum phi value:                                                ',f8.3)
9850    FORMAT ('    Maximum phi value:                                                ',f8.3)
9860    FORMAT ('    Number of gamma value sets read                                   ',i8)
9870    FORMAT ('    Minimum unit weight value (',A8,'):              gamt')
9871    FORMAT ('    Minimum unit weight values (',A8,'):              gamps   gams ')
9872    FORMAT ('    Minimum unit weight value (',A8,'):              gams ')
9875    FORMAT ('    Minimum unit weight value:                           gamt')
9876    FORMAT ('    Minimum unit weight values:                           gamps   gams')
9877    FORMAT ('    Minimum unit weight value:                           gams')
9880    FORMAT ('                                                   ',f8.3)
9881    FORMAT ('                                                   ',2(f8.3,' '))
9890    FORMAT ('    Maximum unit weight value  (',A8,'):              gamt')
9891    FORMAT ('    Maximum unit weight values (',A8,'):              gamps   gams')
9892    FORMAT ('    Maximum unit weight value  (',A8,'):              gams')
9895    FORMAT ('    Maximum unit weight value:                          gamt')
9896    FORMAT ('    Maximum unit weight values:                          gamps   gams')
9897    FORMAT ('    Maximum unit weight value:                          gams')

9900    FORMAT ('    Linearly interpolate between vertical values:                          ',A3)

9905    FORMAT ('DEFINED FAILURE SURFACE')
9910    FORMAT ('Input file name for failure surface: ',A)
9920    FORMAT ('    Number of non-null cells in failure surface grid:                 ',i8)
9930    FORMAT ('    Minimum failure surface elevation (',A2,'):                         ',es10.3)
9935    FORMAT ('    Minimum failure surface elevation:                              ',es10.3)
9940    FORMAT ('    Maximum failure surface elevation (',A2,'):                         ',es10.3)
9945    FORMAT ('    Maximum failure surface elevation:                              ',es10.3)

9950    FORMAT ('EARTHQUAKE LOADING')
10000   FORMAT ('Horizontal pseudo-acceleration coefficient (dimensionless)(eq):       ',f8.3)
10100   FORMAT ('OUTPUT FILTERS ')
10200   FORMAT ('Eliminate surfaces with absolute value of m-alpha')
10250   FORMAT ('     less than (absminma)                                           ',es10.3)
10400   FORMAT ('ADDITIONAL OUTPUT FILES AND PARAMETERS')
10410   FORMAT ('Create new DEM file (remove):                                                ',A)
10420   FORMAT ('   (no surfaces removed)')
10430   FORMAT ('   (all surfaces with F<foscut removed)')
10440   FORMAT ('   (surface with minimum F<foscut removed)')
10450   FORMAT ('   (largest surface with F<foscut removed)')
10460   FORMAT ('F cutoff for removing material from new DEM (foscut):               ',es10.3)
10500   FORMAT ('isqout (search quality files):                                               ',i1)
10510   FORMAT ('irelfos (relative F file):                                                   ',i1)
10515   FORMAT ('icritlattice (3D search lattice for critical nodes):                         ',i1)
10520   FORMAT ('ilattice (3D search lattice):                                                ',i1)
10530   FORMAT ('isubsurf (3D subsurface factor of safety):                                   ',i1)      
10540   FORMAT ('zfrac:                                                                ',f8.3)
10555   FORMAT ('WARNING - isubsurf is not available for a single search')

10900   FORMAT ('# of retrogressions:                                                      ',i4)
13000   FORMAT ('Potential failure masses with F < ',f6.3)
13100   FORMAT ('     xcen            ycen           zcen      radius    angle      volume     area     F_',A4)
13110   FORMAT ('      xcen_',A1,'         ycen_',A1,'         zcen_',A1,'    radius_',A1,'    angle   volume_',&
               & A1,'^3  area_',A1,'^2  F_',A4)
13120   FORMAT ('      xcen_',A2,'       ycen_',A2,'         zcen_',A2,'   radius_',A2,'   angle  volume_',A2,&
              &'^3 area_',A2,'^2  F_',A4)                        	 

        END SUBROUTINE readin
     
     
