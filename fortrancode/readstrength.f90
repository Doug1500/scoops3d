        SUBROUTINE Readstrength(sfilin,strcoords,strnum)
        
!        This subroutine reads and processes the ASCII data file tha
        
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!        This subroutine reads and processes the ASCII data file that 
!        contains 3D strength and unit weight information.  
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!       Called in Readin
!
!       VARIABLES
!
!       cflag -- indicates whether 3-d cohesion data are provided; 0=yes,1=no
!       ceect -- total number of cohesion values read
!       cohes(nx,ny,strmaxk) -- array of 3-d cohesion data
!       coh1 -- cohesion read from input file
!       coords -- flag for coordinates of strength file
!                 1 = ijk, 2 = ijz, 3 = xyz
!       dcohdz(nx,ny,strmaxk) -- array of vertical cohesion gradients
!       ddepth -- depth between strength data elevations
!       delxy -- DEM grid spacing
!       depth -- depth from DEM to strength data elevations
!       depth1,depth2 -- thickness between strength data elevations above and 
!              below piezometric surface
!       dflag -- indicates whether 3-d unit weight data are provided; 0=yes,1=no
!       dfricdz(nx,ny,strmaxk) -- array of vertical friction gradients
!       duwtdz(nx,ny,strmaxk) -- array of vertical unit weight gradients 
!       dzstr -- z spacing of 3-d strength data
!       error -- file allocation flag
!       fflag -- indicates whether 3-d friction data are provided; 0=yes,1=no
!       filein -- strength file name
!       fric(nx*ny*strmaxk) -- friction read from input file
!       fricct -- total number of friction angles read
!       gamsame -- flag for identical partial sat and sat unit weights
!       heading -- character variable for reading input
!       headlong -- length of heading string
!       i,j -- DEM grid array location
!       ilast,jlast -- used to determine z order of data, last x and y values
!       ijorder -- flag for whether strength data is given as all z data for 
!          each i,j point (=1), or all i,j data given for each z location (=0)
!       inull -- integer null value used throughout program
!       invert -- flag indicating whether to invert z data so lowest z is 
!          associated with lowest k.  1=invert
!       ios -- status of read calls
!       ios2 -- status of file opening calls
!       iwater -- flag indicating method for modeling water pressure.
!          0=no water pressures, 1=ru approximation, 2= piezometric surface, 
!          3=input 3-d pressure file, 4=variably saturated conditions
!       ktop -- max k value at particular i,j
!       l -- moist array parameter
!       linterp -- linear interpolation flag; 1=use linear interpolation, other-use
!          nearest node for strength data.
!       maxi - maximum i value from strength file
!       maxj - maximum j value from strength file
!       maxk -- keeps track of maximum k value for each i,j
!       maxstrk(nx,ny) -- array of highest k value of strength data at each DEM cell
!       mincee,maxcee -- overall min and max cohesion
!       minfric,maxfric -- overall min and max friction angles
!       minuwt(3),maxuwt(3) -- overall min and max unit weights for each saturation
!       moist -- local value for moist
!       moist -- determines which element of soil unit weight array to use: 1-dry weight
!          2-partially saturated, 3-fully saturated.
!       numstr -- max number of strength values allowed by dimensions specified
!       nx -- number of DEM cells in x direction
!       ny -- number of DEM cells in y direction
!       piezo(nx,ny) -- array of piezometric elevations at each DEM cell
!       pzsurf -- flag for whether piezometric surface was found; 1=yes,0=no
!       rnull -- real null value used throughout program 
!       sfilin -- root name for strength files
!       strgrid -- flag for strength data grid; 1=regular, 2=irregular grid
!       strcount -- total number of 3-d strength data points
!       strnum -- number of horizontal cells with strength data given
!       strmaxk -- number of 3-d strength values at each cell
!       strz(nx,ny,strmaxk) -- array of z locations of irregular strength data
!       tempc,tempf,tempd,tempz -- temporary cohesion, friction, unit weight, and 
!          elevation arrays used to invert order of data
!       tfric(nx,ny,strmaxk) -- array of 3-d friction data
!       uwt(nx,ny,strmaxk,3) -- array of unit weights
!       uwt1(3) -- unit weights read from input file
!       uwtct -- total number of unit weights read
!       uwt3d(nx,ny,strmaxk) -- unit weight at each depth (depth weighted average)
!          with correct moist for water conditions
!       uwtsum -- sum of unit weights weighted for depth
!       x,y,z -- location parameters read from input file
!       xll,yll -- x and y origin of DEM grid, read in header lines
!       zdem(nx,ny) -- DEM elevations
!       zlast -- elevation of previous data point
!       zloc -- location of z values 1='top',2='center' or 3='bottom'
!       zstrmin -- minimum z elevation of strength data
!
!       INPUT FILES
!         unit filename
!           13 'inputfilename.#' -- file of 3D strength for each DEM cell
!
!       OUTPUT FILES
!           20 'inputfilename_out.txt' -- echo of input and results containing
!                overall minimum slip surface data, written in subroutines
!                Readin, Readpiezo, Readpressh, Readsearch, Readstrat,
!                Readstrength, and Writeout.
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        USE GridData, ONLY: nx,ny,pi,delxy,zdem,xll,yll
        USE MaterialData, ONLY: gsameall,gsameeach
        USE StrengthData
        USE WaterData, ONLY: iwater,pzsurf,piezo,moist
        
        IMPLICIT NONE
        
        INTEGER, INTENT(out) :: strnum
        INTEGER :: ilast,jlast,ijorder,maxi,maxj,ktop,i,j,k,invert,headlong,mmoist,kcount
        INTEGER :: l,strcount,error,gamsame,numstr,coords,ios,fileerror,zloc,headlinecount
        REAL(pr) :: x,y,z,xlast,ylast,depth,depth1,depth2,uwtsum
        REAL(pr) :: coh1, uwt1(3),zlast,ddepth,delzpold
        
        REAL(pr), ALLOCATABLE :: tempc(:,:,:),tempf(:,:,:),tempd(:,:,:,:) 
        REAL(pr), ALLOCATABLE :: tempz(:,:,:),fric(:)

        CHARACTER*220, INTENT(in) :: sfilin
        CHARACTER*60 :: heading
        CHARACTER*3, INTENT(out) :: strcoords 
        CHARACTER*70 :: problemtype
        CHARACTER*120 :: errmessage,solution
      
        errmessage = ' '
        solution = ' '
        problemtype = 'reading 3D material properties file'  
        maxi = 0
        maxj = 0 
        k = 1
        ktop = 0
        strgrid = 1
        strmaxk = 0
        dzstr = 0
        ilast = 0
        jlast = 0
        x = rnull
        y = rnull
        z = rnull
        headlinecount = 0 
        xlast = rnull
        ylast = rnull       
        zlast = rnull
        ijorder = -1
        zstrmin = 999999.0_pr
        zloc = 0
        invert = -1
        gamsame = 0
        minuwt = -rnull
        maxuwt = rnull
        coords = 0
        ios = 1
        fileerror = 0
        uwt1 = 0.0_pr
        coh1 = 0.0_pr
        strcount = 1
        error = 0
        headlong = 0

        OPEN (13,STATUS = 'old',FILE = sfilin,IOSTAT=ios)
        IF (ios.ne.0) THEN
          errmessage = 'error opening 3D material properties file:'
          solution = 'check for existence of file with specified name and location'
          CLOSE(13)
          Call WriteError(1,errmessage,problemtype,solution,'ch ',0,sfilin)         
        END IF

        DO 
          headlinecount = headlinecount + 1        
          READ (13,1000,IOSTAT=ios) heading
          IF (ios.ne.0) EXIT
          headlong = LEN_TRIM(heading)
          IF (heading(1:1).ne.'#'.and.headlong.gt.0) EXIT
        END DO  

        IF (ios.eq.0) THEN
          DO 
! If coordinate system header line is not after comments, then exit           
            IF (heading(1:4).ne.'coor'.and.heading(1:4).ne.'COOR'&
                .and.heading(1:4).ne.'Coor') EXIT
! Read in coordinate system, z spacing format and maximum number of z values               
            READ (13,*,IOSTAT=ios) strcoords         
            IF (ios.ne.0) EXIT              
            headlinecount = headlinecount + 1                            
            SELECT CASE (strcoords)
              CASE ('ijk','IJK')
                coords = 1
                strcoords = 'ijk'
              CASE ('ijz','IJZ')              
                coords = 2
                strcoords = 'ijz'                
              CASE ('xyz','XYZ')               
                coords = 3
                strcoords = 'xyz' 
              CASE DEFAULT
                EXIT               
            END SELECT              
! If coordinate system is ijk, read in vertical spacing
            IF (coords.eq.1) THEN
              READ (13,1000,IOSTAT=ios) heading
              headlinecount = headlinecount+1              
               IF (ios.ne.0.or.(heading(1:4).ne.'delz'.and.heading(1:4).ne.'DELZ')&
                 .and.heading(1:4).ne.'Delz') EXIT 
               READ (13,*,IOSTAT=ios) dzstr,zstrmin 
               headlinecount = headlinecount+1                                           
            END IF                                        
            EXIT
          END DO
        END IF

        solution = '3D material properties file must be in 3D file format, as described in Scoops3D manual'
        IF (coords.eq.0) THEN
          errmessage = 'error reading 3D material properties file -- header lines must specify coordinates: ijk,ijz or xyz'              
          Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')           
          error = 1
        END IF  
        
        IF (coords.eq.1.and.dzstr.eq.0) THEN
          errmessage = 'missing or incorrect header data in 3D material properties file -- "delz" is required with ijk coordinates'              
          Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')                                      
          error = 1
        END IF                    

        IF (coords.eq.1.and.zstrmin.eq.999999.0_pr) THEN 
          errmessage = 'missing or incorrect header data in 3D material properties file -- "zmin" is required with ijk coordinates'              
          Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')                                        
          error = 1
        END IF 

        IF (error.eq.1.or.ios.ne.0) THEN
          errmessage = 'error reading 3D pressure head file' 
          CLOSE (13)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')             
        END IF                 

! Read in location of z coordinates (top, center or bottom) 
        DO      
          READ (13,1000,IOSTAT=ios) heading
          IF (ios.ne.0) EXIT         
          headlinecount = headlinecount + 1   
          IF (heading(1:4).ne.'zloc'.and.heading(1:4).ne.'ZLOC'&
                 .and.heading(1:4).ne.'Zloc') THEN
            BACKSPACE(13)
            EXIT         
          END IF                     
          READ (13,1000,IOSTAT=ios) heading
          IF (ios.ne.0) EXIT 
          headlinecount = headlinecount + 1              
          SELECT CASE (heading(1:3))
            CASE ('top','TOP')          
              zloc = 1  
            CASE ('cen','CEN')            
              zloc = 2
            CASE ('bot','BOT')            
              zloc = 3
            CASE DEFAULT   
              EXIT  
          END SELECT                                     
          EXIT 
        END DO    
     
        IF (coords.eq.0) THEN
          errmessage = 'missing or incorrect header data in 3D material properties file -- "coords" is required'              
          Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')  
          error = 1
        END IF              
              
        IF (linterp.ne.1) THEN
          IF (zloc.eq.0.and.coords.ne.1) THEN 
            errmessage = 'missing or incorrect header data in 3D material properties file -- "zlocation" is required'
            Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                                         
            error = 1
          END IF      
        ELSE
          IF (zloc.ne.2.or.zloc.ne.0) THEN
            errmessage = '"zlocation" is ignored when using interpolation option'              
            Call WriteError(0,errmessage,problemtype,solution,'no ',0,' ')            
          END IF   
        END IF 
         
        IF (error.eq.1.or.ios.ne.0) THEN
          errmessage = 'error reading 3D material properties file' 
          CLOSE (13)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')          
        END IF  

!    Read  data from file first time to determine spacing, and other parameters.
        DO
! read xyz or ijk - the distinction is not important here, so the variables x,y,z are used
          READ (13,*,IOSTAT=ios)x,y,z
          IF (ios.ne.0) EXIT
! If coordinate system is not ijk, determine vertical spacing, minimum z, etc.            
          IF (coords.ne.1) THEN
            dzstr = z - zlast
! if this is not the first line, check to see if spacing is irregular    
            IF (zlast.ne.rnull) THEN       
              IF ((xlast.eq.x.and.ylast.eq.y).and.dzstr.ne.delzpold) strgrid = 2
            END IF   
            zlast = z           
            delzpold = dzstr
            IF (z.lt.zstrmin) zstrmin = z
          END IF  
          IF (x.eq.xlast.and.y.eq.ylast) THEN
            kcount = kcount+1
            IF (strmaxk.lt.kcount) strmaxk = kcount
          ELSE 
            kcount = 1
            xlast = x
            ylast = y
          END IF   
        END DO     

! check for some conditions that are not allowed with irregular spacing
        IF (linterp.ne.1.and.strgrid.eq.2) THEN
          IF (zloc.eq.2) THEN  
            errmessage = '"zlocation" = center is not valid for irregular spacing when not interpolating' 
            Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')                                         
            error = 1
          END IF
          IF (zloc.eq.1) THEN  
            errmessage = '"zlocation" = top with irregular spacing ' 
            Call WriteError(0,errmessage,problemtype,solution,'no ',0,' ')            
            PRINT *,'Scoops3D will not be able to analyze surfaces below the lowest elevation in each column'     
            WRITE (39,*) 'Scoops3D will not be able to analyze surfaces below the lowest elevation in each column'                
          END IF                      
        END IF  

        IF ((ios.ne.0.and.x.eq.rnull).or.(error.eq.1)) THEN
          errmessage = 'error reading 3D material properties file' 
          CLOSE (13)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')         
        END IF    

! adjust z elevations if using top or middle
        SELECT CASE (zloc)
          CASE (1)           
            IF (linterp.ne.1.and.strgrid.eq.1) zstrmin=zstrmin-dzstr
          CASE (2)            
            IF (linterp.ne.1.and.strgrid.eq.1) zstrmin=zstrmin-0.5_pr*dzstr
        END SELECT      
                              
! rewind file to read a second time        
        REWIND (13)
 
 
! skip through header lines when reading the file the second time      
        DO i=1,headlinecount
          READ (13,1000,IOSTAT=ios) heading
          IF (ios.ne.0) EXIT
        END DO       
        
        IF (ios.ne.0) THEN
          errmessage = 'error reading 3D material properties file' 
          CLOSE (13)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')         
        END IF            
        
        solution = 'reduce memory requirements. See "Practical Considerations" chapter of Scoops3D manual'
!     Allocate array of z data locations, and highest k array for each DEM cell.
        ALLOCATE (strz(nx,ny,strmaxk+1),maxstrk(nx,ny),STAT=error)   
        IF (error.ne.0) THEN
          errmessage = '3D material property arrays not allocated successfully' 
          CLOSE (13)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')          
        END IF
        strz = rnull   
        maxstrk = 0   

!    Determine which material parameters are included in file and allocate
!    needed arrays for cohesion, friction angle,and unit weights.
        IF (cflag.eq.0) THEN
          ALLOCATE (cohes(nx,ny,strmaxk+1),STAT=error)
!    If using linear interpolation allocate cohesion gradient array.          
          IF (linterp.eq.1) ALLOCATE (dcohdz(nx,ny,strmaxk+1),STAT=error)
          IF (error.ne.0) THEN
            errmessage = '3D material property array for cohesion not allocated successfully' 
            CLOSE (13)
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')             
          END IF
          cohes = rnull
          IF (ALLOCATED(dcohdz)) dcohdz = 0.0_pr
        END IF

        IF (fflag.eq.0) THEN
          numstr = nx*ny*strmaxk
          ALLOCATE (tfric(nx,ny,strmaxk+1),fric(numstr),STAT=error)
!    If using linear interpolation allocate friction angle gradient array. 
          IF (linterp.eq.1) ALLOCATE (dfricdz(nx,ny,strmaxk+1),STAT=error)
          IF (error.ne.0) THEN
            errmessage = '3D material property array for friction not allocated successfully' 
            CLOSE (13)
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')           
          END IF
          tfric = rnull
          fric = rnull
          IF (ALLOCATED(dfricdz)) dfricdz = 0.0_pr
        END IF

        IF (dflag.eq.0) THEN
          ALLOCATE (uwt(nx,ny,strmaxk+1,3),STAT=error)
          ALLOCATE (uwt3d(nx,ny,strmaxk+1),STAT=error) 
          ALLOCATE (duwtdz(nx,ny,strmaxk+1),STAT=error)
          IF (error.ne.0) THEN
            errmessage = '3D material property array for unit weight not allocated successfully' 
            CLOSE (13)
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')           
          END IF
          uwt = rnull
          uwt3d = rnull
          duwtdz = 0.0_pr
        END IF

!    Read specified material parameters.
        DO
          IF (cflag.eq.0) THEN
            IF (fflag.eq.0) THEN
              IF (dflag.eq.0) THEN 
                SELECT CASE (iwater)
                  CASE(0,1)       ! no water pressure or ru
                    READ (13,*,IOSTAT=ios) x,y,z,coh1,fric(strcount),uwt1(1) 
                    uwt1(2)=uwt1(1)
                    uwt1(3)=uwt1(1)                  
                  CASE(2,3)       ! peizometric surface or 3d pressure file
                    READ (13,*,IOSTAT=ios) x,y,z,coh1,fric(strcount),uwt1(2),uwt1(3)
                    IF (uwt1(3).eq.0.0_pr) uwt1(3)=uwt1(2)                                  
                END SELECT                
              ELSE  
                READ (13,*,IOSTAT=ios) x,y,z,coh1,fric(strcount)
              END IF
            ELSE
              IF (dflag.ne.0) THEN
                READ (13,*,IOSTAT=ios) x,y,z,coh1
              ELSE
                SELECT CASE (iwater)
                  CASE(0,1)       ! no water pressure, ru, partially saturated
                    READ (13,*,IOSTAT=ios) x,y,z,coh1,uwt1(1) 
                    uwt1(2)=uwt1(1)
                    uwt1(3)=uwt1(1)                  
                  CASE(2,3)       ! peizometric surface or 3d pressure file
                    READ (13,*,IOSTAT=ios) x,y,z,coh1,uwt1(2),uwt1(3)
                    IF (uwt1(3).eq.0.0_pr) uwt1(3)=uwt1(2)                                  
                END SELECT          
              END IF
            END IF
          ELSE
            IF (fflag.eq.0) THEN
              IF (dflag.eq.0) THEN
                SELECT CASE (iwater)
                  CASE(0,1)       ! no water pressure, ru, partially saturated
                    READ (13,*,IOSTAT=ios) x,y,z,fric(strcount),uwt1(1)
                    uwt1(2)=uwt1(1)
                    uwt1(3)=uwt1(1)                  
                  CASE(2,3)       ! peizometric surface or 3d pressure file
                    READ (13,*,IOSTAT=ios) x,y,z,fric(strcount),uwt1(2),uwt1(3)
                    IF (uwt1(3).eq.0.0_pr) uwt1(3)=uwt1(2)                                  
                END SELECT                
              ELSE
                READ (13,*,IOSTAT=ios) x,y,z,fric(strcount)
              END IF
            ELSE
              SELECT CASE (iwater)
                CASE(0,1)       ! no water pressure, ru, partially saturated
                  READ (13,*,IOSTAT=ios) x,y,z,uwt1(1)
                  uwt1(2)=uwt1(1)
                  uwt1(3)=uwt1(1)                  
                CASE(2,3)       ! peizometric surface or 3d pressure file
                  READ (13,*,IOSTAT=ios) x,y,z,uwt1(2),uwt1(3)
                  IF (uwt1(3).eq.0.0_pr) uwt1(3)=uwt1(2)                                  
                END SELECT
            END IF
          END IF

          IF (ios.ne.0) EXIT

          solution = '3D material properties file must be in 3D file format, as described in Scoops3D manual'
!    If necessary, calculate i,j,k, and z for each line of data.                   
          IF (coords.eq.1) THEN
            i = NINT(x)
            j = NINT(y)
            k = NINT(z)
            IF (x.ne.REAL(i).or.y.ne.REAL(j).or.k.ne.REAL(z)) THEN
                errmessage = 'Scoops3D found non-integer values for ijk coordinates' 
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')
                PRINT *,'            i = ',x,'j = ',y,'k = ',z                 
                WRITE (39,*) '            i = ',x,'j = ',y,'k = ',z 
                CLOSE (13)                
                CLOSE(20)
                CLOSE(39)
                STOP               
            END IF
            z = REAL(k-1,pr)*dzstr + zstrmin
          ELSE 
            IF (coords.eq.2) THEN
              i = NINT(x)
              j = NINT(y)
              IF (x.ne.REAL(i).or.y.ne.REAL(j)) THEN
                errmessage = 'Scoops3D found non-integer values for ij coordinates' 
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')  
                PRINT *,'            i = ',x,'j = ',y
                WRITE (39,*)  '            i = ',x,'j = ',y             
                CLOSE (13)                
                CLOSE(20)
                CLOSE(39)
                STOP               
              END IF              
            ELSE
              i = INT((x-xll)/delxy) + 1
              j = INT((y-yll)/delxy) + 1              
            END IF
            IF (strgrid.eq.1) THEN
              IF (linterp.ne.1) THEN
                IF (zloc.eq.1) z = z - dzstr
                IF (zloc.eq.2) z = z - 0.5_pr*dzstr
              END IF
              k = NINT((z-zstrmin)/dzstr) + 1 
            END IF              
          END IF

          IF (i.lt.1.or.j.lt.1.or.k.lt.1) THEN
            errmessage = 'Scoops3D found negative values for i,j,or k in 3D material properties file' 
            Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')          
            PRINT *, '         ijk = ',i,j,k
            WRITE (39,*) '        ijk = ',i,j,k  
            CLOSE (13)                            
            CLOSE (20)
            CLOSE (39)
            STOP
          END IF
          
!     If irregular grid, determine order of data, whether all z given for
!     each i and j, or all i,j for each z.      
          IF (strgrid.eq.2) THEN
            IF (ilast.ne.0.and.jlast.ne.0) THEN ! If not first line of input
              IF (ijorder.eq.-1) THEN  
!    Compare first and second line of input to determine order of data.              
                IF (i.eq.ilast.and.j.eq.jlast) THEN
!    All z values are entered at once for each i and j.
                  ijorder = 1
                  IF (z.lt.zlast) THEN
!    Data must be inverted so min z is associated with min k.                  
                    invert = 1
                  ELSE
                    invert = 0
                  END IF
                ELSE
!    All i and j are entered for each z value.
                  ijorder = 0
                END IF
              END IF
!    Determine when to increment k or restart at k=1 at new i,j cell.              
              IF (ijorder.eq.1) THEN
                IF (i.eq.ilast.and.j.eq.jlast) THEN
                  k = k + 1
                ELSE
                  k = 1
                END IF
              ELSE
                IF (z.ne.zlast) THEN
                  k = k + 1
                  IF (invert.eq.-1) THEN
                    IF (z.lt.zlast) THEN
                      invert = 1
                    ELSE
                      invert = 0
                    END IF
                  END IF
                END IF
              END IF
            END IF
          END IF
                    
          IF (maxi.lt.i) maxi = i
          IF (maxj.lt.j) maxj = j
      
          ilast = i
          jlast = j
          zlast = z     
 
          solution = 'check coordinates of 3D material properties file -- horizontal dimensions should agree with DEM'
          IF (i.gt.nx) THEN
            errmessage = 'number of 3D blocks in x-direction for 3D material properties file does not agree with DEM' 
            Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')           
            PRINT *,'maximum value for i, in 3D material properties file = ',i,'. ncols in DEM  = ',nx
            WRITE (39,*) 'maximum value for i, in 3D material properties file = ',i,'. ncols in DEM  = ',nx
            CLOSE (13)                            
            CLOSE (20)
            CLOSE (39)
            STOP
          END IF
          IF (j.gt.ny) THEN
            errmessage = 'number of 3D blocks in y-direction for 3D material properties file does not agree with DEM' 
            Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')           
            PRINT *,'maximum value for j, in 3D material properties file = ',j,'. nrows in DEM  = ',ny          
            WRITE (39,*) 'maximum value for j, in 3D material properties file = ',j,'. nrows in DEM  = ',ny 
            CLOSE (13)                
            CLOSE (20)
            CLOSE (39)
            STOP
          END IF            
           
          maxstrk(i,j) = maxstrk(i,j) + 1  
          strz(i,j,k) = z          

          IF (cflag.eq.0) THEN
             IF (cohes(i,j,k).eq.rnull) THEN
                cohes(i,j,k) = coh1
             ELSE
               errmessage = 'cohesion was previously assigned to 3D block' 
               Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')  
               PRINT *, '     i = ',i,'  j = ', j,'   k  ',k
               WRITE (39,*) '     i = ',i,'  j = ', j,'   k  ',k               
               error = 2                
             END IF
          END IF         
          IF (fflag.eq.0) THEN
             IF (tfric(i,j,k).eq.rnull) THEN
               tfric(i,j,k) = TAN(fric(strcount)*pi/180.0_pr)
             ELSE
               errmessage = 'friction was previously assigned to 3D block' 
               Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')              
               PRINT *, '     i = ',i,'  j = ', j,'   k  ',k
               WRITE (39,*) '     i = ',i,'  j = ', j,'   k  ',k
               error = 2                
             END IF            
          END IF
          IF (dflag.eq.0) THEN                 
              DO l = 1,3
               IF (uwt(i,j,k,l).eq.rnull) THEN   
                uwt(i,j,k,l) = uwt1(l)
                IF (minuwt(l).gt.uwt1(l)) minuwt(l) = uwt1(l)
                IF (maxuwt(l).lt.uwt1(l)) maxuwt(l) = uwt1(l)
               ELSE
                 errmessage = 'unit weight was previously assigned to 3D block' 
                 Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')                
                 PRINT *, '     i = ',i,'  j = ', j,'   k  ',k
                 WRITE (39,*) '     i = ',i,'  j = ', j,'   k  ',k
                 error = 2                
               END IF                       
              END DO  
           END IF           

         solution = '3D material properties file must be in 3D file format, as described in Scoops3D manual'
         IF (error.eq.2) THEN
           errmessage = 'error in 3D material properties file' 
           CLOSE (13)
           Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')          
         END IF  
         
          IF (strgrid.eq.2) THEN 
            IF (z.lt.zstrmin) zstrmin = z
          END IF
          
          strcount=strcount+1
                        
        END DO

        IF ((x.eq.rnull.or.y.eq.rnull.or.z.eq.rnull).and.ios.ne.0) THEN
          errmessage = 'error in 3D material properties file - 3D value was not defined for all three dimensions' 
          CLOSE (13)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')          
       END IF            
        
!    Find value count and min and max value of each strength parameter
        IF (cflag.eq.0) THEN
          mincee = MINVAL(cohes,MASK=cohes.ne.rnull)
          maxcee = MAXVAL(cohes,MASK=cohes.ne.rnull)
          ceect = COUNT(Mask=cohes.ne.rnull)
        END IF
        IF (fflag.eq.0) THEN
          minfric = MINVAL(fric,MASK=fric.ne.rnull)
          maxfric = MAXVAL(fric,MASK=fric.ne.rnull)
          fricct = COUNT(Mask=fric.ne.rnull) 
          DEALLOCATE (fric)   
        END IF    
        IF (dflag.eq.0) THEN
          uwtct = COUNT(MASK=uwt.ne.rnull)
        END IF
       
        solution = 'reduce memory requirements. See "Practical Considerations" chapter of Scoops3D manual'
!    If grid has low z associated with high k, invert k element of arrays.        
        IF (invert.eq.1) THEN
          ALLOCATE (tempz(nx,ny,strmaxk+1),STAT=error)
          IF (cflag.eq.0) ALLOCATE (tempc(nx,ny,strmaxk+1),STAT=error)          
          IF (fflag.eq.0) ALLOCATE (tempf(nx,ny,strmaxk+1),STAT=error)
          IF (dflag.eq.0) ALLOCATE (tempd(nx,ny,strmaxk+1,3),STAT=error)
          IF (error.ne.0) THEN
            errmessage = 'arrays not allocated successfully' 
            CLOSE(13)
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')           
          END IF
          DO j=1,ny
            DO i=1,nx
              IF (zdem(i,j).eq.rnull) CYCLE              
              DO k=1,maxstrk(i,j)    
                tempz(i,j,k) = strz(i,j,maxstrk(i,j)-(k-1))         
                IF (cflag.eq.0) THEN
                  tempc(i,j,k) = cohes(i,j,maxstrk(i,j)-(k-1))
                END IF  
                IF (fflag.eq.0) THEN
                  tempf(i,j,k) = tfric(i,j,maxstrk(i,j)-(k-1))
                END IF
                IF (dflag.eq.0) THEN
                  tempd(i,j,k,1) = uwt(i,j,maxstrk(i,j)-(k-1),1)
                  tempd(i,j,k,2) = uwt(i,j,maxstrk(i,j)-(k-1),2)
                  tempd(i,j,k,3) = uwt(i,j,maxstrk(i,j)-(k-1),3)
                END IF
              END DO
            END DO
          END DO       
          strz = tempz
          DEALLOCATE (tempz)
          IF (cflag.eq.0) cohes = tempc
          IF (fflag.eq.0) tfric = tempf
          IF (dflag.eq.0) uwt = tempd        
          IF (cflag.eq.0) DEALLOCATE (tempc)
          IF (fflag.eq.0) DEALLOCATE (tempf)
          IF (dflag.eq.0) DEALLOCATE (tempd)
        END IF              
        
!    If irregular grid, not interpolating, and zloc=top, adjust strength values to
!    next block down and remove top z values. Bottom strength values will not be used. 
        IF (strgrid.eq.2.and.linterp.ne.1.and.zloc.eq.1) THEN 
          IF (cflag.eq.0) ALLOCATE (tempc(nx,ny,strmaxk+1),STAT=error)          
          IF (fflag.eq.0) ALLOCATE (tempf(nx,ny,strmaxk+1),STAT=error)
          IF (dflag.eq.0) ALLOCATE (tempd(nx,ny,strmaxk+1,3),STAT=error)
          IF (error.ne.0) THEN
            errmessage = 'temporary arrays not allocated successfully' 
            CLOSE(13)
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')            
          END IF
          DO j=ny,1,-1
            DO i=nx,1,-1
              IF (zdem(i,j).eq.rnull) CYCLE    
              IF (maxstrk(i,j).eq.0) THEN
                errmessage = 'missing material properties data for 3D block ' 
                Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')                
                PRINT *,'       i = ',i,' j = ',j
                WRITE (39,*) '       i = ',i,' j = ',j
                CLOSE(13)
                CLOSE (20)
                CLOSE (39)
                STOP
              END IF

              maxstrk(i,j) = maxstrk(i,j)-1
              DO k = maxstrk(i,j),1,-1
                IF (cflag.eq.0) tempc(i,j,k) = cohes(i,j,k+1) 
                IF (fflag.eq.0) tempf(i,j,k) = tfric(i,j,k+1)
                IF (dflag.eq.0) THEN
                  tempd(i,j,k,1) = uwt(i,j,k+1,1)
                  tempd(i,j,k,2) = uwt(i,j,k+1,2)
                  tempd(i,j,k,3) = uwt(i,j,k+1,3)
                END IF
              END DO
            END DO
          END DO
          IF (cflag.eq.0) cohes = tempc
          IF (fflag.eq.0) tfric = tempf
          IF (dflag.eq.0) uwt = tempd        
          IF (cflag.eq.0) DEALLOCATE (tempc)
          IF (fflag.eq.0) DEALLOCATE (tempf)
          IF (dflag.eq.0) DEALLOCATE (tempd)
        END IF
           
!     Calculate depth-weighted unit weights at each strength data location.
        DO j=ny,1,-1
          DO i=nx,1,-1
            IF (zdem(i,j).eq.rnull) CYCLE    
            IF (maxstrk(i,j).eq.0) THEN
              errmessage = 'missing material properties data for column ' 
              Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')                
              PRINT *,'       i = ',i,' j = ',j
              WRITE (39,*) '       i = ',i,' j = ',j
              CLOSE(13)              
              CLOSE (20)
              CLOSE (39)
              STOP
            END IF 
            uwtsum = 0.0_pr
            ktop = maxstrk(i,j)  
! if the elevation of the topmost strength value is below the dem, 
! but within one 3D block from the top (if regular grid), add one more value
            IF (strz(i,j,ktop).lt.zdem(i,j)) THEN
              IF (strgrid.eq.1) THEN
                IF (strz(i,j,ktop).lt.zdem(i,j)-dzstr) THEN
                  errmessage = 'maximum elevation of 3D data must be less than one vertical spacing (delz) below DEM elevation' 
                  Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')
                  PRINT *, '     i =',i,'j = ',j,'DEM =  ',zdem(i,j),'max elevation of 3D data = ',strz(i,j,ktop)
                  WRITE (39,*) '       i = ',i,'j = ',j,'DEM =  ',zdem(i,j),'max elevation of 3D data = ',strz(i,j,ktop)
                  CLOSE (13)
                  CLOSE (20)
                  CLOSE (39)
                  STOP
                END IF
              END IF
              ktop = ktop + 1
              maxstrk(i,j) = maxstrk(i,j) + 1
              strz(i,j,ktop) = zdem(i,j)
!    Assign values of top 3D block               
              IF (cflag.eq.0)cohes(i,j,ktop) = cohes(i,j,ktop-1)
              IF (fflag.eq.0) tfric(i,j,ktop) = tfric(i,j,ktop-1)
              IF (dflag.eq.0) THEN
                uwt(i,j,ktop,1) = uwt(i,j,ktop-1,1)
                uwt(i,j,ktop,2) = uwt(i,j,ktop-1,2)
                uwt(i,j,ktop,3) = uwt(i,j,ktop-1,3) 
              END IF 
            END IF 

            mmoist = moist                             
            DO k = maxstrk(i,j),1,-1
!     If strength data are given for z value above DEM, set z equal to DEM.
!     Change the maximum k value to refer to the first one to fall equal to or above the DEM.            
              IF (strz(i,j,k).ge.zdem(i,j)) THEN
                IF (k.gt.1) THEN
                  IF (strz(i,j,k-1).gt.zdem(i,j)) CYCLE
                END IF
                strz(i,j,k) = zdem(i,j)
                ktop = k
              END IF
              IF (k.lt.ktop) THEN
                ddepth = strz(i,j,k+1) - strz(i,j,k)
                depth = zdem(i,j) - strz(i,j,k) 
              END IF
              IF (dflag.eq.0) THEN 
                IF (uwt(i,j,k,2).eq.uwt(i,j,k,3)) gamsame = 1
!    If no piezometric surface or all 3 unit weights are same
                IF ((pzsurf.eq.0).or.(gamsame.eq.1)) THEN
!    Find depth-weighted unit weights  
                  IF (k.eq.ktop) THEN
                    uwt3d(i,j,k) = uwt(i,j,k,mmoist)
                    uwtsum = 0
                  ELSE
                    IF (linterp.ne.1) THEN    
                      uwtsum = uwtsum + uwt(i,j,k,mmoist)*ddepth 
                    ELSE  !  linear interpolation, use average for weights.
                      uwtsum = uwtsum + ((uwt(i,j,k+1,mmoist)+uwt(i,j,k,mmoist))/2.0_pr)*ddepth
                    END IF
                    uwt3d(i,j,k) = uwtsum/depth
                  END IF
!    If piezometric surface determine when to use partially saturated or 
!    saturated unit weights based on location of water table.                 
                ELSE  
                  IF (pzsurf.eq.1) THEN             
                    IF (k.eq.ktop) THEN
                      uwtsum = 0.0_pr
!    If piezometric surface is above top strength elevation.                
                      IF (piezo(i,j).ge.strz(i,j,ktop)) THEN
                        uwt3d(i,j,k) = uwt(i,j,k,3)
                        mmoist = 3                    
                      ELSE ! Piez. surface below top strength elevation, use partial saturation.
                        uwt3d(i,j,k) = uwt(i,j,k,2)
                        mmoist = 2
                      END IF
                    ELSE 
!    If piezometric surface falls between next two strength elevations                   
                      IF (piezo(i,j).gt.strz(i,j,k).and.mmoist.eq.2) THEN
                        mmoist = 3
                        depth1 = strz(i,j,k+1) - piezo(i,j)
                        depth2 = piezo(i,j)-strz(i,j,k)
                        IF (linterp.ne.1) THEN    
                          uwtsum = uwtsum + depth1*uwt(i,j,k,2)+depth2*uwt(i,j,k,3)  
                        ELSE  !  linear interpolation, distribute average wet and dry unit weights.
                          uwtsum = uwtsum + ((uwt(i,j,k+1,2)+uwt(i,j,k,2))/2.0_pr)*depth1&
                                   + ((uwt(i,j,k+1,3)+uwt(i,j,k,3))/2.0_pr)*depth2
                        END IF
                        uwt3d(i,j,k) = uwtsum/depth 
                      ELSE  !  If piezometric surface above or below elevations of interest
                        IF (linterp.ne.1) THEN    
                          uwtsum = uwtsum + uwt(i,j,k,mmoist)*ddepth 
                        ELSE  !  linear interpolation, use average for weights.
                          uwtsum = uwtsum + ((uwt(i,j,k+1,mmoist)+uwt(i,j,k,mmoist))/2.0_pr)*ddepth
                        END IF
                        uwt3d(i,j,k) = uwtsum/depth   
                      END IF
                    END IF  !  IF (k.eq.ktop)                       
                  END IF  !  IF (pzsurf.eq.1)
                END IF  !  IF (pzsurf.eq.0) 
              END IF  !  IF (dflag.eq.0)
              IF (k.lt.ktop.and.ddepth.ne.0.0_pr) THEN
                IF (dflag.eq.0) duwtdz(i,j,k) = (uwt3d(i,j,k+1)-uwt3d(i,j,k))/ddepth
              END iF
!    If using linear interpolation find cohesion and friction gradients.              
              IF (linterp.eq.1.and.k.lt.ktop) THEN
                IF (ddepth.ne.0.0_pr) THEN
                  IF (cflag.eq.0) dcohdz(i,j,k) = (cohes(i,j,k+1)-cohes(i,j,k))/ddepth
                  IF (fflag.eq.0) dfricdz(i,j,k) = (tfric(i,j,k+1)-tfric(i,j,k))/ddepth
                END IF
              END IF
            END DO
            maxstrk(i,j) = ktop
          END DO
        END DO


        IF (ALLOCATED(uwt)) DEALLOCATE (uwt)
        
        strnum = maxi*maxj
        
        CLOSE(13)

        RETURN
        
1000    FORMAT(A)

        END
        
