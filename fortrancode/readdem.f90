        SUBROUTINE Readdem (demfile)
        
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!        This subroutine reads an ASCII data file that contains
!        DEM information
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!        Called by Readin
!
!        VARIABLES
!
!        delxy -- DEM grid resolution (delta x, delta y)
!        demfile -- DEM file name
!        error -- file allocation flag
!        heading -- character variable for reading input
!        i,j -- DEM grid array location
!        ios -- status of read calls
!        nullsame -- 'true' if the file value for no data is the same as value of
!           rnull, set in CommonData module
!        nx -- number of DEM cells in x direction 
!        ny -- number of DEM cells in y direction
!        ozb(nx,ny) -- final array of DEM after scoops < FOS cut off removed.
!        readline -- counter for number of header lines read
!        rnull -- real null value used throughout program 
!        stuff -- character variable for reading input
!        xll,yll -- x and y origin of DEM grid, read in header lines
!        xllcorner,yllcorner -- character version of these header lines
!        zdem(nx,ny) -- DEM cell elevations
!        zdemnodes(nx+1,ny+1) -- DEM node elevations (average of 4 surrounding cells)
!        zmin -- minimum elevation value
!        zmax -- maximum elevation value
!
!        INPUT FILES READ
!         unit filename
!           13 'inputfilename.#' -- file of DEM input information
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        USE CommonData
        USE GridData, ONLY: nx,ny,delxy,zdem,zmin,zmax,nbdy,cellarea,&
                            xllcorner,yllcorner,xll,yll,zdemnodes
        USE FOSData, ONLY: ozb
        
        IMPLICIT NONE

        INTEGER :: ios,readline,i,j,error,numdem,zcount
        REAL(pr) :: nodata,zd(4)
        
        LOGICAL nullsame

        CHARACTER*220, INTENT(in) :: demfile
        CHARACTER*60 :: heading
        CHARACTER*12 :: stuff
        CHARACTER*70 :: problemtype
        CHARACTER*120 :: errmessage,solution
                      
        errmessage = ' '
        solution = ' '
        problemtype = 'reading DEM file'                  
        zmin = 99999.0_pr
        zmax = 0.0_pr
        nx = 0
        ny = 0
        delxy = 0.0_pr
        nodata = 0.0_pr 
        nullsame = .TRUE.
        readline = 0
        error=0
        numdem = 0

        OPEN (13,STATUS = 'old',FILE = demfile,IOSTAT=ios)
        IF (ios.ne.0) THEN
          errmessage = 'error opening DEM file:'
          solution = 'check for existence of file with specified name and location'
          Call WriteError(1,errmessage,problemtype,solution,'ch ',0,demfile)            
        END IF
          	
        DO 
          READ (13,1000,IOSTAT=ios) heading
          IF (ios.ne.0.or.readline.eq.6) EXIT

          solution = 'DEM file must be in Esri ASCII raster format, as described in Scoops3D manual'
          SELECT CASE (heading(1:5))
        
          CASE ('ncols','NCOLS')   
            BACKSPACE 13     
            READ (13,*,IOSTAT=ios) stuff,nx
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "ncols" header line of DEM file '              
              CLOSE(13)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
            END IF

          CASE ('nrows','NROWS')   
            BACKSPACE 13     
            READ (13,*,IOSTAT=ios) stuff,ny
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "nrows" header line of DEM file '              
              CLOSE(13)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
            END IF

          CASE ('xllco','XLLCO')
            xllcorner = heading  
            BACKSPACE 13     
            READ (13,*,IOSTAT=ios) stuff,xll
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "xllcorner" header line of DEM file '              
              CLOSE(13)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
            END IF
                           
          CASE ('yllco','YLLCO')
            yllcorner = heading   
            BACKSPACE 13     
            READ (13,*,IOSTAT=ios) stuff,yll
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "yllcorner" header line of DEM file '              
              CLOSE(13)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
            END IF      

          CASE ('cells','CELLS')   
            BACKSPACE 13     
            READ (13,*,IOSTAT=ios) stuff,delxy
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "cellsize" header line of DEM file '              
              CLOSE(13)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
            END IF 
            
          CASE ('NODAT','nodat')   
            BACKSPACE 13     
            READ (13,*,IOSTAT=ios) stuff,nodata
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "nodata" header line of DEM file '              
              CLOSE(13)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
            END IF
            EXIT

          END SELECT
        END DO          
        
        IF (readline.ne.6) THEN
          errmessage = 'error reading DEM file -- missing or incorrect header data'              
          CLOSE(13)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')           
        END IF
        
        
!     Set DEM array size and initialize
        ALLOCATE (zdem(nx,ny),zdemnodes(nx+1,ny+1),STAT=error)
        IF (error.ne.0) THEN                
          errmessage = 'DEM data arrays not allocated successfully'                   
          solution = 'reduce memory requirements; see "Practical Considerations" chapter of Scoops3D manual' 
          CLOSE(13)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')          
        END IF    
        ALLOCATE (ozb(nx,ny),nbdy(nx+1,ny+1),STAT=error)
        IF (error.ne.0) THEN
          errmessage = 'data arrays not allocated successfully'                   
          solution = 'reduce memory requirements; see "Practical Considerations" chapter of Scoops3D manual' 
          CLOSE(13)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')          
        END IF    
        
        zdem = rnull
        zdemnodes = rnull
        nbdy = 1

        IF (nodata.ne.rnull) nullsame = .FALSE.

        solution = 'check format of input grid data; see "Program Input" chapter of Scoops3D manual'        
!     Read DEM values
        DO j = ny,1,-1
          READ (13,*,IOSTAT=ios) (zdem(i,j),i=1,nx)
          IF (ios.ne.0.and.i.le.nx) THEN
            errmessage = 'error reading DEM file on line number '              
            CLOSE(13)
            Call WriteError(1,errmessage,problemtype,solution,'int',ny-j+1,' ')            
          END IF
          numdem = numdem + nx
        END DO
        
        IF (numdem.ne.(nx*ny)) THEN        
          solution = 'check dimensions of input data, number of data records should equal ncols*nrows'               
          errmessage = 'DEM file dimensions are incorrect'              
          CLOSE(13)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')         
        END IF
        
!    If DEM null is not same as program null, reset DEM nulls to rnull.        
        IF (.NOT.nullsame) WHERE (zdem.eq.nodata) zdem=rnull

!     Average 4 surrounding cells to find node elevation for continuous 
!     surface approximation      
        DO j = ny+1,1,-1
          DO i = 1,nx+1
            zcount = 0
            zd = 0
            IF (i.gt.1) THEN
              IF (j.gt.1) THEN
                IF (zdem(i-1,j-1).ne.rnull) THEN
                  zd(1) = zdem(i-1,j-1)
                  zcount = zcount + 1
                END IF
              END IF
              IF (j.lt.ny+1) THEN
                IF (zdem(i-1,j).ne.rnull) THEN
                  zd(4) = zdem(i-1,j)
                  zcount = zcount + 1
                END IF
              END IF
            END IF
            IF (i.lt.nx+1) THEN
              IF (j.gt.1) THEN
                IF (zdem(i,j-1).ne.rnull) THEN
                  zd(2) = zdem(i,j-1)
                  zcount = zcount + 1
                END IF
              END IF
              IF (j.lt.ny+1) THEN
                IF (zdem(i,j).ne.rnull) THEN
                  zd(3) = zdem(i,j)
                  zcount = zcount + 1
                END IF
              END IF
            END IF
            IF (zcount.gt.0) zdemnodes(i,j) = SUM(zd)/REAL(zcount,pr)
!     If all 4 surrounding cells are nonnull, flag node as not on the boundary
            IF (zcount.eq.4) nbdy(i,j) = 0          
          END DO
        END DO   

!     Find minimum and maximum DEM elevations.        
        zmin = MINVAL(zdem,MASK=zdem.ne.rnull)
        zmax = MAXVAL(zdem,MASK=zdem.ne.rnull)
        ozb = zdem
        
        cellarea = delxy*delxy

        CLOSE(13) 

        RETURN
        
1000    FORMAT (A)


        END SUBROUTINE Readdem
